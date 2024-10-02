from EMDButils import get_emdb_entries, is_spa_entry, check_entry_level, get_subsets
import os
import configparser
import subprocess
import concurrent.futures
import mysql.connector
from datetime import datetime
from random import sample
from time import sleep
import sys
import argparse
import re

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.yaml'))
EMDB_entries_path = config['EMDB'].get('ENTRIES_PATH')
log_folder = config['EMDB'].get('LOG_PATH')
scipionProjects_path = config['SCIPION'].get('SCIPIONPROJECTS_PATH')
scipion_launcher = config['SCIPION'].get('SCIPION_LAUNCHER')
validation_server_launcher = config['EM-VALIDATION'].get('VALIDATION_SERVER_LAUNCHER')
cleanOriginalData = config['INTERMEDIATE_DATA'].getboolean('CLEAN_ORIGINAL_DATA')

def connect_to_ddbb():
    connection = mysql.connector.connect(host='localhost', user='vrs', password='', database='vrs')
    return connection

def create_ddbb_data():
    print("Creating database...")
    connection = connect_to_ddbb()
    cursor = connection.cursor()
    print("Creating table...")
    cursor.execute('''CREATE TABLE IF NOT EXISTS launch (
                        entry_id VARCHAR(11) NOT NULL,
                        version INT NOT NULL,
                        levels TEXT NOT NULL, 
                        started INT NOT NULL CHECK(started IN (0, 1)),
                        start_date INT DEFAULT NULL, 
                        finished INT NOT NULL CHECK(finished IN (0, 1)),
                        finish_date INT DEFAULT NULL,
                        failed INT NOT NULL CHECK(failed IN (0, 1)),
                        report_path TEXT DEFAULT NULL,
                        fail_reason TEXT DEFAULT NULL,
                        PRIMARY KEY(entry_id, version)
                    );''')
    connection.commit()
    connection.close()

def launcher(entry, cmd, log_file, levels):
    try:
        print("Launching", entry)
        connection = connect_to_ddbb()
        cursor = connection.cursor()

        cursor.execute("SELECT COUNT(*) FROM launch WHERE entry_id = '%s'" % entry)
        n_launchs = cursor.fetchone()[0]

        data = (entry, n_launchs+1, levels, 1, int(datetime.now().timestamp()), 0, None, 0, None, None)
        cursor.execute(
            'INSERT INTO launch (entry_id, version, levels, started, start_date, finished, finish_date, failed, report_path, fail_reason) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)',
            data)
        connection.commit()
        connection.close()

        log_file = log_file + '_' + str(n_launchs+1) + '.log'
        with open(log_file, 'w') as log_file:
            process = subprocess.Popen(cmd.split(), stdout=log_file, stderr=subprocess.PIPE, text=True)
            stderr = process.communicate()[1]
            if stderr:
                log_file.write(stderr)

        reportPath = os.path.join(scipionProjects_path, entry, 'validationReport', 'report.pdf')
        data = (1, int(datetime.now().timestamp()), 0 if process.returncode == 0 and os.path.exists(reportPath) else 1,
                reportPath if process.returncode == 0 and os.path.exists(reportPath) else None,
                stderr if process.returncode != 0 else None, entry, n_launchs+1)
        connection = connect_to_ddbb()
        cursor = connection.cursor()
        cursor.execute(
            'UPDATE launch SET finished = %s, finish_date = %s, failed = %s, report_path = %s, fail_reason = %s WHERE entry_id = %s AND version = %s',
            data)
        connection.commit()
        connection.close()
        # remove scipion project
        if cleanOriginalData and process.returncode == 0 and os.path.exists(reportPath):
            cmd = 'rm -rf %s' % os.path.join(scipionProjects_path, entry)
            subprocess.run(cmd, shell=True)
    except Exception as e:
        print("Exception:", e)

def get_EMDB_entry_subsets():
    print('Getting EMDB entries...')
    # Get all EMDB entries
    get_emdb_entries(os.path.join(EMDB_entries_path, 'EMDB_all_entries.txt'))

    with open(os.path.join(EMDB_entries_path, 'EMDB_all_entries.txt'), 'r') as input:
        all_emdb_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]

    # Keep just the spa ones
    spa_emdb_entries = []
    for entry in all_emdb_entries:
        if is_spa_entry(entry):
            spa_emdb_entries.append(entry)

    with open(os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries.txt'), 'a') as output:
        output.writelines('%s\n' % entry for entry in spa_emdb_entries)

    # Check if entries have half-maps and atomic models associated
    check_entry_level(os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries.txt'),
                      os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_level1.txt'),
                      os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_levelA.txt'),
                      os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_failures.txt'))

    # Create subsets of entries
    get_subsets(os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries.txt'),
                os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_level1.txt'),
                os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_levelA.txt'),
                os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_level0notAnot1.txt'),
                os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_level0Anot1.txt'),
                os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_level01notA.txt'),
                os.path.join(EMDB_entries_path, 'EMDB_all_spa_entries_level01A.txt'))

def get_fails():
    connection = connect_to_ddbb()
    cursor = connection.cursor()
    # get entries whose last launch failed
    cursor.execute("SELECT launch.entry_id, launch.levels FROM launch as launch INNER JOIN (SELECT entry_id, MAX(version) AS last_version FROM launch GROUP BY entry_id) AS last_launchs ON launch.entry_id = last_launchs.entry_id AND launch.version = last_launchs.last_version WHERE launch.finished = 1 AND launch.failed = 1;")
    failed_entries = cursor.fetchall()
    connection.close()
    return failed_entries

def launch(levels, n_entries, start_entry=1, random=False):
    print('Launching validations over EMDB entries...')
    # Create database for keeping track of the validations launched
    create_ddbb_data()

    if 'A' in levels:
        if '1' not in levels:
            doLevels = '0,A'
            subset = 'EMDB_all_spa_entries_level0Anot1.txt'
        else:
            doLevels = '0,1,A'
            subset = 'EMDB_all_spa_entries_level01A.txt'
    elif '1' in levels:
        doLevels = '0,1'
        subset = 'EMDB_all_spa_entries_level01notA.txt'
    else:
        doLevels = '0'
        subset = 'EMDB_all_spa_entries_level0notAnot1.txt'

    # Launch the validations over the subset in a concurrent way
    with open(os.path.join(EMDB_entries_path, subset), 'r') as input:
        emdb_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]

    if start_entry:
        emdb_entries = emdb_entries[start_entry-1:start_entry+n_entries-1]
    elif random:
        emdb_entries = sample(emdb_entries, k=n_entries)

    cmd = '%s python %s EMDBid=%s doLevels=%s'
    cmds = []
    output_files = []

    for entry in emdb_entries:
        cmds.append(cmd % (scipion_launcher, validation_server_launcher, entry, doLevels))
        output_files.append(os.path.join(log_folder, entry))

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for cmd, output_file, entry in zip(cmds, output_files, emdb_entries):
            executor.submit(launcher, entry, cmd, output_file, doLevels)
            sleep(60)

def launch_fails(exceptions=[]):
    print('Launching validations again over previous fails...')
    cmd = '%s python %s EMDBid=%s doLevels=%s'
    cmds = []
    output_files = []
    emdb_entries = []
    doLevels = []

    fails = get_fails()
    for entry, level in fails:
        if entry not in exceptions:
            emdb_entries.append(entry)
            doLevels.append(level)
            cmds.append(cmd % (scipion_launcher, validation_server_launcher, entry, level))
            output_files.append(os.path.join(log_folder, entry))

    if emdb_entries:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for cmd, output_file, entry, level in zip(cmds, output_files, emdb_entries, doLevels):
                executor.submit(launcher, entry, cmd, output_file, level)
                sleep(60)
    else:
        print('There are no failed entries to relaunch')


def launch_list(input_list, doLevels):
    # Create database for keeping track of the validations launched
    create_ddbb_data()

    if 'A' in doLevels:
        if '1' not in doLevels:
            doLevels = '0,A'
        else:
            doLevels = '0,1,A'
    elif '1' in doLevels:
        doLevels = '0,1'
    else:
        doLevels = '0'

    with open(input_list, 'r') as file:
        emdb_entries = [line.strip() for line in file.readlines()]
        print(f"Launching validations over input list: {file.__str__()}...")
        print(f"Processing the following entries: {emdb_entries}.")

        cmd = '%s python %s EMDBid=%s doLevels=%s'
        cmds = []
        output_files = []

        for entry in emdb_entries:
            cmds.append(cmd % (scipion_launcher, validation_server_launcher, entry, doLevels))
            output_files.append(os.path.join(log_folder, entry))

        with concurrent.futures.ThreadPoolExecutor() as executor:
            for cmd, output_file, entry in zip(cmds, output_files, emdb_entries):
                executor.submit(launcher, entry, cmd, output_file, doLevels)
                sleep(60)

def EMDB_pattern_validator(emdb_list):
    emdb_list = [emdb_entry for emdb_entry in emdb_list.replace(' ', '').split(',')]
    for emdb_entry in emdb_list:
        if not re.match(r'^EMD-\d{4,}$', emdb_entry):
            raise argparse.ArgumentTypeError('The EMDB entry must follow the pattern: EMD-xxxx')
    return emdb_list

def main(argv):
    parser = argparse.ArgumentParser(description='Launch validations over EMDB entries')
    main_group = parser.add_mutually_exclusive_group()
    subgroup = parser.add_mutually_exclusive_group()

    # Get EMDB entries subsets
    main_group.add_argument('--getEMDBentries', '-g', help='retrieve EMDB entries', action='store_true')

    # Launch validation over all EMDB entries
    main_group.add_argument('--launchAll', '-la', help='launch validations over all EMDB entries', action='store_true')
    parser.add_argument('--level', '-l', help='when --launchAll or launchList: which level launch', choices=['0', '0,A', '0,1', 'O,A,1'])
    parser.add_argument('--nEntries', '-n', type=int, help='when --launchAll: how many EMDB entries (i.e: 100)')
    subgroup.add_argument('--startEntry', '-start', type=int, help='when --launchAll: starting EMDB position entry from list (i.e:1)')
    subgroup.add_argument('--random', '-r', help='when --launchAll: select nEntries random entries from list', action='store_true')

    # Launch validation on failed entries
    main_group.add_argument('--launchFails', '-lf', help='repeat validations over failed EMDB entries', action='store_true')
    parser.add_argument('--exceptions', '-e', type=EMDB_pattern_validator, help='when --launchFails: list of comma separated EMDB entries you want to avoid launch (i.e: EMD-1111, EMD-4374)')

    # Launch validation over a specific list of EMDB entries
    main_group.add_argument('--launchList', '-ll', help='launch validations over a specific list of EMDB entries', action='store_true')
    parser.add_argument('--inputList', '-i', help='path to the input file containing the list of EMDB entries (it must have one entry per line). \
                        Bear in mind that all entries in the list must have the same level since only one can be entered as part of the flag --level \
                        and that will be the level specified in the database.')

    args = parser.parse_args()

    if args.getEMDBentries:
        get_EMDB_entry_subsets()

    elif args.launchAll:
        level = args.level
        n_entries = args.nEntries
        start_entry = args.startEntry
        random = args.random

        if not level:
            parser.error('--launchAll requires --level')
        if not n_entries:
            parser.error('--launchAll requires --nEntries')
        if not (start_entry or random):
            parser.error('--launchAll requires either --startEntry or --random')

        launch(level, n_entries, start_entry=start_entry, random=random)

    elif args.launchFails:
        exceptions = args.exceptions
        launch_fails(exceptions=exceptions if exceptions else [])

    elif args.launchList:
        input_list = args.inputList
        level = args.level

        if not input_list:
            parser.error('--launchList requires --inputList')
        launch_list(input_list, level)

    else:
        print('You must use a valid option. Use -h or --help to see the help.')

if __name__ == '__main__':
    main(sys.argv[1:])
