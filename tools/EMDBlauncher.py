from EMDButils import get_emdb_entries, is_spa_entry, check_entry_level, get_subsets
import os
import configparser
import subprocess
import concurrent.futures
import sqlite3
from datetime import datetime
from random import sample

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.yaml'))
EMDB_entries_path = config['EMDB'].get('ENTRIES_PATH')
ddbb_path = config['EMDB'].get('SQLITE_DDBB_PATH')
log_folder = config['EMDB'].get('LOG_PATH')
scipionProjects_path = config['SCIPION'].get('SCIPIONPROJECTS_PATH')
scipion_launcher = config['SCIPION'].get('SCIPION_LAUNCHER')
validation_server_launcher = config['EM-VALIDATION'].get('VALIDATION_SERVER_LAUNCHER')

def create_ddbb_data():
    print("Creating database...")
    connection = sqlite3.connect(ddbb_path)
    cursor = connection.cursor()
    print("Creating table...")
    cursor.execute('''CREATE TABLE IF NOT EXISTS launch (
                        entry_id TEXT NOT NULL,
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
                    )''')
    connection.commit()
    connection.close()

def launcher(entry, cmd, log_file, levels):
    print("Launching", entry)
    connection = sqlite3.connect(ddbb_path)
    cursor = connection.cursor()

    cursor.execute("SELECT COUNT(*) FROM launch WHERE entry_id = '%s'" % entry)
    n_launchs = cursor.fetchone()[0]

    data = (entry, n_launchs+1, levels, 1, int(datetime.now().timestamp()), 0, None, 0, None, None)
    cursor.execute(
        'INSERT INTO launch (entry_id, version, levels, started, start_date, finished, finish_date, failed, report_path, fail_reason) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
        data)
    connection.commit()

    with open(log_file, 'w') as log_file:
        process = subprocess.Popen(cmd.split(), stdout=log_file, stderr=subprocess.PIPE, text=True)
        stderr = process.communicate()[1]
        if stderr:
            log_file.write(stderr)

    data = (1, int(datetime.now().timestamp()), 0 if process.returncode == 0 else 1,
            os.path.join(scipionProjects_path, entry, 'validationReport', 'report.pdf') if process.returncode == 0 else None,
            stderr if process.returncode != 0 else None, entry)
    cursor.execute(
        'UPDATE launch SET finished = ?, finish_date = ?, failed = ?, report_path = ?, fail_reason = ? WHERE entry_id = ?',
        data)
    connection.commit()
    connection.close()


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

# Create database for keeping track of the validations launched
create_ddbb_data()

# For every subset of entries, launch the validations over them in a concurrent way
subsets = ['EMDB_all_spa_entries_level0notAnot1.txt', 'EMDB_all_spa_entries_level0Anot1.txt', 'EMDB_all_spa_entries_level01notA.txt', 'EMDB_all_spa_entries_level01A.txt']
doLevels = ['0', '0,A', '0,1', '0,1,A']
start_entry = 1  # start entry
random = False  # random selection of entries
n = 5  # amount of test entries

i = 0
for subset in subsets:
    with open(os.path.join(EMDB_entries_path, subset), 'r') as input:
        emdb_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]

    if start_entry:
        emdb_entries = emdb_entries[start_entry-1:start_entry+n-1]
    elif random:
        emdb_entries = sample(emdb_entries, k=n)

    cmd = '%s python %s EMDBid=%s doLevels=%s'
    cmds = []
    output_files = []

    for entry in emdb_entries:
        cmds.append(cmd % (scipion_launcher, validation_server_launcher, entry, doLevels[i]))
        output_files.append(os.path.join(log_folder, entry + '.log'))

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for cmd, output_file, entry in zip(cmds, output_files, emdb_entries):
            executor.submit(launcher, entry, cmd, output_file, doLevels[i])
    i+=1