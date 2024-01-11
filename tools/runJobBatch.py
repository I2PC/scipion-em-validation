#!/usr/bin/env python
import os
import sys
import argparse
import csv
from utils import *

# runJob_EMD-34897.sh
FILENAME_PATTERN = 'runJob_EMD-*.sh'

def main(argv):

    parser = argparse.ArgumentParser(
        description="Run VRS jobs, batch mode")
    parser.add_argument(
        "-i", "--inputDir", help="input data directory, this shoud be the 'uploads' dir after running the script 'prepareJobBatch.py'", required=True)
    parser.add_argument(
        "-p", "--process", help="number of entries to process from the input directory (1 by default)", required=False)
    parser.add_argument(
        "-l", "--logFile", help="log file. By default 'runJobBatch.log' in a dedicated 'logs' folder.", required=False)
    parser.add_argument(
        "-t", "--test", help="perform a trial run with no changes made", required=False, action='store_true')
    args = parser.parse_args()
    if args.logFile:
        logFile = args.logFile
        logsDir = os.path.dirname(logFile)
        logFilename = os.path.basename(logFile)
    else:
        logsDir = os.path.join(PATH_TOOLS_DIR, DIR_TOOLS_LOGS)
        logFilename = os.path.join(
            logsDir, os.path.splitext(Path(__file__).name)[0] + '.log')
    logger = logSetup(__name__, logsDir, logFilename)

    input_dir = args.inputDir
    process_num = int(args.process) if args.process else 1
    test_only = args.test
    if test_only:
        logger.warning('Performing a trial run. No permanent changes will be made')

    # get list of prepared jobs from input dir
    logger.info('Processing %g entries from %s' % (process_num, input_dir,))
    entries = []

    p = Path(input_dir)
    if not p.is_dir() or not p.exists():
        # print("Input path is not a directory. Can't continue.")
        logger.error("Input path not found or is not a directory. Can't continue.")
        exit(1)

    entries = getFilesInPath(input_dir, FILENAME_PATTERN, recursive= True)
    num_entries = len(entries)
    logger.info('Found %g entries' % num_entries)
    if num_entries < process_num:
        logger.warning('Found %g entries (requested %g)' %
                       (num_entries, process_num))
        process_num = num_entries
    elif num_entries >= process_num:
        logger.warning('Found %s entries. Will process first %g' %
                       (num_entries, process_num))

    # run num='process' jobs sequencially
    for i in range(0, process_num):

        entry_path = entries[i]
        logger.info('Processing entry %g - %s' % (i, entry_path))
        if not test_only:
            # - mark entry as running
            # -- write 'flag' file
            flag_path = entry_path.parent / 'flag.running'
            try:
                # fail if flag file already exists, job still running
                flag_path.touch(exist_ok=False)
            except Exception as ex:
                logger.warning("Found %s. This entry is still running, can't continue." % flag_path)
                break
        
            # run prepared job
            command = 'bash %s' % entry_path
            logger.info('Run entry %g ->  %s' % (i, command))
            return_value = os.system(command)
            # - mark entry as finished
            # -- delete 'running flag' file
            flag_path.unlink()

            if return_value == 0:
                # -- write 'finished flag' file
                flag_path = entry_path.parent / 'flag.success'
                flag_path.touch(exist_ok=True)
            else:
                # -- write 'error flag' file
                flag_path = entry_path.parent / 'flag.error'
                flag_path.touch(exist_ok=True)


if __name__ == '__main__':
    main(sys.argv[1:])
