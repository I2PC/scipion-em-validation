#!/usr/bin/env python
import os
import sys
import argparse
import logging
import csv
from utils import *

LOG_FILENAME = 'prepareJobBatch.log'
logging.basicConfig(filename=LOG_FILENAME,
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(levelname)s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

def main(argv):

    parser = argparse.ArgumentParser(description="Prepare VRS jobs, batch mode")
    parser.add_argument(
        "-i", "--inputFile", help="input data file, a comma-separated values (csv) file with EMDB & PDB IDs to be processed", required=True)
    parser.add_argument(
        "-p", "--process", help="number of entries to process from an input file (1 by default)", required=False)
    parser.add_argument(
        "-t", "--test", help="perform a trial run with no changes made", required=False, action='store_true')
    args = parser.parse_args()
    input_file = args.inputFile
    process_num = int(args.process) if args.process else 1
    test_only = args.test
    if test_only:
        logger.info('Performing a trial run. No permanent changes will be made')

    # read input file
    logger.info('Processing %g entries from %s' % (process_num, input_file,))
    entries = []
    with open(input_file, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            for row in reader:
                # EMD-0001	REL	2018-07-04	2020-11-25	Cryo-EM structure of bacterial RNA polymerase-sigma54 holoenzyme transcription open complex
                emdb_id = row[0]
                status = row[1]
                sattus_date = row[2]
                sattus_date2 = row[3]
                title = row[4]
                if status != 'REL':
                    logger.info(
                        'Skipped entry, not in RELEASE status: %s, %s' % (status, row) )
                else:
                    # logger.info('Added entry %s' % row)
                    entries.append(emdb_id)
        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(
                input_file, reader.line_num, e))
    num_entries = len(entries)
    logger.info('Found %g entries' % num_entries)

    if num_entries < process_num:
        logger.warning('Found %g entries (requested %g)' %
                       (num_entries, process_num))
        process_num = num_entries
    elif num_entries >= process_num:
        logger.warning('Found %s entries. Will process first %g' %
                       (num_entries, process_num))

    # prepare job to run
    for i in range(0, process_num):
        emdb_id = entries[i]
        logger.info('Processing entry %g - %s' % (i, emdb_id))
        if not test_only:
            os.system('./prepareJob.py -m %s' % emdb_id)

if __name__ == '__main__':
    main(sys.argv[1:])