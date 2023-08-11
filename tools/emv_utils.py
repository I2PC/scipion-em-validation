import os
import sys
import argparse
import time
from utils import *


def get_chains_data(data_file):
    logger.info('Process data file %s ' % data_file)
    chains_data = []
    with open(data_file) as f:
        lines_data = f.readlines()

        chain_data = {}
        current_chain = ""
        chain_data["seqData"] = []
        current_residue = 0
        for line in lines_data:
            if (line.startswith('ATOM') or line.startswith('HETATM')):
                # read fields
                # 0....v...10....v...20....v...30....v...40....v...50....v...60....v...70....v...80
                # 0....v....|....v....|....v....|....v....|....v....|....v....|....v....|....v....|
                # ATOM      1  N   ALA A  27     -30.825  45.721 -14.623  1.00  0.00           N
                # ATOM      2  CA  ALA A  27     -31.968  44.974 -14.115  1.00  0.00           C
                # ATOM      3  C   ALA A  27     -31.582  43.527 -13.827  1.00  0.00           C
                # ATOM      4  O   ALA A  27     -31.564  42.694 -14.734  1.00  0.00           O
                # ATOM      5  CB  ALA A  27     -33.115  45.027 -15.108  1.00  0.00           C
                # ATOM      6  N   TYR A  28     -31.270  43.238 -12.563  1.00  0.00           N
                label_comp_id = line[17:20] # residue name
                label_asym_id = line[21:23] # chain name (provided by the author)
                label_seq_id = line[23:26] # residue seq number (provided by the author)
                score = line[60:66] # score mean value for the whole residue

                # skip data for atoms in the same residue
                if current_residue == label_seq_id:
                    continue

                residue_data = {
                    "resSeqName": label_comp_id,
                    "resSeqNumber": label_seq_id,
                }
                residue_data["scoreValue"] = score
                chain_data["seqData"].append(residue_data)                    
                current_residue = label_seq_id

                if current_chain != label_asym_id or line == lines_data[-1]:
                    if current_chain:
                        # save current chain data
                        chains_data.append(chain_data)
                    # start new chain data set
                    chain_data = {}
                    chain_data["name"] = label_asym_id
                    chain_data["seqData"] = []
                    current_chain = label_asym_id

    return chains_data


def get_emv_data_header(emdbId, pdbId, method, proc_date):

    logger.info('Get EMV %s data header %s %s' % (method, emdbId, pdbId))
    emv_data = {}
    emv_data["resource"] = "EMV-%s-Scores" % method
    entry_data = {
        "volume_map": emdbId,
        "atomic_model": pdbId,
        "date": proc_date
    }
    source_data = {
        "method": "%s-Score" % method,
        "citation": "Sorzano C.O.S., Vilas J.L., Ramírez-Aportela E., et al. Image processing tools for the validation of CryoEM maps .Faraday Discuss., 2022,240, 210-227.",
        "doi": "doi.org/10.1039/D2FD00059H",
        "source":
        "https://biocomp.cnb.csic.es/EMValidationService/",
    }
    entry_data["source"] = source_data
    emv_data["entry"] = entry_data
    return emv_data


def convert_2_json(emdb_Id, pdbdb_Id, method, input_file, output_file):
    if method.lower() == 'mapq':
        method = 'MapQ'
    elif method.lower() == 'daq':
        method = 'DAQ'
    else:
        print('Invalid method %s. Must be "mapq" or "daq"' % method)
        return
    print('Convert a EMV PDB-like file %s file  to a JSON file %s' % (input_file, output_file))
    
    proc_date = time.strftime('%Y-%m-%d',
                              time.gmtime(os.path.getmtime(input_file)))
    emv_data = get_emv_data_header(emdb_Id, pdbdb_Id, method, proc_date)
    emv_data["chains"] = []
    chains_data = get_chains_data(input_file)
    emv_data["chains"] = chains_data

    path = os.path.dirname(output_file)
    filename = os.path.basename(output_file)
    out_file = save_json(emv_data, path, filename)

    return out_file