import time
import tempfile
import argparse
import logging
import numpy as np
import pandas as pd
import sys, os, subprocess
from datetime import datetime
import xml.etree.ElementTree as ET
from multiprocessing import Pool, cpu_count

logging.basicConfig(level=logging.INFO)

def create_symlink_to_output_dir(output_dir):
    temp_dir = tempfile.mkdtemp()
    symlink_path = os.path.join(temp_dir, 'output_symlink')
    
    if os.path.islink(symlink_path) or os.path.exists(symlink_path):
        os.unlink(symlink_path)
                    
    os.symlink(output_dir, symlink_path)
    
    return symlink_path

def remove_symlink(symlink_path):
    if os.path.islink(symlink_path):
        os.unlink(symlink_path)
    os.rmdir(os.path.dirname(symlink_path))
    
def meta_parser(meta_file):
    with open(meta_file, 'r') as f:
        meta = f.readlines()
    meta_dict = {}
    for line in meta:
        key, value = line.strip().split(":", 1)
        meta_dict[key.strip()] = value.strip()
    return meta_dict
    
def fold_with_pulsarx(meta_dict, output_dir, cand_file):
    # nsamples = int(meta_dict['Nsamples'])
    subint_length = int(meta_dict['SubintLength'])
    nsubband = int(meta_dict['Nsubband'])
    clfd_q_value = float(meta_dict['ClfdQValue'])
    Telescope = meta_dict['Telescope']
    template = meta_dict['Template']
    beam_name = meta_dict['BeamName']
    beam_id = meta_dict['BeamId']
    utc_beam = meta_dict['UTCBeam']
    filterbank_file = meta_dict['FilterbankFile']
    tstart = meta_dict['StartTime']
    fast_nbins = meta_dict['Fast_nbins']
    slow_nbins = meta_dict['Slow_nbins']
    pulsarx_threads = meta_dict['Threads']
    nbins_string = "-b {} --nbinplan 0.1 {}".format(fast_nbins, slow_nbins)
    
    if 'ifbf' in beam_name:
        beam_tag = "--incoherent -i {}".format(int(beam_name.strip("ifbf")))
    elif 'cfbf' in beam_name:
        beam_tag = "-i {}".format(int(beam_name.strip("cfbf")))
    else:
        beam_tag = ""

    cand_file_name = os.path.basename(cand_file)
    logging.info(f"folding {cand_file_name}")
    output_rootname = os.path.join(output_dir, f"{cand_file_name.split('_')[0]}_{cand_file_name.split('_')[-1].split('.')[0]}")
    logging.info(f"Output path: {output_rootname}")
    os.makedirs(output_rootname, exist_ok=True)
    #print output with the cand_file
    print("Processing cand_file: ", cand_file)
    
    script = "psrfold_fil --plotx -v -t {} --candfile {} -n {} {} {} --template {} --clfd {} -L {} -f {} -o {} --pepoch {}".format(
        pulsarx_threads, cand_file, nsubband, nbins_string, beam_tag, template, clfd_q_value, subint_length, filterbank_file, output_rootname, tstart)
    print(script)
    subprocess.check_output(script, shell=True)

def main():
    parser = argparse.ArgumentParser(description='Fold all candidates from Peasoup xml file')
    parser.add_argument('-o', '--output_path', help='Output path to save results',  default=os.getcwd(), type=str)
    parser.add_argument('-meta', help='Meta file', type=str, required=True)
    parser.add_argument('-cands', help='Candidate file', type=str, required=True)
    args = parser.parse_args()

    if not args.meta:
        print("You need to provide a meta file")
        sys.exit()
        
    if not args.cands:
        print("You need to provide a candidate file")
        sys.exit()
    
    start_time = time.time()
    
    meta_file = args.meta
    cand_file = args.cands
    output_dir = create_symlink_to_output_dir(args.output_path)
    logging.info(f"symlink path : {output_dir}")
    meta_dict = meta_parser(meta_file)
    fold_with_pulsarx(meta_dict, output_dir, cand_file)
    
    logging.info(f"Total time taken: {time.time() - start_time} seconds")
    
    remove_symlink(output_dir)

if __name__ == "__main__":
    main()
