import re
import time
import tempfile
import argparse
import logging
import numpy as np
import pandas as pd
import sys, os, subprocess
from datetime import datetime
import xml.etree.ElementTree as ET

logging.basicConfig(level=logging.INFO)

MASKS = {
    1: 'kadane 4 8 zdot zap 1430 1490 zap 1550 1560 zap 1620 1630 zap 1750 1860',
    2: 'kadane 4 8 zdot zap 1980 2000 zap 2100 2150 zap 2190 2200 zap 2238 2242 zap 2250 2280 zap 2300 2335 zap 2400 2480 zap 2500 2600',
    3: 'kadane 4 8 zdot zap 2600 3000 zap 3930 4010 zap 3650 3660 zap 3680 3700',
    4: 'kadane 4 8 zdot zap 4240 4350 zap 4230 4240 zap 5220 5250',
    5: 'kadane 4 8 zdot'
}

def define_mask(band):
    """Define the mask based on the band number."""
    if band:
        return MASKS.get(band)
    else:
        return 'zdot'

def generate_pulsarX_cand_file(df, cands_dir, tsamp, fft_size, beam_name, GC, dm_name, cands_per_node=96):
    cand_dms = df['dm'].values
    cand_accs = df['acc'].values
    cand_period = df['period'].values
    pdot = a_to_pdot(cand_period, cand_accs)
    # Prepfold is used here on purpose. Period modified to beginning of tobs with epoch pointing to tstart
    cand_mod_period_beginning_tobs = period_correction_for_prepfold(cand_period, pdot, tsamp, fft_size)
    cand_mod_frequencies = 1/cand_mod_period_beginning_tobs
    cand_snrs = df['snr'].values
    cand_file_name = 'pulsarx_candfile.cands' 
    os.makedirs(cands_dir, exist_ok=True)
    cand_file_path = os.path.join(cands_dir, cand_file_name)
    number_of_candidates = len(cand_mod_frequencies)
    with open(cand_file_path, 'w') as f:
        f.write("#id DM accel F0 F1 S/N\n")
        for i in range(number_of_candidates):
            f.write("%d %f %f %f 0 %f\n" % (i, cand_dms[i], cand_accs[i], cand_mod_frequencies[i], cand_snrs[i]))
    
    cand_files = split_cand_file(cand_file_path, number_of_candidates, cands_per_node, cands_dir, beam_name, GC, dm_name)
    logging.info(f"{len(cand_files)} Candidate files created")
    return cand_files
    
def split_cand_file(cand_file_path, number_of_candidates, cands_per_node, cands_dir, beam_name, GC, dm_name):
    # Split the candidates into multiple candfiles
    cand_files = []
    num_files = number_of_candidates // cands_per_node + 1
    candidates_df = pd.read_csv(cand_file_path, sep=' ', header=0)
        
    for i in range(num_files):
        cand_file = f"{GC}_{beam_name}_{dm_name}_{i+1}.candfile"
        cand_file_path = os.path.join(cands_dir, cand_file)
        cand_files.append(cand_file_path)
        start_idx = i * cands_per_node
        end_idx = (i + 1) * cands_per_node
        cands = candidates_df.iloc[start_idx:end_idx].reset_index(drop=True)
        logging.info(f"Writing candidates {start_idx} to {end_idx} to a new candfile")
        try:
            with open(cand_file_path, 'w') as f:
                f.write('#id DM accel F0 F1 S/N\n')
                for index, row in cands.iterrows():
                    f.write(f"{index} {row['DM']} {row['accel']} {row['F0']} 0 {row['S/N']}\n")
        
        except IOError as e:
            logging.error(f"Error writing candidate files {e}")
            return None

    return cand_files
    
def period_correction_for_pulsarx(p0, pdot, no_of_samples, tsamp, fft_size):
    return p0 - pdot * float(fft_size - no_of_samples) * tsamp / 2

def period_correction_for_prepfold(p0,pdot,tsamp,fft_size):
    return p0 - pdot*float(fft_size)*tsamp/2

def a_to_pdot(P_s, acc_ms2):
    LIGHT_SPEED = 2.99792458e8                 # Speed of Light in SI
    return P_s * acc_ms2 /LIGHT_SPEED

def main():
    parser = argparse.ArgumentParser(description='Parse XML file to candfile')
    parser.add_argument('-i', '--input_file', help='Name of the input xml file', type=str)
    parser.add_argument('-dm_name', help='Name of the dm file', type=str)
    parser.add_argument('-o', '--output_path', help='Output path to save results',  default=os.getcwd(), type=str)
    parser.add_argument('-n', '--nh', help='Filter candidates with nh value', type=int, default=0)
    parser.add_argument('-f', '--fold_technique', help='Technique to use for folding (presto or pulsarx)', type=str, default='pulsarx')
    parser.add_argument('-sub', '--subint_length', help='Subint length (s). Default is tobs/64', type=int, default=None)
    parser.add_argument('-cpn', '--cands_per_node', help='Number of candidates to fold per node', type=int, default=96)
    parser.add_argument('-nsub', '--nsubband', help='Number of subbands', type=int, default=64)
    parser.add_argument('-clfd', '--clfd_q_value', help='CLFD Q value', type=float, default=2.0)
    parser.add_argument('-fn', '--fast_nbins', help='High profile bin limit for slow-spinning pulsars', type=int, default=128)
    parser.add_argument('-sn', '--slow_nbins', help='Low profile bin limit for fast-spinning pulsars', type=int, default=64)
    parser.add_argument('-c', '--chan_mask', help='Peasoup Channel mask file to be passed onto pulsarx', type=str, default='')
    args = parser.parse_args()

    if not args.input_file:
        print("You need to provide an xml file to read")
        sys.exit()
    
    start_time = time.time()
    
    xml_file = args.input_file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    peasoup_params = root[0]
    header_params = root[1]
    search_params = root[2]
    candidates = root[6]
    utc_beam = str(peasoup_params.find("utc_datetime").text + ":00")
    filterbank_file = str(search_params.find("infilename").text)
    GC = filterbank_file.split('/')[-1].split('_')[0]
    beam_name = ("Band" + re.sub(r'.*_Band_(\d+).*', r'\1', filterbank_file))
    band = int(re.sub(r'.*_Band_(\d+).*', r'\1', filterbank_file))
    rfi_filter = define_mask(band)
    tsamp = float(header_params.find("tsamp").text)
    fft_size = int(search_params.find("size").text)
    nsamples = int(root.find("header_parameters/nsamples").text)
    tstart = float(header_params.find("tstart").text)
    ignored_entries = ['candidate', 'opt_period', 'folded_snr', 'byte_offset', 'is_adjacent', 'is_physical', 'ddm_count_ratio', 'ddm_snr_ratio']
    rows = []
    for candidate in candidates:
        cand_dict = {}
        for cand_entry in candidate.iter():
            if not cand_entry.tag in ignored_entries:
                cand_dict[cand_entry.tag] = cand_entry.text
        cand_dict['cand_id_in_file'] = candidate.attrib.get("id")
        rows.append(cand_dict)

    df = pd.DataFrame(rows)
    df = df.astype({"snr": float, "dm": float, "period": float, "nh":int, "acc": float, "nassoc": int})
    df = df[df['nh'] >= args.nh]
    
    try:
        if args.fold_technique == 'pulsarx':
            if args.subint_length is None:
                subint_length = max(1, int(nsamples * tsamp / 64))
            else:
                subint_length = args.subint_length
            
            cand_files = generate_pulsarX_cand_file(df, args.output_path, tsamp, fft_size, beam_name, GC, args.dm_name, args.cands_per_node)
            # Create meta file
            meta_file = f"{args.output_path}/{GC}_{beam_name}_{args.dm_name}_meta.txt"
            with open(meta_file, 'w') as f:
                f.write(f"Ncandidates: {len(df)}\n")
                f.write(f"Nsamples: {nsamples}\n")
                f.write(f"SubintLength: {subint_length}\n")
                f.write(f"Nsubband: {args.nsubband}\n")
                f.write(f"ClfdQValue: {args.clfd_q_value}\n")
                f.write(f"RFIFilter: {rfi_filter}\n")
                f.write(f"BeamName: {beam_name}\n")
                f.write(f"Band: {band}\n")
                f.write(f"UTCBeam: {utc_beam}\n")
                f.write(f"Candidates_per_node: {args.cands_per_node}\n")
                f.write(f"FilterbankFile: {filterbank_file}\n")
                f.write(f"StartTime: {tstart}\n")
                f.write(f"Numberofharmonics: {args.nh}\n")
                f.write(f"Fast_nbins: {args.fast_nbins}\n")
                f.write(f"Slow_nbins: {args.slow_nbins}\n")
                # add fold template. 
            logging.info(f"Done")
        else:
            subint_length = args.subint_length
            logging.info(f"Prepfold is not available, exiting.")
            sys.exit()
    finally:
        end_time = time.time()
        print(f"Total time taken: {end_time-start_time} seconds")

if __name__ == "__main__":
    main()
    