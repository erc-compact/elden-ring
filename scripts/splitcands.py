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

def generate_pulsarX_cand_file(df, cands_dir, tsamp, fft_size, chunk_id, beam_name, source_name_prefix, cands_per_node=96):
    cand_dms = df['dm'].values
    cand_accs = df['acc'].values
    cand_period = df['period'].values
    pdot = a_to_pdot(cand_period, cand_accs)
    # Prepfold is used here on purpose. Period modified to beginning of tobs with epoch pointing to tstart
    cand_mod_period_beginning_tobs = period_correction_for_prepfold(cand_period, pdot, tsamp, fft_size)
    cand_mod_frequencies = 1/cand_mod_period_beginning_tobs
    cand_snrs = df['snr'].values
    cand_file_name = f"{source_name_prefix}_{beam_name}_ck{chunk_id}_allCands.txt"
    os.makedirs(cands_dir, exist_ok=True)
    cand_file_path = os.path.join(cands_dir, cand_file_name)
    number_of_candidates = len(cand_mod_frequencies)
    with open(cand_file_path, 'w') as f:
        f.write("#id DM accel F0 F1 S/N\n")
        for i in range(number_of_candidates):
            f.write("%d %f %f %f 0 %f\n" % (i, cand_dms[i], cand_accs[i], cand_mod_frequencies[i], cand_snrs[i]))
    
    cand_files = split_cand_file(cand_file_path, number_of_candidates, cands_per_node, cands_dir, chunk_id, beam_name, source_name_prefix)
    logging.info(f"{len(cand_files)} Candidate files created")
    return cand_files
    
def split_cand_file(cand_file_path, number_of_candidates, cands_per_node, cands_dir, chunk_id, beam_name, source_name_prefix):
    # Split the candidates into multiple candfiles
    cand_files = []
    num_files = number_of_candidates // cands_per_node + 1
    candidates_df = pd.read_csv(cand_file_path, sep=' ', header=0)
        
    for i in range(num_files):
        cand_file = f"{source_name_prefix}_{beam_name}_ck{chunk_id}_{i+1}.candfile"
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
    # parser = argparse.ArgumentParser(description='Parse XML file to candfile')
    # parser.add_argument('-i', '--input_file', help='Name of the input xml file', type=str)
    # parser.add_argument('-dm_name', help='Name of the dm file', type=str)
    # parser.add_argument('-o', '--output_path', help='Output path to save results',  default=os.getcwd(), type=str)
    # parser.add_argument('-n', '--nh', help='Filter candidates with nh value', type=int, default=0)
    # parser.add_argument('-f', '--fold_technique', help='Technique to use for folding (presto or pulsarx)', type=str, default='pulsarx')
    # parser.add_argument('-sub', '--subint_length', help='Subint length (s). Default is tobs/64', type=int, default=None)
    # parser.add_argument('-cpn', '--cands_per_node', help='Number of candidates to fold per node', type=int, default=96)
    # parser.add_argument('-nsub', '--nsubband', help='Number of subbands', type=int, default=64)
    # parser.add_argument('-clfd', '--clfd_q_value', help='CLFD Q value', type=float, default=2.0)
    # parser.add_argument('-fn', '--fast_nbins', help='High profile bin limit for slow-spinning pulsars', type=int, default=128)
    # parser.add_argument('-sn', '--slow_nbins', help='Low profile bin limit for fast-spinning pulsars', type=int, default=64)
    # parser.add_argument('-c', '--chan_mask', help='Peasoup Channel mask file to be passed onto pulsarx', type=str, default='')
    # args = parser.parse_args()
    start_time = time.time()
    parser = argparse.ArgumentParser(description='Fold all candidates from Peasoup xml file')
    parser.add_argument('-o', '--output_path', help='Output path to save results',  default=os.getcwd(), type=str)
    parser.add_argument('-i', '--input_file', help='Name of the input csv file', type=str)
    parser.add_argument('-fil', '--input_fil_file', help='Name of the input filterbank file', type=str)
    parser.add_argument('-m', '--mask_file', help='Mask file for prepfold', type=str)
    parser.add_argument('-t', '--fold_technique', help='Technique to use for folding (presto or pulsarx)', type=str, default='presto')
    parser.add_argument('-n', '--nh', help='Filter candidates with nh value', type=int, default=0)
    parser.add_argument('-f', '--fast_nbins', help='High profile bin limit for slow-spinning pulsars', type=int, default=128)
    parser.add_argument('-s', '--slow_nbins', help='Low profile bin limit for fast-spinning pulsars', type=int, default=64)
    parser.add_argument('-sub', '--subint_length', help='Subint length (s). Default is tobs/64', type=int, default=None)
    parser.add_argument('-nsub', '--nsubband', help='Number of subbands', type=int, default=64)
    parser.add_argument('-b', '--beam_name', help='Beam name string', type=str, default='cfbf00000')
    parser.add_argument('--telescope', help='Telescope name', type=str, default='meerkat')
    parser.add_argument('--beam_id', help='Beam ID', type=int, default=0)
    parser.add_argument('-utc', '--utc_beam', help='UTC beam name string', type=str, default='2024-01-01-00:00:00')
    parser.add_argument('-c', '--channel_mask', help='Peasoup Channel mask file to be passed onto pulsarx', type=str, default=None)
    parser.add_argument('-threads', '--pulsarx_threads', help='Number of threads to be used for pulsarx', type=int, default='24')
    parser.add_argument('-pthreads', '--presto_threads', help='Number of threads to be used for prepfold', type=int, default='12')
    parser.add_argument('-p', '--pulsarx_fold_template', help='Fold template pulsarx folder', type=str, default='meerkat_fold.template')
    parser.add_argument('-ncands', '--no_cands_to_fold', help='Number of candidates to fold. Folding starts with the highest SNR candidates', type=int, default=None)
    parser.add_argument('-clfd', '--clfd_q_value', help='CLFD Q value', type=float, default=2.0)
    parser.add_argument('-cpn', '--cands_per_node', help='Number of candidates to fold per node', type=int, default=96)
    parser.add_argument('--snr_min', help='Minimum SNR value to fold', type=float, default=0)
    parser.add_argument('--metafile', help='Metafile in', type=str, default='metafile.meta')

    args = parser.parse_args()

    if not args.input_file:
        print("You need to provide an csv file to read")
        sys.exit()
    prepfold_threads = args.presto_threads
    filterbank_file = args.input_fil_file
    df = pd.read_csv(args.input_file)
    df = df[df["fil_file"] ==  os.path.basename(args.input_fil_file)]
    meta = pd.read_json(args.metafile, typ='series')
    segment_start_sample = meta['segment_start_sample']
    segment_nsamples = meta['segment_nsamples']
    xml_segment_pepoch = meta['xml_segment_pepoch']
    total_nsamples = meta['total_nsamples']
    tstart = meta['tstart']
    fft_size = meta['fft_size']
    tsamp = meta['tsamp']
    source_name_prefix = meta['source_name']
    chunk_id = meta['chunk_id']
    beam_name = args.beam_name
    threads = args.pulsarx_threads

    if args.no_cands_to_fold is not None:
        # Limit the candidate df to N number of candidates
        df = df.nlargest(args.no_cands_to_fold, 'snr')
    if args.snr_min > 0:
        df = df[df['snr'] >= args.snr_min]
        
    if args.telescope == 'meerkat':
        PulsarX_Template = f"{args.pulsarx_fold_template}/meerkat_fold.template"
    elif args.telescope == 'effelsberg':
        PulsarX_Template = f"{args.pulsarx_fold_template}/Effelsberg_{args.beam_id}.template"
        
    start_fraction = round(segment_start_sample / total_nsamples, 3)
    end_fraction = round((segment_start_sample + segment_nsamples) / total_nsamples, 3)
    
    if end_fraction > 1:
        end_fraction = 1
        
    effective_tobs = tsamp * segment_nsamples
    
    try:
        if args.fold_technique == 'pulsarx':
            if args.subint_length is None:
                subint_length = int(effective_tobs / 64)
            else:
                subint_length = args.subint_length
            
            generate_pulsarX_cand_file(df, args.output_path, tsamp, fft_size, chunk_id, beam_name, source_name_prefix, args.cands_per_node)
            # Create meta file
            meta_file = f"{args.output_path}/{source_name_prefix}_{beam_name}_ck{chunk_id}_meta.txt"
            with open(meta_file, 'w') as f:
                f.write(f"Ncandidates: {len(df)}\n")
                f.write(f"start_fraction: {start_fraction}\n")
                f.write(f"end_fraction: {end_fraction}\n")
                f.write(f"fft_size: {fft_size}\n")
                f.write(f"ChunkId: {chunk_id}\n")
                f.write(f"Effective_tobs: {effective_tobs}\n")
                f.write(f"SegmentNSamples: {segment_nsamples}\n")
                f.write(f"SegmentPepoch: {xml_segment_pepoch}\n")
                f.write(f"Nsamples: {total_nsamples}\n")
                f.write(f"SubintLength: {subint_length}\n")
                f.write(f"Nsubband: {args.nsubband}\n")
                f.write(f"ClfdQValue: {args.clfd_q_value}\n")
                f.write(f"Telescope: {args.telescope}\n")
                f.write(f"Template: {PulsarX_Template}\n")
                f.write(f"BeamName: {args.beam_name}\n")
                f.write(f"BeamId: {args.beam_id}\n")
                f.write(f"UTCBeam: {args.utc_beam}\n")
                f.write(f"Candidates_per_node: {args.cands_per_node}\n")
                f.write(f"FilterbankFile: {filterbank_file}\n")
                f.write(f"StartTime: {tstart}\n")
                f.write(f"Numberofharmonics: {args.nh}\n")
                f.write(f"Fast_nbins: {args.fast_nbins}\n")
                f.write(f"Slow_nbins: {args.slow_nbins}\n")
                f.write(f"Threads: {threads}\n")
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
    