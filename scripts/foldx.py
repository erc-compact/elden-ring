import xml.etree.ElementTree as ET
import sys, os, subprocess
import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count

def generate_pulsarX_cand_file(df, tsamp, fft_size, beam_name, GC, dm_name, cands_per_node=96):
    cand_dms = df['dm'].values
    cand_accs = df['acc'].values
    cand_period = df['period'].values
    pdot = a_to_pdot(cand_period, cand_accs)
    # Prepfold is used here on purpose. Period modified to beginning of tobs with epoch pointing to tstart
    cand_mod_period_beginning_tobs = period_correction_for_prepfold(cand_period, pdot, tsamp, fft_size)
    cand_mod_frequencies = 1/cand_mod_period_beginning_tobs
    cand_snrs = df['snr'].values
    cand_file_path = f"{GC}_{beam_name}_{dm_name}_allCands.txt"
    number_of_candidates = len(cand_mod_frequencies)
    with open(cand_file_path, 'w') as f:
        f.write("#id DM accel F0 F1 S/N\n")
        for i in range(number_of_candidates):
            f.write("%d %f %f %f 0 %f\n" % (i, cand_dms[i], cand_accs[i], cand_mod_frequencies[i], cand_snrs[i]))
    
    cand_files = split_cand_file(cand_file_path, number_of_candidates, cands_per_node, beam_name, GC, dm_name)
    logging.info(f"{len(cand_files)} Candidate files created")
    return cand_files
    
def split_cand_file(cand_file_path, number_of_candidates, cands_per_node, beam_name, GC, dm_name):
    # Split the candidates into multiple candfiles
    cand_files = []
    num_files = number_of_candidates // cands_per_node + 1
    candidates_df = pd.read_csv(cand_file_path, sep=' ', header=0)
        
    for i in range(num_files):
        cand_file_path = f"{GC}_{beam_name}_{dm_name}_{i+1}.candfile"
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



def run_prepfold(args):
    row, filterbank_file, tsamp, fft_size, source_name_prefix, rfifind_mask = args
    fold_period, pdot, cand_id, dm = row
    output_filename = source_name_prefix + '_Peasoup_fold_candidate_id_' + str(int(cand_id) + 1)

    # Base command
    cmd = "prepfold -fixchi -noxwin -topo"
    
    # Check for fold_period and add -slow flag if needed
    if fold_period > 0.1:
        cmd += " -slow"
    
    # Add mask value if rfifind_mask is provided
    if rfifind_mask:
        cmd += " -mask %s" % rfifind_mask

    
    # fold parameters
    cmd += " -p %.16f -dm %.2f -pd %.16f -o %s %s" % (fold_period, dm, pdot, output_filename, filterbank_file)

    try:
        subprocess.check_output(cmd, shell=True)
        return (True, cand_id)
    except subprocess.CalledProcessError as e:
        return (False, cand_id, str(e))


def fold_with_presto(df, filterbank_file, tsamp, fft_size, source_name_prefix, prepfold_threads, rfifind_mask=None):
    #num_cores = min(cpu_count(), len(df))  # Use all available cores but no more than the number of rows
    num_cores = min(prepfold_threads, len(df))  # Use all available cores but no more than the number of rows

    period = df['period'].values
    acc = df['acc'].values
    pdot = a_to_pdot(period, acc)
    fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)

    merged_data = np.column_stack((fold_period, pdot, df['cand_id_in_file'].values, df['dm'].values))
    args_list = [(row, filterbank_file, tsamp, fft_size, source_name_prefix, rfifind_mask) for row in merged_data]

    pool = Pool(num_cores)
    results = pool.map(run_prepfold, args_list)
    pool.close()
    pool.join()

    for result in results:
        if not result[0]:  # If success is False
            print("Error with candidate ID %s: %s" % (result[1], result[2]))


def fold_with_pulsarx(df, input_filenames, tsamp, fft_size, source_name_prefix, tstart, fast_nbins, slow_nbins, subint_length, nsubband, utc_beam, beam_name, pulsarx_threads, TEMPLATE, cmask=None):
    cand_dms = df['dm'].values
    cand_accs = df['acc'].values
    cand_period = df['period'].values
    pdot = a_to_pdot(cand_period, cand_accs)
    # Prepfold is used here on purpose. Period modified to beginning of tobs with epoch pointing to tstart
    cand_mod_period_beginning_tobs = period_correction_for_prepfold(cand_period, pdot, tsamp, fft_size)
    cand_mod_frequencies = 1/cand_mod_period_beginning_tobs
    cand_snrs = df['snr'].values
    pulsarx_predictor = generate_pulsarX_cand_file(df, tsamp, fft_size, beam_name, source_name_prefix, args.dm_name, args.cands_per_node)
    nbins_string = "-b {} --nbinplan 0.1 {}".format(fast_nbins, slow_nbins)
    output_rootname = utc_beam + "_" + beam_name

    if 'ifbf' in beam_name:
        beam_tag = "--incoherent"
    elif 'cfbf' in beam_name :
        beam_tag = "-i {}".format(int(beam_name.strip("cfbf")))
    else:
        beam_tag = ""

    zap_string = ""
    if cmask is not None:
        cmask = cmask.strip()
        if cmask:
            try:
                zap_string = " ".join(["--rfi zap {} {}".format(
                    *i.split(":")) for i in cmask.split(",")])
            except Exception as error:
                raise Exception("Unable to parse channel mask: {}".format(
                    str(error)))
    #--rfi zap 1323 1384
    script = "psrfold_fil --plotx -v -t {} --candfile {} -n {} {} {} --template {} --clfd 2.0  -L {} -f {} --rfi zdot {} -o {} --srcname {} --pepoch {}".format(
              pulsarx_threads, pulsarx_predictor, nsubband, nbins_string, beam_tag, TEMPLATE, subint_length, input_filenames, zap_string, output_rootname, source_name_prefix, tstart)
    subprocess.check_output(script, shell=True)
    
def main():
    parser = argparse.ArgumentParser(description='Fold all candidates from Peasoup xml file')
    parser.add_argument('-o', '--output_path', help='Output path to save results',  default=os.getcwd(), type=str)
    parser.add_argument('-i', '--input_file', help='Name of the input csv file', type=str)
    parser.add_argument('-f', '--input_fil_file', help='Name of the input filterbank file', type=str)
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
    nsamples = meta['nsamples']
    tstart = meta['tstart']
    fft_size = meta['fft_size']
    tsamp = meta['tsamp']
    source_name_prefix = meta['source_name']

    if args.no_cands_to_fold is not None:
        # Limit the candidate df to N number of candidates
        df = df.nlargest(args.no_cands_to_fold, 'snr')
    if args.snr_min > 0:
        df = df[df['snr'] >= args.snr_min]

    if args.telescope == 'meerkat':
        PulsarX_Template = f"{args.pulsarx_fold_template}/meerkat_fold.template"
    elif args.telescope == 'effelsberg':
        PulsarX_Template = f"{args.pulsarx_fold_template}/Effelsberg_{args.beam_id}.template"

    if args.fold_technique == 'presto':
        fold_with_presto(df, filterbank_file, tsamp, fft_size, source_name_prefix, prepfold_threads)
    else:
        if args.subint_length is None:
            subint_length = int(nsamples * tsamp / 64)
        else:
            subint_length = args.subint_length

        fold_with_pulsarx(df, filterbank_file, tsamp, fft_size, source_name_prefix, tstart, args.fast_nbins, args.slow_nbins, subint_length, args.nsubband, args.utc_beam, args.beam_name, args.pulsarx_threads, PulsarX_Template,  args.channel_mask)


if __name__ == "__main__":
    main()
