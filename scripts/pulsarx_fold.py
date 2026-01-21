import time
import shlex
import threading
import tempfile
import argparse
import logging
import numpy as np
import pandas as pd
import sys, os, subprocess
import glob
import csv
import re
from datetime import datetime
import xml.etree.ElementTree as ET
from multiprocessing import Pool, cpu_count

logging.basicConfig(level=logging.INFO)


def buffered_stream_output(pipe, logger, log_level=logging.INFO, flush_interval=1.0):
    """
    Reads lines from 'pipe' and logs them every 'flush_interval' seconds.
    All lines are flushed at the end, ensuring none are skipped.
    """
    buffer = []
    last_flush = time.time()

    try:
        for line in iter(pipe.readline, ""):
            # Some commands produce empty lines occasionally; don't skip them
            buffer.append(line)

            # Flush if we've exceeded the interval
            if (time.time() - last_flush) >= flush_interval:
                for msg in buffer:
                    logger.log(log_level, msg.rstrip("\n"))
                buffer.clear()
                last_flush = time.time()

        # Final flush in case anything’s left
        for msg in buffer:
            logger.log(log_level, msg.rstrip("\n"))
        buffer.clear()

    finally:
        pipe.close()


def create_symlink_to_output_dir(output_dir):
    # temp_dir = tempfile.mkdtemp()
    # symlink_path = os.path.join(temp_dir, 'output_symlink')
    symlink_path = "output_symlink"

    if os.path.islink(symlink_path) or os.path.exists(symlink_path):
        os.unlink(symlink_path)

    os.symlink(output_dir, symlink_path)

    return symlink_path


def remove_symlink(symlink_path):
    if os.path.islink(symlink_path):
        os.unlink(symlink_path)
    os.rmdir(os.path.dirname(symlink_path))


def meta_parser(meta_file):
    with open(meta_file, "r") as f:
        meta = f.readlines()
    meta_dict = {}
    for line in meta:
        key, value = line.strip().split(":", 1)
        meta_dict[key.strip()] = value.strip()
    return meta_dict


# ============================================================================
# PERIOD CORRECTION FUNCTIONS (from foldx.py)
# These are essential for correct folding of accelerated pulsars
# ============================================================================

LIGHT_SPEED = 2.99792458e8  # Speed of Light in SI (m/s)


def a_to_pdot(P_s, acc_ms2):
    """
    Convert acceleration to period derivative.

    Args:
        P_s: Period in seconds
        acc_ms2: Acceleration in m/s^2

    Returns:
        Period derivative (dimensionless)
    """
    return P_s * acc_ms2 / LIGHT_SPEED


def period_correction_for_prepfold(p0, pdot, tsamp, fft_size):
    """
    Correct period to beginning of observation for prepfold.

    The period from peasoup/PRESTO is measured at the center of the observation.
    This function corrects it to the beginning of the observation, which is
    what prepfold expects when using -topo flag.

    Args:
        p0: Period in seconds (at center of observation)
        pdot: Period derivative
        tsamp: Sampling time in seconds
        fft_size: FFT size used in search

    Returns:
        Corrected period at beginning of observation
    """
    return p0 - pdot * float(fft_size) * tsamp / 2


def period_correction_for_pulsarx(p0, pdot, no_of_samples, tsamp, fft_size):
    """
    Correct period for PulsarX folding.

    Args:
        p0: Period in seconds
        pdot: Period derivative
        no_of_samples: Number of samples in the observation
        tsamp: Sampling time in seconds
        fft_size: FFT size used in search

    Returns:
        Corrected period for PulsarX
    """
    return p0 - pdot * float(fft_size - no_of_samples) * tsamp / 2


# ============================================================================
# CANDIDATE FILE GENERATION
# ============================================================================

def generate_pulsarx_candfile(df, tsamp, fft_size, source_name, dm_name, output_dir="."):
    """
    Generate a PulsarX-compatible candidate file from a DataFrame.

    The period is corrected to the beginning of observation using
    period_correction_for_prepfold (as per foldx.py - epoch points to tstart).

    Args:
        df: DataFrame with columns: dm, acc, period, snr
        tsamp: Sampling time in seconds
        fft_size: FFT size used in search
        source_name: Source name prefix
        dm_name: DM identifier for filename
        output_dir: Output directory

    Returns:
        Path to the generated candidate file
    """
    cand_dms = df['dm'].values
    cand_accs = df['acc'].values
    cand_period = df['period'].values
    cand_snrs = df['snr'].values

    # Calculate pdot from acceleration
    pdot = a_to_pdot(cand_period, cand_accs)

    # Correct period to beginning of observation (pepoch = tstart)
    # Using prepfold formula as per foldx.py comment:
    # "Prepfold is used here on purpose. Period modified to beginning of tobs with epoch pointing to tstart"
    cand_mod_period = period_correction_for_prepfold(cand_period, pdot, tsamp, fft_size)
    cand_mod_frequencies = 1.0 / cand_mod_period

    cand_file_path = os.path.join(output_dir, f"{source_name}_{dm_name}_allCands.candfile")
    number_of_candidates = len(cand_mod_frequencies)

    with open(cand_file_path, 'w') as f:
        f.write("#id DM accel F0 F1 S/N\n")
        for i in range(number_of_candidates):
            f.write("%d %f %f %.15f 0 %f\n" % (
                i, cand_dms[i], cand_accs[i], cand_mod_frequencies[i], cand_snrs[i]
            ))

    logging.info(f"Created PulsarX candfile: {cand_file_path} with {number_of_candidates} candidates")
    return cand_file_path


def split_candfile(cand_file_path, cands_per_node, source_name, dm_name, output_dir="."):
    """
    Split a candidate file into multiple smaller files for parallel processing.

    Args:
        cand_file_path: Path to the original candidate file
        cands_per_node: Number of candidates per split file
        source_name: Source name prefix
        dm_name: DM identifier
        output_dir: Output directory

    Returns:
        List of paths to split candidate files
    """
    candidates_df = pd.read_csv(cand_file_path, sep=' ', comment='#',
                                 names=['id', 'DM', 'accel', 'F0', 'F1', 'S/N'])
    number_of_candidates = len(candidates_df)

    if cands_per_node <= 0 or number_of_candidates <= cands_per_node:
        return [cand_file_path]

    cand_files = []
    num_files = (number_of_candidates + cands_per_node - 1) // cands_per_node

    for i in range(num_files):
        split_path = os.path.join(output_dir, f"{source_name}_{dm_name}_{i+1}.candfile")
        start_idx = i * cands_per_node
        end_idx = min((i + 1) * cands_per_node, number_of_candidates)
        cands = candidates_df.iloc[start_idx:end_idx].reset_index(drop=True)

        logging.info(f"Writing candidates {start_idx} to {end_idx} to {split_path}")

        with open(split_path, 'w') as f:
            f.write('#id DM accel F0 F1 S/N\n')
            for index, row in cands.iterrows():
                f.write(f"{index} {row['DM']} {row['accel']} {row['F0']:.15f} 0 {row['S/N']}\n")

        cand_files.append(split_path)

    return cand_files


# ============================================================================
# PREPFOLD FOLDING (PRESTO)
# ============================================================================

def run_prepfold_single(args):
    """
    Run prepfold on a single candidate.

    Args:
        args: Tuple of (row, filterbank_file, source_name_prefix, rfifind_mask)
              where row = (fold_period, pdot, cand_id, dm)

    Returns:
        Tuple of (success, cand_id, [error_message])
    """
    row, filterbank_file, source_name_prefix, rfifind_mask = args
    fold_period, pdot, cand_id, dm = row
    output_filename = f"{source_name_prefix}_fold_candidate_id_{int(cand_id) + 1}"

    # Base command with topocentric folding
    cmd = "prepfold -fixchi -noxwin -topo"

    # Add -slow flag for slow pulsars
    if fold_period > 0.1:
        cmd += " -slow"

    # Add mask if provided
    if rfifind_mask:
        cmd += f" -mask {rfifind_mask}"

    # Add fold parameters
    cmd += f" -p {fold_period:.16f} -dm {dm:.2f} -pd {pdot:.16e} -o {output_filename} {filterbank_file}"

    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        return (True, cand_id)
    except subprocess.CalledProcessError as e:
        return (False, cand_id, str(e))


def fold_with_prepfold(meta_dict, output_dir, cand_file, rfifind_mask=None, prepfold_threads=12):
    """
    Fold candidates using PRESTO's prepfold.

    Uses topocentric folding with period corrected to beginning of observation.

    Args:
        meta_dict: Metadata dictionary with filterbank info
        output_dir: Output directory
        cand_file: Path to candidate CSV file
        rfifind_mask: Optional path to rfifind mask file
        prepfold_threads: Number of parallel prepfold processes
    """
    filterbank_file = meta_dict["filterbank_file"]
    tsamp = float(meta_dict.get("tsamp", meta_dict.get("sampling_time", 0.0)))
    fft_size = float(meta_dict["fft_size"])
    source_name_prefix = meta_dict["source_name"]

    # Read candidates
    df = pd.read_csv(cand_file)

    if len(df) == 0:
        logging.warning("No candidates to fold!")
        return

    # Calculate period derivative and corrected period
    period = df['period'].values
    acc = df['acc'].values
    pdot = a_to_pdot(period, acc)
    fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)

    # Get candidate IDs
    if 'cand_id_in_file' in df.columns:
        cand_ids = df['cand_id_in_file'].values
    elif 'cand_id' in df.columns:
        cand_ids = df['cand_id'].values
    else:
        cand_ids = np.arange(len(df))

    # Prepare arguments for parallel processing
    merged_data = np.column_stack((fold_period, pdot, cand_ids, df['dm'].values))
    args_list = [(row, filterbank_file, source_name_prefix, rfifind_mask) for row in merged_data]

    # Run in parallel
    num_cores = min(prepfold_threads, len(df))
    logging.info(f"Folding {len(df)} candidates with prepfold using {num_cores} threads")

    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(output_dir)

    try:
        pool = Pool(num_cores)
        results = pool.map(run_prepfold_single, args_list)
        pool.close()
        pool.join()

        # Report results
        success_count = sum(1 for r in results if r[0])
        fail_count = len(results) - success_count
        logging.info(f"prepfold completed: {success_count} successful, {fail_count} failed")

        for result in results:
            if not result[0]:
                logging.error(f"Error with candidate ID {result[1]}: {result[2]}")
    finally:
        os.chdir(original_dir)


def fold_with_pulsarx(meta_dict, output_dir, cand_file):
    subint_length = float(meta_dict["subint_length"])
    start_fraction = float(meta_dict["start_fraction"])
    end_fraction = float(meta_dict["end_fraction"])
    pepoch = float(meta_dict["xml_segment_pepoch"])
    fft_size = float(meta_dict["fft_size"])
    chunk_id = str(meta_dict["chunk_id"])
    nsubband = int(meta_dict["nsubband"])
    clfd_q_value = float(meta_dict["clfd_q_value"])
    Telescope = meta_dict["telescope"]
    template = meta_dict["template"]
    beam_name = meta_dict["beam_name"]
    beam_id = int(meta_dict["beam_id"])
    utc_beam = meta_dict["utc_beam"]
    filterbank_file = meta_dict["filterbank_file"]
    nbins = meta_dict["nbins"]
    binplan = meta_dict["binplan"]
    nsubband = int(meta_dict["nsubband"])
    pulsarx_threads = meta_dict["threads"]
    nbins_string = "-b {} --nbinplan {}".format(nbins, binplan)
    cmask = str(meta_dict["cmask"])
    rfi_filter = str(meta_dict["rfi_filter"])
    coherent_dm = meta_dict["cdm"]
    source_name_prefix = meta_dict["source_name"]

    if "ifbf" in beam_name:
        beam_tag = "--incoherent -i {}".format(int(beam_name.strip("ifbf")))
    elif "cfbf" in beam_name:
        beam_tag = "-i {}".format(int(beam_name.strip("cfbf")))
    else:
        beam_tag = ""

    zap_string = ""
    if cmask is not None:
        cmask = cmask.strip()
        if cmask:
            try:
                zap_string = " ".join(
                    ["--rfi zap {} {}".format(*i.split(":")) for i in cmask.split(",")]
                )
            except Exception as error:
                raise Exception(f"Unable to parse channel mask: {error}")
    # if rfi_filter is not None:
    if rfi_filter:
        additional_flags = f"--rfi {rfi_filter}"
    else:
        additional_flags = ""

    cand_file_name = os.path.basename(cand_file)
    logging.info(f"folding {cand_file_name}")
    # output_rootname = os.path.join(output_dir, f"{cand_file_name.split('_')[0]}_{cand_file_name.split('_')[-1].split('.')[0]}_ck{chunk_id}")
    output_rootname = f"{cand_file_name.split('_')[0]}_{cand_file_name.split('_')[-1].split('.')[0]}_cdm_{coherent_dm}_ck{chunk_id}"
    # output_rootname = output_dir
    logging.info(f"Output path: {output_rootname}")
    # os.makedirs(output_rootname, exist_ok=True)
    # print output with the cand_file
    print("Processing cand_file: ", cand_file)

    # # Build the base command
    # script1 = (
    #     "psrfold_fil -v --render --plotx --output_width --cdm {} -t {} --candfile {} -n {} {} {} --template {} --clfd {} -L {} -f {} {} -o {} --srcname {} --pepoch {} --frac {} {} {}"
    # ).format(
    #     coherent_dm,
    #     pulsarx_threads,
    #     cand_file,
    #     nsubband,
    #     nbins_string,
    #     beam_tag,
    #     template,
    #     clfd_q_value,
    #     subint_length,
    #     filterbank_file,
    #     zap_string,
    #     output_rootname,
    #     source_name_prefix,
    #     pepoch,
    #     start_fraction,
    #     end_fraction,
    #     additional_flags,
    # )
    #
    # script2 = (
    #     "psrfold_fil2 -v --render --plotx --output_width --cdm {} -t {} --candfile {} -n {} {} {} --template {} --clfd {} -L {} -f {} {} -o {} --srcname {} --pepoch {} --frac {} {} {}"
    # ).format(
    #     coherent_dm,
    #     pulsarx_threads,
    #     cand_file,
    #     nsubband,
    #     nbins_string,
    #     beam_tag,
    #     template,
    #     clfd_q_value,
    #     subint_length,
    #     filterbank_file,
    #     zap_string,
    #     output_rootname,
    #     source_name_prefix,
    #     pepoch,
    #     start_fraction,
    #     end_fraction,
    #     additional_flags,
    # )
    #
    # if beam_id in [1, 2]:
    #     script = script1
    # else:
    #     script = script2
    #
    script = (
        "psrfold_fil2 -v --render --plotx --output_width --cdm {} -t {} --candfile {} -n {} {} {} --template {} --clfd {} -L {} -f {} {} -o {} --srcname {} --pepoch {} --frac {} {} {}"
    ).format(
        coherent_dm,
        pulsarx_threads,
        cand_file,
        nsubband,
        nbins_string,
        beam_tag,
        template,
        clfd_q_value,
        subint_length,
        filterbank_file,
        zap_string,
        output_rootname,
        source_name_prefix,
        pepoch,
        start_fraction,
        end_fraction,
        additional_flags,
    )

    logging.info(f"Running PulsarX command: {script}")

    process = subprocess.Popen(
        shlex.split(script),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,  # line-buffered
    )

    # Two threads, one for stdout (INFO level), one for stderr (WARNING level)
    stdout_thread = threading.Thread(
        target=buffered_stream_output,
        args=(process.stdout, logging.getLogger(), logging.INFO, 10.0),
        daemon=True,
    )
    stderr_thread = threading.Thread(
        target=buffered_stream_output,
        args=(process.stderr, logging.getLogger(), logging.WARNING, 10.0),
        daemon=True,
    )

    stdout_thread.start()
    stderr_thread.start()

    stdout_thread.join()
    stderr_thread.join()

    return_code = process.wait()

    if return_code != 0:
        logging.error(f"psrfold_fil2 returned non-zero exit status {return_code}")
        sys.exit(1)


def fold_candidates_from_csv(meta_dict, output_dir, cand_csv, fold_backend='pulsarx',
                              cands_per_node=350, rfifind_mask=None, prepfold_threads=12):
    """
    Fold candidates from a CSV file using either PulsarX or prepfold.

    This is the main entry point for folding peasoup/PRESTO candidates.
    It handles period correction and uses tstart as pepoch (as per foldx.py).

    Args:
        meta_dict: Metadata dictionary with observation parameters
        output_dir: Output directory for folded products
        cand_csv: Path to CSV file with candidates (columns: dm, acc, period, snr)
        fold_backend: 'pulsarx' or 'presto' (prepfold)
        cands_per_node: Number of candidates per candfile (for PulsarX splitting)
        rfifind_mask: Optional rfifind mask file (for prepfold)
        prepfold_threads: Number of parallel prepfold processes

    The CSV file should have columns:
        - dm: Dispersion measure
        - acc: Acceleration in m/s^2
        - period: Period in seconds (at center of observation)
        - snr: Signal-to-noise ratio
    """
    # Read observation parameters from meta file
    filterbank_file = meta_dict["filterbank_file"]
    tsamp = float(meta_dict.get("tsamp", meta_dict.get("sampling_time", 0.0)))
    fft_size = float(meta_dict["fft_size"])
    source_name = meta_dict["source_name"]

    # For pepoch: use xml_segment_pepoch if available (for existing pipeline),
    # otherwise fall back to tstart (for PRESTO/peasoup CSV candidates)
    if "xml_segment_pepoch" in meta_dict:
        pepoch = float(meta_dict["xml_segment_pepoch"])
        logging.info(f"Using pepoch from xml_segment_pepoch = {pepoch:.15f}")
    elif "tstart" in meta_dict:
        pepoch = float(meta_dict["tstart"])
        logging.info(f"Using pepoch from tstart = {pepoch:.15f}")
    else:
        raise ValueError("No pepoch available. meta file must have xml_segment_pepoch or tstart.")

    # Read candidates
    df = pd.read_csv(cand_csv)

    if len(df) == 0:
        logging.warning("No candidates to fold!")
        return

    logging.info(f"Folding {len(df)} candidates using {fold_backend}")

    if fold_backend == 'presto':
        # Use prepfold for folding
        fold_with_prepfold(meta_dict, output_dir, cand_csv, rfifind_mask, prepfold_threads)
    else:
        # Use PulsarX for folding
        # Generate candidate file with corrected periods/frequencies
        dm_name = "allDMs"
        candfile = generate_pulsarx_candfile(df, tsamp, fft_size, source_name, dm_name, output_dir)

        # Split if needed
        if cands_per_node > 0 and len(df) > cands_per_node:
            candfiles = split_candfile(candfile, cands_per_node, source_name, dm_name, output_dir)
        else:
            candfiles = [candfile]

        # Update meta_dict with pepoch
        meta_dict_updated = meta_dict.copy()
        meta_dict_updated['xml_segment_pepoch'] = str(pepoch)

        # Fold each candfile
        for cf in candfiles:
            logging.info(f"Folding candfile: {cf}")
            fold_with_pulsarx(meta_dict_updated, output_dir, cf)


def main():
    parser = argparse.ArgumentParser(
        description="Fold candidates using PulsarX or prepfold (supports Peasoup candfiles and CSV candidates)"
    )
    parser.add_argument(
        "-o",
        "--output_path",
        help="Output path to save results",
        default=os.getcwd(),
        type=str,
    )
    parser.add_argument("-meta", help="Meta file", type=str, required=True)
    parser.add_argument("-cands", help="Candidate file (candfile or CSV)", type=str, required=True)

    # Folding backend selection
    parser.add_argument(
        "--fold_backend",
        help="Folding backend: 'pulsarx' (default) or 'presto' (prepfold)",
        type=str,
        default='pulsarx',
        choices=['pulsarx', 'presto']
    )

    # CSV mode for peasoup/PRESTO candidates
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Candidate file is a CSV (not a PulsarX candfile). "
             "CSV must have columns: dm, acc, period, snr. "
             "Period will be corrected and pepoch=tstart will be used."
    )

    # Additional options
    parser.add_argument(
        "--cands_per_node",
        help="Number of candidates per candfile for parallel processing (default: 350)",
        type=int,
        default=350
    )
    parser.add_argument(
        "--rfifind_mask",
        help="Path to rfifind mask file (for prepfold)",
        type=str,
        default=None
    )
    parser.add_argument(
        "--prepfold_threads",
        help="Number of parallel prepfold processes (default: 12)",
        type=int,
        default=12
    )

    args = parser.parse_args()

    if not args.meta:
        print("You need to provide a meta file")
        sys.exit(1)

    if not args.cands:
        print("You need to provide a candidate file")
        sys.exit(1)

    start_time = time.time()

    meta_file = args.meta
    cand_file = args.cands
    output_dir = args.output_path

    logging.info(f"Output path: {output_dir}")
    logging.info(f"Fold backend: {args.fold_backend}")

    meta_dict = meta_parser(meta_file)

    if args.csv:
        # CSV mode: fold candidates from CSV file with period correction
        # Uses tstart as pepoch (correct approach from foldx.py)
        logging.info("Using CSV candidate folding mode (period correction applied, pepoch=tstart)")
        fold_candidates_from_csv(
            meta_dict,
            output_dir,
            cand_file,
            fold_backend=args.fold_backend,
            cands_per_node=args.cands_per_node,
            rfifind_mask=args.rfifind_mask,
            prepfold_threads=args.prepfold_threads
        )
    else:
        # Standard candfile mode: fold using existing PulsarX candfile
        # Assumes pepoch is already set correctly in meta file
        if args.fold_backend == 'presto':
            logging.warning("prepfold requires --csv mode with period/dm/acc values. "
                           "Falling back to PulsarX for candfile folding.")
        logging.info("Using PulsarX candfile folding mode")
        fold_with_pulsarx(meta_dict, output_dir, cand_file)

    logging.info(f"Total time taken: {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
