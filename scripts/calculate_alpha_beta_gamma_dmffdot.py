#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
import pandas as pd
import numpy as np
from astropy.io import fits
from multiprocessing import Pool, cpu_count

def setup_logging(verbose=False):
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')

def set_archive_to_zero_dm(archive_filename, output_filename):
    hdul = fits.open(archive_filename)
    hdul['PRIMARY'].header['CHAN_DM'] = 0
    hdul['SUBINT'].header['DM'] = 0
    hdul.writeto(output_filename, overwrite=True)

def run_dm_ffdot(archive_filename):
    # Example: You might want to limit concurrency of external calls to avoid heavy I/O load
    cmds = f"dmffdot --nosearch --saveimage -f {archive_filename}"
    subprocess.check_output(cmds, shell=True)
    # Replace the extension of the archive file with .px
    px_filename = archive_filename.replace('.ar', '.px')
    return px_filename

def extract_snr_from_px(px_filename):
    hdul = fits.open(px_filename)
    try:
        # Example: The S/N is inside extension #15, row #11, column #4
        snr = float(hdul[15].data[11][4].split('=')[1])
    except:
        logging.error(f"Could not extract S/N from {px_filename}. Assigning 0.0.")
        snr = 0.0
    return snr

def process_one_candidate(task):
    """
    task: (idx, archive_filename, snr_candidate_dm)

    Returns a tuple (idx, zero_dm_filename, zero_dm_snr, alpha).
    If something fails or snr_candidate_dm <= 0, fill them with NaNs/None.
    logging module is not multiprocess-safe, so use print for logging.
    """
    idx, archive_filename, snr_candidate_dm = task

    # If snr_candidate_dm <= 0, skip the computation, but still return an entry so we can mark it.
    if snr_candidate_dm <= 0:
        print(f"Row {idx} has non-positive S/N_new={snr_candidate_dm}; alpha=NaN.")
        return (idx, np.nan, np.nan, np.nan)

    if not os.path.exists(archive_filename):
        print(f"File {archive_filename} does not exist. Skipping row {idx}.")
        return (idx, np.nan, np.nan, np.nan)
    # 1) Construct name for DM=0 file
    output_filename = f"DM0_{archive_filename}"

    try:
        # 2) Set DM=0
        set_archive_to_zero_dm(archive_filename, output_filename)

        # 3) Run dmffdot
        px_filename = run_dm_ffdot(output_filename)

        # 4) Extract SNR
        snr_zero_dm = extract_snr_from_px(px_filename)

        # 5) alpha = ratio of zero_dm_snr to candidate_dm_snr
        alpha_val = snr_zero_dm / snr_candidate_dm 

        #Clean up intermediate files
        os.remove(output_filename)
        os.remove(px_filename)

    except Exception as e:
        # In case of any crash
        print(f"Failed processing row {idx}, file: {archive_filename}, reason: {e}")
        # Clean up partial files if needed
        if os.path.exists(output_filename):
            os.remove(output_filename)
        if os.path.exists(px_filename):
            os.remove(px_filename)
        return (idx, np.nan, np.nan, np.nan)

    
    return (idx, output_filename, snr_zero_dm, alpha_val)

def main():
   
    parser = argparse.ArgumentParser(
        description="Correct DM of folded archive to zero, "
                    "then compute alpha/beta/gamma, in parallel."
    )

    parser.add_argument(
        "-i", "--input_csv", required=True,
        help="Path to the search_fold_merged.csv file."
    )

    parser.add_argument(
        "-o", "--output_csv", default="search_fold_alpha_beta_gamma_merged.csv",
        help="Final output CSV after zero-DM refolding."
    )

    parser.add_argument(
        "-t", "--threshold", type=float, default=0.0,
        help="S/N_new threshold for refolding. "
            "If <=0, all rows are processed. Default=0.0 (all rows)."
    )

    parser.add_argument(
        "-n", "--num_workers", type=int, default=0,
        help="Number of parallel processes to use. "
            "If 0, uses all available cores."
    )

    # parser.add_argument('-g', '--generate_uuid', action='store_true', 
    #                     help='Generate a unique UUID string for each candidate in each model.')

    parser.add_argument(
        "-c", "--create_symlinks", action='store_true',
        help="Create symlinks to publishdir."
    )

    parser.add_argument(
        "-s", "--sourcedir", default=os.getcwd(),
        help="Source directory for symlinks. Default is the current working directory."
    )

    parser.add_argument(
        "-p", "--publishdir", default="publishdir",
        help="Directory to create symlinks in. Default is 'publishdir'."
    )

    parser.add_argument(
        "-v", "--verbose", action='store_true',
        help="Enable verbose logging."
    )

    args = parser.parse_args()
    setup_logging(verbose=args.verbose)

    ###########################################################################
    # 1) Load original merged CSV
    ###########################################################################
    logging.info(f"Reading input CSV: {args.input_csv}")
    df = pd.read_csv(args.input_csv)

    # Compute beta and gamma if possible
    if {'dm_old', 'dm_new', 'dm_err'}.issubset(df.columns):
        df['beta'] = (df['dm_old'] - df['dm_new']).abs() / df['dm_old']
        df['gamma'] = df['dm_err'] / df['dm_new']
    else:
        df['beta'] = np.nan
        df['gamma'] = np.nan

    # Create columns for DM=0 outputs (NaN by default)
    df['zero_dm_filename'] = np.nan
    df['zero_dm_snr'] = np.nan
    df['alpha'] = np.nan

    ###########################################################################
    # 2) Determine which rows to process
    ###########################################################################
    thresh = args.threshold
    if thresh <= 0:
        rows_to_process = df.index
        logging.info(f"Threshold <= 0. Processing all {len(rows_to_process)} rows.")
    else:
        rows_to_process = df.index[df['S/N_new'] > thresh]
        logging.info(f"Threshold = {thresh}. Processing {len(rows_to_process)} rows above threshold.")

    if len(rows_to_process) == 0:
        logging.info("No rows to process. zero_dm_* columns remain NaN.")
        df.to_csv(args.output_csv, index=False)
        logging.info(f"Wrote {args.output_csv}.")
        sys.exit(0)

    ###########################################################################
    # 3) Prepare tasks for multiprocessing
    ###########################################################################
    tasks = []
    for idx in rows_to_process:
        archive_filename = df.at[idx, 'fold_cands_filename']
        snr_candidate_dm = df.at[idx, 'S/N_new']
        tasks.append((idx, archive_filename, snr_candidate_dm))

    ###########################################################################
    # 4) Use a Pool to process tasks in parallel
    ###########################################################################
    nprocs = args.num_workers if args.num_workers > 0 else cpu_count()
    logging.info(f"Spawning {nprocs} parallel workers.")

    with Pool(processes=nprocs) as pool:
        # results is a list of (idx, zero_dm_filename, zero_dm_snr, alpha_val)
        results = pool.map(process_one_candidate, tasks)

    ###########################################################################
    # 5) Assign results back into df
    ###########################################################################
    processed_indices = []
    zero_dm_filenames = []
    zero_dm_snrs = []
    alphas = []
    logging.info("Assigning results back to DataFrame.")
    for (idx, zfn, zsnr, a) in results:
        processed_indices.append(idx)
        zero_dm_filenames.append(zfn)
        zero_dm_snrs.append(zsnr)
        alphas.append(a)

    df.loc[processed_indices, 'zero_dm_filename'] = zero_dm_filenames
    df.loc[processed_indices, 'zero_dm_snr'] = zero_dm_snrs
    df.loc[processed_indices, 'alpha'] = alphas
    logging.info("Results assigned back to DataFrame.")

    # if args.generate_uuid:
    #     logging.info("Generating UUIDs for each candidate.")
    #     uuid_col_name = "cand_tracker_database_uuid_alpha"
    #     df[uuid_col_name] = UUIDUtility.generate_uuid_list(len(df))
    #     uuid_col_name = "cand_tracker_database_uuid_beta"
    #     df[uuid_col_name] = UUIDUtility.generate_uuid_list(len(df))
    #     uuid_col_name = "cand_tracker_database_uuid_gamma"
    #     df[uuid_col_name] = UUIDUtility.generate_uuid_list(len(df))
    #     logging.info("UUIDs generated.")

    ###########################################################################
    # 6) Write final CSV
    ###########################################################################
    df.to_csv(args.output_csv, index=False)
    logging.info(f"Wrote zero-DM results to {args.output_csv}.")
    ###########################################################################
    # 7) Create symlinks if requested
    ###########################################################################

    if args.create_symlinks:
        logging.info("Creating symlinks to publishdir.")
        try:
            os.makedirs(args.publishdir, exist_ok=True)

            # Define thresholds and labels for symlinks
            thresholds = {
                'alpha': [0.25, 0.5, 1.0],
                'beta': [0.20],
                'gamma': [0.5, 1.0]
            }

            for key, threshold_list in thresholds.items():
                for thresh in threshold_list:
                    symlink_dir = os.path.join(args.publishdir, f"{key}_below_{thresh}")
                    os.makedirs(symlink_dir, exist_ok=True)
                    
                    # Filter dataframe based on threshold
                    shortlist = df.loc[df[key] <= thresh]
                    filenames = shortlist['fold_cands_filename'].values
                    filenames = [os.path.splitext(f)[0] + '.png' for f in filenames]
                    
                    # Create symlinks
                    for f in filenames:
                        try:
                            source = os.path.join(args.sourcedir, f)
                            destination = os.path.join(symlink_dir, f)
                            os.symlink(source, destination)
                        except Exception as e:
                            logging.error(f"Symlink creation failed for {f} at {key}_below_{thresh}. Reason: {e}")
                            continue

        except Exception as e:
            logging.error(f"Overall symlink creation failed. Reason: {e}")



if __name__ == "__main__":
    main()
