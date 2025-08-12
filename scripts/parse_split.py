import xml.etree.ElementTree as ET
import sys, os, subprocess
import argparse
import pandas as pd
import numpy as np
import logging
import time
import shlex
import threading
from multiprocessing import Pool, cpu_count
import re
import json
import glob

###############################################################################
# Logging Setup
###############################################################################


def setup_logging(verbose=False):
    """
    Set up logging configuration. DEBUG level if verbose, else INFO.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


################################################################################ Read and Parse XML files
###############################################################################


def read_candidate_files(
    files,
    cdm,
    chunk_id,
    ouput_dir,
    pulsarx_threads=20,
    pepoch_override=None,
    start_frac=None,
    end_frac=None,
    config_file=None,
    verbose=False,
):
    """
    Reads XML candidate files and aggregates the candidates in a single pandas
        DataFrame.

    Args:
        files (list): List of file paths to read.
        verbose (bool): If True, prints additional debug information.

    Returns:
        df_candidates (DataFrame): DataFrame containing all candidates from the files.
        obs_meta_data (dict): Dictionary containing metadata.
    """

    # Number of files found
    if verbose:
        print(f"{len(files)} candidates files found.")

    # Initialize list to store all rows and counters
    all_rows = []
    file_index = 0
    corrupted_list = []
    obs_meta_data = {}

    # Loop over all files
    for file in files:
        file = file.replace(",", "")
        try:
            if verbose:
                print("Reading file: {}".format(file))
            # Parse XML file
            tree = ET.parse(file)
        except:
            # If file is corrupted, skip it and increment corrupted count
            if verbose:
                print("{} is corrupted, skipping file".format(file))
            corrupted_list.append(file)
            continue
        root = tree.getroot()
        header_params = root[1]
        search_params = root[2]
        segment_params = root[3]
        candidates = root[7]

        # Read the fil file from the XML file
        filterbank_file = str(search_params.find("infilename").text)

        ignored_entries = [
            "candidate",
            "opt_period",
            "folded_snr",
            "byte_offset",
            "is_adjacent",
            "is_physical",
            "ddm_count_ratio",
            "ddm_snr_ratio",
        ]

        # Get metadata from the XML file
        tsamp = float(header_params.find("tsamp").text)
        total_nsamples = float(header_params.find("nsamples").text)
        source_name = str(header_params.find("source_name").text).strip()
        segment_start_sample = int(segment_params.find("segment_start_sample").text)
        segment_nsamples = int(segment_params.find("segment_nsamples").text)
        xml_segment_pepoch = float(segment_params.find("segment_pepoch").text)

        # Decide which pepoch to use
        if pepoch_override is not None:
            segment_pepoch = pepoch_override
            logging.info(f"Using user-provided pepoch = {segment_pepoch}")
        else:
            segment_pepoch = xml_segment_pepoch
            logging.info(f"Using pepoch from XML = {segment_pepoch}")

        # If user hasn't supplied explicit start_frac / end_frac, compute from XML data
        if start_frac is not None:
            user_start_fraction = start_frac
            logging.info(f"Using user-provided start_frac = {user_start_fraction}")
        else:
            user_start_fraction = round(segment_start_sample / total_nsamples, 3)
            logging.info(f"Using start_frac derived from XML = {user_start_fraction}")

        if end_frac is not None:
            user_end_fraction = end_frac
            logging.info(f"Using user-provided end_frac = {user_end_fraction}")
        else:
            user_end_fraction = round(
                (segment_start_sample + segment_nsamples) / total_nsamples, 3
            )
            logging.info(f"Using end_frac derived from XML = {user_end_fraction}")

        if user_end_fraction > 1:
            user_end_fraction = 1
            logging.info(f"End fraction > 1, setting to 1.0")

        # Create a row for each candidate and add it to the total list
        # all_rows.extend(create_row(root, candidates, file, file_index, filterbank_file))
        for candidate in candidates:
            cand_dict = {}
            for cand_entry in candidate.iter():
                if cand_entry.tag not in ignored_entries:
                    cand_dict[cand_entry.tag] = cand_entry.text
            cand_dict["cand_id_in_file"] = candidate.attrib.get("id")
            cand_dict["file"] = file
            cand_dict["filterbank_file"] = filterbank_file
            cand_dict["file_index"] = file_index
            cand_dict["pepoch"] = segment_params.find("segment_pepoch").text
            all_rows.append(cand_dict)

        # Grab needed meta data
        if file_index == 0:
            effective_tobs = tsamp * segment_nsamples
            tstart = float(root[1].find("tstart").text)
            fft_size = int(search_params.find("size").text)
            obs_meta_data = {
                "tsamp": tsamp,
                "segment_start_sample": segment_start_sample,
                "segment_nsamples": segment_nsamples,
                "xml_segment_pepoch": segment_pepoch,
                "chunk_id": chunk_id,
                "total_nsamples": total_nsamples,
                "obs_length": effective_tobs,
                "fft_size": fft_size,
                "tstart": tstart,
                "source_name": source_name,
                "threads": pulsarx_threads,
                "start_fraction": user_start_fraction,
                "end_fraction": user_end_fraction,
                "filterbank_file": filterbank_file,
            }
        file_index += 1

    if verbose:
        corrupted_count = len(corrupted_list)
        print(
            "{} inputs were corrupted ({}%).".format(
                corrupted_count, 100.0 * corrupted_count / len(files)
            )
        )
    obs_meta_data["xml_files"] = [x for x in files if x not in corrupted_list]
    obs_meta_data["corrupted_xml_files"] = corrupted_list

    if len(all_rows) == 0:
        # If all files were corrupted, return an empty DataFrame
        return pd.DataFrame(all_rows), obs_meta_data

    # Create a DataFrame from the list of rows
    df_candidates = pd.DataFrame(all_rows)

    # Cast columns to the appropriate data types
    df_candidates = df_candidates.astype(
        {
            "snr": float,
            "dm": float,
            "period": float,
            "nh": int,
            "acc": float,
            "nassoc": int,
            "cand_id_in_file": int,
        }
    )

    if verbose:
        print(f"{len(df_candidates)} candidates read.")

    logging.info(f"Dumping the candidates to unfiltered_df_for_folding.csv")
    df_candidates.to_csv(
        f"{ouput_dir}/unfiltered_for_folding_cdm_{cdm}.csv",
        index=False,
        float_format="%.18f",
    )

    # If a config file is provided, filter the dataframe accordingly
    if config_file:
        logging.info(
            f"XML file contains {len(df_candidates)} candidates before filtering."
        )
        logging.info(f"Applying folding configuration from {config_file}")
        df_candidates = apply_folding_configuration(df_candidates, config_file)
        logging.info(
            f"After filtering, {len(df_candidates)} candidates remain for folding."
        )
    else:
        logging.info(
            f"No configuration file provided, folding all {len(df_candidates)} candidates."
        )

    # Sort DataFrame by snr
    df_candidates.sort_values("snr", inplace=True, ascending=False)
    df_candidates.reset_index(inplace=True, drop=True)

    # Dump the candidates selected for folding to a CSV
    logging.info(f"Dumping the selected candidates to filtered_df_for_folding.csv")
    df_candidates.to_csv(
        f"{ouput_dir}/filtered_candidates_file_cdm_{cdm}.csv",
        index=False,
        float_format="%.18f",
    )

    return df_candidates, obs_meta_data


###############################################################################
# Pre-Select Candidates based on JSON Configuration
###############################################################################


def apply_folding_configuration(df: pd.DataFrame, config_file: str) -> pd.DataFrame:
    """
    Filter candidates from the dataframe based on the folding configuration
    specified in the JSON config file.

    The configuration file is expected to have the following structure:
    {
        "first_run": [
            {
                "spin_period": {"min": X, "max": Y},
                "dm": {"min": X, "max": Y},
                "fft_snr": {"min": X, "max": Y},
                "nh": {"min": X, "max": Y},
                "acc": {"min": X, "max": Y},
                "total_cands_limit": N
            },
            ...
        ]
    }

    For each filter, rows in `df` are selected if:
      - period is between spin_period.min and spin_period.max,
      - dm is between dm.min and dm.max,
      - snr is between fft_snr.min and fft_snr.max,
      - nh is between nh.min and nh.max,
      - acc is between acc.min and acc.max.

    Only up to total_cands_limit rows are selected per filter.
    The filtered candidates from all filters are then concatenated and duplicates dropped.

    Args:
        df: DataFrame containing candidate data.
        config_file: Path to JSON configuration file.

    Returns:
        A DataFrame with pre-selected candidates.
    """
    with open(config_file, "r") as f:
        config = json.load(f)

    # Assuming one group key; e.g., "first_run"
    group_key = list(config.keys())[0]
    filters = config[group_key]

    # filtered_dfs = []
    # for filter_def in filters:
    #     # Apply filter conditions
    #     sel = df.loc[
    #         (df['period'].between(filter_def['spin_period']['min'], filter_def['spin_period']['max'])) &
    #         (df['dm'].between(filter_def['dm']['min'], filter_def['dm']['max'])) &
    #         (df['snr'].between(filter_def['fft_snr']['min'], filter_def['fft_snr']['max'])) &
    #         (df['nh'].between(filter_def['nh']['min'], filter_def['nh']['max'])) &
    #         (df['acc'].between(filter_def['acc']['min'], filter_def['acc']['max']))
    #     ]

    filtered_dfs = []
    for filter_def in filters:
        conditions = []
        # Always apply period, snr, and nh bounds
        conditions.append(
            df["period"].between(
                filter_def["spin_period"]["min"], filter_def["spin_period"]["max"]
            )
        )
        conditions.append(
            df["snr"].between(
                filter_def["fft_snr"]["min"], filter_def["fft_snr"]["max"]
            )
        )
        conditions.append(
            df["nh"].between(filter_def["nh"]["min"], filter_def["nh"]["max"])
        )

        # Apply dm lower and upper bounds only if not null
        dm_min = filter_def["dm"].get("min")
        dm_max = filter_def["dm"].get("max")
        if dm_min is not None and dm_max is not None:
            conditions.append(df["dm"].between(dm_min, dm_max))
        elif dm_min is not None:
            conditions.append(df["dm"] >= dm_min)
        elif dm_max is not None:
            conditions.append(df["dm"] <= dm_max)
        # If both are None, do not filter dm

        # Apply acc lower and upper bounds only if not null
        acc_min = filter_def["acc"].get("min")
        acc_max = filter_def["acc"].get("max")
        if acc_min is not None and acc_max is not None:
            conditions.append(df["acc"].between(acc_min, acc_max))
        elif acc_min is not None:
            conditions.append(df["acc"] >= acc_min)
        elif acc_max is not None:
            conditions.append(df["acc"] <= acc_max)
        # If both are None, do not filter acc

        # Combine all conditions
        sel = df
        for cond in conditions:
            sel = sel[cond]

        logging.info(f"Filter {filter_def} selected {len(sel)} candidates.")
        # Limit to the desired number of candidates
        limit = filter_def.get("total_cands_limit", None)
        if limit is not None and len(sel) > limit:
            logging.info(
                f"Filter {filter_def} returned more candidates than the limit {limit}. Truncating to {limit} for folding."
            )
            sel = sel.head(limit)

        filtered_dfs.append(sel)

    if filtered_dfs:
        # Concatenate and drop duplicate candidates (assuming cand_id_in_file is unique)
        df_filtered = pd.concat(filtered_dfs)
        logging.info(f"Total candidates after Concatenation: {len(df_filtered)}")
        df_filtered = df_filtered.drop_duplicates(subset="cand_id_in_file")
        df_filtered = df_filtered.sort_values(by="cand_id_in_file")
        logging.info(f"Total candidates after dropping duplicates: {len(df_filtered)}")
        return df_filtered
    else:
        return df


###############################################################################
# PulsarX Candidate File
###############################################################################

# def generate_pulsarX_cand_file(cand_mod_frequencies, cand_dms, cand_accs, cand_snrs):
#     cand_file_path = 'pulsarx.candfile'
#     with open(cand_file_path, 'w') as f:
#         f.write("#id DM accel F0 F1 F2 S/N\n")
#         for i in range(len(cand_mod_frequencies)):
#             f.write("%d %f %f %f 0 0 %f\n" % (i, cand_dms[i], cand_accs[i], cand_mod_frequencies[i], cand_snrs[i]))

#     return cand_file_path


def generate_pulsarX_cand_file(df, meta, cdm, cands_dir, beam_name, cands_per_node=96):
    ## get the meta data from the XML file
    chunk_id = meta["chunk_id"]
    source_name_prefix = meta["source_name"]

    cand_dms = df["dm"].values
    cand_accs = df["acc"].values
    cand_period = df["period"].values
    cand_freq = 1 / cand_period
    cand_snrs = df["snr"].values

    cand_file_name = (
        f"{source_name_prefix}_{beam_name}_cdm_{cdm}_ck{chunk_id}_allCands.txt"
    )
    os.makedirs(cands_dir, exist_ok=True)
    cand_file_path = os.path.join(cands_dir, cand_file_name)
    number_of_candidates = len(cand_freq)
    with open(cand_file_path, "w") as f:
        f.write("#id DM accel F0 F1 F2 S/N\n")
        for i in range(number_of_candidates):
            f.write(
                "%d %f %f %f 0 0 %f\n"
                % (i, cand_dms[i], cand_accs[i], cand_freq[i], cand_snrs[i])
            )

    cand_files = split_cand_file(
        cand_file_path,
        number_of_candidates,
        cands_per_node,
        cdm,
        cands_dir,
        chunk_id,
        beam_name,
        source_name_prefix,
    )
    logging.info(f"{len(cand_files)} Candidate files created")
    return cand_files


def split_cand_file(
    cand_file_path,
    number_of_candidates,
    cands_per_node,
    cdm,
    cands_dir,
    chunk_id,
    beam_name,
    source_name_prefix,
):
    # Split the candidates into multiple candfiles
    cand_files = []
    candidates_df = pd.read_csv(cand_file_path, sep=" ", header=0)

    # if cands_per_node is not zero
    if cands_per_node == 0:
        cands_per_node = number_of_candidates
        num_files = 1
    else:
        num_files = number_of_candidates // cands_per_node + 1

    for i in range(num_files):
        cand_file = (
            f"{source_name_prefix}_{beam_name}_cdm_{cdm}_ck{chunk_id}_{i + 1}.candfile"
        )
        cand_file_path = os.path.join(cands_dir, cand_file)
        cand_files.append(cand_file_path)
        start_idx = i * cands_per_node
        end_idx = (i + 1) * cands_per_node
        cands = candidates_df.iloc[start_idx:end_idx].reset_index(drop=True)
        logging.info(f"Writing candidates {start_idx} to {end_idx} to a new candfile")
        try:
            with open(cand_file_path, "w") as f:
                f.write("#id DM accel F0 F1 F2 S/N\n")
                for index, row in cands.iterrows():
                    f.write(
                        f"{index} {row['DM']} {row['accel']} {row['F0']} 0 0 {row['S/N']}\n"
                    )

        except IOError as e:
            logging.error(f"Error writing candidate files {e}")
            return None

    return cand_files


###############################################################################
# Period Correction for PulsarX
###############################################################################


def period_correction_for_pulsarx(p0, pdot, no_of_samples, tsamp, fft_size):
    return p0 - pdot * float(fft_size - no_of_samples) * tsamp / 2


def period_correction_for_prepfold(p0, pdot, tsamp, fft_size):
    return p0 - pdot * float(fft_size) * tsamp / 2


def a_to_pdot(P_s, acc_ms2):
    LIGHT_SPEED = 2.99792458e8
    return P_s * acc_ms2 / LIGHT_SPEED


###############################################################################
# Folding with Presto
###############################################################################


def run_prepfold(args):
    """
    args is a tuple:
      (row, filterbank_file, tsamp, fft_size, source_name_prefix, rfifind_mask, extra_args)
    """
    (
        row,
        filterbank_file,
        tsamp,
        fft_size,
        source_name_prefix,
        rfifind_mask,
        extra_args,
    ) = args
    fold_period, pdot, cand_id, dm = row
    output_filename = (
        source_name_prefix + "_Peasoup_fold_candidate_id_" + str(int(cand_id) + 1)
    )

    cmd = "prepfold -fixchi -noxwin -topo"

    # Slow search if period > 0.1s
    if fold_period > 0.1:
        cmd += " -slow"

    if rfifind_mask:
        cmd += f" -mask {rfifind_mask}"

    # Add the fundamental folding parameters
    cmd += " -p %.16f -dm %.2f -pd %.16f -o %s %s" % (
        fold_period,
        dm,
        pdot,
        output_filename,
        filterbank_file,
    )

    # Append extra_args if provided
    if extra_args:
        cmd += f" {extra_args}"

    try:
        logging.debug(f"Running Presto command: {cmd}")
        subprocess.check_output(cmd, shell=True)
        return (True, cand_id)
    except subprocess.CalledProcessError as e:
        return (False, cand_id, str(e))


def fold_with_presto(
    df,
    filterbank_file,
    tsamp,
    fft_size,
    source_name_prefix,
    prepfold_threads,
    rfifind_mask=None,
    extra_args=None,
):
    num_cores = min(prepfold_threads, len(df))

    period = df["period"].values
    acc = df["acc"].values
    pdot = a_to_pdot(period, acc)
    fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)

    merged_data = np.column_stack(
        (fold_period, pdot, df["cand_id_in_file"].values, df["dm"].values)
    )
    args_list = [
        (
            row,
            filterbank_file,
            tsamp,
            fft_size,
            source_name_prefix,
            rfifind_mask,
            extra_args,
        )
        for row in merged_data
    ]

    pool = Pool(num_cores)
    results = pool.map(run_prepfold, args_list)
    pool.close()
    pool.join()

    for result in results:
        if not result[0]:  # If success is False
            logging.error(f"Error with candidate ID {result[1]}: {result[2]}")


###############################################################################
# Folding with PulsarX
###############################################################################


def fold_with_pulsarx(
    df,
    segment_start_sample,
    segment_nsamples,
    pepoch,
    total_nsamples,
    input_filenames,
    source_name_prefix,
    nbins,
    binplan,
    subint_length,
    nsubband,
    utc_beam,
    beam_name,
    pulsarx_threads,
    TEMPLATE,
    clfd_q_value,
    rfi_filter,
    cmask=None,
    start_fraction=None,
    end_fraction=None,
    extra_args=None,
    output_rootname=None,
    coherent_dm=0.0,
):
    """
    Fold candidates with pulsarx (psrfold_fil).
    'pepoch', 'start_fraction', and 'end_fraction' are either user-provided or derived.
    """

    cand_dms = df["dm"].values
    cand_accs = df["acc"].values
    cand_period = df["period"].values
    cand_freq = 1 / cand_period
    cand_snrs = df["snr"].values

    pulsarx_predictor = generate_pulsarX_cand_file(
        cand_freq, cand_dms, cand_accs, cand_snrs
    )

    nbins_string = "-b {} --nbinplan {}".format(nbins, binplan)

    if output_rootname is None:
        output_rootname = utc_beam

    if "ifbf" in beam_name:
        beam_tag = "--incoherent"
    elif "cfbf" in beam_name:
        beam_tag = "-i {}".format(int(beam_name.strip("cfbf")))
    else:
        # pulsarx does not take strings as beam names
        numeric_part = re.search(r"\d+$", beam_name).group()
        beam_tag = "-i {}".format(numeric_part)

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

    if rfi_filter:
        additional_flags = f"--rfi {rfi_filter}"
    else:
        additional_flags = ""

    # Build the base command
    script = (
        "psrfold_fil2 -v --render --output_width --cdm {} -t {} --candfile {} -n {} {} {} --template {} "
        "--clfd {} -L {} -f {} {} -o {} --srcname {} --pepoch {} --frac {} {} {}"
    ).format(
        coherent_dm,
        pulsarx_threads,
        pulsarx_predictor,
        nsubband,
        nbins_string,
        beam_tag,
        TEMPLATE,
        clfd_q_value,
        subint_length,
        input_filenames,
        zap_string,
        output_rootname,
        source_name_prefix,
        pepoch,
        start_fraction,
        end_fraction,
        additional_flags,
    )

    # Append extra_args if provided
    if extra_args:
        script += f" {extra_args}"

    logging.debug(f"Running PulsarX command: {script}")

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


###############################################################################
# Main
###############################################################################


def main():
    parser = argparse.ArgumentParser(
        description="Fold all candidates from Peasoup xml file"
    )
    parser.add_argument(
        "-o",
        "--output_path",
        help="Output path to save results",
        default=os.getcwd(),
        type=str,
    )
    parser.add_argument(
        "-i", "--input_file", help="Name of the input xml file", type=str, nargs="+"
    )
    parser.add_argument(
        "-ck", "--chunk_id", help="Chunk ID for the input file", type=str, default="10"
    )
    parser.add_argument("-m", "--mask_file", help="Mask file for prepfold", type=str)
    parser.add_argument(
        "-t",
        "--fold_technique",
        help="Technique to use for folding (presto or pulsarx)",
        type=str,
        default="pulsarx",
    )
    parser.add_argument(
        "-u", "--nbins_default", help="Default number of bins", type=int, default=32
    )
    parser.add_argument(
        "-l",
        "--binplan",
        help="Binplan for pulsarx.",
        type=str,
        default="0.005 32 0.01 64 0.1 128",
    )
    parser.add_argument(
        "-sub",
        "--subint_length",
        help="Subint length (s). Default is tobs/64",
        type=int,
        default=None,
    )
    parser.add_argument(
        "-nsub", "--nsubband", help="Number of subbands", type=int, default=64
    )
    parser.add_argument(
        "-clfd", "--clfd_q_value", help="CLFD Q value", type=float, default=2.0
    )
    parser.add_argument(
        "-rfi", "--rfi_filter", help="RFI filter value", type=str, default=""
    )
    parser.add_argument(
        "--source_name",
        help="Source name to be used for folding",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-b", "--beam_name", help="Beam name string", type=str, default="cfbf00000"
    )
    parser.add_argument("-b_id", "--beam_id", help="Beam ID", type=str, default="0")
    parser.add_argument(
        "-utc",
        "--utc_beam",
        help="UTC beam name string",
        type=str,
        default="2024-01-01-00:00:00",
    )
    parser.add_argument(
        "-c",
        "--chan_mask",
        help="Peasoup Channel mask file to be passed onto pulsarx",
        type=str,
        default="",
    )
    parser.add_argument(
        "-threads",
        "--pulsarx_threads",
        help="Number of threads to be used for pulsarx",
        type=int,
        default=20,
    )
    parser.add_argument(
        "-pthreads",
        "--presto_threads",
        help="Number of threads to be used for prepfold",
        type=int,
        default=12,
    )
    parser.add_argument(
        "-p", "--pulsarx_fold_template", help="Fold template pulsarx", type=str
    )
    parser.add_argument(
        "--template_dir",
        help="Template directory with meerkat and effelsberg templates",
        type=str,
    )
    parser.add_argument(
        "--telescope",
        help="Telescope name (e.g., meerkat, effelsberg)",
        type=str,
        default="meerkat",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    parser.add_argument(
        "-f",
        "--filterbank_publish_dir",
        help="Optional Path to filterbank publish directory. If None, use dir in XML file",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--config_file",
        type=str,
        help="Path to JSON configuration file to pre-select candidates to fold.",
        default=None,
    )
    parser.add_argument(
        "--filtered_candidates_file",
        type=str,
        default="filtered_df_for_folding.csv",
        help="Path to CSV file for pre-selected candidates for folding.",
    )
    # New arguments for pepoch, start_frac, end_frac overrides, and extra arguments
    parser.add_argument(
        "--pepoch_override",
        type=float,
        default=None,
        help="Override Pepoch value. If not provided, read from XML.",
    )
    parser.add_argument(
        "--start_frac",
        type=float,
        default=None,
        help="Override start fraction [0..1]. If not provided, computed from XML data.",
    )
    parser.add_argument(
        "--end_frac",
        type=float,
        default=None,
        help="Override end fraction [0..1]. If not provided, computed from XML data.",
    )
    parser.add_argument(
        "--extra_args",
        type=str,
        default=None,
        help="Extra arguments to pass to prepfold or psrfold_fil.",
    )
    parser.add_argument(
        "--cdm",
        type=float,
        default=0.0,
        help="Coherent DM to use for folding. Default is 0.0.",
    )
    parser.add_argument(
        "--cands_per_node",
        type=int,
        default=0,
        help="Number of candidates per node for splitting cand files.",
    )
    args = parser.parse_args()
    setup_logging(verbose=args.verbose)
    start_time = time.time()

    if not args.input_file:
        logging.error("You need to provide an XML file to read.")
        sys.exit(1)

    cands, meta = read_candidate_files(
        args.input_file,
        args.cdm,
        args.chunk_id,
        args.output_path,
        args.pulsarx_threads,
        args.pepoch_override,
        args.start_frac,
        args.end_frac,
        args.config_file,
        args.verbose,
    )

    if args.pulsarx_fold_template:
        PulsarX_Template = args.pulsarx_fold_template
    else:
        if args.telescope is None:
            logging.error("You need to provide a telescope name or a fold template.")
            sys.exit(1)
        elif args.telescope == "meerkat":
            PulsarX_Template = f"{args.template_dir}/meerkat_fold.template"
        elif args.telescope == "effelsberg":
            PulsarX_Template = f"{args.template_dir}/Effelsberg_{args.beam_id}.template"

    try:
        if args.fold_technique == "pulsarx":
            if args.subint_length is None:
                subint_length = int(meta["obs_length"] / 64)
                print(
                    f"Subint length not provided, using default value of {subint_length} seconds"
                )
            else:
                subint_length = args.subint_length
            generate_pulsarX_cand_file(
                cands,
                meta,
                args.cdm,
                args.output_path,
                args.beam_name,
                args.cands_per_node,
            )

            meta["subint_length"] = subint_length
            meta["nsubband"] = args.nsubband
            meta["clfd_q_value"] = args.clfd_q_value
            meta["telescope"] = args.telescope
            meta["template"] = PulsarX_Template
            meta["beam_name"] = args.beam_name
            meta["beam_id"] = args.beam_id
            meta["utc_beam"] = args.utc_beam
            meta["cands_per_node"] = args.cands_per_node
            meta["nbins"] = args.nbins_default
            meta["binplan"] = args.binplan
            meta["threads"] = args.pulsarx_threads
            meta["cdm"] = args.cdm
            meta["cmask"] = args.chan_mask
            meta["rfi_filter"] = args.rfi_filter

            # Create meta file
            meta_file = f"{args.output_path}/{meta['source_name']}_{args.beam_name}_ck{args.chunk_id}.meta"

            # dump the meta data to a file
            with open(meta_file, "w") as f:
                for key, value in meta.items():
                    f.write(f"{key}: {value}\n")
            logging.info(f"Meta data written to {meta_file}")

        else:
            subint_length = args.subint_length
            logging.info(f"Prepfold is not available, exiting.")
            sys.exit()
    finally:
        end_time = time.time()
        print(f"Total time taken: {end_time - start_time} seconds")

        # fold_with_pulsarx(
        #     df,
        #     segment_start_sample,
        #     segment_nsamples,
        #     segment_pepoch,
        #     total_nsamples,
        #     filterbank_file,
        #     source_name_prefix,
        #     args.nbins_high,
        #     args.nbins_low,
        #     subint_length,
        #     args.nsubband,
        #     args.utc_beam,
        #     args.beam_name,
        #     args.pulsarx_threads,
        #     PulsarX_Template,
        #     args.clfd_q_value,
        #     args.rfi_filter,
        #     cmask=args.chan_mask,
        #     start_fraction=user_start_fraction,
        #     end_fraction=user_end_fraction,
        #     extra_args=args.extra_args,
        #     output_rootname=args.output_rootname,
        #     coherent_dm=args.cdm
        # )


if __name__ == "__main__":
    main()
