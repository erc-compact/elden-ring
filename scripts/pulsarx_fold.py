import time
import shlex
import threading
import tempfile
import argparse
import logging
import numpy as np
import sys, os, subprocess
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
    parser.add_argument("-meta", help="Meta file", type=str, required=True)
    parser.add_argument("-cands", help="Candidate file", type=str, required=True)
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
    # output_dir = create_symlink_to_output_dir(args.output_path)
    output_dir = args.output_path
    logging.info(f"symlink path : {output_dir}")
    meta_dict = meta_parser(meta_file)
    fold_with_pulsarx(meta_dict, output_dir, cand_file)

    logging.info(f"Total time taken: {time.time() - start_time} seconds")

    # remove_symlink(output_dir)


if __name__ == "__main__":
    main()
