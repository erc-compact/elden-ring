#!/usr/bin/env python3
"""
Unified state saving script for PRESTO pipeline.

Replaces 5 separate inline scripts:
- save_presto_rfi_state
- save_presto_birdies_state
- save_presto_dedisperse_state
- save_presto_search_state
- save_presto_sift_fold_state

Usage:
    python3 presto_save_state.py --stage rfi \
        --input-file input.fil --rfi-mask mask.mask --rfi-stats rfifind.stats --rfi-inf rfifind.inf

    python3 presto_save_state.py --stage birdies \
        --input-file input.fil --rfi-mask mask.mask --rfi-stats rfifind.stats \
        --rfi-inf rfifind.inf --zerodm-inf zerodm.inf --zaplist birdies.zaplist

    python3 presto_save_state.py --stage dedisperse \
        --input-file input.fil --rfi-mask mask.mask --rfi-stats rfifind.stats --zaplist birdies.zaplist

    python3 presto_save_state.py --stage search \
        --input-file input.fil --rfi-mask mask.mask --rfi-stats rfifind.stats

    python3 presto_save_state.py --stage sift_fold \
        --input-file input.fil --sifted-csv sifted.csv --provenance-csv provenance.csv
"""

import argparse
import json
import os
import glob
import sys


# Stage to next workflow mapping
NEXT_WORKFLOW = {
    "rfi": "presto_birdies",
    "birdies": "presto_dedisperse",
    "dedisperse": "presto_search",
    "search": "presto_sift_fold",
    "sift_fold": "presto_postprocess",
}


def get_abs_path(path):
    """Get absolute path, return None for non-existent optional files."""
    if not path or path == "NO_ZAPLIST" or path == "NO_MASK":
        return None
    if not os.path.exists(path):
        return None
    return os.path.abspath(path)


def save_rfi_state(args):
    """Save RFI detection state."""
    state = {
        "stage": "rfi",
        "input_file": get_abs_path(args.input_file),
        "rfi_mask": get_abs_path(args.rfi_mask),
        "rfi_stats": get_abs_path(args.rfi_stats),
        "rfi_inf": get_abs_path(args.rfi_inf),
        "next_workflow": NEXT_WORKFLOW["rfi"],
    }
    return state, "rfi_state.json"


def save_birdies_state(args):
    """Save birdies detection state."""
    state = {
        "stage": "birdies",
        "input_file": get_abs_path(args.input_file),
        "rfi_mask": get_abs_path(args.rfi_mask),
        "rfi_stats": get_abs_path(args.rfi_stats),
        "rfi_inf": get_abs_path(args.rfi_inf),
        "zerodm_inf": get_abs_path(args.zerodm_inf),
        "zaplist": get_abs_path(args.zaplist),
        "next_workflow": NEXT_WORKFLOW["birdies"],
    }
    return state, "birdies_state.json"


def save_dedisperse_state(args):
    """Save dedispersion state with dat/inf file lists."""
    # Collect dat and inf files from current directory
    dat_files = sorted([os.path.abspath(f) for f in glob.glob("*.dat")])
    inf_files = sorted([os.path.abspath(f) for f in glob.glob("*.inf")])

    state = {
        "stage": "dedisperse",
        "input_file": get_abs_path(args.input_file),
        "rfi_mask": get_abs_path(args.rfi_mask),
        "rfi_stats": get_abs_path(args.rfi_stats),
        "zaplist": get_abs_path(args.zaplist),
        "dat_files": dat_files,
        "inf_files": inf_files,
        "next_workflow": NEXT_WORKFLOW["dedisperse"],
    }
    return state, "dedisperse_state.json"


def save_search_state(args):
    """Save acceleration search state with accel/cand file lists."""
    # Collect ACCEL files from current directory
    accel_all = sorted([os.path.abspath(f) for f in glob.glob("*_ACCEL_*")])
    # Main ACCEL outputs (exclude sidecars like .cand/.txtcand/.inf)
    accel_files = [f for f in accel_all if not any(f.endswith(ext) for ext in [".cand", ".txtcand", ".inf"])]
    # Retain cand files for compatibility
    cand_files = [f for f in accel_all if f.endswith(".cand")]

    state = {
        "stage": "search",
        "input_file": get_abs_path(args.input_file),
        "rfi_mask": get_abs_path(args.rfi_mask),
        "rfi_stats": get_abs_path(args.rfi_stats),
        "accel_files": accel_files,
        "cand_files": cand_files,
        "next_workflow": NEXT_WORKFLOW["search"],
    }
    return state, "search_state.json"


def save_sift_fold_state(args):
    """Save sift/fold state with pfd file list."""
    # Collect PFD files from current directory
    pfd_files = sorted([os.path.abspath(f) for f in glob.glob("*.pfd")])

    state = {
        "stage": "sift_fold",
        "input_file": get_abs_path(args.input_file),
        "sifted_csv": get_abs_path(args.sifted_csv),
        "provenance_csv": get_abs_path(args.provenance_csv),
        "pfd_files": pfd_files,
        "next_workflow": NEXT_WORKFLOW["sift_fold"],
    }
    return state, "sift_fold_state.json"


def main():
    parser = argparse.ArgumentParser(description="Save PRESTO pipeline state to JSON file")
    parser.add_argument("--stage", required=True, choices=["rfi", "birdies", "dedisperse", "search", "sift_fold"],
                        help="Pipeline stage name")

    # Common arguments
    parser.add_argument("--input-file", required=True, help="Input filterbank file")

    # RFI stage arguments
    parser.add_argument("--rfi-mask", help="RFI mask file")
    parser.add_argument("--rfi-stats", help="RFI stats file")
    parser.add_argument("--rfi-inf", help="RFI info file")

    # Birdies stage arguments
    parser.add_argument("--zerodm-inf", help="Zero-DM info file")
    parser.add_argument("--zaplist", help="Birdies zaplist file")

    # Sift/fold stage arguments
    parser.add_argument("--sifted-csv", help="Sifted candidates CSV")
    parser.add_argument("--provenance-csv", help="Provenance CSV")

    args = parser.parse_args()

    # Dispatch to appropriate handler
    handlers = {
        "rfi": save_rfi_state,
        "birdies": save_birdies_state,
        "dedisperse": save_dedisperse_state,
        "search": save_search_state,
        "sift_fold": save_sift_fold_state,
    }

    state, output_file = handlers[args.stage](args)

    # Write state file
    with open(output_file, "w") as f:
        json.dump(state, f, indent=2)

    # Print summary
    if args.stage == "dedisperse":
        print(f"Dedisperse state saved with {len(state['dat_files'])} dat files")
    elif args.stage == "search":
        print(f"Search state saved with {len(state['accel_files'])} ACCEL files and {len(state['cand_files'])} cand files")
    elif args.stage == "sift_fold":
        print(f"Sift/fold state saved with {len(state['pfd_files'])} PFD files")
    else:
        print(f"{args.stage.upper()} state saved: {output_file}")


if __name__ == "__main__":
    main()
