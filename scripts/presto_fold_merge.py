#!/usr/bin/env python3
"""
Unified fold merge script for PRESTO pipeline.

Replaces two separate inline scripts:
- presto_fold_merge (PRESTO backend with .pfd/.bestprof files)
- presto_fold_merge_pulsarx (PulsarX backend with PNG only)

Usage:
    python3 presto_fold_merge.py --backend presto \
        --sifted-csv sifted.csv --provenance-csv provenance.csv \
        --pointing P123 --cluster NGC6544 --beam-name cfbf00000 \
        --beam-id 0 --utc-start 2024-01-01 --ra 12:00:00 --dec -30:00:00 \
        --cdm 30.0 --filterbank-path /path/to/file.fil

    python3 presto_fold_merge.py --backend pulsarx \
        --sifted-csv sifted.csv --provenance-csv provenance.csv \
        --meta-file observation.meta \
        --pointing P123 ...
"""

import argparse
import csv
import glob
import os
import re
import shutil
import sys


def parse_bestprof(path):
    """Parse PRESTO .bestprof file to extract pepoch, dm_opt, sn_fold."""
    info = {}
    try:
        with open(path, 'r') as f:
            for line in f:
                if line.startswith("# Epoch_bary (MJD)"):
                    val = line.split("=", 1)[1].strip()
                    try:
                        info["pepoch"] = float(val)
                        info["mjd_start"] = float(val)
                    except ValueError:
                        pass
                elif line.startswith("# Epoch_topo"):
                    val = line.split("=", 1)[1].strip()
                    try:
                        if "pepoch" not in info:
                            info["pepoch"] = float(val)
                    except ValueError:
                        pass
                elif line.startswith("# Best DM"):
                    val = line.split("=", 1)[1].strip()
                    try:
                        info["dm_opt"] = float(val)
                    except ValueError:
                        pass
                elif line.startswith("# Prob(Noise)"):
                    match = re.search(r"([0-9.]+)\s*sigma", line)
                    if match:
                        try:
                            info["sn_fold"] = float(match.group(1))
                        except ValueError:
                            pass
    except IOError:
        pass
    return info


def parse_meta_file(path):
    """Parse meta file to extract pepoch (PulsarX backend)."""
    pepoch = ""
    try:
        with open(path, "r") as mf:
            for line in mf:
                if line.startswith("xml_segment_pepoch:"):
                    pepoch = line.split(":", 1)[1].strip()
                    break
                if line.startswith("tstart:"):
                    pepoch = line.split(":", 1)[1].strip()
        if pepoch:
            try:
                pepoch = float(pepoch)
            except Exception:
                pass
    except Exception:
        pepoch = ""
    return pepoch


def build_meta_defaults(args):
    """Build metadata defaults dict from command-line args."""
    meta = {
        "pointing_id": args.pointing or "",
        "beam_id": args.beam_id or "",
        "beam_name": args.beam_name or "",
        "utc_start": args.utc_start or "",
        "ra": args.ra or "",
        "dec": args.dec or "",
        "cdm": args.cdm or "",
        "filterbank_path": args.filterbank_path or "",
        "cluster": args.cluster or "",
        "segment_id": "0",
        "source_name": "",
    }
    if meta["utc_start"]:
        meta["metafile_path"] = "meta/%s.meta" % meta["utc_start"]
    else:
        meta["metafile_path"] = ""
    return meta


def read_sifted_csv(path):
    """Read sifted candidates CSV and build lookup dict."""
    sifted_data = {}
    sifted_rows = []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sifted_rows.append(row)
            cand_key = row.get('cand_id') or row.get('id') or row.get('cand_num', '')
            if cand_key:
                sifted_data[str(cand_key)] = row
            cand_num = row.get('cand_num', '')
            if cand_num:
                sifted_data[str(cand_num)] = row
    return sifted_data, sifted_rows


def read_provenance_csv(path):
    """Read provenance CSV and build lookup dict."""
    prov_data = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            prov_id = row.get('id', '')
            if prov_id:
                prov_data[str(prov_id)] = row
    return prov_data


def apply_provenance(result, cand_id, prov_data):
    """Apply provenance data to result dict."""
    if cand_id in prov_data:
        result["accel_file"] = prov_data[cand_id].get("accel_file", "")
        result["cand_num"] = prov_data[cand_id].get("cand_num", "")
        result["sn_fft"] = prov_data[cand_id].get("sn_fft", prov_data[cand_id].get("sigma", ""))
        if result.get("sn_fft") not in (None, ""):
            try:
                result["sn_fft"] = float(result["sn_fft"])
            except Exception:
                pass


def apply_field_normalization(result):
    """Normalize dm/f0/f1 field names."""
    # sn_fold fallback
    if result.get("sn_fold") in (None, ""):
        sigma_val = result.get("sigma") or result.get("snr") or result.get("sn_fft")
        if sigma_val not in (None, ""):
            try:
                result["sn_fold"] = float(sigma_val)
            except Exception:
                result["sn_fold"] = sigma_val

    # dm normalization
    if "dm" in result and "dm_user" not in result:
        result["dm_user"] = result.get("dm")
    if "dm_opt" not in result and "dm" in result:
        result["dm_opt"] = result.get("dm")

    # f0/f1 normalization
    if "f0" in result and "f0_user" not in result:
        result["f0_user"] = result.get("f0")
    if "f1" in result and "f1_user" not in result:
        result["f1_user"] = result.get("f1")


def merge_presto_backend(args, meta_defaults, sifted_data, prov_data):
    """Merge results for PRESTO backend (with .pfd/.bestprof)."""
    os.makedirs("pfd_files", exist_ok=True)
    os.makedirs("png_files", exist_ok=True)

    # Copy pfd and png files
    for pfd in glob.glob("*.pfd"):
        shutil.copy(pfd, "pfd_files/")
    for png in glob.glob("*.png"):
        shutil.copy(png, "png_files/")

    results = []

    for bp in glob.glob("*.bestprof"):
        raw_cand_id = bp.split('_cand_')[1].replace('.pfd.bestprof', '') if '_cand_' in bp else ''
        cand_id = raw_cand_id
        if raw_cand_id:
            cand_id = raw_cand_id.split('_')[0]
        basename = os.path.basename(bp).replace('.pfd.bestprof', '')

        png_file = ''
        if os.path.exists(basename + '.png'):
            png_file = basename + '.png'
        elif os.path.exists(basename + '.pfd.png'):
            png_file = basename + '.pfd.png'

        result = {
            'basename': basename,
            'cand_id': cand_id,
            'cand_label': raw_cand_id,
            'pfd_file': basename + '.pfd',
            'png_file': png_file
        }

        result.update(meta_defaults)

        if cand_id in sifted_data:
            result.update(sifted_data[cand_id])

        apply_provenance(result, cand_id, prov_data)

        bp_info = parse_bestprof(bp)
        result.update(bp_info)

        if (not result.get("utc_start")) and ("mjd_start" in result):
            result["utc_start"] = str(result.get("mjd_start"))
        if "pepoch" not in result and "mjd_start" in result:
            result["pepoch"] = result.get("mjd_start")

        apply_field_normalization(result)
        results.append(result)

    return results


def merge_pulsarx_backend(args, meta_defaults, sifted_data, sifted_rows, prov_data):
    """Merge results for PulsarX backend (PNG only, no .pfd)."""
    os.makedirs("png_files", exist_ok=True)

    # Copy png files
    for png in glob.glob("*.png"):
        shutil.copy(png, "png_files/")

    if not glob.glob("png_files/*.png"):
        open("png_files/NO_PNG", "w").close()

    # Get pepoch from meta file
    pepoch = ""
    if args.meta_file:
        pepoch = parse_meta_file(args.meta_file)

    pngs = sorted(glob.glob("*.png"))
    results = []

    for idx, row in enumerate(sifted_rows):
        cand_id = row.get('cand_id') or row.get('id') or row.get('cand_num', '')
        basename = f"cand_{cand_id}" if cand_id else f"cand_{idx}"
        png_file = ""
        if idx < len(pngs):
            png_file = os.path.basename(pngs[idx])
            basename = os.path.splitext(png_file)[0]

        result = {
            'basename': basename,
            'cand_id': cand_id,
            'pfd_file': '',
            'png_file': png_file,
        }

        result.update(meta_defaults)
        if pepoch != "":
            result["pepoch"] = pepoch
        result.update(row)

        apply_provenance(result, cand_id, prov_data)
        apply_field_normalization(result)
        results.append(result)

    return results


def write_results(results, output_file="merged_results.csv"):
    """Write results to CSV."""
    if results:
        fieldnames = list(results[0].keys())
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
    else:
        with open(output_file, 'w') as f:
            f.write("basename,cand_id,pfd_file,png_file\n")


def main():
    parser = argparse.ArgumentParser(description="Merge PRESTO fold results into CSV")
    parser.add_argument("--backend", required=True, choices=["presto", "pulsarx"],
                        help="Folding backend used")
    parser.add_argument("--sifted-csv", required=True, help="Sifted candidates CSV")
    parser.add_argument("--provenance-csv", required=True, help="Provenance CSV")
    parser.add_argument("--meta-file", help="Meta file (PulsarX backend)")

    # Metadata arguments
    parser.add_argument("--pointing", default="", help="Pointing ID")
    parser.add_argument("--cluster", default="", help="Cluster name")
    parser.add_argument("--beam-name", default="", help="Beam name")
    parser.add_argument("--beam-id", default="", help="Beam ID")
    parser.add_argument("--utc-start", default="", help="UTC start time")
    parser.add_argument("--ra", default="", help="Right ascension")
    parser.add_argument("--dec", default="", help="Declination")
    parser.add_argument("--cdm", default="", help="Central DM")
    parser.add_argument("--filterbank-path", default="", help="Filterbank file path")

    parser.add_argument("--output", default="merged_results.csv", help="Output CSV file")

    args = parser.parse_args()

    meta_defaults = build_meta_defaults(args)
    sifted_data, sifted_rows = read_sifted_csv(args.sifted_csv)
    prov_data = read_provenance_csv(args.provenance_csv)

    if args.backend == "presto":
        results = merge_presto_backend(args, meta_defaults, sifted_data, prov_data)
    else:
        results = merge_pulsarx_backend(args, meta_defaults, sifted_data, sifted_rows, prov_data)

    write_results(results, args.output)
    print(f"Merged {len(results)} candidates to {args.output}")


if __name__ == "__main__":
    main()
