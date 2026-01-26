#!/usr/bin/env python3
"""
Batch prepfold script for PRESTO pipeline.

Folds multiple candidates in parallel using prepfold.

Usage:
    python3 presto_prepfold_batch.py \
        --candidate-csv sifted.csv \
        --input-file input.fil \
        --basename input \
        --npart 64 \
        --max-cands 100 \
        --mask-opt "-mask rfifind.mask" \
        --start-frac 0.0 \
        --end-frac 1.0
"""

import argparse
import csv
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def fold_candidate(args):
    """Fold a single candidate with prepfold."""
    f0, f1, f2, dm, cand_id, basename, input_file, mask_opt, npart, start_frac, end_frac, extra_flags = args
    outname = f"{basename}_cand_{cand_id}"
    cmd = (
        f"prepfold -noxwin "
        f"-f {f0:.15f} -fd {f1:.15f} -fdd {f2:.15f} "
        f"-dm {dm:.6f} -npart {npart} {mask_opt} "
        f"-start {start_frac} -end {end_frac} "
        f"{extra_flags} "
        f"-o {outname} {input_file}"
    )

    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        return True, cand_id
    except Exception as e:
        return False, str(e)


def main():
    parser = argparse.ArgumentParser(description="Batch prepfold for PRESTO candidates")
    parser.add_argument("--candidate-csv", required=True, help="Sifted candidates CSV")
    parser.add_argument("--input-file", required=True, help="Input filterbank file")
    parser.add_argument("--basename", required=True, help="Output file basename")
    parser.add_argument("--npart", type=int, default=64, help="Number of parts for folding")
    parser.add_argument("--max-cands", type=int, default=100, help="Maximum candidates to fold")
    parser.add_argument("--mask-opt", default="", help="RFI mask option string")
    parser.add_argument("--start-frac", type=float, default=0.0, help="Start fraction of data")
    parser.add_argument("--end-frac", type=float, default=1.0, help="End fraction of data")
    parser.add_argument("--extra-flags", default="", help="Extra prepfold flags")
    parser.add_argument("--max-workers", type=int, default=8, help="Max parallel workers")

    args = parser.parse_args()

    # Read candidates from CSV
    # Expected columns: f0,f1,f2,dm (from sifter output)
    candidates = []
    with open(args.candidate_csv, 'r') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= args.max_cands:
                break

            try:
                f0 = float(row.get('f0', 0.0))
                f1 = float(row.get('f1', 0.0))
                f2 = float(row.get('f2', 0.0))
                dm = float(row.get('dm', row.get('DM', 0)))
                cand_id = row.get('cand_id', row.get('id', i))

                candidates.append((
                    f0,
                    f1,
                    f2,
                    dm,
                    cand_id,
                    args.basename,
                    args.input_file,
                    args.mask_opt,
                    args.npart,
                    args.start_frac,
                    args.end_frac,
                    args.extra_flags
                ))
            except (KeyError, ValueError) as e:
                print(f"Warning: Could not parse candidate row: {e}")
                continue

    print(f"Folding {len(candidates)} candidates with prepfold")
    if not candidates:
        raise SystemExit(0)

    # Fold in parallel
    num_workers = min(args.max_workers, len(candidates))
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(fold_candidate, c) for c in candidates]
        for future in as_completed(futures):
            success, result = future.result()
            if not success:
                print(f"Warning: Folding failed: {result}")


if __name__ == "__main__":
    main()
