#!/usr/bin/env python3
"""
Create CandyJar-compatible tarballs from PRESTO pipeline results.
Similar to create_candyjar_tarball.py but handles PRESTO pfd files.
"""

import argparse
import csv
import glob
import os
import sys
import tarfile
import shutil
from datetime import datetime


def parse_pfd_info(pfd_file):
    """
    Extract basic info from a pfd filename.
    Expected format: beam_DM{dm}_{cand_num}.pfd
    """
    info = {
        'filename': os.path.basename(pfd_file),
        'dm': 0.0,
        'cand_num': 0
    }

    base = os.path.basename(pfd_file).replace('.pfd', '')
    parts = base.split('_')

    for i, part in enumerate(parts):
        if part.startswith('DM'):
            try:
                info['dm'] = float(part[2:])
            except ValueError:
                pass
        if part.isdigit():
            info['cand_num'] = int(part)

    return info


def read_presto_csv(csv_file):
    """
    Read the PRESTO candidates CSV file.
    """
    candidates = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)

    return candidates


def read_classification_scores(score_file):
    """
    Read classification scores (if available).
    """
    scores = {}

    if not os.path.exists(score_file):
        return scores

    with open(score_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cand_id = row.get('id') or row.get('candidate_id')
            if cand_id:
                scores[cand_id] = row

    return scores


def create_tarball(output_name, candidates, png_dir, pfd_dir, output_dir,
                   pics_scores=None, alpha_threshold=1.0, pics_threshold=0.1,
                   snr_threshold=6.0, verbose=False):
    """
    Create CandyJar-compatible tarball from PRESTO results.

    Args:
        output_name: Name for output tarball
        candidates: List of candidate dictionaries
        png_dir: Directory containing PNG files
        pfd_dir: Directory containing PFD files
        output_dir: Output directory for tarball
        pics_scores: Optional PICS classification scores
        alpha_threshold: Alpha threshold for filtering
        pics_threshold: PICS score threshold for filtering
        snr_threshold: SNR threshold for filtering
        verbose: Print verbose output
    """
    os.makedirs(output_dir, exist_ok=True)

    # Create working directory
    work_dir = os.path.join(output_dir, 'tarball_work')
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)

    # Prepare candidate lists
    all_candidates = []
    alpha_below_one = []
    pics_above_threshold = []

    for cand in candidates:
        cand_entry = cand.copy()

        # Get classification scores if available
        cand_id = str(cand.get('id', cand.get('cand_num', '')))
        if pics_scores and cand_id in pics_scores:
            score_data = pics_scores[cand_id]
            cand_entry['pics_score'] = score_data.get('pics_score', score_data.get('score', 0.0))
            cand_entry['alpha'] = score_data.get('alpha', 0.0)
            cand_entry['beta'] = score_data.get('beta', 0.0)
            cand_entry['gamma'] = score_data.get('gamma', 0.0)

        # Find associated PNG and PFD files
        if 'png_file' not in cand_entry:
            # Search for matching PNG
            cand_name = cand.get('filename', '').replace('.pfd', '')
            png_patterns = [
                os.path.join(png_dir, f"*{cand_name}*.png"),
                os.path.join(png_dir, f"*{cand_id}*.png"),
            ]
            for pattern in png_patterns:
                matches = glob.glob(pattern)
                if matches:
                    cand_entry['png_file'] = matches[0]
                    break

        if 'pfd_file' not in cand_entry:
            # Search for matching PFD
            cand_name = cand.get('filename', '').replace('.png', '')
            pfd_patterns = [
                os.path.join(pfd_dir, f"*{cand_name}*.pfd"),
                os.path.join(pfd_dir, f"*{cand_id}*.pfd"),
            ]
            for pattern in pfd_patterns:
                matches = glob.glob(pattern)
                if matches:
                    cand_entry['pfd_file'] = matches[0]
                    break

        all_candidates.append(cand_entry)

        # Filter by alpha
        alpha = float(cand_entry.get('alpha', 2.0))
        if alpha < alpha_threshold:
            alpha_below_one.append(cand_entry)

        # Filter by PICS score
        pics = float(cand_entry.get('pics_score', 0.0))
        snr = float(cand_entry.get('sigma', cand_entry.get('sn_fold', cand_entry.get('snr', 0.0))))
        if pics >= pics_threshold and snr >= snr_threshold:
            pics_above_threshold.append(cand_entry)

    # Create CSV files
    csv_fields = [
        'id', 'dm', 'period_ms', 'freq_hz', 'sigma', 'snr', 'accel',
        'alpha', 'beta', 'gamma', 'pics_score',
        'png_file', 'pfd_file', 'source_file'
    ]

    def write_candidates_csv(cands, filename):
        with open(filename, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction='ignore')
            writer.writeheader()
            for c in cands:
                writer.writerow(c)

    all_csv = os.path.join(work_dir, 'candidates.csv')
    alpha_csv = os.path.join(work_dir, 'candidates_alpha_below_one.csv')
    pics_csv = os.path.join(work_dir, 'candidates_pics_above_threshold.csv')

    write_candidates_csv(all_candidates, all_csv)
    write_candidates_csv(alpha_below_one, alpha_csv)
    write_candidates_csv(pics_above_threshold, pics_csv)

    # Create plots directory and copy PNG files
    plots_dir = os.path.join(work_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    for cand in all_candidates:
        png_file = cand.get('png_file')
        if png_file and os.path.exists(png_file):
            shutil.copy2(png_file, plots_dir)

    # Create archives directory and copy PFD files
    archives_dir = os.path.join(work_dir, 'archives')
    os.makedirs(archives_dir, exist_ok=True)

    for cand in all_candidates:
        pfd_file = cand.get('pfd_file')
        if pfd_file and os.path.exists(pfd_file):
            shutil.copy2(pfd_file, archives_dir)

    # Create metadata file
    metadata = {
        'pipeline': 'ELDEN-RING (PRESTO)',
        'created': datetime.now().isoformat(),
        'total_candidates': len(all_candidates),
        'alpha_below_threshold': len(alpha_below_one),
        'pics_above_threshold': len(pics_above_threshold),
        'alpha_threshold': alpha_threshold,
        'pics_threshold': pics_threshold,
        'snr_threshold': snr_threshold
    }

    meta_file = os.path.join(work_dir, 'metadata.txt')
    with open(meta_file, 'w') as f:
        for key, value in metadata.items():
            f.write(f"{key}: {value}\n")

    # Create tarball
    tarball_path = os.path.join(output_dir, f"{output_name}.tar.gz")
    with tarfile.open(tarball_path, 'w:gz') as tar:
        tar.add(work_dir, arcname=output_name)

    # Cleanup
    shutil.rmtree(work_dir)

    if verbose:
        print(f"Created tarball: {tarball_path}")
        print(f"  Total candidates: {len(all_candidates)}")
        print(f"  Alpha < {alpha_threshold}: {len(alpha_below_one)}")
        print(f"  PICS >= {pics_threshold}: {len(pics_above_threshold)}")

    # Also copy individual CSV files to output directory
    final_all_csv = os.path.join(output_dir, 'candidates.csv')
    final_alpha_csv = os.path.join(output_dir, 'candidates_alpha_below_one.csv')
    final_pics_csv = os.path.join(output_dir, 'candidates_pics_above_threshold.csv')

    # Re-write CSVs to output dir (since we deleted work_dir)
    write_candidates_csv(all_candidates, final_all_csv)
    write_candidates_csv(alpha_below_one, final_alpha_csv)
    write_candidates_csv(pics_above_threshold, final_pics_csv)

    return tarball_path, all_csv, alpha_csv, pics_csv


def main():
    parser = argparse.ArgumentParser(
        description='Create CandyJar-compatible tarball from PRESTO results'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input candidates CSV file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output tarball name (without .tar.gz extension)'
    )
    parser.add_argument(
        '-d', '--output-dir',
        default='.',
        help='Output directory for tarball (default: current directory)'
    )
    parser.add_argument(
        '--png-dir',
        default='.',
        help='Directory containing PNG files'
    )
    parser.add_argument(
        '--pfd-dir',
        default='.',
        help='Directory containing PFD files'
    )
    parser.add_argument(
        '--scores',
        default=None,
        help='Classification scores CSV file (optional)'
    )
    parser.add_argument(
        '--alpha-threshold',
        type=float,
        default=1.0,
        help='Alpha threshold for filtering (default: 1.0)'
    )
    parser.add_argument(
        '--pics-threshold',
        type=float,
        default=0.1,
        help='PICS score threshold for filtering (default: 0.1)'
    )
    parser.add_argument(
        '--snr-threshold',
        type=float,
        default=6.0,
        help='SNR threshold for filtering (default: 6.0)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Read candidates
    candidates = read_presto_csv(args.input)

    if not candidates:
        print(f"ERROR: No candidates found in {args.input}")
        sys.exit(1)

    if args.verbose:
        print(f"Read {len(candidates)} candidates from {args.input}")

    # Read scores if available
    scores = None
    if args.scores:
        scores = read_classification_scores(args.scores)
        if args.verbose:
            print(f"Read {len(scores)} classification scores")

    # Create tarball
    create_tarball(
        output_name=args.output,
        candidates=candidates,
        png_dir=args.png_dir,
        pfd_dir=args.pfd_dir,
        output_dir=args.output_dir,
        pics_scores=scores,
        alpha_threshold=args.alpha_threshold,
        pics_threshold=args.pics_threshold,
        snr_threshold=args.snr_threshold,
        verbose=args.verbose
    )


if __name__ == '__main__':
    main()
