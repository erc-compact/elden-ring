#!/usr/bin/env python3
"""
Sift PRESTO acceleration search candidates.
Based on PRESTO's ACCEL_sift.py with modifications for ELDEN pipeline integration.
"""

import os
import sys
import glob
import argparse
from collections import defaultdict


def parse_accel_cand_file(filename):
    """
    Parse a PRESTO ACCEL candidate file.

    Returns a list of dictionaries with candidate info.
    """
    candidates = []

    if not os.path.exists(filename):
        return candidates

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Skip header lines
    in_header = True
    for line in lines:
        line = line.strip()

        # Empty line or comment
        if not line or line.startswith('#'):
            continue

        # End of header section
        if 'Cand' in line and 'Sigma' in line:
            in_header = False
            continue

        if in_header:
            continue

        # Parse candidate line
        parts = line.split()
        if len(parts) < 8:
            continue

        try:
            cand = {
                'cand_num': int(parts[0]),
                'sigma': float(parts[1]),
                'summed_power': float(parts[2]),
                'coherent_power': float(parts[3]),
                'num_harm': int(parts[4]),
                'period_ms': float(parts[5]),
                'freq_hz': float(parts[6]),
                'fft_r': float(parts[7]),
                'z': float(parts[8]) if len(parts) > 8 else 0.0,
                'accel': float(parts[9]) if len(parts) > 9 else 0.0,
                'filename': filename
            }
            candidates.append(cand)
        except (ValueError, IndexError) as e:
            continue

    return candidates


def remove_duplicates(candidates, period_tol=1e-6, dm_tol=0.5):
    """
    Remove duplicate candidates based on period and DM tolerance.
    """
    if not candidates:
        return []

    # Sort by sigma (descending)
    sorted_cands = sorted(candidates, key=lambda x: x['sigma'], reverse=True)

    unique = []
    for cand in sorted_cands:
        is_dup = False
        for uniq in unique:
            # Check if periods are similar (within tolerance)
            period_diff = abs(cand['period_ms'] - uniq['period_ms']) / uniq['period_ms']
            if period_diff < period_tol:
                is_dup = True
                break
        if not is_dup:
            unique.append(cand)

    return unique


def remove_harmonics(candidates, harm_tol=0.001):
    """
    Remove harmonically related candidates, keeping the fundamental.
    """
    if not candidates:
        return []

    # Sort by frequency
    sorted_cands = sorted(candidates, key=lambda x: x['freq_hz'])

    unique = []
    for cand in sorted_cands:
        is_harmonic = False
        for uniq in unique:
            # Check if this is a harmonic of an existing candidate
            freq_ratio = cand['freq_hz'] / uniq['freq_hz']
            # Check for integer ratio (within tolerance)
            nearest_int = round(freq_ratio)
            if nearest_int > 0:
                ratio_diff = abs(freq_ratio - nearest_int) / nearest_int
                if ratio_diff < harm_tol:
                    is_harmonic = True
                    break
        if not is_harmonic:
            unique.append(cand)

    return unique


def filter_by_period(candidates, min_period=0.001, max_period=15.0):
    """
    Filter candidates by period range.
    """
    return [c for c in candidates if min_period <= c['period_ms'] / 1000.0 <= max_period]


def filter_by_sigma(candidates, min_sigma=5.0):
    """
    Filter candidates by minimum sigma (SNR).
    """
    return [c for c in candidates if c['sigma'] >= min_sigma]


def sift_candidates(accel_files, dm_values=None, min_sigma=5.0,
                    min_period=0.001, max_period=15.0,
                    remove_dups=True, remove_harms=True):
    """
    Main sifting function for PRESTO acceleration search results.

    Args:
        accel_files: List of ACCEL candidate files
        dm_values: Dictionary mapping filenames to DM values
        min_sigma: Minimum sigma threshold
        min_period: Minimum period in seconds
        max_period: Maximum period in seconds
        remove_dups: Whether to remove duplicates
        remove_harms: Whether to remove harmonics

    Returns:
        List of sifted candidate dictionaries
    """
    all_candidates = []

    for accel_file in accel_files:
        cands = parse_accel_cand_file(accel_file)

        # Add DM info if available
        if dm_values and accel_file in dm_values:
            dm = dm_values[accel_file]
            for c in cands:
                c['dm'] = dm
        else:
            # Try to extract DM from filename
            try:
                base = os.path.basename(accel_file)
                # Look for DM pattern like "DM100.00"
                for part in base.split('_'):
                    if part.startswith('DM'):
                        dm = float(part[2:])
                        for c in cands:
                            c['dm'] = dm
                        break
            except:
                pass

        all_candidates.extend(cands)

    # Apply filters
    filtered = filter_by_sigma(all_candidates, min_sigma)
    filtered = filter_by_period(filtered, min_period, max_period)

    if remove_dups:
        filtered = remove_duplicates(filtered)

    if remove_harms:
        filtered = remove_harmonics(filtered)

    # Sort by sigma (descending)
    filtered.sort(key=lambda x: x['sigma'], reverse=True)

    return filtered


def write_candfile(candidates, output_file, dm=None):
    """
    Write candidates to a candfile format compatible with prepfold/psrfold.
    """
    with open(output_file, 'w') as f:
        f.write('#id dm acc F0 F1 F2 S/N\n')
        for i, cand in enumerate(candidates, 1):
            cand_dm = cand.get('dm', dm if dm is not None else 0.0)
            # F0 in Hz, F1 in Hz/s (approximated from acceleration)
            f0 = cand['freq_hz']
            # Approximate F1 from acceleration: a = c * f1 / f0
            # f1 = a * f0 / c (speed of light)
            accel = cand.get('accel', 0.0)
            f1 = accel * f0 / 299792458.0 if accel != 0 else 0.0

            f.write(f"{i} {cand_dm:.3f} {accel:.6f} {f0:.10f} {f1:.15e} 0 {cand['sigma']:.2f}\n")


def write_csv(candidates, output_file):
    """
    Write candidates to CSV format.
    """
    with open(output_file, 'w') as f:
        f.write('id,dm,period_ms,freq_hz,sigma,accel,z,num_harm,source_file\n')
        for i, cand in enumerate(candidates, 1):
            dm = cand.get('dm', 0.0)
            f.write(f"{i},{dm:.3f},{cand['period_ms']:.10f},{cand['freq_hz']:.10f},"
                    f"{cand['sigma']:.3f},{cand.get('accel', 0.0):.6f},"
                    f"{cand.get('z', 0.0):.2f},{cand['num_harm']},"
                    f"{os.path.basename(cand['filename'])}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Sift PRESTO acceleration search candidates'
    )
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input ACCEL candidate files (can use glob patterns)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output filename (will create .candfile and .csv)'
    )
    parser.add_argument(
        '--min-sigma',
        type=float,
        default=5.0,
        help='Minimum sigma threshold (default: 5.0)'
    )
    parser.add_argument(
        '--min-period',
        type=float,
        default=0.001,
        help='Minimum period in seconds (default: 0.001)'
    )
    parser.add_argument(
        '--max-period',
        type=float,
        default=15.0,
        help='Maximum period in seconds (default: 15.0)'
    )
    parser.add_argument(
        '--no-duplicates',
        action='store_true',
        help='Do not remove duplicates'
    )
    parser.add_argument(
        '--no-harmonics',
        action='store_true',
        help='Do not remove harmonics'
    )
    parser.add_argument(
        '--dm',
        type=float,
        default=None,
        help='Default DM value for candidates'
    )
    parser.add_argument(
        '--max-cands',
        type=int,
        default=500,
        help='Maximum number of candidates to output (default: 500)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Expand glob patterns
    accel_files = []
    for pattern in args.input:
        if '*' in pattern or '?' in pattern:
            accel_files.extend(glob.glob(pattern))
        else:
            accel_files.append(pattern)

    if not accel_files:
        print("ERROR: No ACCEL files found")
        sys.exit(1)

    if args.verbose:
        print(f"Processing {len(accel_files)} ACCEL files...")

    # Sift candidates
    candidates = sift_candidates(
        accel_files,
        min_sigma=args.min_sigma,
        min_period=args.min_period,
        max_period=args.max_period,
        remove_dups=not args.no_duplicates,
        remove_harms=not args.no_harmonics
    )

    if args.verbose:
        print(f"Found {len(candidates)} candidates after sifting")

    # Limit number of candidates
    if args.max_cands and len(candidates) > args.max_cands:
        candidates = candidates[:args.max_cands]
        if args.verbose:
            print(f"Limited to top {args.max_cands} candidates")

    # Write outputs
    base_output = args.output.replace('.candfile', '').replace('.csv', '')
    candfile = f"{base_output}.candfile"
    csvfile = f"{base_output}.csv"

    write_candfile(candidates, candfile, dm=args.dm)
    write_csv(candidates, csvfile)

    print(f"Wrote {len(candidates)} candidates to:")
    print(f"  {candfile}")
    print(f"  {csvfile}")


if __name__ == '__main__':
    main()
