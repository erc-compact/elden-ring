#!/usr/bin/env python3
"""
Sift pulsar candidates using PRESTO's sifting module (NEXTO-style).
"""
import sys
import argparse
from operator import attrgetter
import glob
import os
import sifting

SOL = 299792458.0

def main():
    parser = argparse.ArgumentParser(
        description='Sift pulsar candidates using PRESTO sifting module'
    )
    parser.add_argument('accel_files', nargs='+',
                        help='ACCEL candidate files to sift')
    parser.add_argument('--min-period', type=float, default=0.0005,
                        help='Minimum period in seconds (default: 0.0005)')
    parser.add_argument('--max-period', type=float, default=15.0,
                        help='Maximum period in seconds (default: 15.0)')
    parser.add_argument('--sigma-threshold', type=float, default=6.0,
                        help='Minimum sigma threshold (default: 6.0)')
    parser.add_argument('--cpow-threshold', type=float, default=100.0,
                        help='Minimum coherent power threshold (default: 100.0)')
    parser.add_argument('--min-num-dms', type=int, default=1,
                        help='Minimum number of DMs for detection (default: 1)')
    parser.add_argument('--low-dm-cutoff', type=float, default=2.0,
                        help='Lowest DM to consider as real pulsar (default: 2.0)')
    parser.add_argument('--remove-duplicates', action='store_true',
                        help='Remove duplicate candidates')
    parser.add_argument('--remove-harmonics', action='store_true',
                        help='Remove harmonic candidates')
    parser.add_argument('--max-cands-to-fold', type=int, default=100,
                        help='Maximum candidates to fold (default: 100)')
    parser.add_argument('--output', default='sifted_candidates',
                        help='Output file prefix (no extension)')
    parser.add_argument('--tobs', type=float, default=None,
                        help='Observation time in seconds (auto-detected if not provided)')

    args = parser.parse_args()

    # Expand globs
    accel_files = []
    for pattern in args.accel_files:
        if '*' in pattern or '?' in pattern:
            accel_files.extend(glob.glob(pattern))
        else:
            accel_files.append(pattern)

    if not accel_files:
        print("No candidates found!")
        open(f"{args.output}.candfile", 'w').close()
        open(f"{args.output}.csv", 'w').close()
        return

    sifting.sigma_threshold = args.sigma_threshold
    sifting.c_pow_threshold = args.cpow_threshold
    sifting.short_period = args.min_period
    sifting.long_period = args.max_period
    sifting.r_err = 1.1
    sifting.harm_pow_cutoff = 8.0
    sifting.known_birds_p = []
    sifting.known_birds_f = []

    cands = sifting.read_candidates(accel_files)
    if not cands:
        print("No candidates found after reading ACCEL files")
        open(f"{args.output}.candfile", 'w').close()
        open(f"{args.output}.csv", 'w').close()
        return

    dm_set = {cand.DM for cand in cands}
    dmstrs = ["%.2f" % dm for dm in sorted(dm_set)]

    if args.remove_duplicates and len(cands):
        cands = sifting.remove_duplicate_candidates(cands)

    if len(cands) and args.min_num_dms > 1:
        cands = sifting.remove_DM_problems(cands, args.min_num_dms, dmstrs, args.low_dm_cutoff)

    if args.remove_harmonics and len(cands):
        cands = sifting.remove_harmonics(cands)

    if not cands:
        open(f"{args.output}.candfile", 'w').close()
        open(f"{args.output}.csv", 'w').close()
        return

    cands.sort(key=attrgetter('sigma'), reverse=True)

    if args.tobs is not None:
        tobs = args.tobs
    else:
        tobs = cands[0].T if hasattr(cands[0], 'T') else 300.0

    # Limit output
    if args.max_cands_to_fold and len(cands) > args.max_cands_to_fold:
        cands = cands[:args.max_cands_to_fold]

    candfile = f"{args.output}.candfile"
    csvfile = f"{args.output}.csv"

    with open(candfile, 'w') as cf, open(csvfile, 'w') as csvf:
        cf.write("#id dm acc F0 F1 F2 S/N\n")
        csvf.write("id,dm,f0,f1,f2,snr,sigma,accel,z,w,period_ms,file\n")
        for k, cand in enumerate(cands, 1):
            z0 = cand.z - 0.5 * cand.w
            r0 = cand.r - 0.5 * z0 - cand.w / 6.0
            f = r0 / cand.T
            fd = z0 / (cand.T * cand.T)
            fdd = cand.w / (cand.T * cand.T * cand.T)
            period_ms = (1.0 / f) * 1000.0 if f != 0 else 0.0
            cf.write("%d\t%.3f\t%.15f\t%.15f\t%.15f\t%.15f\t%.2f\n" % (k, cand.DM, 0., f, fd, fdd, cand.snr))
            csvf.write(f"{k},{cand.DM:.3f},{f:.15f},{fd:.15e},{fdd:.15e},{cand.snr:.2f},{cand.sigma:.2f},{0.0:.6f},{cand.z:.2f},{cand.w:.2f},{period_ms:.3f},{os.path.basename(cand.filename)}\n")

if __name__ == '__main__':
    main()
