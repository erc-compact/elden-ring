#!/usr/bin/env python3
"""
Generate meta file for CSV-based folding (PRESTO/peasoup candidates).

Reads filterbank header and creates a meta file with observation parameters.

Usage:
    python3 presto_generate_fold_meta.py \
        --fil-file input.fil \
        --output fold_meta.txt \
        --fft-size 16777216 \
        --telescope meerkat \
        --source-name NGC6544 \
        --beam-name cfbf00000 \
        --beam-id 0 \
        --utc-start 2024-01-01 \
        --ra 12:00:00 \
        --dec -30:00:00 \
        --cdm 30.0
"""

import argparse
import struct


def read_fil_header(filepath):
    """Read filterbank header to extract tstart, tsamp, nsamples."""
    header = {}
    try:
        with open(filepath, 'rb') as f:
            # Skip to find header keywords
            f.seek(0)
            content = f.read(4096)  # Read first 4KB

            # Parse common keywords
            if b'tstart' in content:
                idx = content.find(b'tstart')
                f.seek(idx + 6)
                header['tstart'] = struct.unpack('d', f.read(8))[0]

            if b'tsamp' in content:
                idx = content.find(b'tsamp')
                f.seek(idx + 5)
                header['tsamp'] = struct.unpack('d', f.read(8))[0]

            if b'nsamples' in content:
                idx = content.find(b'nsamples')
                f.seek(idx + 8)
                header['nsamples'] = struct.unpack('i', f.read(4))[0]
    except Exception:
        pass

    return header


def main():
    parser = argparse.ArgumentParser(description="Generate fold meta file from filterbank")
    parser.add_argument("--fil-file", required=True, help="Input filterbank file")
    parser.add_argument("--output", default="fold_meta.txt", help="Output meta file")
    parser.add_argument("--fft-size", type=int, default=0, help="FFT size")
    parser.add_argument("--telescope", default="meerkat", help="Telescope name")
    parser.add_argument("--source-name", default="", help="Source/cluster name")
    parser.add_argument("--beam-name", default="", help="Beam name")
    parser.add_argument("--beam-id", default="", help="Beam ID")
    parser.add_argument("--utc-start", default="", help="UTC start time")
    parser.add_argument("--ra", default="", help="Right ascension")
    parser.add_argument("--dec", default="", help="Declination")
    parser.add_argument("--cdm", type=float, default=0.0, help="Coherent DM")
    parser.add_argument("--subint-length", type=float, default=10.0, help="Subintegration length")
    parser.add_argument("--nsubband", type=int, default=64, help="Number of subbands")
    parser.add_argument("--clfd", type=float, default=2.0, help="CLFD Q value")
    parser.add_argument("--nbins", type=int, default=64, help="Number of bins")
    parser.add_argument("--binplan", default="0.005 32 0.01 64 0.1 128", help="Bin plan")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    parser.add_argument("--template", default="", help="Template file path")

    args = parser.parse_args()

    # Try to read header, use defaults if fails
    header = read_fil_header(args.fil_file)
    tstart = header.get('tstart', 0.0)
    tsamp = header.get('tsamp', 64e-6)
    nsamples = header.get('nsamples', 0)

    # Write meta file
    with open(args.output, 'w') as f:
        f.write(f"filterbank_file: {args.fil_file}\n")
        f.write(f"fft_size: {args.fft_size}\n")
        f.write(f"tsamp: {tsamp}\n")
        f.write(f"tstart: {tstart}\n")
        f.write(f"nsamples: {nsamples}\n")
        f.write(f"source_name: {args.source_name}\n")
        f.write(f"telescope: {args.telescope}\n")
        f.write(f"beam_name: {args.beam_name}\n")
        f.write(f"beam_id: {args.beam_id}\n")
        f.write(f"utc_beam: {args.utc_start}\n")
        f.write(f"ra: {args.ra}\n")
        f.write(f"dec: {args.dec}\n")
        f.write(f"cdm: {args.cdm}\n")
        f.write(f"subint_length: {args.subint_length}\n")
        f.write(f"nsubband: {args.nsubband}\n")
        f.write(f"clfd_q_value: {args.clfd}\n")
        f.write(f"nbins: {args.nbins}\n")
        f.write(f"binplan: {args.binplan}\n")
        f.write(f"threads: {args.threads}\n")
        f.write(f"template: {args.template}\n")
        f.write("start_fraction: 0.0\n")
        f.write("end_fraction: 1.0\n")
        f.write("chunk_id: 0\n")
        f.write("cmask: \n")
        f.write("rfi_filter: \n")

    print(f"Meta file created with tstart={tstart}, tsamp={tsamp}")


if __name__ == "__main__":
    main()
