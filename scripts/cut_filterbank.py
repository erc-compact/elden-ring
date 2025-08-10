#!/usr/bin/env python3
import os
import argparse
import numpy as np
from sigpyproc.readers import FilReader


def split_filterbank_plan(
    infile: str, cutoff_mhz: float, out_low: str, out_high: str, gulp: int = 16384
):
    """
    Split a filterbank file into two at a specified cutoff frequency.
    Written to split the lowest baseband file from UBB into high and low bands.
    """
    # Open input and read header
    fb = FilReader(infile)
    hdr = fb.header
    fch1, foff, nchans = hdr.fch1, hdr.foff, hdr.nchans
    total = hdr.nsamples

    # Compute channel split index
    if foff >= 0:
        raw_idx = (cutoff_mhz - fch1) / foff
        idx = int(np.floor(raw_idx))
    else:
        raw_idx = (cutoff_mhz - fch1) / foff
        idx = int(np.ceil(raw_idx))
    idx = max(0, min(nchans - 1, idx))

    # Determine counts for low/high channels
    if foff >= 0:
        n_low = idx + 1
        n_high = nchans - n_low
    else:
        n_low = nchans - idx
        n_high = nchans - n_low

    # Compute starting frequency for outputs
    low_fch1 = fch1 if foff >= 0 else fch1 + n_high * foff
    high_fch1 = fch1 + n_low * foff if foff >= 0 else fch1

    # Prepare output headers and files
    hdr_low = hdr.new_header({"nchans": n_low, "fch1": low_fch1})
    hdr_high = hdr.new_header({"nchans": n_high, "fch1": high_fch1})
    w_low = hdr_low.prep_outfile(out_low)
    w_high = hdr_high.prep_outfile(out_high)

    # Process data in gulps
    for nsamps_r, _, data in fb.read_plan(start=0, nsamps=total, gulp=gulp):
        # Reshape into (samples, channels) to match sample-major order
        samp_chan = np.asarray(data).reshape(nsamps_r, nchans)

        # Slice frequency channels
        if foff >= 0:
            low_mat = samp_chan[:, :n_low]
            high_mat = samp_chan[:, n_low:]
        else:
            low_mat = samp_chan[:, n_high:]
            high_mat = samp_chan[:, :n_high]

        # Write back in sample-major order
        w_low.cwrite(low_mat.flatten("C"))
        w_high.cwrite(high_mat.flatten("C"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a filterbank file by frequency."
    )
    parser.add_argument("-i", "--infile", required=True, help="Input .fil file path")
    parser.add_argument(
        "-c", "--cutoff", type=float, required=True, help="Cutoff frequency in MHz"
    )
    parser.add_argument(
        "-l", "--out_low", required=True, help="Output file for lower band"
    )
    parser.add_argument(
        "-u", "--out_high", required=True, help="Output file for upper band"
    )
    parser.add_argument(
        "-g", "--gulp", type=int, default=16384, help="Number of samples per gulp"
    )
    args = parser.parse_args()

    # Remove existing outputs if present
    for fn in (args.out_low, args.out_high):
        if os.path.exists(fn):
            os.remove(fn)

    split_filterbank_plan(
        args.infile, args.cutoff, args.out_low, args.out_high, args.gulp
    )
    print("Success: Files split correctly.")
