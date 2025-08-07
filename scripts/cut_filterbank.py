#!/usr/bin/env python3
import os
import argparse
import numpy as np
from sigpyproc.readers import FilReader

def split_filterbank_plan(
    infile: str,
    cutoff_mhz: float,
    out_low: str,
    out_high: str,
    gulp: int = 16384
):
    """
    Split a filterbank file into two at a specified cutoff frequency.
    Written to split the lowest baseband file from UBB into high and low 1300MHz
    
    Author: Fazal Kareem
    Date: 18.05.2025
    """
    fb = FilReader(infile)
    hdr = fb.header
    fch1, foff, nchans = hdr.fch1, hdr.foff, hdr.nchans
    total = hdr.nsamples

    # Compute channel split index
    if foff >= 0:
        idx = int(np.floor((cutoff_mhz - fch1) / foff))
    else:
        idx = int(np.ceil((cutoff_mhz - fch1) / foff))
    idx = max(0, min(nchans - 1, idx))

    # Determine low/high channels based on frequency order
    n_low = idx + 1 if foff >= 0 else nchans - idx
    n_high = nchans - n_low

    # Adjust fch1 for output headers
    low_fch1 = fch1 if foff >= 0 else fch1 + n_high * foff
    high_fch1 = (fch1 + n_low * foff) if foff >= 0 else fch1

    # Prepare output files
    hdr_low = hdr.new_header({"nchans": n_low, "fch1": low_fch1})
    hdr_high = hdr.new_header({"nchans": n_high, "fch1": high_fch1})
    w_low = hdr_low.prep_outfile(out_low)
    w_high = hdr_high.prep_outfile(out_high)

    # Process data in chunks (read_plan handles partial final gulp)
    for nsamps_r, _, data in fb.read_plan(start=0, nsamps=total, gulp=gulp):
        block = np.asarray(data).reshape(nchans, nsamps_r)
        if foff >= 0:
            low_block = block[:n_low, :]
            high_block = block[n_low:, :]
        else:
            low_block = block[n_high:, :]
            high_block = block[:n_high, :]

        w_low.cwrite(low_block.flatten("C"))
        w_high.cwrite(high_block.flatten("C"))

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Split a filterbank file by frequency.")
    p.add_argument("-i", "--infile", required=True)
    p.add_argument("-c", "--cutoff", type=float, required=True)
    p.add_argument("-l", "--out_low", required=True)
    p.add_argument("-u", "--out_high", required=True)
    p.add_argument("-g", "--gulp", type=int, default=16384)
    args = p.parse_args()

    for fn in (args.out_low, args.out_high):
        if os.path.exists(fn):
            os.remove(fn)

    split_filterbank_plan(args.infile, args.cutoff, args.out_low, args.out_high, args.gulp)
    print("Success: Files split without overflow errors.")