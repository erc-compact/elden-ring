#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
from sigpyproc.readers import FilReader


def freq_stack(infiles, outfile, block_size=1024):
    """Stack multiple .fil files in frequency, time-aligned with proper flipping."""
    # 1) Open all inputs
    readers = [FilReader(f) for f in infiles]
    headers = [r.header for r in readers]

    # 2) Validate consistent channel width (absolute value)
    foff_abs_vals = [abs(h.foff) for h in headers]
    if len(set(foff_abs_vals)) > 1:
        sys.exit("ERROR: All files must have the same channel width (|foff|)")
    channel_width = foff_abs_vals[0]

    # Use first file's foff sign as reference
    global_foff = headers[0].foff

    # 3) Time alignment
    tstarts = [h.tstart for h in headers]
    tsamps = [h.tsamp for h in headers]
    ends = [
        t0 + (h.nsamples - 1) * dt / 86400.0
        for h, t0, dt in zip(headers, tstarts, tsamps)
    ]
    mjd_start = max(tstarts)
    mjd_end = min(ends)
    if mjd_end <= mjd_start:
        sys.exit("ERROR: No overlapping time window!")

    # 4) Compute sample skips & common length
    skips = [
        int(round((mjd_start - t0) * 86400.0 / dt)) for t0, dt in zip(tstarts, tsamps)
    ]
    avails = [h.nsamples - sk for h, sk in zip(headers, skips)]
    total = min(avails)

    # 5) Build global frequency axis
    all_freqs = np.hstack([h.chan_freqs for h in headers])
    fmin, fmax = all_freqs.min(), all_freqs.max()
    nchan_total = int(round(abs((fmax - fmin) / channel_width))) + 1
    is_ascending = global_foff > 0
    new_fch1 = fmin if is_ascending else fmax

    # 6) Prepare output
    new_hdr = headers[0].new_header(
        {
            "nchans": nchan_total,
            "foff": global_foff,
            "fch1": new_fch1,
            "nsamples": total,
        }
    )

    writer = new_hdr.prep_outfile(outfile)

    # 7) Determine flipping needs and calculate FINAL offsets
    global_freqs = new_hdr.chan_freqs
    file_props = []
    for r in readers:
        # Determine if flipping is needed based on foff sign comparison
        needs_flip = (r.header.foff > 0) != (global_foff > 0)

        if needs_flip:
            start_freq = r.header.chan_freqs[-1]  # Last becomes first after flip
        else:
            start_freq = r.header.fch1

        offset = int(round(abs(start_freq - new_fch1) / channel_width))
        file_props.append(
            {
                "reader": r,
                "offset": offset,
                "needs_flip": needs_flip,
                "nchans": r.header.nchans,
            }
        )
        print(
            f"[INFO] File {os.path.basename(r.filename)}: "
            f"offset={offset}, flip={'YES' if needs_flip else 'NO'}"
        )

    # 8) Process data with correct offsets
    written = 0
    while written < total:
        nread = min(block_size, total - written)
        merged = np.zeros((nread, nchan_total), dtype=np.uint8)

        for props in file_props:
            r = props["reader"]
            arr = r.read_block(skips[readers.index(r)] + written, nread).data.T

            if props["needs_flip"]:
                arr = np.flip(arr, axis=1)

            start = props["offset"]
            end = start + props["nchans"]

            # Handle potential overflow
            if end > nchan_total:
                actual_end = nchan_total
                arr = arr[:, : actual_end - start]
            else:
                actual_end = end

            merged[:, start:actual_end] = arr

        writer.cwrite(merged.flatten("C"))
        written += nread
        print(
            f"[PROGRESS] Written {written}/{total} samples ({written / total:.1%})",
            end="\r",
        )

    writer.close()
    print(f"\nSUCCESS: Created {outfile} ({total} samples, {nchan_total} channels)")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Frequency-stack filterbank files")
    p.add_argument("infiles", nargs="+", help="Input .fil files")
    p.add_argument("-o", "--output", required=True, help="Output file")
    p.add_argument(
        "-B", "--block-size", type=int, default=1024, help="Samples per block"
    )
    args = p.parse_args()

    freq_stack(args.infiles, args.output, args.block_size)

