#!/usr/bin/env python3
"""
Convert PRESTO .pfd files to PNG plots.
This script reads pfd files and generates diagnostic PNG plots similar to prepfold output.
"""

import argparse
import os
import sys
import subprocess
import glob


def convert_pfd_to_png(pfd_file, output_dir=None):
    """
    Convert a single pfd file to PNG using PRESTO's show_pfd.

    Args:
        pfd_file: Path to the .pfd file
        output_dir: Output directory for PNG files (default: same as pfd)

    Returns:
        Path to the generated PNG file or None on failure
    """
    if not os.path.exists(pfd_file):
        print(f"ERROR: PFD file not found: {pfd_file}")
        return None

    base_name = os.path.basename(pfd_file).replace('.pfd', '')

    if output_dir is None:
        output_dir = os.path.dirname(pfd_file) or '.'

    os.makedirs(output_dir, exist_ok=True)

    # Generate PNG using show_pfd (non-interactive mode)
    # show_pfd generates a .ps file, then we convert to PNG
    ps_file = os.path.join(output_dir, f"{base_name}.pfd.ps")
    png_file = os.path.join(output_dir, f"{base_name}.png")

    try:
        # Run show_pfd to generate PostScript file
        cmd = ['show_pfd', '-noxwin', pfd_file]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)

        if result.returncode != 0:
            print(f"WARNING: show_pfd failed for {pfd_file}: {result.stderr}")
            # Try alternative method using prepfold's built-in plotting
            return None

        # Find the generated .ps file
        expected_ps = pfd_file + '.ps'
        if not os.path.exists(expected_ps):
            # Try looking in output_dir
            ps_files = glob.glob(os.path.join(output_dir, '*.ps'))
            if ps_files:
                expected_ps = ps_files[0]
            else:
                print(f"WARNING: No PS file generated for {pfd_file}")
                return None

        # Convert PS to PNG using ghostscript or convert
        try:
            # Try ghostscript first
            gs_cmd = [
                'gs', '-dSAFER', '-dBATCH', '-dNOPAUSE',
                '-sDEVICE=png16m', '-r150',
                f'-sOutputFile={png_file}',
                expected_ps
            ]
            result = subprocess.run(gs_cmd, capture_output=True, text=True)

            if result.returncode != 0:
                # Try ImageMagick convert as fallback
                convert_cmd = ['convert', '-density', '150', expected_ps, png_file]
                result = subprocess.run(convert_cmd, capture_output=True, text=True)

                if result.returncode != 0:
                    print(f"WARNING: Could not convert PS to PNG for {pfd_file}")
                    return expected_ps  # Return PS file path at least
        except FileNotFoundError:
            # Try convert directly
            try:
                convert_cmd = ['convert', '-density', '150', expected_ps, png_file]
                result = subprocess.run(convert_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    return expected_ps
            except FileNotFoundError:
                print(f"WARNING: Neither ghostscript nor ImageMagick found")
                return expected_ps

        # Clean up PS file if PNG was created
        if os.path.exists(png_file):
            if os.path.exists(expected_ps):
                os.remove(expected_ps)
            return png_file

        return expected_ps

    except Exception as e:
        print(f"ERROR converting {pfd_file}: {e}")
        return None


def batch_convert(pfd_files, output_dir=None, parallel=1):
    """
    Convert multiple pfd files to PNG.

    Args:
        pfd_files: List of pfd file paths
        output_dir: Output directory for PNG files
        parallel: Number of parallel conversions

    Returns:
        List of (pfd_path, png_path) tuples
    """
    results = []

    if parallel > 1:
        from multiprocessing import Pool

        def convert_wrapper(pfd_file):
            return (pfd_file, convert_pfd_to_png(pfd_file, output_dir))

        with Pool(parallel) as pool:
            results = pool.map(convert_wrapper, pfd_files)
    else:
        for pfd_file in pfd_files:
            png_path = convert_pfd_to_png(pfd_file, output_dir)
            results.append((pfd_file, png_path))

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Convert PRESTO .pfd files to PNG plots'
    )
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input .pfd file(s) or glob pattern'
    )
    parser.add_argument(
        '-o', '--output-dir',
        default=None,
        help='Output directory for PNG files (default: same as input)'
    )
    parser.add_argument(
        '-p', '--parallel',
        type=int,
        default=1,
        help='Number of parallel conversions (default: 1)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Expand glob patterns
    pfd_files = []
    for pattern in args.input:
        if '*' in pattern or '?' in pattern:
            pfd_files.extend(glob.glob(pattern))
        else:
            pfd_files.append(pattern)

    if not pfd_files:
        print("ERROR: No .pfd files found")
        sys.exit(1)

    if args.verbose:
        print(f"Converting {len(pfd_files)} PFD files...")

    results = batch_convert(pfd_files, args.output_dir, args.parallel)

    # Print summary
    success = sum(1 for _, png in results if png is not None)
    print(f"Converted {success}/{len(results)} files successfully")

    if args.verbose:
        for pfd, png in results:
            status = "OK" if png else "FAILED"
            print(f"  {os.path.basename(pfd)}: {status}")

    # Return non-zero if any failures
    sys.exit(0 if success == len(results) else 1)


if __name__ == '__main__':
    main()
