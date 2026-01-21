#!/usr/bin/env python3
"""
Convert PRESTO PostScript (.ps) files to PNG plots.

prepfold outputs .pfd.ps files directly, so this script simply converts
the existing PS files to PNG format using ghostscript or ImageMagick.
"""

import argparse
import os
import sys
import subprocess
import glob
from multiprocessing import Pool


def convert_ps_to_png(ps_file, output_dir=None):
    """
    Convert a single PostScript file to PNG.

    Args:
        ps_file: Path to the .ps file
        output_dir: Output directory for PNG files (default: same as ps)

    Returns:
        Path to the generated PNG file or None on failure
    """
    if not os.path.exists(ps_file):
        print(f"ERROR: PS file not found: {ps_file}")
        return None

    # Generate PNG filename - remove .pfd.ps or .ps suffix
    base_name = os.path.basename(ps_file)
    if base_name.endswith('.pfd.ps'):
        base_name = base_name[:-7]  # Remove .pfd.ps
    elif base_name.endswith('.ps'):
        base_name = base_name[:-3]  # Remove .ps

    if output_dir is None:
        output_dir = os.path.dirname(ps_file) or '.'

    os.makedirs(output_dir, exist_ok=True)
    png_file = os.path.join(output_dir, f"{base_name}.png")

    try:
        # Try ghostscript first (preferred for quality)
        gs_cmd = [
            'gs', '-dSAFER', '-dBATCH', '-dNOPAUSE',
            '-sDEVICE=png16m', '-r150',
            f'-sOutputFile={png_file}',
            ps_file
        ]
        result = subprocess.run(gs_cmd, capture_output=True, text=True)

        if result.returncode == 0 and os.path.exists(png_file):
            return png_file

        # Fallback to ImageMagick convert
        convert_cmd = ['convert', '-density', '150', ps_file, png_file]
        result = subprocess.run(convert_cmd, capture_output=True, text=True)

        if result.returncode == 0 and os.path.exists(png_file):
            return png_file

        print(f"WARNING: Could not convert {ps_file} to PNG: {result.stderr}")
        return None

    except FileNotFoundError as e:
        print(f"ERROR: Required converter not found: {e}")
        return None
    except Exception as e:
        print(f"ERROR converting {ps_file}: {e}")
        return None


def batch_convert(ps_files, output_dir=None, workers=1):
    """
    Convert multiple PS files to PNG.

    Args:
        ps_files: List of PS file paths
        output_dir: Output directory for PNG files
        workers: Number of parallel workers

    Returns:
        List of (ps_path, png_path) tuples
    """
    if workers > 1:
        def convert_wrapper(ps_file):
            return (ps_file, convert_ps_to_png(ps_file, output_dir))

        with Pool(workers) as pool:
            return pool.map(convert_wrapper, ps_files)
    else:
        results = []
        for ps_file in ps_files:
            png_path = convert_ps_to_png(ps_file, output_dir)
            results.append((ps_file, png_path))
        return results


def main():
    parser = argparse.ArgumentParser(
        description='Convert PRESTO PostScript (.ps) files to PNG plots'
    )
    parser.add_argument(
        '--input_dir',
        default='.',
        help='Input directory containing .ps files (default: current dir)'
    )
    parser.add_argument(
        '--output_dir',
        default='.',
        help='Output directory for PNG files (default: current dir)'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=4,
        help='Number of parallel workers (default: 4)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find all PS files - look for .pfd.ps specifically (prepfold output)
    ps_files = glob.glob(os.path.join(args.input_dir, '*.pfd.ps'))
    if not ps_files:
        # Fallback to any .ps files
        ps_files = glob.glob(os.path.join(args.input_dir, '*.ps'))

    if not ps_files:
        print("WARNING: No .ps files found in input directory")
        # Create empty placeholder to satisfy Nextflow output requirements
        with open(os.path.join(args.output_dir, 'no_plots.txt'), 'w') as f:
            f.write("No PostScript files found to convert\n")
        sys.exit(0)

    if args.verbose:
        print(f"Converting {len(ps_files)} PS files to PNG...")

    results = batch_convert(ps_files, args.output_dir, args.workers)

    # Print summary
    success = sum(1 for _, png in results if png is not None)
    print(f"Converted {success}/{len(results)} files successfully")

    if args.verbose:
        for ps, png in results:
            status = "OK" if png else "FAILED"
            print(f"  {os.path.basename(ps)}: {status}")

    # Return non-zero only if ALL conversions failed
    sys.exit(0 if success > 0 else 1)


if __name__ == '__main__':
    main()
