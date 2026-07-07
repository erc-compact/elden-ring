#!/usr/bin/env python3
"""
Split a CandyJar candidate tarball into N tarballs, one per person, for
distributed manual classification.

The input is an existing CandyJar tarball (as produced by
create_candyjar_tarball.py) containing:
    candidates.csv          top-level merged candidate CSV
    plots/*.png             one folded-candidate plot per row (row.png_path)
    metafiles/*.meta        one metafile per utc_start (row.metafile_path)

Candidates are distributed using STRATIFIED ROUND-ROBIN on sn_fold: rows are
sorted by SNR and dealt one-by-one across the N sets, so every person receives
a near-identical SNR distribution (nobody gets only high- or only low-SNR
candidates). Distribution is deterministic for a given --seed.

Overlap:
    --overlap        Each candidate is seen by exactly 2 people (assigned to 2
                     sets instead of 1), for inter-rater cross-checking.
    (default)        Disjoint partition -- each candidate assigned to 1 person.

Each output set is repackaged as its own CandyJar tarball with only the
candidates.csv rows, plots and metafiles it needs.

Usage:
    python split_candidates_for_classification.py \
        -i candidates.tar.gz -n 4 -o /path/to/output/ \
        [--overlap] [--min-snr 8] [--meta new.meta] [--seed 42]

Output:
    <output_dir>/<prefix>_set_1.tar.gz ... _set_N.tar.gz
"""

import argparse
import logging
import os
import sys
import tarfile
import tempfile
import shutil

import pandas as pd


def setup_logging(verbose: bool) -> logging.Logger:
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if verbose else logging.INFO,
    )
    return logging.getLogger(__name__)


def extract_tarball(tarball_path: str, dest: str, logger: logging.Logger) -> str:
    """Extract the input tarball into dest. Returns the directory that holds
    candidates.csv (handles tarballs with or without a top-level folder)."""
    logger.info("Extracting input tarball: %s", tarball_path)
    with tarfile.open(tarball_path, "r:*") as tar:
        tar.extractall(dest)

    # candidates.csv may be at dest/ or dest/<somedir>/
    for root, _dirs, files in os.walk(dest):
        if "candidates.csv" in files:
            logger.debug("Found candidates.csv in %s", root)
            return root
    logger.error("candidates.csv not found inside %s", tarball_path)
    sys.exit(1)


def assign_sets(df: pd.DataFrame, n_sets: int, overlap: bool,
                logger: logging.Logger) -> dict:
    """Stratified round-robin assignment of rows to sets.

    Returns {set_index (1-based): [row positional indices]}.
    With overlap=True each row is dealt to two consecutive sets, so every
    candidate is seen by exactly 2 people.
    """
    if "sn_fold" not in df.columns:
        logger.error("Input CSV has no 'sn_fold' column; cannot stratify by SNR.")
        sys.exit(1)

    # Sort by SNR so the round-robin deal is stratified across the SNR range.
    order = df["sn_fold"].sort_values(kind="mergesort").index.tolist()

    sets = {i: [] for i in range(1, n_sets + 1)}
    copies = 2 if overlap else 1
    if overlap and n_sets < 2:
        logger.error("--overlap requires -n >= 2 (need 2 people per candidate).")
        sys.exit(1)

    for deal_pos, row_idx in enumerate(order):
        for c in range(copies):
            # Offset the second copy by 1 so the pair lands on different people,
            # while still rotating through all sets evenly.
            target = ((deal_pos + c) % n_sets) + 1
            sets[target].append(row_idx)

    for i, rows in sets.items():
        logger.info("Set %d: %d candidates", i, len(rows))
    return sets


def collect_member_files(subset: pd.DataFrame, src_root: str, meta_override: str,
                         logger: logging.Logger) -> (list, list):
    """Resolve the plots and metafiles referenced by a subset, relative to the
    extracted tarball root.

    Returns (png_paths, meta_members) where meta_members is a list of
    (arcname_basename, source_path) pairs -- arcname_basename is the name the
    metafile must keep inside the tarball (so CandyJar's metafile_path resolves),
    source_path is where its content is read from.

    --meta override: replace the metafile content with the provided file.
      - single metafile in this set -> use the provided content under that name.
      - multiple metafiles         -> replace only the one whose name matches the
                                       provided file's basename (utc match); others
                                       keep their original content.
    """
    png_paths = []
    for rel in subset["png_path"].dropna().unique():
        p = os.path.join(src_root, rel)
        if os.path.isfile(p):
            png_paths.append(p)
        else:
            logger.warning("Plot missing in source tarball: %s", rel)

    # Collect referenced metafile arcnames and their original source paths.
    referenced = []  # list of (arcname_basename, original_source_path_or_None)
    if "metafile_path" in subset.columns:
        for rel in subset["metafile_path"].dropna().unique():
            arc = os.path.basename(rel)
            src = os.path.join(src_root, rel)
            referenced.append((arc, src if os.path.isfile(src) else None))

    meta_members = []
    override_base = os.path.basename(meta_override) if meta_override else None
    single = len(referenced) == 1

    for arc, orig in referenced:
        if meta_override and (single or arc == override_base):
            # Use provided content, keep the existing arcname.
            meta_members.append((arc, meta_override))
            logger.info("Replacing metafile '%s' with provided --meta content.", arc)
        elif orig is not None:
            meta_members.append((arc, orig))
        else:
            logger.warning("Metafile missing in source tarball: %s", arc)

    if meta_override and not single and override_base not in {a for a, _ in referenced}:
        logger.warning(
            "--meta basename '%s' matched no metafile in this set (%d metafiles); "
            "left originals unchanged.", override_base, len(referenced))

    return png_paths, meta_members


def create_set_tarball(subset: pd.DataFrame, png_paths: list, meta_members: list,
                       tarball_path: str, logger: logging.Logger) -> None:
    """Write one CandyJar tarball: candidates.csv + plots/ + metafiles/.

    meta_members is a list of (arcname_basename, source_path) -- content comes
    from source_path but the name inside metafiles/ stays arcname_basename.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as tmp:
        subset.to_csv(tmp.name, index=False)
        csv_tmp = tmp.name
    try:
        with tarfile.open(tarball_path, "w:gz", dereference=True) as tar:
            tar.add(csv_tmp, arcname="candidates.csv")
            for png in png_paths:
                tar.add(png, arcname=os.path.join("plots", os.path.basename(png)))
            for arc, src in meta_members:
                tar.add(src, arcname=os.path.join("metafiles", arc))
        logger.info("Created %s (%d cands, %d plots, %d metafiles)",
                    tarball_path, len(subset), len(png_paths), len(meta_members))
    finally:
        os.unlink(csv_tmp)


def parse_arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Split a CandyJar tarball into N tarballs for distributed classification."
    )
    p.add_argument("-i", "--input_tarball", required=True,
                   help="Path to the existing CandyJar .tar.gz to split.")
    p.add_argument("-n", "--nsets", type=int, required=True,
                   help="Number of people / output tarballs.")
    p.add_argument("-o", "--output_path", default=".",
                   help="Directory to write the output tarballs (default: cwd).")
    p.add_argument("--overlap", action="store_true",
                   help="Assign every candidate to exactly 2 people for cross-checking.")
    p.add_argument("--min-snr", type=float, default=None, dest="min_snr",
                   help="Drop candidates with sn_fold below this value before splitting.")
    p.add_argument("--meta", default=None,
                   help="Replace metafile content with this file, keeping the existing "
                        "metafile name(s). Single metafile per set -> content copied; "
                        "multiple -> matched to the metafile whose name equals this "
                        "file's basename (utc match).")
    p.add_argument("--seed", type=int, default=42,
                   help="Seed for deterministic assignment (default: 42).")
    p.add_argument("--prefix", default=None,
                   help="Output tarball name prefix (default: input tarball basename).")
    p.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return p.parse_args()


def main():
    args = parse_arguments()
    logger = setup_logging(args.verbose)

    if args.nsets < 1:
        logger.error("-n must be >= 1")
        sys.exit(1)
    if not os.path.isfile(args.input_tarball):
        logger.error("Input tarball not found: %s", args.input_tarball)
        sys.exit(1)
    if args.meta is not None and not os.path.isfile(args.meta):
        logger.error("--meta file not found: %s", args.meta)
        sys.exit(1)

    os.makedirs(args.output_path, exist_ok=True)
    prefix = args.prefix or os.path.basename(args.input_tarball)
    for ext in (".tar.gz", ".tgz", ".tar"):
        if prefix.endswith(ext):
            prefix = prefix[: -len(ext)]
            break

    workdir = tempfile.mkdtemp(prefix="split_cands_")
    try:
        src_root = extract_tarball(args.input_tarball, workdir, logger)
        df = pd.read_csv(os.path.join(src_root, "candidates.csv"))
        logger.info("Loaded %d candidates from tarball.", len(df))
        if df.empty:
            logger.error("candidates.csv is empty; nothing to split.")
            sys.exit(1)

        if args.min_snr is not None:
            if "sn_fold" not in df.columns:
                logger.error("--min-snr given but no 'sn_fold' column in candidates.csv.")
                sys.exit(1)
            before = len(df)
            df = df[df["sn_fold"] >= args.min_snr].reset_index(drop=True)
            logger.info("Applied --min-snr %.2f: kept %d of %d candidates.",
                        args.min_snr, len(df), before)
            if df.empty:
                logger.error("No candidates remain after --min-snr filter; nothing to split.")
                sys.exit(1)

        # Deterministic ordering for reproducibility, then stratified deal.
        df = df.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
        sets = assign_sets(df, args.nsets, args.overlap, logger)

        for i, row_idx in sets.items():
            subset = df.loc[row_idx].reset_index(drop=True)
            png_paths, meta_members = collect_member_files(subset, src_root, args.meta, logger)
            out = os.path.join(args.output_path, f"{prefix}_set_{i}.tar.gz")
            create_set_tarball(subset, png_paths, meta_members, out, logger)

        mode = "overlapping (2 viewers/candidate)" if args.overlap else "disjoint"
        logger.info("Done. Split %d candidates into %d %s sets.",
                    len(df), args.nsets, mode)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
