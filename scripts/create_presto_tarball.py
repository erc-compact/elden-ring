#!/usr/bin/env python3
"""
Create CandyJar-compatible tarballs from PRESTO pipeline results.

This script creates tarballs that match the exact structure expected by CandyJar,
including the same CSV fields and folder structure as create_candyjar_tarball.py.

Tarball Structure:
    {tarball_name}/
    ├── candidates.csv
    ├── candidates_pics_above_threshold.csv
    ├── metafiles/
    │   └── {utc_start}.meta
    └── plots/
        └── *.png

Note: PRESTO doesn't produce alpha/beta/gamma values (those come from pulsarx fold
optimization), so candidates_alpha_below_one.csv is not created for PRESTO.

Author: ELDEN-RING Pipeline
"""

import argparse
import logging
import os
import sys
import tarfile
import shutil
from datetime import datetime

import pandas as pd
import numpy as np

# Optional imports for galactic coordinate calculation
try:
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.time import Time
    import pygedm
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False


def setup_logging(verbose: bool) -> logging.Logger:
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=log_level,
    )
    return logging.getLogger(__name__)


def add_galactic_info(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Add galactic coordinates, MJD, and max DM from YMW16 model.
    Matches the processing done in create_candyjar_tarball.py.
    """
    if not HAS_ASTROPY:
        logger.warning("Astropy/pygedm not available. Skipping galactic coordinate calculation.")
        df["gl"] = 0.0
        df["gb"] = 0.0
        df["mjd_start"] = 0.0
        df["maxdm_ymw16"] = 0.0
        df["dist_ymw16"] = 0.0
        return df

    logger.info("Adding galactic coordinate info.")

    # Get unique RA/Dec combinations
    if "ra" in df.columns and "dec" in df.columns:
        unique = df[["ra", "dec"]].drop_duplicates().copy()

        # Handle different RA/Dec formats
        try:
            # Try hourangle format first (HH:MM:SS)
            coords = SkyCoord(ra=unique["ra"], dec=unique["dec"], unit=(u.hourangle, u.deg))
        except Exception:
            try:
                # Try degree format
                coords = SkyCoord(ra=unique["ra"], dec=unique["dec"], unit=(u.deg, u.deg))
            except Exception as e:
                logger.warning(f"Could not parse coordinates: {e}")
                df["gl"] = 0.0
                df["gb"] = 0.0
                df["maxdm_ymw16"] = 0.0
                df["dist_ymw16"] = 0.0
                return df

        unique["gl"] = coords.galactic.l.deg
        unique["gb"] = coords.galactic.b.deg

        # Calculate max DM using YMW16 model
        max_dm = []
        for idx, row in unique.iterrows():
            try:
                dm, _ = pygedm.dist_to_dm(row["gl"], row["gb"], 50000, method="ymw16")
                max_dm.append(dm.value)
            except Exception as err:
                logger.debug(f"DM calc failed at index {idx}: {err}")
                max_dm.append(float("nan"))
        unique["maxdm_ymw16"] = max_dm

        # Merge back
        df = df.merge(unique[["ra", "dec", "gl", "gb", "maxdm_ymw16"]], on=["ra", "dec"], how="left")
    else:
        df["gl"] = 0.0
        df["gb"] = 0.0
        df["maxdm_ymw16"] = 0.0

    # Add MJD start if utc_start is available
    if "utc_start" in df.columns:
        try:
            # Try different datetime formats
            for fmt in ["%Y-%m-%dT%H:%M:%S", "%Y-%m-%d-%H:%M:%S", "%Y/%m/%d %H:%M:%S"]:
                try:
                    df["mjd_start"] = df["utc_start"].apply(
                        lambda x: Time(datetime.strptime(str(x), fmt)).mjd if pd.notna(x) else 0.0
                    )
                    break
                except Exception:
                    continue
            else:
                df["mjd_start"] = 0.0
        except Exception as e:
            logger.warning(f"Could not calculate MJD: {e}")
            df["mjd_start"] = 0.0
    else:
        df["mjd_start"] = 0.0

    # Add dist_ymw16 placeholder (would need DM to calculate properly)
    df["dist_ymw16"] = 0.0

    logger.debug("Galactic info added.")
    return df


def normalize_presto_columns(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Normalize PRESTO output columns to match CandyJar expected format.
    """
    logger.info("Normalizing PRESTO columns to CandyJar format.")

    # Column mapping from PRESTO output to CandyJar format
    rename_map = {
        # Period/frequency columns
        "period": "period_s",
        "p0": "period_s",
        "frequency": "f0_user",
        "freq": "f0_user",
        "f0": "f0_user",

        # DM columns
        "dm": "dm_user",
        "DM": "dm_user",

        # Acceleration columns
        "acc": "acc_user",
        "accel": "acc_user",
        "acceleration": "acc_user",

        # S/N columns
        "snr": "sn_fft",
        "sigma": "sn_fft",
        "S/N": "sn_fft",
        "sn": "sn_fft",

        # Candidate ID
        "cand_id": "id",
        "cand_num": "id",
        "candidate_id": "id",
        "#id": "id",

        # File paths
        "pfd_file": "archive_path",
        "png_file": "png_path",
        "fil_file": "filterbank_path",
        "filterbank_file": "filterbank_path",
    }

    # Apply renames (only for columns that exist)
    for old_name, new_name in rename_map.items():
        if old_name in df.columns and new_name not in df.columns:
            df.rename(columns={old_name: new_name}, inplace=True)
            logger.debug(f"Renamed column: {old_name} -> {new_name}")

    # Calculate frequency from period if needed
    if "f0_user" not in df.columns and "period_s" in df.columns:
        df["f0_user"] = 1.0 / df["period_s"].replace(0, np.nan)
        logger.debug("Calculated f0_user from period_s")

    # Add missing columns with default values
    default_columns = {
        "pointing_id": 0,
        "beam_id": "00",
        "beam_name": "",
        "source_name": "",
        "segment_id": 0,
        "ra": "",
        "dec": "",
        "gl": 0.0,
        "gb": 0.0,
        "mjd_start": 0.0,
        "utc_start": "",
        "f0_user": 0.0,
        "f0_opt": 0.0,
        "f0_opt_err": 0.0,
        "f1_user": 0.0,
        "f1_opt": 0.0,
        "f1_opt_err": 0.0,
        "acc_user": 0.0,
        "acc_opt": 0.0,
        "acc_opt_err": 0.0,
        "dm_user": 0.0,
        "dm_opt": 0.0,
        "dm_opt_err": 0.0,
        "sn_fft": 0.0,
        "sn_fold": 0.0,
        "maxdm_ymw16": 0.0,
        "dist_ymw16": 0.0,
        "cdm": 0.0,
        "png_path": "",
        "metafile_path": "",
        "filterbank_path": "",
        "candidate_tarball_path": "",
    }

    for col, default_val in default_columns.items():
        if col not in df.columns:
            df[col] = default_val
            logger.debug(f"Added missing column: {col} = {default_val}")

    # Copy f0_user to f0_opt if not set (PRESTO doesn't optimize)
    if df["f0_opt"].eq(0).all() and "f0_user" in df.columns:
        df["f0_opt"] = df["f0_user"]

    # Copy dm_user to dm_opt if not set
    if df["dm_opt"].eq(0).all() and "dm_user" in df.columns:
        df["dm_opt"] = df["dm_user"]

    # Copy acc_user to acc_opt if not set
    if df["acc_opt"].eq(0).all() and "acc_user" in df.columns:
        df["acc_opt"] = df["acc_user"]

    # Use sn_fft as sn_fold if sn_fold not available
    if df["sn_fold"].eq(0).all() and "sn_fft" in df.columns:
        df["sn_fold"] = df["sn_fft"]

    return df


def add_pics_columns(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Ensure PICS score columns exist. They should be populated from the PICS classifier.
    """
    pics_columns = [
        "pics_trapum_ter5",
        "pics_palfa",
        "pics_palfa_meerkat_l_sband_best_fscore",
        "pics_meerkat_l_sband_combined_best_recall",
    ]

    # Map from raw PICS output names to CandyJar names
    pics_rename = {
        "clfl2_trapum_Ter5": "pics_trapum_ter5",
        "PALFA_MeerKAT_L_SBAND_Best_Fscore": "pics_palfa_meerkat_l_sband_best_fscore",
        "clfl2_PALFA": "pics_palfa",
        "MeerKAT_L_SBAND_COMBINED_Best_Recall": "pics_meerkat_l_sband_combined_best_recall",
    }

    # Apply renames
    for old_name, new_name in pics_rename.items():
        if old_name in df.columns and new_name not in df.columns:
            df.rename(columns={old_name: new_name}, inplace=True)
            logger.debug(f"Renamed PICS column: {old_name} -> {new_name}")

    # Add missing PICS columns with default 0.0
    for col in pics_columns:
        if col not in df.columns:
            df[col] = 0.0
            logger.debug(f"Added missing PICS column: {col}")

    return df


def create_tarball(
    candidates_df: pd.DataFrame,
    png_dir: str,
    metafile_dir: str,
    output_tarball: str,
    tarball_prefix: str,
    pics_threshold: float,
    snr_threshold: float,
    logger: logging.Logger,
) -> str:
    """
    Create CandyJar-compatible tarball from PRESTO results.
    """
    logger.info(f"Creating tarball: {output_tarball}")

    # Create working directory
    work_dir = os.path.join(os.path.dirname(output_tarball) or ".", "tarball_work")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)

    # Determine metafile subdir (prefer existing metafile_path prefix)
    meta_subdir = "metafiles"
    if "metafile_path" in candidates_df.columns:
        try:
            if candidates_df["metafile_path"].astype(str).str.startswith("meta/").any():
                meta_subdir = "meta"
        except Exception:
            pass

    # Create subdirectories
    plots_dir = os.path.join(work_dir, "plots")
    meta_dir = os.path.join(work_dir, meta_subdir)
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    # Set paths in dataframe
    candidates_df["candidate_tarball_path"] = output_tarball

    # Update png_path to be relative to tarball
    candidates_df["png_path"] = candidates_df.apply(
        lambda row: os.path.join("plots", os.path.basename(str(row.get("png_path", ""))) or f"cand_{row.name}.png"),
        axis=1,
    )

    # Update metafile_path only if missing/blank
    if "utc_start" in candidates_df.columns:
        if "metafile_path" not in candidates_df.columns:
            candidates_df["metafile_path"] = ""
        candidates_df["metafile_path"] = candidates_df.apply(
            lambda row: row.get("metafile_path")
            if str(row.get("metafile_path", "")).strip()
            else (f"{meta_subdir}/{row['utc_start']}.meta" if pd.notna(row.get("utc_start")) and str(row.get("utc_start")) else f"{meta_subdir}/unknown.meta"),
            axis=1,
        )

    # Define required columns in order (matching CandyJar)
    required_cols = [
        "pointing_id",
        "beam_id",
        "beam_name",
        "source_name",
        "segment_id",
        "ra",
        "dec",
        "gl",
        "gb",
        "mjd_start",
        "utc_start",
        "f0_user",
        "f0_opt",
        "f0_opt_err",
        "f1_user",
        "f1_opt",
        "f1_opt_err",
        "acc_user",
        "acc_opt",
        "acc_opt_err",
        "dm_user",
        "dm_opt",
        "dm_opt_err",
        "sn_fft",
        "sn_fold",
        "maxdm_ymw16",
        "dist_ymw16",
        "cdm",
        "pics_trapum_ter5",
        "pics_palfa",
        "pics_palfa_meerkat_l_sband_best_fscore",
        "pics_meerkat_l_sband_combined_best_recall",
        "png_path",
        "metafile_path",
        "filterbank_path",
        "candidate_tarball_path",
    ]

    # Get extra columns not in required list
    extra_cols = [col for col in candidates_df.columns if col not in required_cols]

    # Reorder columns
    final_cols = [col for col in required_cols if col in candidates_df.columns] + extra_cols
    final_df = candidates_df[final_cols].copy()

    # Filter by SNR threshold
    if snr_threshold > 0:
        original_count = len(final_df)
        final_df = final_df[final_df["sn_fold"] >= snr_threshold]
        logger.info(f"SNR filter: {original_count} -> {len(final_df)} candidates (threshold: {snr_threshold})")

    # Save main candidates CSV
    candidates_csv = os.path.join(work_dir, "candidates.csv")
    final_df.to_csv(candidates_csv, index=False)
    logger.info(f"Saved {len(final_df)} candidates to candidates.csv")

    # Filter by PICS threshold and save
    pics_cols = [
        "pics_trapum_ter5",
        "pics_palfa",
        "pics_palfa_meerkat_l_sband_best_fscore",
        "pics_meerkat_l_sband_combined_best_recall",
    ]
    pics_condition = (
        (final_df["pics_trapum_ter5"] >= pics_threshold)
        | (final_df["pics_palfa"] >= pics_threshold)
        | (final_df["pics_palfa_meerkat_l_sband_best_fscore"] >= pics_threshold)
        | (final_df["pics_meerkat_l_sband_combined_best_recall"] >= pics_threshold)
    )
    pics_df = final_df[pics_condition]
    pics_csv = os.path.join(work_dir, "candidates_pics_above_threshold.csv")
    pics_df.to_csv(pics_csv, index=False)
    logger.info(f"Saved {len(pics_df)} PICS-filtered candidates (threshold: {pics_threshold})")

    # Emit CSVs to current working directory for workflow outputs
    final_df.to_csv(os.path.join(os.getcwd(), "candidates.csv"), index=False)
    pics_df.to_csv(os.path.join(os.getcwd(), "candidates_pics_above_threshold.csv"), index=False)

    # Copy PNG files
    png_count = 0
    if png_dir and os.path.isdir(png_dir):
        for png_file in os.listdir(png_dir):
            if png_file.endswith(".png"):
                src = os.path.join(png_dir, png_file)
                dst = os.path.join(plots_dir, png_file)
                shutil.copy2(src, dst)
                png_count += 1
    logger.info(f"Copied {png_count} PNG files")

    # Copy metafiles
    meta_count = 0
    if metafile_dir and os.path.isdir(metafile_dir):
        unique_utc = final_df["utc_start"].dropna().unique()
        for utc in unique_utc:
            meta_file = os.path.join(metafile_dir, f"{utc}.meta")
            if os.path.isfile(meta_file):
                shutil.copy2(meta_file, meta_dir)
                meta_count += 1
    logger.info(f"Copied {meta_count} metafiles")

    # Create tarball
    tarball_name = f"{tarball_prefix}_presto"
    with tarfile.open(output_tarball, "w:gz") as tar:
        tar.add(work_dir, arcname=tarball_name)
    logger.info(f"Created tarball: {output_tarball}")

    # Cleanup
    shutil.rmtree(work_dir)

    return output_tarball


def main():
    parser = argparse.ArgumentParser(
        description="Create CandyJar-compatible tarball from PRESTO results"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input candidates CSV file (from PRESTO pipeline)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output tarball path (e.g., output_presto.tar.gz)"
    )
    parser.add_argument(
        "--png-dir",
        default=".",
        help="Directory containing PNG files"
    )
    parser.add_argument(
        "--metafile-dir",
        default=None,
        help="Directory containing .meta files"
    )
    parser.add_argument(
        "--tarball-prefix",
        default="presto",
        help="Prefix for tarball internal directory name"
    )
    parser.add_argument(
        "--pics-threshold",
        type=float,
        default=0.1,
        help="PICS score threshold for filtering (default: 0.1)"
    )
    parser.add_argument(
        "--snr-threshold",
        type=float,
        default=6.0,
        help="SNR threshold for filtering (default: 6.0)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    logger.info("Starting PRESTO tarball creation")

    # Read input CSV
    logger.info(f"Reading input file: {args.input}")
    try:
        df = pd.read_csv(args.input)
        logger.info(f"Loaded {len(df)} candidates")
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        sys.exit(1)

    if len(df) == 0:
        logger.warning("No candidates found in input file")
        sys.exit(0)

    # Normalize columns to CandyJar format
    df = normalize_presto_columns(df, logger)

    # Add galactic info
    df = add_galactic_info(df, logger)

    # Add PICS columns
    df = add_pics_columns(df, logger)

    # Create tarball
    create_tarball(
        candidates_df=df,
        png_dir=args.png_dir,
        metafile_dir=args.metafile_dir,
        output_tarball=args.output,
        tarball_prefix=args.tarball_prefix,
        pics_threshold=args.pics_threshold,
        snr_threshold=args.snr_threshold,
        logger=logger,
    )

    logger.info("PRESTO tarball creation completed successfully")


if __name__ == "__main__":
    main()
