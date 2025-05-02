#!/usr/bin/env python3
"""
A script to process candidate CSV files, merge related data,
and package outputs into tarballs. Optionally, you can split the output
into multiple tarballs, each containing a specified number of unique pointings.
If no splitting is requested (npointings=0), a separate tarball for candidates
with alpha < 1.0 will also be created.

Usage:
    python pipeline.py -i input.csv -o output.tar.gz -d /path/to/output/ -m /path/to/metafiles/ [--threshold 0.1] [--npointings N] [--verbose]

Examples:
    - Create one tarball with all candidates and an extra tarball for alpha < 1.0:
        python pipeline.py -i input.csv -o output.tar.gz
    - Create tarballs with 2 unique pointings each:
        python pipeline.py -i input.csv -o output.tar.gz --npointings 2

Author: Vishnu
Date: 03-03-2025

Updated: 03-03-2025
By: Fazal
"""

import argparse
import logging
import os
import sys
import tarfile
import multiprocessing
from datetime import datetime

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
import pygedm


def setup_logging(verbose: bool) -> logging.Logger:
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=log_level,
    )
    return logging.getLogger(__name__)


class TarballCreator:
    @staticmethod
    def create_tarball(
        output_csv: str,
        meta_files: list,
        png_files: list,
        tarball_name: str,
        additional_files: list,
        logger: logging.Logger = None,
    ) -> None:
        # If a logger is provided (i.e. in main process), log events; otherwise, skip logging.
        if logger:
            logger.info("Packaging files into tarball: %s", tarball_name)
        try:
            with tarfile.open(tarball_name, "w:gz", dereference=True) as tar:
                # Add the main CSV file.
                tar.add(output_csv, arcname="candidates.csv")
                if logger:
                    logger.debug("Added main CSV: %s", output_csv)
                # Add metafiles into the 'metafiles' folder.
                for metafile in meta_files:
                    tar.add(
                        metafile,
                        arcname=os.path.join("metafiles", os.path.basename(metafile)),
                    )
                    if logger:
                        logger.debug("Added metafile: %s", metafile)
                # Add PNG images into the 'plots' folder.
                for png in png_files:
                    # check if the file exists
                    if os.path.isfile(png):
                        tar.add(
                            png, arcname=os.path.join("plots", os.path.basename(png))
                        )
                        if logger:
                            logger.debug("Added PNG: %s", png)
                    else:
                        if logger:
                            # Add a warning if the file does not exist, but continue
                            logger.warning("PNG file does not exist: %s", png)

                # Add any extra files at the top level.
                for file in additional_files:
                    tar.add(file, arcname=os.path.basename(file))
                    if logger:
                        logger.debug("Added extra file: %s", file)
            if logger:
                logger.info("Tarball '%s' created successfully.", tarball_name)
        except Exception as err:
            if logger:
                logger.error("Failed to create tarball: %s", err)
            sys.exit(1)


class CandidateProcessor:
    """
    Processes the candidate CSV file by cleaning data, adding galactic information,
    merging candidate-related files, and saving outputs.
    """

    def __init__(
        self,
        input_file: str,
        output_tarball: str,
        output_tarball_path: str,
        metafile_source_path: str,
        threshold: float,
        logger: logging.Logger,
    ):
        self.input_file = input_file
        self.output_tarball = output_tarball
        self.output_tarball_path = output_tarball_path
        self.metafile_source_path = metafile_source_path
        self.threshold = threshold
        self.logger = logger

        # Intermediate CSV filenames.
        self.candidates_csv = "candidates.csv"
        self.alpha_csv = "candidates_alpha_below_one.csv"
        self.pics_csv = "candidates_pics_above_threshold.csv"

    def load_input_dataframe(self) -> pd.DataFrame:
        self.logger.info("Reading input file: %s", self.input_file)
        try:
            df = pd.read_csv(self.input_file)
            self.logger.debug("Loaded %d rows.", len(df))
            return df
        except Exception as err:
            self.logger.error("Error reading input file: %s", err)
            sys.exit(1)

    def preprocess_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        self.logger.info("Preprocessing candidate data.")
        df["beam_id"] = df["beam_id"].apply(lambda x: str(x).zfill(2))

        # Combine segment and segment_id into a two-digit segment_id
        df["segment_id"] = df["segments"].astype(str) + df["segment_id"].astype(str)

        # Convert the combined column to integer for numerical consistency
        df["segment_id"] = df["segment_id"].astype(int)

        # Drop the original segment column as it's now merged
        df = df.drop(columns=["segments"])

        df = df.sort_values(by=["beam_id", "segment_id"]).reset_index(drop=True)
        pointing_ids = {
            pointing: idx for idx, pointing in enumerate(df["pointing"].unique())
        }
        df["pointing_id"] = df["pointing"].map(pointing_ids)
        # df["beam_id"] = df["beam"]

        # df["utc_start"] = df["pointing"].apply(
        #     lambda x: datetime.strptime(x, "%Y-%m-%d-%H:%M:%S").strftime("%Y-%m-%dT%H:%M:%S")
        # )

        df["metafile_path"] = df["utc_start"].apply(lambda x: f"metafiles/{x}.meta")
        self.logger.debug("Preprocessing done. Rows: %d", len(df))
        return df

    def add_galactic_info(self, df: pd.DataFrame) -> pd.DataFrame:
        self.logger.info("Adding galactic coordinate info.")
        unique = df[["utc_start", "beam", "ra", "dec"]].drop_duplicates().copy()
        coords = SkyCoord(ra=unique["ra"], dec=unique["dec"], unit=(u.hourangle, u.deg))
        unique["gl"] = coords.galactic.l.deg
        unique["gb"] = coords.galactic.b.deg
        unique["mjd_start"] = Time(
            unique["utc_start"].tolist(), format="isot", scale="utc"
        ).mjd

        max_dm = []
        for idx, row in unique.iterrows():
            try:
                dm, _ = pygedm.dist_to_dm(row["gl"], row["gb"], 50000, method="ymw16")
                max_dm.append(dm.value)
            except Exception as err:
                self.logger.error("DM calc failed at index %d: %s", idx, err)
                max_dm.append(float("nan"))
        unique["maxdm_ymw16"] = max_dm

        df = df.merge(
            unique[["utc_start", "beam", "gl", "gb", "mjd_start", "maxdm_ymw16"]],
            on=["utc_start", "beam"],
            how="left",
        )
        self.logger.debug("Galactic info added.")
        return df

    def merge_candidate_files(self, df: pd.DataFrame) -> pd.DataFrame:
        self.logger.info("Merging candidate files.")
        merged_results = []
        total_rows = 0

        for _, row in df.iterrows():
            alpha_beta_file = row["alpha_beta_file"]
            pics_file = row["pics_file"]
            try:
                df_alpha = pd.read_csv(alpha_beta_file)
                df_pics = pd.read_csv(pics_file)
                merged = df_alpha.merge(
                    df_pics,
                    left_on="fold_cands_filename",
                    right_on="filename",
                    how="inner",
                    suffixes=("", "_drop"),
                )
                merged = merged.loc[:, ~merged.columns.str.endswith("_drop")]
                if "filename" in merged.columns:
                    merged.drop(columns=["filename"], inplace=True)

                merged["pointing_id"] = row["pointing_id"]
                merged["beam_id"] = row["beam_id"]
                merged["beam_name"] = row["beam"]
                merged["utc_start"] = row["utc_start"]
                merged["source_name"] = row["target"]
                merged["ra"] = row["ra"]
                merged["dec"] = row["dec"]
                merged["gl"] = row["gl"]
                merged["gb"] = row["gb"]
                merged["mjd_start"] = row["mjd_start"]
                merged["maxdm_ymw16"] = row["maxdm_ymw16"]
                merged["metafile_path"] = row["metafile_path"]
                merged["segment_id"] = row["segment_id"]

                total_rows += merged.shape[0]
                merged_results.extend(merged.values.tolist())
            except Exception as err:
                self.logger.error(
                    "Error merging '%s' and '%s': %s", alpha_beta_file, pics_file, err
                )

        if not merged_results:
            self.logger.error("No candidate data merged. Exiting.")
            sys.exit(1)

        final_df = pd.DataFrame(merged_results, columns=merged.columns)
        self.logger.info(
            "Merged candidate files: expected %d rows, got %d rows.",
            total_rows,
            final_df.shape[0],
        )
        return final_df

    def finalize_dataframe(self, df: pd.DataFrame) -> (pd.DataFrame, list, list):
        self.logger.info("Finalizing candidate DataFrame.")
        rename_map = {
            "clfl2_trapum_Ter5": "pics_trapum_ter5",
            "PALFA_MeerKAT_L_SBAND_Best_Fscore": "pics_palfa_meerkat_l_sband_best_fscore",
            "clfl2_PALFA": "pics_palfa",
            "MeerKAT_L_SBAND_COMBINED_Best_Recall": "pics_meerkat_l_sband_combined_best_recall",
            "dm_old": "dm_user",
            "dm_new": "dm_opt",
            "dm_err": "dm_opt_err",
            "f0_err": "f0_opt_err",
            "f1_err": "f1_opt_err",
            "acc_err": "acc_opt_err",
            "f0_old": "f0_user",
            "f0_new": "f0_opt",
            "f1_old": "f1_user",
            "f1_new": "f1_opt",
            "acc_old": "acc_user",
            "acc_new": "acc_opt",
            "S/N": "sn_fft",
            "S/N_new": "sn_fold",
            "#id": "id",
            "filterbank_file": "filterbank_path",
        }
        df.rename(columns=rename_map, inplace=True)

        # Generate PNG file paths.
        df["png_path"] = df.apply(
            lambda row: os.path.join(
                "plots", os.path.splitext(row["fold_cands_filename"])[0] + ".png"
            ),
            axis=1,
        )
        png_files = (
            df.apply(
                lambda row: os.path.join(
                    row["fold_cands_filepath"],
                    os.path.splitext(row["fold_cands_filename"])[0] + ".png",
                ),
                axis=1,
            )
            .unique()
            .tolist()
        )

        df["candidate_tarball_path"] = os.path.join(
            self.output_tarball_path, self.output_tarball
        )

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
            "pics_trapum_ter5",
            "pics_palfa",
            "pics_palfa_meerkat_l_sband_best_fscore",
            "pics_meerkat_l_sband_combined_best_recall",
            "png_path",
            "metafile_path",
            "filterbank_path",
            "candidate_tarball_path",
        ]
        extra_cols = [col for col in df.columns if col not in required_cols]
        final_df = df[required_cols + extra_cols]

        unique_starts = final_df["utc_start"].unique()
        meta_files = [
            os.path.join(self.metafile_source_path, f"{start}.meta")
            for start in unique_starts
        ]

        self.logger.info("Final DataFrame ready with %d rows.", final_df.shape[0])
        return final_df, png_files, meta_files

    def filter_and_save(self, final_df: pd.DataFrame) -> None:
        final_df.to_csv(self.candidates_csv, index=False)
        self.logger.info("Saved candidates to '%s'.", self.candidates_csv)

        alpha_df = final_df.loc[final_df["alpha"] < 1.0]
        alpha_df.to_csv(self.alpha_csv, index=False)
        self.logger.info(
            "Saved %d alpha candidates to '%s'.", alpha_df.shape[0], self.alpha_csv
        )

        condition = (
            (final_df["pics_trapum_ter5"] >= self.threshold)
            | (final_df["pics_palfa"] >= self.threshold)
            | (final_df["pics_palfa_meerkat_l_sband_best_fscore"] >= self.threshold)
            | (final_df["pics_meerkat_l_sband_combined_best_recall"] >= self.threshold)
        )
        pics_df = final_df.loc[condition]
        pics_df.to_csv(self.pics_csv, index=False)
        self.logger.info(
            "Saved %d pics candidates to '%s'.", pics_df.shape[0], self.pics_csv
        )

    def process(self) -> (list, list):
        df = self.load_input_dataframe()
        df = self.preprocess_dataframe(df)
        df = self.add_galactic_info(df)
        merged_df = self.merge_candidate_files(df)
        final_df, png_files, meta_files = self.finalize_dataframe(merged_df)
        self.filter_and_save(final_df)
        return png_files, meta_files


def chunk_list(lst: list, n: int) -> list:
    return [lst[i : i + n] for i in range(0, len(lst), n)]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process candidate CSV files and package outputs into tarball(s)."
    )
    parser.add_argument(
        "-i", "--input_file", help="Path to the input CSV file", required=True
    )
    parser.add_argument(
        "-o", "--output_file", help="Base name for the output tarball(s)", required=True
    )
    parser.add_argument(
        "-d",
        "--output_path",
        help="Output directory for tarball(s)",
        default="/hercules/scratch/fkareem/CANDIDATE_TARBALLS/",
    )
    parser.add_argument(
        "-m",
        "--metafile_source_path",
        help="Directory for metafiles",
        default="/hercules/scratch/fkareem/elden-ring/include/metafiles/",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Threshold for pics score filtering (default: 0.1)",
    )
    parser.add_argument(
        "--npointings",
        type=int,
        default=0,
        help=(
            "Number of unique pointings per tarball. "
            "If 0, all candidates are packaged in one tarball "
            "and a separate tarball for alpha < 1.0 is created."
        ),
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def process_pointing_group(args_tuple):
    # Multiprocessing worker (logging omitted to avoid inter-process issues).
    (
        group,
        idx,
        candidate_csv,
        output_prefix,
        metafile_source_path,
        threshold,
        output_path,
    ) = args_tuple

    # Each worker reads the candidate CSV from disk and filters for its group.
    subset_df = pd.read_csv(candidate_csv)
    subset_df = subset_df[subset_df["utc_start"].isin(group)]

    # Define subset CSV filenames.
    subset_candidates_csv = f"{output_prefix}_set_{idx}_candidates.csv"
    subset_alpha_csv = f"{output_prefix}_set_{idx}_alpha.csv"
    subset_pics_csv = f"{output_prefix}_set_{idx}_pics.csv"

    subset_df.to_csv(subset_candidates_csv, index=False)
    alpha_df = subset_df[subset_df["alpha"] < 1.0]
    alpha_df.to_csv(subset_alpha_csv, index=False)

    condition = (
        (subset_df["pics_trapum_ter5"] >= threshold)
        | (subset_df["pics_palfa"] >= threshold)
        | (subset_df["pics_palfa_meerkat_l_sband_best_fscore"] >= threshold)
        | (subset_df["pics_meerkat_l_sband_combined_best_recall"] >= threshold)
    )
    pics_df = subset_df.loc[condition]
    pics_df.to_csv(subset_pics_csv, index=False)

    # Prepare metafile and PNG file lists.
    subset_pointings = subset_df["utc_start"].unique()
    subset_meta_files = [
        os.path.join(metafile_source_path, f"{utc}.meta") for utc in subset_pointings
    ]
    subset_png_files = (
        subset_df.apply(
            lambda row: os.path.join(
                row["fold_cands_filepath"],
                os.path.splitext(row["fold_cands_filename"])[0] + ".png",
            ),
            axis=1,
        )
        .unique()
        .tolist()
    )

    # Create main tarball.
    tarball_name = f"{output_prefix}_set_{idx}.tar.gz"
    # additional_files = [subset_alpha_csv, subset_pics_csv]
    # Comment this out if you want the tarball with all the candidates
    # TarballCreator.create_tarball(
    #     output_csv=subset_candidates_csv,
    #     meta_files=subset_meta_files,
    #     png_files=subset_png_files,
    #     tarball_name=tarball_name,
    #     additional_files=additional_files,
    #     logger=None
    # )

    # Create separate alpha tarball.
    alpha_png_files = (
        alpha_df.apply(
            lambda row: os.path.join(
                row["fold_cands_filepath"],
                os.path.splitext(row["fold_cands_filename"])[0] + ".png",
            ),
            axis=1,
        )
        .unique()
        .tolist()
    )
    alpha_pointings = alpha_df["utc_start"].unique()
    alpha_meta_files = [
        os.path.join(metafile_source_path, f"{utc}.meta") for utc in alpha_pointings
    ]
    alpha_output_file = os.path.join(
        output_path, tarball_name.replace(".tar.gz", "_alpha_below_one.tar.gz")
    )
    print(f"Creating tarball: {alpha_output_file}")
    TarballCreator.create_tarball(
        output_csv=subset_alpha_csv,
        meta_files=alpha_meta_files,
        png_files=alpha_png_files,
        tarball_name=alpha_output_file,
        additional_files=[subset_pics_csv],
        logger=None,
    )

    # Create separate pics tarball.
    pics_png_files = (
        pics_df.apply(
            lambda row: os.path.join(
                row["fold_cands_filepath"],
                os.path.splitext(row["fold_cands_filename"])[0] + ".png",
            ),
            axis=1,
        )
        .unique()
        .tolist()
    )
    pics_pointings = pics_df["utc_start"].unique()
    pics_meta_files = [
        os.path.join(metafile_source_path, f"{utc}.meta") for utc in pics_pointings
    ]
    pics_output_file = os.path.join(
        output_path,
        tarball_name.replace(".tar.gz", f"_pics_above_threshold_{threshold}.tar.gz"),
    )
    print(f"Creating tarball: {pics_output_file}")
    TarballCreator.create_tarball(
        output_csv=subset_pics_csv,
        meta_files=pics_meta_files,
        png_files=pics_png_files,
        tarball_name=pics_output_file,
        additional_files=[subset_alpha_csv],
        logger=None,
    )


def create_single_tarball(processor, meta_files, png_files, args, logger):

    # Comment this out if you want the tarball with all the candidates
    # main_tarball = os.path.join(args.output_path, args.output_file)
    # additional_files = [processor.alpha_csv, processor.pics_csv]
    # TarballCreator.create_tarball(
    #     output_csv=processor.candidates_csv,
    #     meta_files=meta_files,
    #     png_files=png_files,
    #     tarball_name=main_tarball,
    #     additional_files=additional_files,
    #     logger=logger
    # )
    # logger.info("Created single tarball with all candidates: %s", args.output_file)

    # Create separate alpha tarball.
    alpha_df = pd.read_csv(processor.alpha_csv)
    alpha_png_files = (
        alpha_df.apply(
            lambda row: os.path.join(
                row["fold_cands_filepath"],
                os.path.splitext(row["fold_cands_filename"])[0] + ".png",
            ),
            axis=1,
        )
        .unique()
        .tolist()
    )
    alpha_pointings = alpha_df["utc_start"].unique()
    alpha_meta_files = [
        os.path.join(args.metafile_source_path, f"{utc}.meta")
        for utc in alpha_pointings
    ]
    alpha_output_file = os.path.join(
        args.output_path, args.output_file.replace(".tar.gz", "_alpha_below_one.tar.gz")
    )
    TarballCreator.create_tarball(
        output_csv=processor.alpha_csv,
        meta_files=alpha_meta_files,
        png_files=alpha_png_files,
        tarball_name=alpha_output_file,
        additional_files=[],
        logger=logger,
    )
    logger.info("Created separate alpha tarball: %s", alpha_output_file)

    # Create separate pics tarball.
    pics_df = pd.read_csv(processor.pics_csv)
    pics_png_files = (
        pics_df.apply(
            lambda row: os.path.join(
                row["fold_cands_filepath"],
                os.path.splitext(row["fold_cands_filename"])[0] + ".png",
            ),
            axis=1,
        )
        .unique()
        .tolist()
    )
    pics_pointings = pics_df["utc_start"].unique()
    pics_meta_files = [
        os.path.join(args.metafile_source_path, f"{utc}.meta") for utc in pics_pointings
    ]
    pics_output_file = os.path.join(
        args.output_path,
        args.output_file.replace(".tar.gz", "_pics_above_threshold.tar.gz"),
    )
    TarballCreator.create_tarball(
        output_csv=processor.pics_csv,
        meta_files=pics_meta_files,
        png_files=pics_png_files,
        tarball_name=pics_output_file,
        additional_files=[],
        logger=logger,
    )
    logger.info("Created separate pics tarball: %s", pics_output_file)


def main():
    args = parse_arguments()
    logger = setup_logging(args.verbose)
    logger.info("Starting candidate processing pipeline.")

    processor = CandidateProcessor(
        input_file=args.input_file,
        output_tarball=args.output_file,
        output_tarball_path=args.output_path,
        metafile_source_path=args.metafile_source_path,
        threshold=args.threshold,
        logger=logger,
    )

    # Process the candidate CSV and generate supporting files.
    png_files, meta_files = processor.process()

    if args.npointings > 0:
        logger.info(
            "Splitting output into multiple tarballs with %d unique pointings each.",
            args.npointings,
        )
        full_df = pd.read_csv(processor.candidates_csv)
        unique_pointings = sorted(full_df["utc_start"].unique())
        pointing_groups = chunk_list(unique_pointings, args.npointings)
        output_prefix = args.output_file.replace(".tar.gz", "")
        # Prepare arguments for each worker without sending the full DataFrame.
        args_list = []
        for idx, group in enumerate(pointing_groups, start=1):
            args_list.append(
                (
                    group,
                    idx,
                    processor.candidates_csv,
                    output_prefix,
                    args.metafile_source_path,
                    args.threshold,
                    args.output_path,
                )
            )
        with multiprocessing.Pool(
            processes=min(len(args_list), multiprocessing.cpu_count())
        ) as pool:
            pool.map(process_pointing_group, args_list)
    else:
        create_single_tarball(processor, meta_files, png_files, args, logger)

    logger.info("Candidate processing pipeline completed successfully.")


if __name__ == "__main__":
    main()
