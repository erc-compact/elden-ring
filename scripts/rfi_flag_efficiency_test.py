#!/usr/bin/env python
import pandas as pd
import argparse

def find_known_pulsar(candidate_period, candidate_dm, df_known, period_tol, dm_tol):
    """
    Check each known pulsar. If both the period and DM differences
    are within the specified tolerances, return the pulsar's name.
    Otherwise, return None.
    """
    for _, row in df_known.iterrows():
        if abs(candidate_period - row['Period']) <= period_tol and abs(candidate_dm - row['DM']) <= dm_tol:
            return row['Pulsar']
    return None

def main():
    parser = argparse.ArgumentParser(
        description="Compare candidate pulsars with known pulsars and generate an annotated CSV and report."
    )
    parser.add_argument("candidates_csv", help="CSV file containing candidates (e.g., candidates.csv).")
    parser.add_argument("known_csv", help="CSV file containing known pulsars (e.g., known_pulsars.csv).")
    parser.add_argument("--rfi_flag", type=str, help='RFI Flag string')
    parser.add_argument("--output", default="annotated_candidates.csv", help="Output CSV file name.")
    parser.add_argument("--report", default="report.txt", help="Output report file name.")
    parser.add_argument("--period_tol", type=float, default=0.1, help="Tolerance for period matching.")
    parser.add_argument("--dm_tol", type=float, default=0.5, help="Tolerance for DM matching.")
    args = parser.parse_args()

    # Load the CSV files
    df_candidates = pd.read_csv(args.candidates_csv)
    df_known = pd.read_csv(args.known_csv)

    # Prepare lists to store the new columns
    match_tags = []
    pulsar_names = []

    # For each candidate, try to match with a known pulsar
    for _, row in df_candidates.iterrows():
        candidate_period = row['period']
        candidate_dm = row['dm']
        match = find_known_pulsar(candidate_period, candidate_dm, df_known, args.period_tol, args.dm_tol)
        if match:
            match_tags.append("KNOWN_PSR")
            pulsar_names.append(match)
        else:
            match_tags.append("UNKNOWN")
            pulsar_names.append("")
    
    # Add the new columns to the candidate DataFrame
    df_candidates['match'] = match_tags
    df_candidates['pulsar_name'] = pulsar_names

    # Write the annotated CSV file
    df_candidates.to_csv(args.output, index=False)
    print(f"Annotated candidate CSV written to {args.output}")

    # Prepare the report information
    known_df = df_candidates[df_candidates['match'] == "KNOWN_PSR"]
    unknown_df = df_candidates[df_candidates['match'] == "UNKNOWN"]
    known_counts = known_df['pulsar_name'].value_counts()
    num_known_distinct = known_counts.size
    num_known_total = known_df.shape[0]
    num_unknown = unknown_df.shape[0]
    total = df_candidates.shape[0]

    # Create the report lines
    report_lines = []
    report_lines.append("Pulsar Candidate Matching Report")
    report_lines.append("--------------------------------")
    report_lines.append(f"Total candidates: {total}")
    report_lines.append(f"KNOWN_PSR candidates: {num_known_total}")
    report_lines.append(f"Distinct known pulsars detected: {num_known_distinct}")
    report_lines.append("")
    report_lines.append("Detections per known pulsar:")
    for pulsar, count in known_counts.items():
        report_lines.append(f"  {pulsar}: {count}")
    report_lines.append("")
    report_lines.append(f"UNKNOWN candidates: {num_unknown}")

    # Write the report file
    with open(args.report, "w") as report_file:
        report_file.write("\n".join(report_lines))
    print(f"Report written to {args.report}")

if __name__ == "__main__":
    main()
