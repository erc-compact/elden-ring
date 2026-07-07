#!/usr/bin/env python3
"""
followup_candfile.py
====================

Candidate followup / confirmation helper for the ELDEN-RING ``run_followup`` workflow.

The idea (modelled on CaReFuL, https://github.com/erc-compact/CaReFuL) is to take the
T1/T2 candidates found in ONE search/observation and try to re-find them in a DIFFERENT
observation.  The candidate CSV supplies *what to look for* (DM, acceleration, F0); the
filterbanks to actually search come from the pipeline's ``inputfile.txt`` (a different
observation) and are handled by Nextflow, not this script.

Two modes:

``prepare``
    Read a CandyJar classification CSV, keep only the T1/T2 candidates and emit:
      * ``input.candfile``  -- pulsarx candfile ``#id DM acc F0 F1 F2 S/N`` listing every
        T1/T2 target (used both to derive the peasoup acceleration range and, later, for
        matching + folding).
      * ``targets.dm``      -- a peasoup ``--dm_file``.  This is the SPARSE UNION of each
        candidate's narrow band ``[dm_opt - dm_tol, dm_opt + dm_tol]`` (stepped by
        ``dm_step``).  It is deliberately NOT a filled continuous range between the
        smallest and largest DM -- well separated candidates leave gaps.
      * (MeerKAT only, ``--per-beam``) per-original-beam candfiles/dm-files so the
        neighbouring-beam logic can join each candidate to its neighbour beams.

``match``
    Given a peasoup ``overview.xml`` and the ``input.candfile`` from ``prepare``, keep the
    XML candidates that fall within ``+/- ptol`` (fractional period) AND
    ``+/- dm_tol`` (absolute DM) of an input candidate.  Optionally reject birdies and top
    up to ``--ncands`` by S/N.  Emit ``output.candfile`` (for folding) and ``matches.csv``
    (the confirmation record: which input candidate matched which XML candidate).
"""

import os
import sys
import argparse
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CandyJar CSV -> candidate targets
# ---------------------------------------------------------------------------
T1T2 = ("T1_CAND", "T2_CAND")


def read_t1t2_candidates(csv_path):
    """Return a DataFrame of T1/T2 candidates with normalised numeric columns."""
    df = pd.read_csv(csv_path)
    if "classification" not in df.columns:
        sys.exit("[followup] input CSV has no 'classification' column: %s" % csv_path)

    df = df[df["classification"].isin(T1T2)].copy()
    if df.empty:
        sys.exit("[followup] no T1_CAND/T2_CAND rows found in %s" % csv_path)

    # f0_opt is the optimised spin frequency (Hz); period = 1/f0
    df["f0"] = df["f0_opt"].astype(float)
    df["f1"] = df.get("f1_opt", 0.0)
    df["f1"] = pd.to_numeric(df["f1"], errors="coerce").fillna(0.0)
    df["dm"] = df["dm_opt"].astype(float)
    df["acc"] = df["acc_opt"].astype(float)
    df["snr"] = pd.to_numeric(df.get("sn_fold", 0.0), errors="coerce").fillna(0.0)
    df["period"] = 1.0 / df["f0"]
    # original beam of the candidate (needed for the MeerKAT neighbour-beam branch)
    if "beam_name" in df.columns:
        df["orig_beam"] = df["beam_name"].astype(str).str.strip()
    else:
        df["orig_beam"] = ""
    df = df.reset_index(drop=True)
    return df


def write_candfile(path, df):
    """Write a pulsarx candfile: ``#id DM acc F0 F1 F2 S/N``."""
    with open(path, "w") as f:
        f.write("#id DM acc F0 F1 F2 S/N\n")
        for i, row in enumerate(df.itertuples(index=False)):
            f.write("%d %f %f %f %g 0 %f\n" % (i, row.dm, row.acc, row.f0, row.f1, row.snr))
    return path


def build_dm_file(path, dm_values, dm_tol, dm_step, use_zero_dm=False):
    """Sparse union of each candidate's ``[dm - dm_tol, dm + dm_tol]`` band."""
    dm_set = set()
    for dm in dm_values:
        band = np.round(np.arange(dm - dm_tol, dm + dm_tol + dm_step, dm_step), 3)
        band = band[band >= 0.0]
        dm_set.update(band.tolist())
    if use_zero_dm:
        dm_set.add(0.0)
    dms = np.array(sorted(dm_set))
    np.savetxt(path, dms, fmt="%f")
    return path


# ---------------------------------------------------------------------------
# prepare mode
# ---------------------------------------------------------------------------
def assign_tiers(df, acc_floor, mid_cap):
    """Split candidates into acceleration tiers and give each an acc search span.

    Rationale: a few extreme-acceleration candidates would otherwise force the
    whole search to a huge acc range. Instead we search each group only over its
    own candidates' DMs, at an acc span sized to that group:

      * Tier 1 (bulk): |acc| <= cap, where cap = max(median|acc|, acc_floor).
                       search span = +/- cap.
      * Tier 2 (mid):  cap < |acc| <= mid_cap.   search span = +/- mid_cap.
      * Tier 3 (high): |acc| > mid_cap.          search span = +/- max(|acc|) of tier.

    Each tier searches ONLY its own candidates' DM bands (the sparse union of
    those bands), so the expensive high-acc span never runs over the full DM set.
    Returns a list of dicts: {name, df, acc_span}.
    """
    ab = df["acc"].abs()
    cap = max(float(np.median(ab)), float(acc_floor))
    print("[followup] median|acc|=%.1f  cap=max(median,%.0f)=%.1f  mid_cap=%.0f"
          % (np.median(ab), acc_floor, cap, mid_cap))

    tiers = []
    t1 = df[ab <= cap]
    if len(t1):
        tiers.append({"name": "t1_bulk", "df": t1, "acc_span": cap})
    t2 = df[(ab > cap) & (ab <= mid_cap)]
    if len(t2):
        tiers.append({"name": "t2_mid", "df": t2, "acc_span": float(mid_cap)})
    t3 = df[ab > mid_cap]
    if len(t3):
        tiers.append({"name": "t3_high", "df": t3, "acc_span": float(t3["acc"].abs().max())})

    for t in tiers:
        print("[followup]   %-8s: %2d cands, acc=+/-%.0f"
              % (t["name"], len(t["df"]), t["acc_span"]))
    return tiers


def mode_prepare(args):
    df = read_t1t2_candidates(args.input_csv)
    print("[followup] %d T1/T2 candidates read from %s" % (len(df), args.input_csv))

    # Full input.candfile is kept for reference/matching top-up; matching itself
    # uses the per-tier candfile so a candidate is only matched to its own tier.
    write_candfile(os.path.join(args.out_dir, "input.candfile"), df)

    tiers = assign_tiers(df, args.acc_floor, args.mid_cap)

    # Emit one candfile + dm file per tier, and a manifest Nextflow fans out over.
    manifest = os.path.join(args.out_dir, "tiers.csv")
    with open(manifest, "w") as mf:
        mf.write("tier,candfile,dmfile,acc_start,acc_end\n")
        for t in tiers:
            cf = "%s.candfile" % t["name"]
            dmf = "%s.dm" % t["name"]
            write_candfile(os.path.join(args.out_dir, cf), t["df"])
            build_dm_file(
                os.path.join(args.out_dir, dmf),
                t["df"]["dm"].values, args.dm_tol, args.dm_step, args.use_zero_dm,
            )
            mf.write("%s,%s,%s,%.6f,%.6f\n"
                     % (t["name"], cf, dmf, -t["acc_span"], t["acc_span"]))
    print("[followup] wrote %d tier(s) and manifest %s" % (len(tiers), manifest))


# ---------------------------------------------------------------------------
# match mode
# ---------------------------------------------------------------------------
IGNORED_XML_TAGS = {
    "candidate", "opt_period", "folded_snr", "byte_offset", "is_adjacent",
    "is_physical", "ddm_count_ratio", "ddm_snr_ratio",
}


def xml_to_df(xml_path):
    """Parse a peasoup overview.xml candidate list into a DataFrame."""
    root = ET.parse(xml_path).getroot()
    candidates = root[7]
    rows = []
    for cand in candidates:
        d = {}
        for entry in cand.iter():
            if entry.tag not in IGNORED_XML_TAGS:
                d[entry.tag] = entry.text
        d["cand_id_in_file"] = cand.attrib.get("id")
        rows.append(d)
    xml_df = pd.DataFrame(rows)
    xml_df = xml_df.astype({
        "snr": float, "dm": float, "period": float,
        "acc": float, "cand_id_in_file": int,
    })
    return xml_df


def reject_birdies(xml_df, birdies):
    if not birdies:
        return xml_df
    tol = 1e-6
    for b in (float(x) for x in birdies.split(",") if x.strip()):
        xml_df = xml_df[~xml_df["period"].between(b - tol, b + tol)]
    return xml_df


def shortlist(input_df, xml_df, ptol, dm_tol, ncands):
    """Keep XML cands within +/-ptol period AND +/-dm_tol DM of an input candidate."""
    matched_frames = []
    match_records = []
    for _, row in input_df.iterrows():
        f0 = float(row["F0"])
        period = 1.0 / f0
        dm = float(row["DM"])
        pmin, pmax = period - period * ptol, period + period * ptol
        dmin, dmax = dm - dm_tol, dm + dm_tol
        hits = xml_df[
            (xml_df["period"] >= pmin) & (xml_df["period"] <= pmax) &
            (xml_df["dm"] >= dmin) & (xml_df["dm"] <= dmax)
        ]
        matched_frames.append(hits)
        for h in hits.itertuples(index=False):
            match_records.append({
                "input_id": int(row["#id"]) if "#id" in row else -1,
                "input_dm": dm, "input_period": period, "input_f0": f0,
                "input_acc": float(row["acc"]),
                "match_cand_id": int(h.cand_id_in_file),
                "match_dm": float(h.dm), "match_period": float(h.period),
                "match_acc": float(h.acc), "match_snr": float(h.snr),
            })

    if matched_frames:
        matched = pd.concat(matched_frames, ignore_index=True)
        matched = matched.drop_duplicates(subset=["cand_id_in_file"])
    else:
        matched = xml_df.iloc[0:0].copy()

    print("[followup] %d XML candidates matched input targets" % len(matched))

    # top up to ncands by S/N from the remaining (unmatched) XML candidates
    num_to_add = ncands - len(matched)
    if num_to_add > 0:
        remaining = xml_df[~xml_df["cand_id_in_file"].isin(matched["cand_id_in_file"])]
        remaining = remaining.sort_values("snr", ascending=False).head(num_to_add)
        matched = pd.concat([matched, remaining], ignore_index=True)

    matched = matched.sort_values("snr", ascending=False)
    matches_df = pd.DataFrame(match_records)
    return matched, matches_df


def write_output_candfile(path, df):
    with open(path, "w") as f:
        f.write("#id DM accel F0 F1 F2 S/N\n")
        for i, row in enumerate(df.itertuples(index=False)):
            f0 = 1.0 / float(row.period)
            f.write("%d %f %f %f 0 0 %f\n" % (i, float(row.dm), float(row.acc), f0, float(row.snr)))
    return path


def read_input_candfile(path):
    df = pd.read_csv(path, sep=r"\s+")
    # header is "#id DM acc F0 F1 F2 S/N" -> first column name is "#id"
    return df


def mode_match(args):
    input_df = read_input_candfile(args.input_candfile)
    xml_df = xml_to_df(args.xml)
    xml_df = reject_birdies(xml_df, args.birdies)

    matched, matches_df = shortlist(input_df, xml_df, args.ptol, args.dm_tol, args.ncands)

    write_output_candfile(os.path.join(args.out_dir, "output.candfile"), matched)
    matches_path = os.path.join(args.out_dir, "matches.csv")
    if matches_df.empty:
        # still emit a header-only file so downstream publishing has something
        pd.DataFrame(columns=[
            "input_id", "input_dm", "input_period", "input_f0", "input_acc",
            "match_cand_id", "match_dm", "match_period", "match_acc", "match_snr",
        ]).to_csv(matches_path, index=False)
    else:
        matches_df.to_csv(matches_path, index=False)
    print("[followup] wrote output.candfile (%d cands) and matches.csv (%d matches)"
          % (len(matched), len(matches_df)))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="ELDEN-RING followup candfile helper")
    p.add_argument("--mode", required=True, choices=["prepare", "match"])
    p.add_argument("--out_dir", default=".", help="Output directory")

    # prepare
    p.add_argument("--input_csv", help="CandyJar classification CSV (prepare mode)")
    p.add_argument("--dm_tol", type=float, default=10.0,
                   help="Half-width of each candidate's DM band (pc/cc)")
    p.add_argument("--dm_step", type=float, default=0.1, help="DM step (pc/cc)")
    p.add_argument("--use_zero_dm", action="store_true", help="Include DM=0 in dm file")
    # acceleration tiering
    p.add_argument("--acc_floor", type=float, default=50.0,
                   help="Tier-1 acc cap = max(median|acc|, acc_floor); "
                        "candidates within +/-cap searched together (default: 50)")
    p.add_argument("--mid_cap", type=float, default=150.0,
                   help="Tier-2 upper acc bound; cands in (cap, mid_cap] searched at "
                        "+/-mid_cap, above at +/-max|acc| of that group (default: 150)")

    # match
    p.add_argument("--input_candfile", help="input.candfile from prepare (match mode)")
    p.add_argument("--xml", help="peasoup overview.xml (match mode)")
    p.add_argument("--ptol", type=float, default=0.01,
                   help="Fractional period tolerance for matching")
    p.add_argument("--ncands", type=int, default=1000,
                   help="Cap on candidates to fold; top-up by S/N after matches")
    p.add_argument("--birdies", default=None,
                   help="Comma-separated periods (s) to reject as birdies")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    if args.mode == "prepare":
        if not args.input_csv:
            sys.exit("[followup] --input_csv is required in prepare mode")
        mode_prepare(args)
    else:
        if not (args.input_candfile and args.xml):
            sys.exit("[followup] --input_candfile and --xml are required in match mode")
        mode_match(args)


if __name__ == "__main__":
    main()
