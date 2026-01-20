# Output Structure

This page describes the directory structure and files produced by ELDEN-RING.

## Overview

ELDEN-RING creates two main output areas:

1. **Shared Cache** (`shared_cache/`) - Reusable intermediate products
2. **Run Directory** (`<runID>/`) - Run-specific outputs

```
basedir/
├── shared_cache/           # Persistent, reusable files
│   └── <cluster>/
│       ├── FITS/           # DADA→FITS conversions
│       └── <beam_name>/
│           ├── RFIFILTER/  # RFI diagnostic plots
│           └── CLEANEDFIL/ # Cleaned filterbanks
│
├── <runID>/                # Run-specific outputs
│   ├── <beam_name>/
│   │   └── segment_<N>/
│   │       └── <seg_id>/
│   │           ├── BIRDIES/
│   │           ├── SEARCH/
│   │           ├── PARSEXML/
│   │           ├── FOLDING/
│   │           ├── ABG/
│   │           └── CLASSIFICATION/
│   │
│   ├── TARBALL_CSV/
│   ├── CANDIDATE_TARBALLS/
│   ├── DMFILES/
│   └── pipeline_summary_*.txt
│
└── .cumulative_runtime_*.txt
```

---

## Shared Cache (`shared_cache/`)

Contains reusable intermediate files that persist across runs.

### FITS Directory

DADA to FITS conversions (when using `run_dada_search`):

```
shared_cache/<cluster>/FITS/
├── <beam_name>_cdm47.0.sf
├── <beam_name>_cdm101.0.sf
└── ...
```

### RFI Filter Directory

RFI diagnostic outputs:

```
shared_cache/<cluster>/<beam_name>/RFIFILTER/
├── <beam_name>_rfi_report.png      # RFI diagnostic plot
├── <beam_name>_bandpass.png        # Bandpass plot
├── <beam_name>_kurtosis.png        # Spectral kurtosis plot
├── <beam_name>_filter_params.txt   # Generated filter parameters
└── <beam_name>_masked_channels.txt # List of masked channels
```

### Cleaned Filterbank Directory

RFI-cleaned filterbank files:

```
shared_cache/<cluster>/<beam_name>/CLEANEDFIL/
├── <beam_name>_cdm47.0_cleaned.fil
├── <beam_name>_cdm101.0_cleaned.fil
└── ...
```

---

## Run Directory (`<runID>/`)

Contains all outputs specific to a particular search run.

### Beam/Segment Structure

Each beam's outputs are organized by segmentation:

```
<runID>/<beam_name>/segment_<N>/<seg_id>/
```

Where:
- `<N>` = segment factor (1=full, 2=half, 4=quarter)
- `<seg_id>` = segment identifier

### BIRDIES Directory

RFI birdie detection outputs:

```
BIRDIES/
├── birdies.xml              # Detected birdies in XML format
├── birdies_list.txt         # Parsed birdie frequencies
└── birdie_candidates.png    # Diagnostic plot
```

### SEARCH Directory

Peasoup search outputs:

```
SEARCH/
├── overview.xml             # Search overview
├── candidates.xml           # Raw candidate list
├── search_log.txt           # Search log
└── search_params.txt        # Parameters used
```

### PARSEXML Directory

Parsed and filtered candidates:

```
PARSEXML/
├── XML/
│   ├── filtered_candidates.xml   # Filtered candidates
│   └── candidate_summary.csv     # Summary statistics
├── candidates_above_threshold.csv
└── duplicate_removed.csv
```

### FOLDING Directory

PulsarX folding outputs:

```
FOLDING/
├── PNG/                     # Diagnostic plots
│   ├── candidate_001.png
│   ├── candidate_002.png
│   └── ...
├── AR/                      # Archive files
│   ├── candidate_001.ar
│   ├── candidate_002.ar
│   └── ...
├── CANDS/                   # Candidate metadata
│   ├── candidate_001.cands
│   └── ...
├── CSV/                     # Merged CSV files
│   └── all_candidates.csv
└── PROVENANCE/              # Tracking files
    └── fold_provenance.json
```

### ABG Directory

Alpha-Beta-Gamma scoring:

```
ABG/
├── abg_scores.csv           # Candidate scores
└── abg_ranking.csv          # Ranked candidates
```

### CLASSIFICATION Directory

ML classification outputs:

```
CLASSIFICATION/
├── pics_scores.csv          # PICS classifier scores
├── ranked_candidates.csv    # Final ranked list
└── classification_summary.txt
```

---

## Final Outputs

### CANDIDATE_TARBALLS Directory

Final CandyJar-compatible tarballs:

```
<runID>/CANDIDATE_TARBALLS/
├── <cluster>_<beam_name>_segment1_candidates.tar
├── <cluster>_<beam_name>_segment2_candidates.tar
└── ...
```

Each tarball contains:
```
candidates.tar/
├── candidates.csv           # Candidate metadata
├── plots/                   # PNG diagnostic plots
├── archives/                # Archive files
└── metadata.json            # Run metadata
```

### TARBALL_CSV Directory

CSV files used for tarball creation:

```
<runID>/TARBALL_CSV/
├── segment1_candidates.csv
├── segment2_candidates.csv
└── merged_all_segments.csv
```

### DMFILES Directory

DM search parameter files:

```
<runID>/DMFILES/
├── dm_range_1.txt
├── dm_range_2.txt
└── ...
```

---

## Summary Files

### Pipeline Summary

Generated after each run:

```
<runID>/pipeline_summary_<timestamp>.txt
```

Contains:
- Run parameters
- Processing statistics
- Candidate counts per stage
- Timing information
- Error summary

### Cumulative Runtime

Tracks total processing time across resume runs:

```
basedir/.cumulative_runtime_<runID>.txt
```

---

## Output File Formats

### Candidate CSV Format

```csv
id,dm,period,snr,acc,pics_score,alpha,beta,gamma,rank
1,47.123,0.00567,12.5,0.0,0.85,0.92,0.78,0.65,1
2,47.456,0.01234,10.2,-5.2,0.72,0.88,0.71,0.58,2
```

| Column | Description |
|--------|-------------|
| `id` | Candidate ID |
| `dm` | Dispersion measure (pc/cm³) |
| `period` | Period (seconds) |
| `snr` | Signal-to-noise ratio |
| `acc` | Acceleration (m/s²) |
| `pics_score` | ML classification score (0-1) |
| `alpha` | Alpha ranking metric |
| `beta` | Beta ranking metric |
| `gamma` | Gamma ranking metric |
| `rank` | Final ranking |

### Archive File Format

PulsarX archives (`.ar` files) are PSRFITS-compatible and can be viewed with:

```bash
# Using PSRCHIVE
pav -G your_candidate.ar

# Using PulsarX
pxview your_candidate.ar
```

---

## Cleaning Up

### Remove Run-Specific Outputs

```bash
rm -rf basedir/<runID>/
```

### Clean Shared Cache

Use the cleanup workflow:

```bash
nextflow run elden.nf -entry cleanup_cache --basedir /path/to/basedir
```

Or manually:

```bash
rm -rf basedir/shared_cache/<cluster>/<beam_name>/
```

### Clean Nextflow Work Directory

```bash
rm -rf work/
nextflow clean -f
```
