# ELDEN-RING

**E**ffelsberg **L**arge-scale **D**ata **E**xploration with **N**extflow for **R**obust **I**dentification of **N**ew **G**lobular cluster pulsars.

![elden-ring-transformed](https://github.com/user-attachments/assets/3e1a1c35-d055-4266-9ff7-5380ab1d463f)

A GPU-accelerated Nextflow pipeline for pulsar candidate detection featuring RFI mitigation, periodicity searches with [peasoup](https://github.com/ewanbarr/peasoup) or [PRESTO](https://github.com/scottransom/presto), candidate folding with [PulsarX](https://github.com/ypmen/PulsarX) or prepfold, and machine learning classification.

## Table of Contents

- [Features](#features)
- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input File Formats](#input-file-formats)
- [Available Workflows](#available-workflows)
- [Configuration](#configuration)
- [Output Structure](#output-structure)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Features

- **Multiple Search Backends**: Choose between peasoup (GPU), PRESTO, or Riptide FFA for periodicity searches
- **GPU-Accelerated Search**: Fast periodicity searches using peasoup on NVIDIA GPUs
- **PRESTO Pipeline**: Full PRESTO support including rfifind, prepsubband, accelsearch, and prepfold
- **Riptide FFA Search**: Fast Folding Algorithm search for long-period pulsars (complementary to FFT-based searches)
- **Hybrid Mode**: Dump peasoup time series for subsequent PRESTO accelsearch or Riptide FFA processing
- **RFI Mitigation**: Automated RFI detection and filtering with spectral kurtosis
- **Multi-Beam Support**: Process multiple beams in parallel
- **Coherent Dedispersion**: Support for DADA baseband data with digifits conversion
- **Filterbank Stacking**: Stack multiple beams by coherent DM for improved sensitivity
- **Segmented Searches**: Search full observation and sub-segments for accelerated pulsars
- **Flexible Folding**: Fold candidates with PulsarX or PRESTO prepfold
- **ML Classification**: PICS-based candidate scoring
- **Alpha-Beta-Gamma Scoring**: Additional candidate ranking metrics
- **Resume Support**: Automatic caching and resume capability via Nextflow
- **Cumulative Runtime Tracking**: Track total processing time across resumed runs
- **Email Notifications**: Optional notifications on completion or failure
- **Input Validation**: Pre-flight checks for parameters and input files

## Pipeline Overview

### Peasoup Pipeline (GPU-accelerated FFT search)
```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Peasoup Pipeline (GPU)                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   DADA Files ──► digifits ──┐                                               │
│                             │                                               │
│   FITS/Filterbanks ─────────┼──► RFI Filter ──► filtool ──► Segmentation   │
│                                                                             │
│   Segmentation ──► birdies ──► peasoup (GPU) ──► XML Parse ──► PulsarX     │
│                                    │                                        │
│                                    └──► [dump_timeseries] ──► PRESTO/FFA   │
│                                                                             │
│   PulsarX ──► Merge Folds ──► PICS Classifier ──► Alpha-Beta-Gamma         │
│                                                                             │
│   Final Output: CandyJar tarball with ranked candidates                     │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### PRESTO Pipeline (CPU-based acceleration search)
```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           PRESTO Pipeline                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Filterbank ──► rfifind ──► prepdata (zero-DM) ──► accelsearch (birdies)  │
│                    │                                                        │
│                    └──► prepsubband ──► accelsearch ──► sift candidates    │
│                            (DM range)      (z/w)           │                │
│                                                            │                │
│                    ┌───────────────────────────────────────┘                │
│                    │                                                        │
│                    ▼                                                        │
│              [fold_backend]                                                 │
│                    │                                                        │
│        ┌──────────┴──────────┐                                             │
│        ▼                     ▼                                              │
│    prepfold              PulsarX ──► PNG plots                             │
│        │                                                                    │
│        └──► show_pfd ──► PNG plots                                         │
│                                                                             │
│   PNG + CSV ──► PICS Classifier ──► CandyJar tarball                       │
│                                                                             │
│   State files saved at each stage for resumption                           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Riptide FFA Pipeline (long-period pulsar search)
```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Riptide FFA Pipeline                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Filterbank ──► [PRESTO or Peasoup dedispersion] ──► .dat/.inf files      │
│                                                                             │
│   Time series ──► rffa (Fast Folding Algorithm) ──► candidates.csv         │
│                                                                             │
│   Output: candidates.csv + peaks.csv + clusters.csv + plots                │
│                                                                             │
│   Can run standalone or alongside Peasoup/PRESTO pipelines                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

![mermaid-diagram-2025-04-24-115529](https://github.com/user-attachments/assets/92b86bb0-0b40-4050-8526-5f380136cf62)

## Requirements

### Software

- **Nextflow** >= 21.10.0
- **Singularity** >= 3.0 (or Docker)
- **NVIDIA GPU** with CUDA support (for peasoup)

### Container Images

The pipeline uses containerized tools. Required images:

| Tool | Purpose |
|------|---------|
| `pulsarx_image` | Candidate folding (PulsarX) |
| `peasoup_image` | GPU periodicity search |
| `presto_image` | Filterbank utilities (readfile) |
| `presto5_image` | PRESTO 5 with pdot support |
| `prestozl_image` | GPU-accelerated PRESTO (PrestoZL) |
| `rfi_mitigation_image` | RFI analysis and filtering |
| `pics_classifier_image` | ML candidate classification |
| `edd_pulsar_image` | DADA to FITS conversion (digifits) |
| `filtools_sig_image` | Filterbank tools with signal injection |
| `rusty_candypicker` | Candidate filtering (Rust implementation) |
| `riptide_image` | Riptide FFA periodicity search |

## Installation

### Option 1: Clone the repository

```bash
git clone https://github.com/erc-compact/elden-ring.git
cd elden-ring
```

### Option 2: Use Nextflow's built-in pull (recommended)

```bash
nextflow pull erc-compact/elden-ring
```

## Quick Start

### 1. Initialize a new project

```bash
nextflow run elden.nf -entry setup_basedir --basedir /path/to/my_project
```

This creates:
```
/path/to/my_project/
├── params.config           # Main configuration (edit this)
├── inputfile.txt           # Input data CSV (edit this)
├── generate_inputfile.sh   # Helper script for input generation
├── meta/                   # Pipeline metadata
└── shared_cache/           # Reusable cached files
```

### 2. Generate your input file

For filterbank/FITS files:
```bash
cd /path/to/my_project
bash generate_inputfile.sh \
    --cluster NGC6544 \
    --ra "18:07:20.5" \
    --dec "-24:59:51" \
    --utc "2024-01-15T10:00:00" \
    --cdm "60.0" \
    /path/to/data/*.fil
```

For DADA baseband directories:
```bash
bash generate_inputfile.sh \
    --dada \
    --cluster 2MASS-GC02 \
    --ra "18:09:36.51" \
    --dec "+20:46:43.99" \
    --utc "2025-12-06T13:08:08" \
    --cdm "156.0 428.0 700.0" \
    /path/to/baseband3 /path/to/baseband4 /path/to/baseband5
```

### 3. Edit configuration

```bash
vim params.config
```

Key parameters to review:
- `basedir` - Output directory (auto-set by setup)
- `runID` - Unique identifier for this search
- `telescope` - Your telescope (effelsberg, meerkat, etc.)
- `ddplan.*` - DM search range
- `peasoup.*` - Search parameters (acceleration, segments, SNR threshold)

### 4. Run the pipeline

```bash
nextflow run elden.nf \
    -entry full \
    -profile hercules \
    -c params.config \
    --runID my_search_v1 \
    -resume
```

## Input File Formats

### Standard Input (inputfile.txt)

CSV format for filterbank/FITS files:

```csv
pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm
0,NGC6544,cfbf00001,1,2024-01-15T10:00:00,18:07:20.5,-24:59:51,/path/to/beam1.fil,60.0
0,NGC6544,cfbf00002,2,2024-01-15T10:00:00,18:07:20.5,-24:59:51,/path/to/beam2.fil,60.0
```

| Column | Description |
|--------|-------------|
| `pointing` | Pointing index (integer) |
| `cluster` | Target name / cluster identifier |
| `beam_name` | Beam identifier (e.g., cfbf00001) |
| `beam_id` | Numeric beam ID |
| `utc_start` | Observation start time (ISO format) |
| `ra` | Right ascension (HH:MM:SS.ss) |
| `dec` | Declination (DD:MM:SS.ss) |
| `fits_files` | Full path to filterbank/FITS file |
| `cdm` | Coherent dedispersion DM |

### DADA Input (dada_files.csv)

CSV format for DADA baseband directories:

```csv
pointing,dada_files,cluster,beam_name,beam_id,utc_start,ra,dec,cdm_list
0,/path/to/baseband3/*dada,2MASS-GC02,cfbf00003,3,2025-12-06T13:08:08,18:09:36.51,+20:46:43.99,156.0 428.0 700.0
0,/path/to/baseband4/*dada,2MASS-GC02,cfbf00004,4,2025-12-06T13:08:08,18:09:36.51,+20:46:43.99,156.0 428.0 700.0
```

Note: `cdm_list` contains space-separated coherent DM values. The pipeline will process each CDM independently.

## Available Workflows

Select a workflow with the `-entry` flag:

### Main Processing Pipelines (Peasoup)

| Workflow | Description |
|----------|-------------|
| `full` | Complete pipeline: intake → RFI → clean → search → fold → classify |
| `run_search_fold` | Search & fold on pre-cleaned filterbanks |
| `run_rfi_clean` | RFI cleaning only (intake → filter → clean) |
| `generate_rfi_filter` | Generate RFI diagnostic plots only |
| `search_pipeline` | Auto-select backend based on `params.search_backend` |

### PRESTO Pipelines

| Workflow | Description |
|----------|-------------|
| `presto_pipeline` | Full PRESTO pipeline for single filterbank file (recommended) |
| `presto_full_entry` | Full PRESTO pipeline from CSV input file |
| `presto_search_fold` | Resume from state file: search → fold → postprocess |
| `run_accelsearch_on_timeseries` | Run accelsearch on pre-dumped .dat/.inf files from peasoup |
| `peasoup_with_presto_search` | Hybrid: peasoup search + optional PRESTO accelsearch on dumped time series |

**PRESTO Pipeline Stages:**
1. **RFI Detection** (`presto_rfi`) - rfifind to create RFI mask
2. **Birdie Detection** (`presto_birdies`) - Zero-DM accelsearch for persistent signals
3. **Dedispersion** (`presto_dedisperse`) - prepsubband across DM range
4. **Acceleration Search** (`presto_search`) - accelsearch with GPU support
5. **Sift and Fold** (`presto_sift_fold`) - Candidate sifting + folding (PulsarX or prepfold)
6. **Post-processing** (`presto_postprocess`) - PNG conversion + tarball creation

Each stage saves a state file to `PRESTO_STATE/` for resumption with `presto_search_fold`.

### Riptide FFA Pipelines

| Workflow | Description |
|----------|-------------|
| `run_riptide` | Standalone Riptide FFA: dedisperse (PRESTO) → FFA search |
| `run_riptide_on_timeseries` | Run Riptide FFA on pre-existing .dat/.inf files |

**Note:** Riptide FFA can also be enabled alongside other pipelines using `--riptide.run_ffa_search true`

### DADA Processing Pipelines

| Workflow | Description |
|----------|-------------|
| `run_dada_search` | Full pipeline starting from DADA baseband files |
| `run_digifits` | Convert DADA to FITS/filterbank only |
| `run_dada_clean_stack` | DADA → FITS → clean → stack (no search) |

### Specialized Workflows

| Workflow | Description |
|----------|-------------|
| `fold_par` | Fold data using a known pulsar ephemeris (.par file) |
| `candypolice` | Re-fold candidates from an existing CandyJar CSV |

### Utility Workflows

| Workflow | Description |
|----------|-------------|
| `help` | Display detailed usage information |
| `setup_basedir` | Initialize a new project directory |
| `validate_inputs` | Validate input files and parameters |
| `cleanup_cache` | Find orphaned files in shared cache |

## Configuration

### Key Parameters

```groovy
// Required
params.basedir = "/path/to/project"
params.runID = "search_v1"
params.files_list = "inputfile.txt"
params.telescope = "effelsberg"

// Search Backend Selection
params.search_backend = "peasoup"  // Options: 'peasoup' (default) or 'presto'

// DM Search Range
params.ddplan.dm_start = -10    // Relative to coherent DM
params.ddplan.dm_end = 10
params.ddplan.dm_step = 0.1

// Peasoup Search
params.peasoup.segments = [1, 2, 4]   // Full, half, quarter segments
params.peasoup.acc_start = -50        // Acceleration range (m/s²)
params.peasoup.acc_end = 50
params.peasoup.min_snr = 8.0
params.peasoup.dump_timeseries = false  // Set true to dump .dat/.inf for PRESTO

// PRESTO Search (when search_backend = 'presto')
params.presto.zmax = 200              // Max z (acceleration) for accelsearch
params.presto.wmax = 0                // Max w (jerk) for accelsearch
params.presto.numharm = 8             // Number of harmonics to sum
params.presto.fold_backend = 'pulsarx'  // Options: 'pulsarx' or 'presto'
params.presto.dm_ranges = [
    [dm_low: 0.0, dm_high: 100.0, dm_step: 0.5, downsamp: 1]
]

// Processing Options
params.filtool.run_filtool = true
params.generateRfiFilter.run_rfi_filter = true
params.stack_by_cdm = false
params.split_fil = false

// Riptide FFA Search (complementary to FFT-based searches)
params.riptide.run_ffa_search = false         // Enable FFA search alongside accelsearch
params.riptide.config_file = "riptide_config.yml"  // YAML config for rffa (edit directly)
params.riptide.backend = 'presto'             // Backend for standalone: 'presto' or 'peasoup'

// Notifications (optional)
params.notification.enabled = true
params.notification.email = "user@example.com"
params.notification.on_complete = true
params.notification.on_fail = true
```

### Cluster Profiles

Select a profile with `-profile`:

| Profile | Description |
|---------|-------------|
| `local` | Local execution (testing) |
| `hercules` | SLURM cluster with GPU nodes |
| `edgar` | Edgar cluster configuration |
| `contra` | Contra cluster configuration |
| `condor` | HTCondor submission |

Create custom profiles in `conf/profiles/`.

## Output Structure

### Peasoup Pipeline Output
```
basedir/
├── shared_cache/                    # Reusable cached files
│   └── <cluster>/
│       ├── FITS/                    # Converted FITS files (from DADA)
│       └── <beam_name>/
│           ├── RFIFILTER/           # RFI diagnostic plots
│           └── CLEANEDFIL/          # Cleaned filterbanks
│
├── <runID>/                         # Run-specific outputs
│   ├── <beam_name>/
│   │   └── segment_<N>/
│   │       └── <seg_id>/
│   │           ├── BIRDIES/         # Birdie detection files
│   │           ├── SEARCH/          # Peasoup XML results
│   │           ├── PARSEXML/        # Parsed candidates
│   │           │   └── XML/         # Filtered XML files
│   │           ├── FOLDING/         # PulsarX outputs
│   │           │   ├── PNG/         # Diagnostic plots
│   │           │   ├── AR/          # Archive files
│   │           │   ├── CANDS/       # .cands files
│   │           │   ├── CSV/         # Merged CSVs
│   │           │   └── PROVENANCE/  # Tracking files
│   │           ├── ABG/             # Alpha-beta-gamma scores
│   │           ├── ZERODM/          # Zero-DM plots
│   │           └── CLASSIFICATION/  # PICS scores
│   │
│   ├── TARBALL_CSV/                 # CSV files for tarball
│   ├── CANDIDATE_TARBALLS/          # Final candidate packages
│   ├── DMFILES/                     # DM search files
│   └── pipeline_summary_*.txt       # Run summary
│
└── .cumulative_runtime_*.txt        # Runtime tracking
```

### PRESTO Pipeline Output
```
basedir/
├── <runID>/
│   ├── sharedcache/presto/          # Shareable/reusable artifacts
│   │   ├── rfi/                     # RFI masks (*.mask)
│   │   ├── subbands/<dm_range>/     # Dedispersed time series
│   │   └── timeseries/<dm_range>/   # Peasoup time series dumps
│   │
│   ├── PRESTO_RFI/                  # RFI stats, inf, out files
│   ├── PRESTO_SEARCH/               # Accelsearch outputs
│   │   ├── full/                    # Full observation search
│   │   └── segmented/               # Segmented search (if enabled)
│   ├── PRESTO_SIFTED/               # Sifted candidate files
│   │   ├── sifted_candidates.csv
│   │   ├── sifted_candidates.candfile
│   │   └── sifted_candidates.provenance.csv
│   ├── PRESTO_FOLDING/
│   │   ├── PFD/                     # .pfd files (prepfold backend)
│   │   ├── BESTPROF/                # .bestprof files
│   │   ├── PNG/                     # Diagnostic plots
│   │   ├── MERGED/                  # Merged results CSV
│   │   └── PULSARX/                 # PulsarX fold outputs
│   ├── PRESTO_CLASSIFICATION/       # PICS classification results
│   ├── PRESTO_TARBALLS/             # Final CandyJar-compatible tarballs
│   └── PRESTO_STATE/                # State files for stage resume
│       ├── rfi_state.json
│       ├── birdies_state.json
│       ├── dedisperse_state.json
│       ├── search_state.json
│       └── sift_fold_state.json
```

### Riptide FFA Output
```
basedir/
├── <runID>/
│   └── RIPTIDE_SEARCH/
│       ├── candidates.csv      # Final candidates (filtered, harmonics removed)
│       ├── peaks.csv           # All detected periodogram peaks
│       ├── clusters.csv        # Peaks grouped by frequency proximity
│       ├── *.json              # JSON file per candidate (loadable with riptide.load_json)
│       └── *.png               # Diagnostic plot per candidate
```

## Advanced Usage

### Resume a Failed Run

```bash
nextflow run elden.nf -entry full -profile hercules -c params.config -resume
```

### Validate Inputs Before Running

```bash
nextflow run elden.nf -entry validate_inputs -c params.config
```

### Clean Up Orphaned Cache Files

```bash
# Dry run (shows what would be deleted)
nextflow run elden.nf -entry cleanup_cache --basedir /path/to/project

# Actually delete orphaned files
bash scripts/cleanup_shared_cache.sh /path/to/project false
```

### Fold with Known Pulsar Ephemeris

```bash
nextflow run elden.nf -entry fold_par \
    -c params.config \
    --parfold.parfile_path /path/to/pulsar.par
```

### Re-fold Candidates from CandyJar

```bash
nextflow run elden.nf -entry candypolice \
    -c params.config \
    --candypolice.input_csv /path/to/candyjar.csv
```

### Run PRESTO Pipeline

```bash
# Full PRESTO pipeline with single filterbank file
nextflow run elden.nf -entry presto_pipeline \
    -profile hercules \
    --input_fil /path/to/file.fil \
    --presto.dm_ranges '[{"dm_low": 0, "dm_high": 100, "dm_step": 0.5, "downsamp": 1}]' \
    --presto.fold_backend pulsarx

# Full PRESTO pipeline with CSV input
nextflow run elden.nf -entry presto_full_entry \
    -profile hercules \
    -c params.config \
    --presto.fold_backend pulsarx

# Resume from dedispersion state file (search + fold + postprocess)
nextflow run elden.nf -entry presto_search_fold \
    -profile hercules \
    --state_file /path/to/PRESTO_STATE/dedisperse_state.json \
    --presto.fold_backend pulsarx

# Run individual stages for debugging (each saves state files)
# Example: Run only RFI detection
nextflow run elden.nf -entry presto_rfi \
    -profile hercules \
    --input_fil /path/to/file.fil
```

### Hybrid Peasoup + PRESTO Search

Run peasoup with time series dumping, then process with PRESTO accelsearch:

```bash
# Step 1: Run peasoup with dump_timeseries enabled
nextflow run elden.nf -entry full \
    -profile hercules \
    -c params.config \
    --peasoup.dump_timeseries true

# Step 2: Run accelsearch on dumped time series
nextflow run elden.nf -entry run_accelsearch_on_timeseries \
    -profile hercules \
    -c params.config
```

### Choose Folding Backend for PRESTO Candidates

PRESTO candidates can be folded with either PulsarX or prepfold:

```bash
# Fold with PulsarX (recommended for better plots)
nextflow run elden.nf -entry presto_full \
    -c params.config \
    --presto.fold_backend pulsarx

# Fold with PRESTO prepfold (traditional .pfd files)
nextflow run elden.nf -entry presto_full \
    -c params.config \
    --presto.fold_backend presto
```

### Run Riptide FFA Search

Riptide performs Fast Folding Algorithm (FFA) searches, which are sensitive to long-period pulsars that may be missed by FFT-based searches. The FFA search can be run standalone or alongside existing pipelines.

**Configuration:** Edit `riptide_config.yml` in your basedir (created during `setup_basedir`). See [riptide documentation](https://riptide-ffa.readthedocs.io/en/latest/pipeline.html) for parameter details.

```bash
# Add FFA search to peasoup pipeline (runs in parallel with peasoup)
nextflow run elden.nf -entry full \
    -profile hercules \
    -c params.config \
    --riptide.run_ffa_search true \
    --peasoup.dump_timeseries true

# Add FFA search to PRESTO pipeline (runs after dedispersion)
nextflow run elden.nf -entry presto_pipeline \
    -profile hercules \
    --input_fil /path/to/file.fil \
    --riptide.run_ffa_search true

# Add FFA search to accelsearch on existing time series
nextflow run elden.nf -entry run_accelsearch_on_timeseries \
    -profile hercules \
    -c params.config \
    --riptide.run_ffa_search true

# Standalone Riptide with PRESTO dedispersion
nextflow run elden.nf -entry run_riptide \
    -profile hercules \
    --input_fil /path/to/file.fil \
    --riptide.backend presto

# Standalone Riptide on pre-existing time series
nextflow run elden.nf -entry run_riptide_on_timeseries \
    -profile hercules \
    --timeseries_input_dir /path/to/TIMESERIES
```

### Copy Data from Remote Cluster

Enable in params.config:
```groovy
params.copy_from_tape.run_copy = true
params.copy_from_tape.remoteUser = "username"
params.copy_from_tape.remoteHost = "remote.cluster.edu"
```

## Troubleshooting

### Check Nextflow Logs

```bash
# View recent log
cat .nextflow.log

# View execution history
nextflow log

# View specific run
nextflow log <run_name> -f name,status,exit,duration
```

### Common Issues

**GPU not detected**
- Ensure CUDA drivers are installed
- Check Singularity GPU bindings: `singularity exec --nv`

**Out of memory**
- Reduce `params.peasoup.segments` to fewer segments
- Adjust SLURM memory requests in profile

**Missing input files**
- Run `validate_inputs` workflow to check paths
- Verify CSV file format matches expected columns

**Cache corruption**
- Delete `work/` directory and re-run with `-resume`
- Clean shared_cache if needed

### Get Help

```bash
nextflow run elden.nf -entry help
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This project is part of the ERC COMPACT project.

## Contact

- Open an issue: [https://github.com/erc-compact/elden-ring/issues](https://github.com/erc-compact/elden-ring/issues)
- Email: fkareem[at]mpifr-bonn.mpg.de
