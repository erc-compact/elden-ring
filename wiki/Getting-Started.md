# Getting Started

This guide will help you set up and run your first ELDEN-RING pipeline.

## Prerequisites

### Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| Nextflow | 21.10.0+ | Pipeline orchestration |
| Singularity/Apptainer | 3.0+ | Container runtime |
| NVIDIA Driver | Compatible with CUDA | GPU acceleration |

### Hardware Requirements

- **CPU**: Multi-core processor (recommended: 16+ cores)
- **RAM**: 64GB+ recommended for large datasets
- **GPU**: NVIDIA GPU with CUDA support (for peasoup searches)
- **Storage**: Sufficient space for input data and outputs

## Installation

### 1. Install Nextflow

```bash
# Using curl
curl -s https://get.nextflow.io | bash

# Or using conda
conda install -c bioconda nextflow
```

### 2. Clone the Repository

```bash
git clone https://github.com/your-org/elden-ring.git
cd elden-ring
```

### 3. Verify Container Access

Ensure Singularity containers are accessible. The pipeline uses pre-built containers for:
- `peasoup` - GPU periodicity search
- `pulsarx` - Candidate folding
- `presto` - Pulsar utilities
- `pics_classifier` - ML classification
- `rfi_mitigation` - RFI filtering
- `filtools` - Filterbank tools

## Quick Start

### Step 1: Setup Project Directory

```bash
nextflow run elden.nf -entry setup_basedir --basedir /path/to/your/project
```

This creates the necessary directory structure:
```
/path/to/your/project/
├── params.config
├── inputfile.txt
├── dada_input.csv (if applicable)
└── known_pulsars.csv
```

### Step 2: Generate Input File

Use the helper script to generate the input CSV:

```bash
bash generate_inputfile.sh \
    --cluster NGC6401 \
    --ra 17:38:36.83 \
    --dec -23:54:33.9 \
    --utc 2025-09-14T16:01:33 \
    --cdm "47.0 101.0" \
    /path/to/filterbank/data/
```

Or manually create `inputfile.txt`:

```csv
pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm
0,NGC6401,cfbf00003,3,2025-09-14T16:01:33.9,17:38:36.53,-23:54:34.6,/data/beam3.fil,47.0
```

### Step 3: Configure Parameters

Edit `params.config` in your project directory:

```groovy
params {
    basedir = "/path/to/your/project"
    runID = "search_v1"
    files_list = "inputfile.txt"
    telescope = "effelsberg"

    ddplan {
        dm_start = 85
        dm_end = 95
        dm_step = 0.1
    }

    peasoup {
        segments = [1, 2, 4]
        acc_start = -50
        acc_end = 50
        min_snr = 7.0
    }
}
```

### Step 4: Validate Inputs

Before running the full pipeline, validate your configuration:

```bash
nextflow run elden.nf -entry validate_inputs \
    -c params.config \
    --basedir /path/to/your/project
```

### Step 5: Run the Pipeline

```bash
# Full pipeline with local profile
nextflow run elden.nf -entry full \
    -profile local \
    -c params.config \
    --runID search_v1

# Or with SLURM cluster profile
nextflow run elden.nf -entry full \
    -profile hercules \
    -c params.config \
    --runID search_v1 \
    -resume
```

## Common Entry Points

| Entry Point | Use Case |
|-------------|----------|
| `-entry full` | Complete pipeline from raw data to classified candidates |
| `-entry run_dada_search` | Process DADA baseband files |
| `-entry run_rfi_clean` | RFI cleaning only |
| `-entry run_search_fold` | Search and fold pre-cleaned files |
| `-entry fold_par` | Fold with known pulsar ephemeris |

## Resuming Interrupted Runs

Nextflow supports resuming interrupted pipelines:

```bash
nextflow run elden.nf -entry full \
    -profile hercules \
    -c params.config \
    -resume
```

The `-resume` flag uses cached results from previous runs.

## Next Steps

- [Configuration](Configuration) - Detailed parameter reference
- [Workflows](Workflows) - Understanding available workflows
- [Input Formats](Input-Formats) - Input file specifications
