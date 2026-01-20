# Workflows

ELDEN-RING provides multiple entry points (workflows) for different use cases. This page describes each workflow and when to use it.

## Main Workflows

### `full` - Complete Pipeline

The primary workflow that runs the entire pulsar search pipeline.

```bash
nextflow run elden.nf -entry full -c params.config
```

**Pipeline stages:**
1. Intake - Load filterbank files from CSV
2. RFI Filter - Generate RFI diagnostic masks
3. RFI Clean - Apply filtool filtering
4. Stack (optional) - Combine beams by coherent DM
5. Segmentation - Split into temporal segments
6. Search - GPU-accelerated peasoup search
7. Parse XML - Filter and parse candidates
8. Fold - PulsarX candidate folding
9. Classify - PICS ML classification
10. Package - Create CandyJar tarballs

**Input**: `inputfile.txt` (filterbank CSV)
**Output**: CandyJar tarballs with ranked candidates

---

### `run_dada_search` - DADA Baseband Processing

Processes DADA baseband files through the complete pipeline.

```bash
nextflow run elden.nf -entry run_dada_search -c params.config
```

**Additional stages:**
- DADA to FITS conversion using digifits
- Coherent dedispersion at specified DM values

**Input**: `dada_input.csv` (DADA file paths)
**Output**: CandyJar tarballs

---

## Specialized Workflows

### `run_rfi_clean` - RFI Cleaning Only

Performs only RFI filtering without searching.

```bash
nextflow run elden.nf -entry run_rfi_clean -c params.config
```

**Use case**: Prepare cleaned filterbanks for later processing or inspection.

**Output**: Cleaned filterbank files in `shared_cache/`

---

### `run_search_fold` - Search Pre-Cleaned Data

Runs search and folding on already-cleaned filterbanks.

```bash
nextflow run elden.nf -entry run_search_fold -c params.config
```

**Use case**: Re-run searches with different parameters on previously cleaned data.

---

### `generate_rfi_filter` - RFI Diagnostics

Generates RFI diagnostic plots without cleaning or searching.

```bash
nextflow run elden.nf -entry generate_rfi_filter -c params.config
```

**Use case**: Analyze RFI characteristics before running the full pipeline.

**Output**: RFI diagnostic plots and filter parameters

---

### `fold_par` - Fold with Ephemeris

Folds data using known pulsar ephemeris files (.par).

```bash
nextflow run elden.nf -entry fold_par \
    -c params.config \
    --par_file /path/to/pulsar.par
```

**Use case**: Confirm known pulsars or timing observations.

**Input**: Filterbank files + ephemeris file
**Output**: Folded archives and diagnostic plots

---

### `candypolice` - Re-fold Candidates

Re-folds candidates from an existing CandyJar CSV file.

```bash
nextflow run elden.nf -entry candypolice \
    -c params.config \
    --candyjar_csv /path/to/candidates.csv
```

**Use case**: Re-process interesting candidates with different folding parameters.

---

### `run_digifits` - DADA to FITS Only

Converts DADA baseband files to FITS format without further processing.

```bash
nextflow run elden.nf -entry run_digifits -c params.config
```

**Use case**: Prepare FITS files for external tools or archival.

---

### `run_dada_clean_stack` - DADA to Stacked Filterbanks

Converts DADA files, cleans, and stacks without searching.

```bash
nextflow run elden.nf -entry run_dada_clean_stack -c params.config
```

**Use case**: Prepare stacked filterbanks for custom analysis.

---

## Utility Workflows

### `setup_basedir` - Initialize Project

Creates the standard project directory structure.

```bash
nextflow run elden.nf -entry setup_basedir --basedir /path/to/project
```

**Creates:**
- `params.config` - Configuration template
- `inputfile.txt` - Input file template
- `dada_input.csv` - DADA input template
- `known_pulsars.csv` - Known pulsar list template

---

### `validate_inputs` - Pre-flight Validation

Validates input files and configuration before running.

```bash
nextflow run elden.nf -entry validate_inputs -c params.config
```

**Checks:**
- Input file existence and format
- Parameter validity
- Container availability
- Directory permissions

---

### `cleanup_cache` - Cache Management

Analyzes and optionally cleans orphaned cache files.

```bash
nextflow run elden.nf -entry cleanup_cache --basedir /path/to/project
```

**Use case**: Reclaim disk space from interrupted runs.

---

### `help` - Display Help

Shows usage information and available workflows.

```bash
nextflow run elden.nf -entry help
```

---

## Workflow Selection Guide

| Scenario | Recommended Workflow |
|----------|---------------------|
| First-time processing of raw filterbank data | `full` |
| Processing DADA baseband files | `run_dada_search` |
| Inspecting data quality before processing | `generate_rfi_filter` |
| Re-running search with new parameters | `run_search_fold` |
| Processing known pulsars | `fold_par` |
| Re-examining specific candidates | `candypolice` |
| Preparing data for external tools | `run_rfi_clean` or `run_digifits` |

## Workflow Dependencies

```
                    ┌─────────────────┐
                    │   setup_basedir │
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
                    │ validate_inputs │
                    └────────┬────────┘
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
┌───────▼───────┐   ┌────────▼────────┐   ┌───────▼───────┐
│ run_rfi_clean │   │      full       │   │run_dada_search│
└───────┬───────┘   └────────┬────────┘   └───────┬───────┘
        │                    │                    │
        └──────────┬─────────┴─────────┬──────────┘
                   │                   │
          ┌────────▼────────┐ ┌────────▼────────┐
          │ run_search_fold │ │   candypolice   │
          └─────────────────┘ └─────────────────┘
```

## Combining Workflows

You can chain workflows by using the outputs of one as inputs to another:

```bash
# Step 1: Clean data
nextflow run elden.nf -entry run_rfi_clean -c params.config

# Step 2: Search with different parameters
nextflow run elden.nf -entry run_search_fold \
    -c params.config \
    --peasoup.acc_end 100 \
    --runID search_high_acc
```
