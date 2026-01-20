# Configuration

This page provides a comprehensive reference for all ELDEN-RING configuration parameters.

## Configuration Files

### Hierarchy

Configuration is loaded in the following order (later overrides earlier):

1. `nextflow.config` - Default pipeline settings
2. `conf/profiles/*.config` - Cluster-specific profiles
3. `params.config` - User parameters (in basedir)
4. Command-line parameters (`--param value`)

### Example Configuration

See `example/params.config.example` for a complete template.

## Core Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `params.basedir` | String | Base directory for all outputs |
| `params.runID` | String | Unique identifier for this run |
| `params.files_list` | String | Path to input CSV file |
| `params.telescope` | String | Telescope name (`effelsberg` or `meerkat`) |

```groovy
params {
    basedir = "/data/NGC6401_search"
    runID = "search_v1"
    files_list = "inputfile.txt"
    telescope = "effelsberg"
}
```

## Dispersion Measure (DM) Parameters

```groovy
params.ddplan {
    dm_start = 85       // Starting DM (relative to coherent DM)
    dm_end = 95         // Ending DM
    dm_step = 0.1       // DM step size
}
```

**Note**: DM values are typically specified relative to the coherent DM of the observation.

## Peasoup Search Parameters

```groovy
params.peasoup {
    // Segmentation
    segments = [1, 2, 4]           // 1=full, 2=half, 4=quarter segments

    // Acceleration search range (m/s²)
    acc_start = -50
    acc_end = 50

    // Signal-to-noise thresholds
    min_snr = 7.0                  // Minimum S/N for candidates
    birdies_min_snr = 11.0         // S/N for birdie detection

    // Harmonic summing
    nharmonics = 4                 // Number of harmonics to sum

    // Additional options
    extra_args = ""                // Extra peasoup arguments
}
```

## RFI Filtering Parameters

### filtool Settings

```groovy
params.filtool {
    run_filtool = true             // Enable RFI filtering
    rfi_filter = ""                // Path to RFI filter file (or auto-generate)
    extra_args = ""                // Extra filtool arguments
}
```

### RFI Filter Generation

```groovy
params.generateRfiFilter {
    run_rfi_filter = true          // Generate dynamic RFI masks
    bandpass_sigma = 3.0           // Sigma threshold for bandpass
    spectral_kurtosis_sigma = 3.0  // Sigma for spectral kurtosis
    time_sigma = 3.0               // Sigma for time-domain RFI
}
```

## Stacking and Splitting

```groovy
params {
    // Coherent DM stacking
    stack_by_cdm = true            // Stack beams by coherent DM value

    // Frequency splitting
    split_fil = true               // Split filterbank by frequency
    split_freq = 1200.0            // Split frequency (MHz)
}
```

## Folding Parameters

```groovy
params.pulsarx {
    nsubint = 64                   // Number of sub-integrations
    nbin = 128                     // Number of phase bins
    extra_args = ""                // Extra PulsarX arguments
}
```

## Classification Parameters

```groovy
params.classification {
    pics_model = "MeerKAT_L_SBAND_COMBINED_Best_Recall.pkl"
    min_pics_score = 0.1           // Minimum PICS score threshold
}
```

## Candidate Filtering

```groovy
params.parsexml {
    // Candidate selection
    min_snr = 7.0                  // Minimum S/N
    max_dm = 1000.0                // Maximum DM
    min_period = 0.0005            // Minimum period (seconds)
    max_period = 30.0              // Maximum period (seconds)

    // Duplicate removal
    dm_tolerance = 0.5             // DM matching tolerance
    period_tolerance = 0.0001      // Period matching tolerance (seconds)
}
```

## Cluster Profiles

### Available Profiles

| Profile | System | Description |
|---------|--------|-------------|
| `local` | Local | Single-machine execution |
| `hercules` | SLURM | Hercules HPC cluster |
| `edgar` | SLURM | Edgar cluster |
| `contra` | SLURM | Contra cluster |
| `condor` | HTCondor | HTCondor batch system |

### Using Profiles

```bash
# Single profile
nextflow run elden.nf -profile hercules

# Multiple profiles
nextflow run elden.nf -profile hercules,gpu
```

### Profile Configuration Example (SLURM)

```groovy
// conf/profiles/hercules.config.example
process {
    executor = 'slurm'
    queue = 'normal'

    withLabel: 'gpu' {
        queue = 'gpu'
        clusterOptions = '--gres=gpu:1'
    }

    withLabel: 'high_memory' {
        memory = '128 GB'
        cpus = 16
    }
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/scratch/singularity_cache'
}
```

## Container Configuration

```groovy
singularity {
    enabled = true
    autoMounts = true
    cacheDir = "${HOME}/.singularity_cache"
}

// Container paths
params.containers {
    peasoup = "/path/to/peasoup.sif"
    pulsarx = "/path/to/pulsarx.sif"
    presto = "/path/to/presto.sif"
    pics_classifier = "/path/to/pics.sif"
    rfi_mitigation = "/path/to/rfi.sif"
    filtools = "/path/to/filtools.sif"
}
```

## Email Notifications

```groovy
params {
    email = "user@example.com"     // Notification email
    email_on_fail = true           // Email on failure
    email_on_complete = true       // Email on completion
}
```

## Advanced Parameters

### Caching and Performance

```groovy
params {
    // Shared cache for reusable outputs
    use_shared_cache = true
    shared_cache_dir = "${params.basedir}/shared_cache"

    // Process-level settings
    max_retries = 3                // Retry failed processes
    error_strategy = 'retry'       // or 'ignore', 'terminate'
}
```

### DADA Processing

```groovy
params.dada {
    dada_files_list = "dada_input.csv"
    cdm_list = "47.0 101.0"        // Coherent DM values to process
}

params.digifits {
    extra_args = ""                // Extra digifits arguments
}
```

## Command-Line Overrides

Any parameter can be overridden from the command line:

```bash
nextflow run elden.nf -entry full \
    --basedir /new/path \
    --runID "test_run" \
    --peasoup.min_snr 8.0 \
    --ddplan.dm_start 50 \
    --ddplan.dm_end 150
```

## Environment Variables

| Variable | Description |
|----------|-------------|
| `NXF_WORK` | Override work directory location |
| `NXF_SINGULARITY_CACHEDIR` | Singularity cache directory |
| `NXF_TEMP` | Temporary directory for Nextflow |

```bash
export NXF_WORK=/scratch/$USER/nextflow_work
export NXF_SINGULARITY_CACHEDIR=/scratch/singularity_cache
```
