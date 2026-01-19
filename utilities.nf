#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ELDEN-RING: Utility Workflows
// ============================================================================
// This file contains utility workflows for setup, validation, cleanup, and help.
// These are separated from the main pipeline workflows for cleaner organization.
// ============================================================================

// ============================================================================
// HELP WORKFLOW - Display usage information
// ============================================================================
workflow help {
    main:
    def helpMessage = """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                           ELDEN-RING Pipeline v1.0                           ║
    ║              GPU-Accelerated Pulsar Periodicity Search Pipeline              ║
    ╚══════════════════════════════════════════════════════════════════════════════╝

    DESCRIPTION:
    ------------
    ELDEN-RING is a Nextflow pipeline for pulsar candidate detection using:
    • RFI mitigation and filterbank cleaning (filtool)
    • GPU-accelerated periodicity searches (peasoup)
    • Candidate folding (pulsarx)
    • Machine learning classification (PICS)
    • Alpha-beta-gamma scoring

    QUICK START:
    ------------
    1. Setup a new project:
       nextflow run elden.nf -entry setup_basedir --basedir /path/to/project

    2. Edit the generated config files:
       - params.config (search parameters)
       - inputfile.txt (input data files)

    3. Run the full pipeline:
       nextflow run elden.nf -entry full -profile hercules -c params.config

    AVAILABLE WORKFLOWS (-entry):
    -----------------------------
    Main Pipelines:
      full                 Full pipeline: intake → RFI → clean → search → fold → classify
      run_search_fold      Search & fold on pre-cleaned filterbanks
      run_rfi_clean        RFI cleaning only (intake → filter → clean)
      generate_rfi_filter  Generate RFI filter plots only

    DADA Pipelines:
      run_dada_search      Full pipeline starting from DADA files
      run_digifits         Convert DADA to FITS only
      run_dada_clean_stack DADA → FITS → clean → stack

    Utility Workflows:
      fold_par             Fold with a known pulsar .par file
      candypolice          Re-fold candidates from a CandyJar CSV
      setup_basedir        Initialize a new project directory
      validate_inputs      Validate input files and parameters before running
      cleanup_cache        Manage shared cache (find orphaned files)
      help                 Show this help message

    REQUIRED PARAMETERS:
    --------------------
      --basedir            Base directory for all outputs
      --runID              Unique identifier for this search run
      --files_list         Path to input CSV file with data files
      --telescope          Telescope name (effelsberg, meerkat, etc.)

    KEY OPTIONAL PARAMETERS:
    ------------------------
    Search Parameters:
      --ddplan.dm_start    DM range start (relative to coherent DM)
      --ddplan.dm_end      DM range end
      --ddplan.dm_step     DM step size
      --peasoup.segments   Segmentation list, e.g., [1,2,4] for full, half, quarter
      --peasoup.acc_start  Acceleration range start (m/s²)
      --peasoup.acc_end    Acceleration range end (m/s²)
      --peasoup.min_snr    Minimum S/N threshold

    Processing Options:
      --filtool.run_filtool           Enable/disable RFI cleaning (true/false)
      --generateRfiFilter.run_rfi_filter  Generate dynamic RFI masks (true/false)
      --stack_by_cdm                  Stack filterbanks by coherent DM (true/false)
      --split_fil                     Split filterbank by frequency (true/false)

    PROFILES (-profile):
    --------------------
      local      Local execution (for testing)
      hercules   SLURM cluster with GPU nodes
      edgar      Edgar cluster configuration
      contra     Contra cluster configuration
      condor     HTCondor submission

    INPUT FILE FORMAT:
    ------------------
    CSV file with columns:
      pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm

    Example:
      pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm
      0,NGC6544,cfbf00001,1,2024-01-15T10:00:00,18:07:20.5,-24:59:51,/path/to/file.fil,60.0

    OUTPUT STRUCTURE:
    -----------------
      basedir/
      ├── shared_cache/          # Reusable cached files (RFI masks, cleaned fils)
      ├── runID/                  # Run-specific outputs
      │   ├── beam_name/
      │   │   ├── RFIFILTER/      # RFI filter plots (symlinks)
      │   │   ├── CLEANEDFIL/     # Cleaned filterbanks (symlinks)
      │   │   └── segment_N/
      │   │       └── NN/
      │   │           ├── BIRDIES/      # Birdie detection
      │   │           ├── SEARCH/       # Peasoup XML results
      │   │           ├── PARSEXML/     # Parsed candidates + candfiles
      │   │           │   └── XML/      # Picked XML files
      │   │           ├── FOLDING/      # Folded candidate outputs
      │   │           │   ├── PNG/      # Diagnostic plots
      │   │           │   ├── AR/       # Archive files
      │   │           │   ├── CANDS/    # PulsarX .cands files
      │   │           │   ├── CSV/      # Merged candidate CSVs
      │   │           │   └── PROVENANCE/ # Candidate tracking files
      │   │           ├── ABG/          # Alpha-beta-gamma scores
      │   │           ├── ZERODM/       # Zero-DM diagnostic plots
      │   │           └── CLASSIFICATION/ # PICS scores
      │   ├── TARBALL_CSV/        # CSV files used for tarball creation
      │   ├── CANDIDATE_TARBALLS/ # Final candidate tarballs
      │   └── DMFILES/            # DM search files

    EXAMPLES:
    ---------
    # Setup new project
    nextflow run elden.nf -entry setup_basedir --basedir /data/NGC6544_search

    # Run full pipeline with custom DM range
    nextflow run elden.nf -entry full -profile hercules \\
        -c params.config \\
        --runID search_v1 \\
        --ddplan.dm_start -10 \\
        --ddplan.dm_end 10

    # Resume a failed/interrupted run
    nextflow run elden.nf -entry full -profile hercules -c params.config -resume

    # Cleanup orphaned cache files
    nextflow run elden.nf -entry cleanup_cache --basedir /data/NGC6544_search

    TIPS:
    -----
    • Always use -resume to avoid re-running completed tasks
    • Use different runID values for different search configurations
    • The shared_cache directory stores reusable RFI cleaning results
    • Check .nextflow.log for detailed execution logs
    • Use 'nextflow log' to see execution history

    DOCUMENTATION:
    --------------
    GitHub: https://github.com/your-repo/elden-ring
    Issues: https://github.com/your-repo/elden-ring/issues

    """.stripIndent()

    println helpMessage
}

// ============================================================================
// SETUP BASEDIR WORKFLOW - Initialize a new project directory
// ============================================================================
workflow setup_basedir {
    main:
    println """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                      ELDEN-RING Project Setup Wizard                         ║
    ╚══════════════════════════════════════════════════════════════════════════════╝
    """

    // Validate basedir parameter
    if (!params.basedir) {
        error """
        ERROR: --basedir parameter is required!

        Usage: nextflow run elden.nf -entry setup_basedir --basedir /path/to/project

        This will create:
          /path/to/project/
          ├── params.config          # Main configuration file
          ├── inputfile.txt          # Template input file
          ├── meta/                  # Metadata directory
          └── shared_cache/          # Cache directory
        """
    }

    def basedir = file(params.basedir)
    def projectDir_nf = projectDir

    println "Setting up project in: ${basedir}"
    println "Pipeline source: ${projectDir_nf}"
    println ""

    // Execute setup script
    def setupScript = """
    #!/bin/bash
    set -e

    BASEDIR="${basedir}"
    PROJECT_DIR="${projectDir_nf}"

    echo "Creating directory structure..."
    mkdir -p "\${BASEDIR}"
    mkdir -p "\${BASEDIR}/meta"
    mkdir -p "\${BASEDIR}/shared_cache"

    echo ""
    echo "Copying configuration templates..."

    # Copy params.config.example
    if [[ -f "\${PROJECT_DIR}/example/params.config.example" ]]; then
        cp "\${PROJECT_DIR}/example/params.config.example" "\${BASEDIR}/params.config"
        # Update basedir in the config
        sed -i "s|basedir.*=.*|basedir = \\"\${BASEDIR}\\"|" "\${BASEDIR}/params.config"
        echo "  ✓ Created params.config"
    else
        echo "  ✗ Warning: params.config.example not found"
    fi

    # Copy inputfile template
    if [[ -f "\${PROJECT_DIR}/example/inputfile.txt" ]]; then
        cp "\${PROJECT_DIR}/example/inputfile.txt" "\${BASEDIR}/inputfile.txt.template"
        # Create empty inputfile with header only
        head -1 "\${PROJECT_DIR}/example/inputfile.txt" > "\${BASEDIR}/inputfile.txt"
        echo "  ✓ Created inputfile.txt (with header)"
        echo "  ✓ Created inputfile.txt.template (with examples)"
    fi

    # Copy known_pulsars.csv if exists
    if [[ -f "\${PROJECT_DIR}/example/known_pulsars.csv" ]]; then
        cp "\${PROJECT_DIR}/example/known_pulsars.csv" "\${BASEDIR}/known_pulsars.csv"
        echo "  ✓ Copied known_pulsars.csv"
    fi

    # Copy rfi_filters.txt if exists
    if [[ -f "\${PROJECT_DIR}/example/rfi_filters.txt" ]]; then
        cp "\${PROJECT_DIR}/example/rfi_filters.txt" "\${BASEDIR}/rfi_filters.txt"
        echo "  ✓ Copied rfi_filters.txt"
    fi

    # Create a README in the project directory
    cat > "\${BASEDIR}/README.txt" << 'READMEEOF'
ELDEN-RING Project Directory
=============================

This directory was created by the ELDEN-RING setup wizard.

Directory Structure:
--------------------
  params.config         - Main configuration file (EDIT THIS)
  inputfile.txt         - Input data files CSV (EDIT THIS)
  inputfile.txt.template - Example input file format
  meta/                 - Metadata files (auto-generated)
  shared_cache/         - Cached intermediate files (auto-managed)
  <runID>/              - Output directories for each search run

Quick Start:
------------
1. Edit params.config to set your search parameters
2. Edit inputfile.txt to list your input data files
3. Run the pipeline:

   nextflow run /path/to/elden.nf -entry full \\
       -profile hercules \\
       -c params.config \\
       --runID my_search_v1

   Or with resume (recommended):
   nextflow run /path/to/elden.nf -entry full \\
       -profile hercules \\
       -c params.config \\
       --runID my_search_v1 \\
       -resume

Input File Format:
------------------
CSV with columns: pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm

See inputfile.txt.template for examples.

Generating Input Files:
-----------------------
To auto-generate inputfile.txt from your data directory:

   bash generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/data

For DADA inputs (directories; files expanded later):

   bash generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/dada_dir1 /path/to/dada_dir2

Tips:
-----
- Use different runID values for different search configurations
- The shared_cache directory stores reusable preprocessing results
- Always use -resume to avoid re-running completed tasks
- Check .nextflow.log for detailed execution logs

READMEEOF

    echo "  ✓ Created README.txt"

    # Create generate_inputfile.sh script
    cat > "\${BASEDIR}/generate_inputfile.sh" << 'GENEOF'
#!/bin/bash
# Generate inputfile.txt or dada_files.csv from data directories.
#
# Usage:
#   bash generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/data
#   bash generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" /path/to/dada_dir1 /path/to/dada_dir2

set -uo pipefail

usage() {
  cat << 'USAGE'
Usage:
  generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" DATA_DIR [DATA_DIR...]
  generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm "47.0 101.0" DADA_DIR [DADA_DIR...]

Options:
  --dada           Generate dada_files.csv using directory/*dada globs.
  --cluster NAME   Cluster name to write in CSV.
  --ra RA          Right ascension (e.g. 18:24:32.89).
  --dec DEC        Declination (e.g. -24:52:11.4).
  --utc UTC        UTC start time (e.g. 2025-05-13T01:41:02).
  --cdm LIST       Space or comma separated CDM values (e.g. "37.0 101.0").
  --output FILE    Output file name (optional).
  -h, --help       Show this help.
USAGE
}

MODE="fits"
OUTPUT=""
CLUSTER=""
RA=""
DEC=""
UTC=""
CDM_LIST=""
DATA_DIRS=()

while [[ \$# -gt 0 ]]; do
  case "\$1" in
    --dada) MODE="dada"; shift ;;
    --cluster) CLUSTER="\${2:-}"; shift 2 ;;
    --ra) RA="\${2:-}"; shift 2 ;;
    --dec) DEC="\${2:-}"; shift 2 ;;
    --utc) UTC="\${2:-}"; shift 2 ;;
    --cdm) CDM_LIST="\${2:-}"; shift 2 ;;
    --output) OUTPUT="\${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) DATA_DIRS+=("\$1"); shift ;;
  esac
done

if [[ -z "\${CLUSTER}" || -z "\${RA}" || -z "\${DEC}" || -z "\${UTC}" || -z "\${CDM_LIST}" || \${#DATA_DIRS[@]} -eq 0 ]]; then
  echo "ERROR: missing required arguments."
  usage
  exit 1
fi

CDM_LIST_CLEAN="\$(echo "\${CDM_LIST}" | tr ',' ' ' | xargs)"
if [[ -z "\${CDM_LIST_CLEAN}" ]]; then
  echo "ERROR: CDM list is empty."
  exit 1
fi
read -r -a CDMS <<< "\${CDM_LIST_CLEAN}"

if [[ -z "\${OUTPUT}" ]]; then
  if [[ "\${MODE}" == "dada" ]]; then
    OUTPUT="dada_files.csv"
  else
    OUTPUT="inputfile.txt"
  fi
fi

echo "Mode: \${MODE}"
echo "Output: \${OUTPUT}"
echo "Cluster: \${CLUSTER}"
echo "RA/DEC: \${RA} \${DEC}"
echo "UTC: \${UTC}"
echo "CDM list: \${CDM_LIST_CLEAN}"
echo "Data dirs: \${DATA_DIRS[*]}"
echo ""

extract_beam() {
  local text="\$1"
  if [[ "\${text}" =~ baseband([0-9]+) ]]; then
    local id="\${BASH_REMATCH[1]}"
    printf "%s %s\n" "\${id}" "\$(printf "cfbf%05d" "\${id}")"
  elif [[ "\${text}" =~ cfbf([0-9]+) ]]; then
    echo "\${BASH_REMATCH[1]} cfbf\${BASH_REMATCH[1]}"
  elif [[ "\${text}" =~ beam([0-9]+) ]]; then
    echo "\${BASH_REMATCH[1]} beam\${BASH_REMATCH[1]}"
  elif [[ "\${text}" =~ Band([0-9]+) ]]; then
    echo "\${BASH_REMATCH[1]} Band\${BASH_REMATCH[1]}"
  else
    echo "0 beam0"
  fi
}

pointing=0

if [[ "\${MODE}" == "dada" ]]; then
  echo "pointing,dada_files,cluster,beam_name,beam_id,utc_start,ra,dec,cdm_list" > "\${OUTPUT}"
  echo "Processing \${#DATA_DIRS[@]} directories..."
  for dir in "\${DATA_DIRS[@]}"; do
    echo "  Processing: \${dir}"
    if [[ ! -d "\${dir}" ]]; then
      echo "  WARNING: \${dir} is not a directory; skipping"
      continue
    fi
    dir="\${dir%/}"
    read -r beam_id beam_name < <(extract_beam "\${dir}")
    dada_glob="\${dir}/*dada"
    echo "\${pointing},\${dada_glob},\${CLUSTER},\${beam_name},\${beam_id},\${UTC},\${RA},\${DEC},\${CDM_LIST_CLEAN}" >> "\${OUTPUT}"
  done
else
  echo "pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm" > "\${OUTPUT}"
  mapfile -t files < <(find "\${DATA_DIRS[@]}" -type f \\( -name "*.fil" -o -name "*.fits" -o -name "*.sf" -o -name "*.rf" \\) | sort)
  if [[ \${#files[@]} -eq 0 ]]; then
    echo "WARNING: no matching files found."
  fi
  for filepath in "\${files[@]}"; do
    filename="\$(basename "\${filepath}")"
    read -r beam_id beam_name < <(extract_beam "\${filename}")
    for cdm in "\${CDMS[@]}"; do
      echo "\${pointing},\${CLUSTER},\${beam_name},\${beam_id},\${UTC},\${RA},\${DEC},\${filepath},\${cdm}" >> "\${OUTPUT}"
      pointing=\$((pointing + 1))
    done
  done
fi

lines=\$(( \$(wc -l < "\${OUTPUT}") - 1 ))
echo "Generated \${OUTPUT} with \${lines} entries"
GENEOF

    chmod +x "\${BASEDIR}/generate_inputfile.sh"
    echo "  ✓ Created generate_inputfile.sh"

    # Create meta/pipeline_info.txt
    cat > "\${BASEDIR}/meta/pipeline_info.txt" << METAEOF
ELDEN-RING Pipeline Info
========================
Setup Date: \$(date)
Pipeline Version: 1.0.0
Pipeline Path: \${PROJECT_DIR}
Project Directory: \${BASEDIR}
METAEOF
    echo "  ✓ Created meta/pipeline_info.txt"

    echo ""
    echo "════════════════════════════════════════════════════════════════════════════"
    echo "Setup complete!"
    echo "════════════════════════════════════════════════════════════════════════════"
    echo ""
    echo "Next steps:"
    echo ""
    echo "  1. Edit params.config to configure your search parameters:"
    echo "     vim ${basedir}/params.config"
    echo ""
    echo "  2. Generate input file from your data directory:"
    echo "     cd ${basedir}"
    echo "     bash generate_inputfile.sh --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm \"47.0 101.0\" /path/to/data"
    echo "     # DADA example:"
    echo "     bash generate_inputfile.sh --dada --cluster CLUSTER --ra RA --dec DEC --utc UTC --cdm \"47.0 101.0\" /path/to/dada_dir1 /path/to/dada_dir2"
    echo ""
    echo "     Or manually edit inputfile.txt"
    echo ""
    echo "  3. Run the pipeline:"
    echo "     nextflow run ${projectDir_nf}/elden.nf -entry full \\\\"
    echo "         -profile hercules \\\\"
    echo "         -c ${basedir}/params.config \\\\"
    echo "         --runID my_first_search \\\\"
    echo "         -resume"
    echo ""
    echo "  4. For help and available options:"
    echo "     nextflow run ${projectDir_nf}/elden.nf -entry help"
    echo ""
    """

    def proc = ['bash', '-c', setupScript].execute()
    proc.waitFor()
    print proc.text
    if (proc.exitValue() != 0) {
        println "STDERR: ${proc.err.text}"
    }
}

// ============================================================================
// CLEANUP WORKFLOW - Manage shared cache
// ============================================================================
workflow cleanup_cache {
    main:
    def scriptPath = "${projectDir}/scripts/cleanup_shared_cache.sh"

    println "================================================================"
    println "Shared Cache Cleanup Utility"
    println "================================================================"
    println "This workflow scans the shared_cache directory and identifies"
    println "files that are not referenced by any runID-specific directories."
    println ""
    println "Usage:"
    println "  DRY RUN (default - shows what would be deleted):"
    println "    nextflow run elden.nf -entry cleanup_cache"
    println ""
    println "  LIVE RUN (actually deletes orphaned files):"
    println "    bash ${scriptPath} ${params.basedir} false"
    println ""
    println "Running dry-run scan now..."
    println "================================================================"

    // Execute the cleanup script in dry-run mode
    def cleanup_cmd = ["bash", scriptPath, params.basedir, "true"]
    def proc = cleanup_cmd.execute()
    proc.waitFor()
    println proc.text
    if (proc.exitValue() != 0) {
        println "ERROR: ${proc.err.text}"
    }
}

// ============================================================================
// INPUT VALIDATION WORKFLOW - Validate input files and parameters
// ============================================================================
workflow validate_inputs {
    main:
    println """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                      ELDEN-RING Input Validation                             ║
    ╚══════════════════════════════════════════════════════════════════════════════╝
    """

    def errors = []
    def warnings = []

    // Check required parameters
    println "Checking required parameters..."

    if (!params.basedir) {
        errors << "ERROR: --basedir is required but not set"
    } else {
        def basedir = file(params.basedir)
        if (!basedir.exists()) {
            warnings << "WARNING: basedir does not exist: ${params.basedir} (will be created)"
        }
        println "  ✓ basedir: ${params.basedir}"
    }

    if (!params.runID) {
        errors << "ERROR: --runID is required but not set"
    } else {
        println "  ✓ runID: ${params.runID}"
    }

    if (!params.files_list) {
        errors << "ERROR: --files_list is required but not set"
    } else {
        def inputFile = file(params.files_list)
        if (!inputFile.exists()) {
            errors << "ERROR: Input file does not exist: ${params.files_list}"
        } else {
            println "  ✓ files_list: ${params.files_list}"
        }
    }

    if (!params.telescope) {
        errors << "ERROR: --telescope is required but not set"
    } else {
        def validTelescopes = ['effelsberg', 'meerkat', 'parkes', 'gbt', 'chime']
        if (!(params.telescope.toLowerCase() in validTelescopes)) {
            warnings << "WARNING: Unrecognized telescope '${params.telescope}'. Expected one of: ${validTelescopes.join(', ')}"
        }
        println "  ✓ telescope: ${params.telescope}"
    }

    println ""
    println "Validating input CSV file..."

    // Validate input CSV structure
    if (params.files_list && file(params.files_list).exists()) {
        def inputCsv = file(params.files_list)
        def lines = inputCsv.readLines()

        if (lines.size() < 2) {
            errors << "ERROR: Input file has no data rows (only header or empty)"
        } else {
            def header = lines[0].split(',')*.trim()
            def requiredCols = ['pointing', 'cluster', 'beam_name', 'beam_id', 'utc_start', 'ra', 'dec', 'fits_files', 'cdm']
            def missingCols = requiredCols.findAll { !(it in header) }

            if (missingCols) {
                errors << "ERROR: Missing required columns in input CSV: ${missingCols.join(', ')}"
            } else {
                println "  ✓ CSV header has all required columns"
            }

            // Validate data rows
            def validRows = 0
            def invalidRows = []
            def missingFiles = []

            lines[1..-1].eachWithIndex { line, idx ->
                if (line.trim()) {
                    def cols = line.split(',')
                    if (cols.size() < requiredCols.size()) {
                        invalidRows << "Row ${idx + 2}: Not enough columns (${cols.size()} < ${requiredCols.size()})"
                    } else {
                        // Check if fits_file exists
                        def fitsIdx = header.findIndexOf { it == 'fits_files' }
                        if (fitsIdx >= 0 && cols.size() > fitsIdx) {
                            def fitsPath = cols[fitsIdx].trim()
                            def fitsFile = file(fitsPath)
                            if (!fitsFile.exists()) {
                                missingFiles << "Row ${idx + 2}: File not found: ${fitsPath}"
                            }
                        }
                        validRows++
                    }
                }
            }

            println "  ✓ Found ${validRows} data rows"

            if (invalidRows) {
                invalidRows.each { errors << "ERROR: ${it}" }
            }

            if (missingFiles.size() > 0 && missingFiles.size() <= 5) {
                missingFiles.each { warnings << "WARNING: ${it}" }
            } else if (missingFiles.size() > 5) {
                warnings << "WARNING: ${missingFiles.size()} input files not found (first 5 shown)"
                missingFiles[0..4].each { warnings << "  ${it}" }
            }
        }
    }

    println ""
    println "Validating search parameters..."

    // Validate DM plan
    if (params.ddplan) {
        if (params.ddplan.dm_start >= params.ddplan.dm_end) {
            errors << "ERROR: ddplan.dm_start (${params.ddplan.dm_start}) must be less than dm_end (${params.ddplan.dm_end})"
        }
        if (params.ddplan.dm_step <= 0) {
            errors << "ERROR: ddplan.dm_step must be positive"
        }
        println "  ✓ DM plan: ${params.ddplan.dm_start} to ${params.ddplan.dm_end} step ${params.ddplan.dm_step}"
    }

    // Validate peasoup parameters
    if (params.peasoup) {
        if (params.peasoup.min_snr < 5) {
            warnings << "WARNING: peasoup.min_snr (${params.peasoup.min_snr}) is very low, may produce many false positives"
        }
        if (params.peasoup.segments) {
            println "  ✓ Segments: ${params.peasoup.segments}"
        }
        println "  ✓ Acceleration range: ${params.peasoup.acc_start} to ${params.peasoup.acc_end}"
    }

    // Check container images exist (if using containers)
    println ""
    println "Checking container images..."

    def images = [
        'pulsarx_image': params.pulsarx_image,
        'peasoup_image': params.peasoup_image,
        'presto_image': params.presto_image,
        'pics_classifier_image': params.pics_classifier_image
    ]

    images.each { name, path ->
        if (path) {
            println "  ✓ ${name}: ${path}"
        } else {
            warnings << "WARNING: ${name} not configured"
        }
    }

    // Print summary
    println ""
    println "════════════════════════════════════════════════════════════════════════════════"
    println "Validation Summary"
    println "════════════════════════════════════════════════════════════════════════════════"

    if (warnings) {
        println ""
        println "WARNINGS (${warnings.size()}):"
        warnings.each { println "  ⚠ ${it}" }
    }

    if (errors) {
        println ""
        println "ERRORS (${errors.size()}):"
        errors.each { println "  ✗ ${it}" }
        println ""
        println "════════════════════════════════════════════════════════════════════════════════"
        println "Validation FAILED - please fix errors before running the pipeline"
        println "════════════════════════════════════════════════════════════════════════════════"
        error "Input validation failed with ${errors.size()} error(s)"
    } else {
        println ""
        println "════════════════════════════════════════════════════════════════════════════════"
        println "Validation PASSED - inputs are valid"
        println "════════════════════════════════════════════════════════════════════════════════"
        println ""
        println "You can now run the pipeline with:"
        println "  nextflow run elden.nf -entry full -profile <profile> -c <config> -resume"
    }
}
