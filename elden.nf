#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ELDEN-RING: Pulsar Search Pipeline
// ============================================================================
// A Nextflow pipeline for pulsar candidate detection using GPU-accelerated
// periodicity searches with peasoup, pulsarx folding, and ML classification.
// ============================================================================

include { filtool } from './processes'
include { generateDMFiles } from './processes'
include { segmented_params } from './processes'
include { peasoup } from './processes'
include { birdies } from './processes'
include { parse_xml } from './processes'
include { search_fold_merge } from './processes'
include { psrfold } from './processes'
include { pics_classifier } from './processes'
include { readfile } from './processes'
include { generateRfiFilter } from './processes'
include { alpha_beta_gamma_test } from './processes'
include { syncFiles } from './processes'
include { create_candyjar_tarball } from './processes'
include { parfold } from './processes'
include { candypolice_pulsarx} from './processes'
include { extract_candidates } from './processes'
include { dada_to_fits } from './processes'
include { merge_filterbanks } from './processes'
include {split_filterbank} from './processes.nf'

// Default params to avoid warnings when running lightweight entries (e.g., help). !!! DO NOT CHANGE THIS
params.basedir = params.basedir ?: '.'
params.runID = params.runID ?: ''
params.notification = params.notification ?: [
    enabled: false,
    email: '',
    on_complete: false,
    on_fail: false,
    on_error: false
]

workflow intake {
    main:
    // Parse the CSV file to get the list of FITS files and parameters
    fits_file_channel_and_meta = Channel.fromPath("${params.files_list}")
        .splitCsv(header : true, sep : ',')
        .map { row -> 
            def pointing = row.pointing.trim()
            def fits_files = row.fits_files.trim()
            def cluster = row.cluster.trim()
            def beam_name = row.beam_name.trim()
            def beam_id = row.beam_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-")
            def ra = row.ra.trim()
            def dec = row.dec.trim()
            def cdm = row.cdm.trim()
            // extract file name from path
            def filename = fits_files.tokenize("/")[-1]
            return tuple(pointing,fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, cdm, filename)
        }

    // Copy from tape? 
    orig_fits_channel = params.copy_from_tape.run_copy
    ? syncFiles(fits_file_channel_and_meta)
    : fits_file_channel_and_meta
    
    orig_fits_channel

    emit:
    orig_fits_channel
}

workflow dada_intake {
    main:
    // Parse the CSV file to get the list of dada files and parameters
    dada_file_channel_and_meta = Channel.fromPath("${params.dada.dada_csv}")
        .splitCsv(header : true, sep : ',')
        .flatMap { row -> 
        def pointing = row.pointing.trim()
        def source = row.dada_files.trim()
        def cluster = row.cluster.trim()
        def beam_name = row.beam_name.trim()
        def beam_id = row.beam_id.trim()
        def utc_start = row.utc_start.trim().replace(" ", "-")
        def ra = row.ra.trim()
        def dec = row.dec.trim()
        def cdms = row.cdm_list.trim().tokenize().collect { it as Double }

        cdms.collect { cdm ->
            tuple(pointing, source, cluster, beam_name, beam_id, utc_start, ra, dec, cdm)
        }
    }

    dada_file_channel_and_meta

    emit:
    dada_file_channel_and_meta
}

// generate_rfi_filter 
workflow rfi_filter {
    take:
    orig_fits_channel
    
    main:
    readfile(orig_fits_channel).set{ rdout }

    if (params.generateRfiFilter.run_rfi_filter) {
        fil_input = generateRfiFilter(rdout).map { p,f,c,bn,bi,u,ra,dec,cdm,rfi,ts,ns,si,png,txt -> 
            return tuple(p,f,c,bn,bi,u,ra,dec,cdm,rfi,ts,ns,si)   
        }
    } else {
        fil_input = rdout.map { p,f,c,bn,bi,u,ra,dec,cdm,tpf,ts,ns,si -> 
            return tuple(
                p,f,c,bn,bi,u,ra,dec,cdm,
                params.filtool.rfi_filter_list[params.telescope],
                ts,ns,si
            )
        }
    }

    emit: 
    fil_input
}

// rfi_clean: read metadata → (opt) generateRfiFilter → filtool
workflow rfi_clean {
    take:
    fil_input

    main:
    new_fil = params.filtool.run_filtool
    ? filtool(fil_input, params.threads, params.telescope)
    : fil_input.map{ p,f,c,bn,bi,u,ra,dec,cdm,rfi,ts,ns,si -> 
        tuple(p,f,c,bn,bi,u,ra,dec,cdm)
    }

    emit:
    new_fil
}

workflow stack_by_cdm {
  take:
  new_fil

  main:
  new_fil
    .groupTuple(by: [0, 2, 5, 6, 7, 8])  // group by pointing, cluster, utc, ra, dec, cdm
    .map { group ->
      // Unpack group values
      def (p, fil_paths, cluster, beam_names, beam_ids, utc, ra, dec, cdm) = group

      // Create beam_id -> file_path map
      def beamIdToFile = [:]
      beam_ids.eachWithIndex { bid, idx -> 
        beamIdToFile[bid as int] = fil_paths[idx] 
      }

      // Generate stacks dynamically
      def results = []
      params.stacks.each { stackName, beamList ->
        // Get files for beams present in current group
        def stackFiles = beamList.findResults { beamId -> 
          beamIdToFile[beamId]  // returns null for missing beams
        }
        // Only include non-empty stacks
        if (stackFiles) {
          results << tuple(p, cluster, utc, ra, dec, cdm, stackName, stackFiles)
        }
      }
      
      return results
    }
    .flatMap { it }  // flatten list of lists
    .set { stacked_group }
  
  merge_filterbanks(stacked_group)
    .set { stacked_fil }

  emit:
  stacked_fil
}

// segmentation: break into segments → format for peasoup
workflow segmentation {
    take:
    new_fil

    main:
    split_params = new_fil
        .flatMap { p,fp,c,bn,bi,u,ra,dec,cdm -> 
            params.peasoup.segments.collect { segments -> 
                tuple(p,fp,c,bn,bi,u,ra,dec,cdm,segments)
            }
        }
    segmented_params(split_params)
        .set{ seg_ch }

    peasoup_input = seg_ch
        .flatMap { p,fp,c,bn,bi,u,ra,dec,cdm,ts,ns,seg,segf -> 
            segf.splitCsv(header : true, sep : ',').collect { row -> 
                tuple(
                    p,fp,c,bn,bi,u,ra,dec,cdm,ts,
                    row.nsamples_per_segment.trim(),
                    seg, row.i.trim(),
                    row.fft_size.trim(), row.start_sample.trim()
                )
            }
        }

    emit:
    peasoup_input
}

workflow dm {
    take:
    bird_out

    main:
    generateDMFiles(bird_out).set{ dm_file }

    emit:
    dm_file
}

// search: birdies → search → grouping
workflow search {
    take:
    segmentation

    main:
    birdies(segmentation)
        .set{ bird_out }
    
    dm(bird_out)
        .flatMap { p,fi,c,bn,bi,u,ra,dec,cdm,ts,ns,seg,seg_id,fft,start,bird,dml ->
            def dm_files = dml instanceof List ? dml : [dml]
            dm_files.collect { dm_file -> 
              tuple(p, fi, c, bn, bi, u, ra, dec, cdm, ts, ns, seg, seg_id, fft, start, bird, dm_file)
            }
        }
        .set{ peasoup_input }

    peasoup(peasoup_input)
        .map { p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,dm_file,fil_file,xml_path,birds,ss,ns ->
            def fil_base = fil_file.getBaseName()
            tuple(p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,dm_file,fil_base,fil_file,xml_path,ss)
        }
        .groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15])
        .map { p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,dm_file,fil_base,fil_file,xml_paths,start_sample ->
            def combined = []
            def fil_file_sorted = fil_file instanceof List ? fil_file.sort() : [fil_file]
            def first_fil_file = (fil_file_sorted instanceof List) ? fil_file_sorted[0] : fil_file_sorted
            (0..<dm_file.size()).each { i -> 
                combined << [dm: dm_file[i], xml: xml_paths[i]]
            }
            def combined_sorted = combined.sort { a, b -> a.dm <=> b.dm }
            def sorted_dm_file = combined_sorted.collect { it.dm }
            def sorted_xml_file = combined_sorted.collect { it.xml}
            tuple(p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,sorted_dm_file,fil_base,first_fil_file,sorted_xml_file,start_sample)
        }
        .set{ search_out }

    emit:
    search_out
}

// xml_parse: parse XML -> split candidates
workflow xml_parse {
    take:
    search_out

    main:
    parse_xml(search_out)
        .flatMap {  it ->
            def (p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,dm_file,fil_base,fil_file,xml_file,start_sample,filtered_candidate_csv,unfiltered_candidates_csv,candfiles,metafile,allCands) = it
            def cList = candfiles instanceof List ? candfiles : [candfiles]
            cList.collect { candfile ->
                tuple(p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,fil_base,fil_file,start_sample,filtered_candidate_csv,candfile,metafile)
            }
        }
        .set{ splitcands_ch }

    emit:
    splitcands_ch
}

// fold : pulsarx -> grouping
workflow fold { 
    take:
    splitcands_ch

    main:
    psrfold(splitcands_ch)
        .groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11])
        .map { p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,fil_base,fil,f_csv,candfile,metatext,pngs,ar,cands ->
            def f_csv_sort = f_csv instanceof List ? f_csv.sort() : [f_csv]
            def first_f_csv = (f_csv_sort instanceof List) ? f_csv_sort[0] : f_csv_sort
            ar = ar instanceof List ? ar : [ar]
            cands = cands instanceof List ? cands : [cands]
            ar = ar.flatten()
            cands = cands.flatten()
            tuple(p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,fil_base,first_f_csv,ar,cands)
        }
        .set{ fold_grouped }

    emit:
    fold_grouped
}
            
// fold_merge: merge pulsarx output
workflow fold_merge {
    take:
    fold_grouped

    main:
    search_fold_merge(fold_grouped)
        .set{ search_fold_merged }

    emit:
    search_fold_merged
}

// classify: alpha_beta_gamma + pics_classifier -> combine -> tarball
workflow classify {
    take:
    search_fold_merged

    main:
    def abg = alpha_beta_gamma_test(search_fold_merged)
    def pics = pics_classifier(search_fold_merged)

    abg.join(pics, by: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11])
        .map { row -> row.join(',') } // Convert each list to a CSV line
        .collectFile(
            name: 'alpha_beta_pics_combined.csv', 
            newLine: true, 
            storeDir: params.basedir
        )
        .set{ abg_pics_combined_csv }

    emit:
    abg_pics_combined_csv
}

workflow candyjar_tarball {
    take:
    abg_pics_combined_csv

    main:
    if (params.alpha_beta_gamma.create_candyjar_tarball) {
        def tar_input = abg_pics_combined_csv.map { f ->
            tuple(f, "${params.alpha_beta_gamma.output_dir_alpha_pics_results}.tar.gz") }

        create_candyjar_tarball(tar_input)
            .set { tarball_output }

        // Log the tarball CSV location
        tarball_output.subscribe { header_csv, candidates_csv, alpha_csv, pics_csv ->
            println "Tarball CSV files created in TARBALL_CSV directory:"
            println "  Header CSV: ${header_csv}"
            println "  Candidates CSV: ${candidates_csv}"
            println "  Alpha below one CSV: ${alpha_csv}"
            println "  PICS above threshold CSV: ${pics_csv}"
        }
    }
}

// ----------------Main workflow ------------------
workflow full {
    main:
    def intake_ch    = intake()
    def rfi_ch       = rfi_filter(intake_ch)
    def cleaned_ch   = rfi_clean(rfi_ch)
    def cut_ch
    if (params.split_fil) {
        cut_ch    = split_filterbank(cleaned_ch)
    } else {
        cut_ch       = cleaned_ch
        }
    def seg_ch
    if (params.stack_by_cdm) {
        def stacked_ch = stack_by_cdm(cut_ch)
        seg_ch         = segmentation(stacked_ch)
    } else {
        seg_ch         = segmentation(cut_ch)
    }

    def search_ch    = search(seg_ch)
    def xml_ch       = xml_parse(search_ch)
    def fold_ch      = fold(xml_ch)
    def merged_ch    = fold_merge(fold_ch)
    def classify_ch  = classify(merged_ch)
    candyjar_tarball(classify_ch)
}

// ----------------DADA-SF-FIL-STACK-SEARCH-FOLD------------ 

workflow run_dada_search {
    main:
    dada_intake()
    def sf = dada_to_fits(dada_intake.out)
    def rfi_ch = rfi_filter(sf)
    def cleaned_ch   = rfi_clean(rfi_ch)

    def cut_ch
    if (params.split_fil) {
        cut_ch    = split_filterbank(cleaned_ch)
    } else {
        cut_ch       = cleaned_ch
        }

    def seg_ch
    if (params.stack_by_cdm) {
        def stacked_ch = stack_by_cdm(cut_ch)
        seg_ch         = segmentation(cut_ch)
    } else {
        seg_ch         = segmentation(cut_ch)
    }
    def search_ch    = search(seg_ch)
    def xml_ch       = xml_parse(search_ch)
    def fold_ch      = fold(xml_ch)
    def merged_ch    = fold_merge(fold_ch)
    def classify_ch  = classify(merged_ch)
    candyjar_tarball(classify_ch)
  }

// -------------DADA TO FITS conversion ----------
workflow run_digifits {
    main:
    dada_intake()
    dada_to_fits(dada_intake.out)
        .set{ digifits_out }
}

// -------------DADA->SF->FIL->STACK--------------
workflow run_dada_clean_stack {
    main:
    dada_intake()
    def sf = dada_to_fits(dada_intake.out)
    def rfi_ch = rfi_filter(sf)
    def cleaned_ch   = rfi_clean(rfi_ch)
    def cut_ch
    if (params.split_fil) {
        cut_ch    = split_filterbank(cleaned_ch)
    } else {
        cut_ch       = cleaned_ch
    }

    if (params.stack_by_cdm) {
        def stack = stack_by_cdm(cut_ch)
      }
}
 
// -------------Generate the rfi plots ----------
workflow generate_rfi_filter {
    intake()
    rfi_filter(intake.out)
}

//-------------Filtool the files ----------------
workflow run_rfi_clean {
    intake()
    rfi_filter(intake.out)
    rfi_clean(rfi_filter.out)
}

// ---------- Run search and fold on filtooled files -----
// run_search assumes rfi_cleaned files inside the files_lists
workflow run_search_fold {
    def cleaned_ch   = intake().map{ p,f,c,bn,bi,u,ra,dec,cdm,fname -> 
        tuple(p,f,c,bn,bi,u,ra,dec,cdm)
    }
    def cut_ch
    if (params.split_fil) {
        cut_ch    = split_filterbank(cleaned_ch)
    } else {
        cut_ch       = cleaned_ch
        }
    def seg_ch
    if (params.stack_by_cdm) {
        def stacked_ch = stack_by_cdm(cut_ch)
        seg_ch         = segmentation(stacked_ch)
    } else {
        seg_ch         = segmentation(cut_ch)
    }

    def search_ch    = search(seg_ch)
    def xml_ch       = xml_parse(search_ch)
    def fold_ch      = fold(xml_ch)
    def merged_ch    = fold_merge(fold_ch)
    def classify_ch  = classify(merged_ch)
    candyjar_tarball(classify_ch)
}

//------------- parfold workflow ---------------
workflow fold_par {
    parfile_ch = Channel.fromPath("${params.parfold.parfile_path}")
    intake()
    rfi_filter(intake.out)
    rfi_clean(rfi_filter.out)
    parfold(rfi_clean.out, parfile_ch)
        .set{ parfold_out }
    
    emit:
        parfold_out
}

workflow candypolice {
    intake()
    readfile(intake.out).map{ p,f,c,bn,bi,u,ra,dec,tpf,ts,ns,si ->
        tuple(p,f,c,bn,bi,u,ra,dec,ts,ns,si)
    }.set{rdout}
    candyjar_csv = Channel.fromPath("${params.candypolice.input_csv}")
    candfiles = extract_candidates(candyjar_csv).flatMap{ cd, cf ->
        def cList = cf instanceof List ? cf : [cf]
        return tuple(cList)
    }
    
    candypolice_pulsarx(rdout, candfiles)
}

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

   bash generate_inputfile.sh /path/to/data

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
# Generate inputfile.txt from data directory
# Usage: bash generate_inputfile.sh /path/to/data [cluster_name] [cdm_value]

DATA_DIR="\${1:-.}"
CLUSTER="\${2:-UNKNOWN}"
CDM="\${3:-0.0}"
OUTPUT="inputfile.txt"

echo "Scanning for data files in: \${DATA_DIR}"
echo "Cluster name: \${CLUSTER}"
echo "Default CDM: \${CDM}"
echo ""

# Write header
echo "pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm" > "\${OUTPUT}"

# Counter for pointing IDs
pointing=0

# Find all supported file types
find "\${DATA_DIR}" -type f \\( -name "*.fil" -o -name "*.fits" -o -name "*.sf" -o -name "*.rf" \\) | sort | while read filepath; do
    filename=\$(basename "\${filepath}")

    # Try to extract beam info from filename (common patterns)
    # Pattern 1: *_cfbfNNNNN_* or *_cfbfN_*
    if [[ "\${filename}" =~ cfbf([0-9]+) ]]; then
        beam_id="\${BASH_REMATCH[1]}"
        beam_name="cfbf\${beam_id}"
    # Pattern 2: *_beamNN_*
    elif [[ "\${filename}" =~ beam([0-9]+) ]]; then
        beam_id="\${BASH_REMATCH[1]}"
        beam_name="beam\${beam_id}"
    # Pattern 3: *_BandN_*
    elif [[ "\${filename}" =~ Band([0-9]+) ]]; then
        beam_id="\${BASH_REMATCH[1]}"
        beam_name="Band\${beam_id}"
    else
        beam_id="0"
        beam_name="beam0"
    fi

    # Try to extract UTC from filename (pattern: YYYY-MM-DDTHH:MM:SS)
    if [[ "\${filename}" =~ ([0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}) ]]; then
        utc_start="\${BASH_REMATCH[1]}"
    elif [[ "\${filename}" =~ ([0-9]{4}-[0-9]{2}-[0-9]{2}) ]]; then
        utc_start="\${BASH_REMATCH[1]}T00:00:00"
    else
        utc_start="1970-01-01T00:00:00"
    fi

    # Try to extract CDM from filename (pattern: cdm_NN.N or cdm_NN)
    if [[ "\${filename}" =~ cdm_([0-9.]+) ]]; then
        file_cdm="\${BASH_REMATCH[1]}"
    else
        file_cdm="\${CDM}"
    fi

    # Default RA/DEC (should be updated manually or from headers)
    ra="00:00:00.0"
    dec="+00:00:00.0"

    # Write entry
    echo "\${pointing},\${CLUSTER},\${beam_name},\${beam_id},\${utc_start},\${ra},\${dec},\${filepath},\${file_cdm}" >> "\${OUTPUT}"

    ((pointing++))
done

echo ""
echo "Generated \${OUTPUT} with \$(( \$(wc -l < \${OUTPUT}) - 1 )) entries"
echo ""
echo "IMPORTANT: Please verify and update the following in \${OUTPUT}:"
echo "  - RA and DEC coordinates (currently set to defaults)"
echo "  - CDM values if not extracted from filenames"
echo "  - Cluster name if different from '\${CLUSTER}'"
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
    echo "     bash generate_inputfile.sh /path/to/your/data CLUSTER_NAME CDM_VALUE"
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
// Cleanup workflow for shared_cache management
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

// ============================================================================
// WORKFLOW COMPLETION HANDLERS - Monitoring and Notifications
// ============================================================================

// Helper function to format duration in human-readable format
def formatDuration(long millis) {
    def totalSeconds = (millis / 1000L) as long
    def seconds = totalSeconds % 60
    def minutes = Math.floorDiv(totalSeconds, 60L) % 60
    def hours = Math.floorDiv(totalSeconds, 3600L) % 24
    def days = Math.floorDiv(totalSeconds, 86400L)

    if (days > 0) {
        return String.format("%dd %dh %dm %ds", days, hours, minutes, seconds)
    } else if (hours > 0) {
        return String.format("%dh %dm %ds", hours, minutes, seconds)
    } else if (minutes > 0) {
        return String.format("%dm %ds", minutes, seconds)
    } else {
        return String.format("%ds", seconds)
    }
}

workflow.onComplete {
    def duration = workflow.duration
    def success = workflow.success
    def exitStatus = workflow.exitStatus

    // ========================================================================
    // Cumulative runtime tracking across resume runs
    // ========================================================================
    def runtimeFile = file("${params.basedir ?: '.'}/.cumulative_runtime_${workflow.sessionId}.txt")
    def cumulativeMillis = 0L
    def runCount = 1

    // Read previous cumulative time if exists
    try {
        if (runtimeFile.exists()) {
            def lines = runtimeFile.readLines()
            if (lines.size() >= 2) {
                cumulativeMillis = lines[0].toLong()
                runCount = lines[1].toInteger() + 1
            }
        }
    } catch (Exception e) {
        println "WARNING: Could not read cumulative runtime file: ${e.message}"
    }

    // Add current run duration
    def currentMillis = duration.toMillis()
    def totalMillis = cumulativeMillis + currentMillis

    // Save updated cumulative time
    try {
        runtimeFile.text = "${totalMillis}\n${runCount}\n"
    } catch (Exception e) {
        println "WARNING: Could not write cumulative runtime file: ${e.message}"
    }

    def cumulativeDuration = formatDuration(totalMillis)
    def currentDuration = formatDuration(currentMillis)

    // Generate completion summary
    def summary = """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                      ELDEN-RING Pipeline Complete                            ║
    ╚══════════════════════════════════════════════════════════════════════════════╝

    Pipeline:      ${workflow.manifest.name ?: 'ELDEN-RING'}
    Version:       ${workflow.manifest.version ?: '1.0.0'}
    Run Name:      ${workflow.runName}
    Session ID:    ${workflow.sessionId}

    Status:        ${success ? 'SUCCESS' : 'FAILED'}
    Exit Code:     ${exitStatus}

    TIMING:
    -------
    This run:              ${currentDuration}
    Total (across resumes): ${cumulativeDuration}
    Number of runs:        ${runCount}

    Completed:     ${workflow.complete}

    Work Dir:      ${workflow.workDir}
    Output Dir:    ${params.basedir ?: 'N/A'}
    Run ID:        ${params.runID ?: 'N/A'}

    Command:       ${workflow.commandLine}
    Profile:       ${workflow.profile}
    Container:     ${workflow.container ?: 'N/A'}
    """.stripIndent()

    println summary

    // Write summary to file
    def summaryFile = file("${params.basedir ?: '.'}/pipeline_summary_${workflow.runName}.txt")
    try {
        summaryFile.text = summary
        println "Pipeline summary written to: ${summaryFile}"
    } catch (Exception e) {
        println "WARNING: Could not write summary file: ${e.message}"
    }

    // Send email notification if enabled
    if (params.notification?.enabled && params.notification?.email) {
        def subject = "[ELDEN-RING] Pipeline ${success ? 'COMPLETED' : 'FAILED'}: ${params.runID ?: workflow.runName}"
        def shouldSend = (success && params.notification?.on_complete) ||
                         (!success && params.notification?.on_fail)

        if (shouldSend) {
            try {
                sendMail(
                    to: params.notification.email,
                    subject: subject,
                    body: summary
                )
                println "Email notification sent to: ${params.notification.email}"
            } catch (Exception e) {
                println "WARNING: Could not send email notification: ${e.message}"
                println "Make sure your system has sendmail configured or SMTP settings are correct."
            }
        }
    }
}

workflow.onError {
    def errorMessage = """
    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║                      ELDEN-RING Pipeline ERROR                               ║
    ╚══════════════════════════════════════════════════════════════════════════════╝

    Error Message: ${workflow.errorMessage}
    Error Report:  ${workflow.errorReport ?: 'N/A'}

    Pipeline:      ${workflow.manifest.name ?: 'ELDEN-RING'}
    Run Name:      ${workflow.runName}
    Work Dir:      ${workflow.workDir}

    To debug, check the .nextflow.log file and the failed task's work directory.
    """.stripIndent()

    println errorMessage

    // Send error notification if enabled
    if (params.notification?.enabled && params.notification?.email && params.notification?.on_error) {
        def subject = "[ELDEN-RING] Pipeline ERROR: ${params.runID ?: workflow.runName}"

        try {
            sendMail(
                to: params.notification.email,
                subject: subject,
                body: errorMessage
            )
            println "Error notification sent to: ${params.notification.email}"
        } catch (Exception e) {
            println "WARNING: Could not send error notification: ${e.message}"
        }
    }
}

// Default to `full` if no --entry is given
workflow {
    full()
}
