#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ELDEN-RING: Pulsar Search Pipeline
// ============================================================================
// A Nextflow pipeline for pulsar candidate detection using GPU-accelerated
// periodicity searches with peasoup, pulsarx folding, and ML classification.
// ============================================================================

// Process includes
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
include { split_filterbank } from './processes.nf'

// PRESTO pipeline processes
include { presto_rfifind } from './processes'
include { presto_prepdata_zerodm } from './processes'
include { presto_accelsearch_zerodm } from './processes'
include { presto_prepsubband } from './processes'
include { presto_accelsearch } from './processes'
include { presto_sift_candidates } from './processes'
include { presto_prepfold_batch } from './processes'
include { presto_pfd_to_png } from './processes'
include { presto_fold_merge } from './processes'
include { presto_create_tarball } from './processes'
include { presto_fold_pulsarx } from './processes'

// Peasoup time series dumping for PRESTO processing
include { peasoup_dump_timeseries } from './processes'
include { generate_fold_meta } from './processes'

// Utility workflow includes
include { help } from './utilities'
include { setup_basedir } from './utilities'
include { cleanup_cache } from './utilities'
include { validate_inputs } from './utilities'

// ============================================================================
// Default parameters to avoid warnings when running lightweight entries (e.g., help)
// These are overridden when a params.config file is provided
// !!! DO NOT CHANGE THIS SECTION !!!
// ============================================================================
params.basedir = null
params.runID = null
params.notification = [
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

    peasoup(peasoup_input).search_results
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
// PRESTO PIPELINE WORKFLOWS
// ============================================================================

/*
 * PRESTO RFI Detection workflow
 */
workflow presto_rfi {
    take:
    input_files

    main:
    presto_rfifind(input_files)
        .set { rfi_out }

    emit:
    rfi_mask = rfi_out.rfi_mask
    rfi_files = rfi_out.rfi_files
    rfi_stats = rfi_out.rfi_stats
}

/*
 * PRESTO Birdie Detection workflow (zero-DM search for RFI identification)
 */
workflow presto_birdies {
    take:
    input_file
    rfi_mask

    main:
    presto_prepdata_zerodm(input_file, rfi_mask)
        .set { zerodm_out }

    presto_accelsearch_zerodm(zerodm_out.dat_file, zerodm_out.inf_file)
        .set { birdie_out }

    emit:
    zaplist = birdie_out.zaplist
}

/*
 * PRESTO Dedispersion workflow
 */
workflow presto_dedisperse {
    take:
    input_file
    rfi_mask
    dm_ranges  // Channel of tuples: (dm_low, dm_high, dm_step, downsamp)

    main:
    presto_prepsubband(input_file, rfi_mask, dm_ranges)
        .set { subband_out }

    emit:
    subband_data = subband_out.subband_data
}

/*
 * PRESTO Acceleration Search workflow
 */
workflow presto_search {
    take:
    dat_files
    inf_files
    zaplist

    main:
    // Handle case where zaplist may not exist
    zaplist_ch = zaplist.ifEmpty(file('NO_ZAPLIST'))

    presto_accelsearch(dat_files, inf_files, zaplist_ch)
        .set { accel_out }

    emit:
    accel_files = accel_out.accel_files
    cand_files = accel_out.cand_files
}

/*
 * PRESTO Sift and Fold workflow
 */
workflow presto_sift_fold {
    take:
    accel_files
    input_file
    rfi_mask

    main:
    presto_sift_candidates(accel_files, input_file)
        .set { sifted_out }

    presto_prepfold_batch(input_file, sifted_out.sifted_csv, rfi_mask)
        .set { fold_out }

    emit:
    pfd_files = fold_out.pfd_files
    sifted_csv = sifted_out.sifted_csv
}

/*
 * PRESTO Post-processing workflow (PNG conversion and tarball creation)
 */
workflow presto_postprocess {
    take:
    pfd_files
    sifted_csv
    input_file

    main:
    presto_pfd_to_png(pfd_files)
        .set { png_out }

    // Collect bestprof files (may be empty)
    bestprof_ch = pfd_files.map { pfd ->
        file("${pfd}.bestprof")
    }.filter { it.exists() }.collect().ifEmpty([])

    presto_fold_merge(
        pfd_files.collect(),
        png_out.png_files.collect(),
        bestprof_ch,
        sifted_csv
    ).set { merged_out }

    presto_create_tarball(
        merged_out.merged_csv,
        merged_out.all_pfd,
        merged_out.all_png,
        input_file
    ).set { tarball_out }

    emit:
    tarball = tarball_out.tarball
    final_csv = tarball_out.final_csv
    png_files = png_out.png_files
}

/*
 * Full PRESTO search pipeline
 * Runs: RFI detection -> Birdie detection -> Dedispersion -> Acceleration search -> Sifting -> Folding -> Post-processing
 *
 * Set params.presto.fold_backend = 'presto' (default) to fold with prepfold
 * Set params.presto.fold_backend = 'pulsarx' to fold with psrfold_fil2
 */
workflow presto_full {
    main:
    // Intake files
    def intake_ch = intake()

    // Extract filterbank files and metadata
    def fil_channel = intake_ch.map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        file(f)
    }

    // Also keep the full metadata for PulsarX folding
    def meta_channel = intake_ch.map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        tuple(p, f, c, bn, bi, u, ra, dec, cdm)
    }

    // Run RFI detection
    presto_rfi(fil_channel)
        .set { rfi_out }

    // Run birdie detection
    presto_birdies(fil_channel, rfi_out.rfi_mask)
        .set { birdie_out }

    // Generate DM ranges from params
    def dm_ranges = Channel.from(params.presto.dm_ranges).map { range ->
        tuple(range.dm_low, range.dm_high, range.dm_step, range.downsamp)
    }

    // Combine input file with rfi_mask for dedispersion
    def dedisperse_input = fil_channel.combine(rfi_out.rfi_mask)

    // Run dedispersion for each DM range
    presto_prepsubband(
        dedisperse_input.map { it[0] },
        dedisperse_input.map { it[1] },
        dm_ranges
    ).set { subband_out }

    // Run acceleration search on all dedispersed data
    presto_accelsearch(
        subband_out.dat_files.flatten(),
        subband_out.inf_files.flatten(),
        birdie_out.zaplist.ifEmpty(file('NO_ZAPLIST'))
    ).set { accel_out }

    // Collect all ACCEL files and sift
    presto_sift_candidates(
        accel_out.accel_files.collect(),
        fil_channel
    ).set { sifted_out }

    // Choose folding backend based on params.presto.fold_backend
    def fold_backend = params.presto?.fold_backend ?: 'presto'

    if (fold_backend == 'pulsarx') {
        // Fold with PulsarX (psrfold_fil2) using pulsarx_fold.py --csv
        // Generate meta file with observation parameters
        def fft_size = params.presto?.fft_size ?: 134217728
        generate_fold_meta(meta_channel, fft_size, fold_backend)
            .set { meta_out }

        // Fold with PulsarX
        presto_fold_pulsarx(
            fil_channel,
            sifted_out.sifted_csv,
            meta_out.meta_file
        ).set { fold_out }

        // Note: PulsarX folding produces .ar and .png files directly
        // No additional tarball creation needed for PulsarX output
    } else {
        // Default: Fold with PRESTO prepfold
        // prepfold doesn't need a meta file - it uses command-line parameters
        // The period correction is handled by the presto_prepfold_batch process
        presto_prepfold_batch(
            fil_channel,
            sifted_out.sifted_csv,
            rfi_out.rfi_mask
        ).set { fold_out }

        // Post-processing: PNG conversion and tarball
        presto_pfd_to_png(fold_out.pfd_files.collect())
            .set { png_out }

        // Create merged results and tarball
        presto_fold_merge(
            fold_out.pfd_files.collect(),
            png_out.png_files.collect(),
            fold_out.bestprof_files.collect().ifEmpty([]),
            sifted_out.sifted_csv
        ).set { merged_out }

        presto_create_tarball(
            merged_out.merged_csv,
            merged_out.all_pfd,
            merged_out.all_png,
            fil_channel
        ).set { tarball_out }
    }
}

/*
 * PRESTO search on pre-cleaned filterbanks
 * Similar to run_search_fold but using PRESTO instead of peasoup+pulsarx
 */
workflow run_presto_search {
    main:
    // Intake pre-cleaned files
    def cleaned_ch = intake().map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        file(f)
    }

    // Run RFI detection
    presto_rfi(cleaned_ch)
        .set { rfi_out }

    // Run birdie detection
    presto_birdies(cleaned_ch, rfi_out.rfi_mask)
        .set { birdie_out }

    // Generate DM ranges
    def dm_ranges = Channel.from(params.presto.dm_ranges).map { range ->
        tuple(range.dm_low, range.dm_high, range.dm_step, range.downsamp)
    }

    // Run dedispersion
    presto_prepsubband(cleaned_ch, rfi_out.rfi_mask, dm_ranges)
        .set { subband_out }

    // Run acceleration search
    presto_accelsearch(
        subband_out.dat_files.flatten(),
        subband_out.inf_files.flatten(),
        birdie_out.zaplist.ifEmpty(file('NO_ZAPLIST'))
    ).set { accel_out }

    // Sift candidates
    presto_sift_candidates(accel_out.accel_files.collect(), cleaned_ch)
        .set { sifted_out }

    // Fold candidates
    presto_prepfold_batch(cleaned_ch, sifted_out.sifted_csv, rfi_out.rfi_mask)
        .set { fold_out }

    // Convert to PNG
    presto_pfd_to_png(fold_out.pfd_files.collect())
        .set { png_out }

    // Merge and create tarball
    presto_fold_merge(
        fold_out.pfd_files.collect(),
        png_out.png_files.collect(),
        fold_out.bestprof_files.collect().ifEmpty([]),
        sifted_out.sifted_csv
    ).set { merged_out }

    presto_create_tarball(
        merged_out.merged_csv,
        merged_out.all_pfd,
        merged_out.all_png,
        cleaned_ch
    ).set { tarball_out }
}

/*
 * Peasoup time series dumping workflow
 * Dumps PRESTO-compatible .dat/.inf files from peasoup for subsequent PRESTO processing
 *
 * Enable with: params.peasoup.dump_timeseries = true
 */
workflow peasoup_timeseries_dump {
    take:
    input_channel  // tuple(pointing, fil_file, cluster, beam_name, beam_id, utc, ra, dec, cdm)
    dm_file
    birdies_file

    main:
    def fil_channel = input_channel.map { p, f, c, bn, bi, u, ra, dec, cdm ->
        file(f)
    }

    peasoup_dump_timeseries(fil_channel, dm_file, birdies_file)
        .set { timeseries_out }

    emit:
    dat_files = timeseries_out.dat_files
    inf_files = timeseries_out.inf_files
    timeseries_data = timeseries_out.timeseries_data
}

/*
 * PRESTO search on peasoup-dumped time series
 *
 * This workflow:
 * 1. Takes peasoup-dumped .dat/.inf files
 * 2. Runs PRESTO accelsearch on them
 * 3. Sifts candidates
 * 4. Folds with either prepfold or PulsarX
 *
 * Usage: nextflow run elden.nf -entry presto_on_peasoup_timeseries
 */
workflow presto_on_peasoup_timeseries {
    take:
    dat_files   // Channel of .dat files from peasoup
    inf_files   // Channel of .inf files from peasoup
    fil_channel // Original filterbank files for folding
    meta_channel // Metadata channel for PulsarX folding

    main:
    // Run acceleration search on peasoup time series
    presto_accelsearch(
        dat_files.flatten(),
        inf_files.flatten(),
        file('NO_ZAPLIST')  // No zaplist for peasoup time series
    ).set { accel_out }

    // Sift candidates
    presto_sift_candidates(accel_out.accel_files.collect(), fil_channel)
        .set { sifted_out }

    // Choose folding backend
    def fold_backend = params.presto?.fold_backend ?: 'pulsarx'

    if (fold_backend == 'pulsarx') {
        // Generate meta file for PulsarX folding
        def fft_size = params.peasoup?.fft_size ?: 134217728
        generate_fold_meta(meta_channel, fft_size, fold_backend)
            .set { meta_out }

        // Fold with PulsarX
        presto_fold_pulsarx(fil_channel, sifted_out.sifted_csv, meta_out.meta_file)
            .set { fold_out }
    } else {
        // Fold with prepfold (no meta file needed)
        presto_prepfold_batch(fil_channel, sifted_out.sifted_csv, file('NO_MASK'))
            .set { fold_out }

        // Post-processing
        presto_pfd_to_png(fold_out.pfd_files.collect())
            .set { png_out }

        presto_fold_merge(
            fold_out.pfd_files.collect(),
            png_out.png_files.collect(),
            fold_out.bestprof_files.collect().ifEmpty([]),
            sifted_out.sifted_csv
        ).set { merged_out }

        presto_create_tarball(
            merged_out.merged_csv,
            merged_out.all_pfd,
            merged_out.all_png,
            fil_channel
        ).set { tarball_out }
    }

    emit:
    sifted_csv = sifted_out.sifted_csv
}

/*
 * Hybrid peasoup+PRESTO workflow
 *
 * 1. Run peasoup search (standard)
 * 2. Optionally dump time series from peasoup
 * 3. Run PRESTO accelsearch on dumped time series
 * 4. Combine candidates from both searches
 *
 * Enable with: params.peasoup.dump_timeseries = true
 */
workflow peasoup_with_presto_search {
    main:
    // Standard intake and cleaning
    def intake_ch = intake()
    def rfi_ch = rfi_filter(intake_ch)
    def cleaned_ch = rfi_clean(rfi_ch)

    // Run standard peasoup search
    def seg_ch = segmentation(cleaned_ch)
    def search_ch = search(seg_ch)

    // Check if time series dumping is enabled
    if (params.peasoup?.dump_timeseries) {
        // Prepare for time series dumping
        def dm_file = generateDMFiles(cleaned_ch.first()).dm_file
        def birdies_file = birdies(cleaned_ch.first()).birdies_file

        // Dump time series
        peasoup_timeseries_dump(cleaned_ch, dm_file, birdies_file)
            .set { timeseries_out }

        // Run PRESTO search on dumped time series
        def fil_ch = cleaned_ch.map { p, f, c, bn, bi, u, ra, dec, cdm ->
            file(f)
        }

        presto_on_peasoup_timeseries(
            timeseries_out.dat_files,
            timeseries_out.inf_files,
            fil_ch,
            cleaned_ch
        ).set { presto_out }
    }

    // Continue with standard peasoup pipeline
    def xml_ch = xml_parse(search_ch)
    def fold_ch = fold(xml_ch)
    def merged_ch = fold_merge(fold_ch)
    def classify_ch = classify(merged_ch)
    candyjar_tarball(classify_ch)
}

/*
 * Main workflow with search backend selection
 * Use params.search_backend = 'presto' or 'peasoup' to select pipeline
 */
workflow search_pipeline {
    main:
    def backend = params.search_backend ?: 'peasoup'

    if (backend == 'presto') {
        presto_full()
    } else {
        full()
    }
}

/*
 * Entry point for running accelsearch on pre-dumped time series from peasoup
 *
 * This workflow is for the second run scenario:
 * 1. First run: normal peasoup search with dump_timeseries=true → tarballs + .dat/.inf files
 * 2. Second run: use this entry point to run accelsearch on the dumped .dat/.inf files
 *
 * Two modes of operation:
 *
 * MODE 1: Multi-file mode (uses same input CSV as first run)
 *   Uses the same params.files_list CSV to find filterbanks and auto-discovers
 *   time series directories based on the output structure from run 1.
 *
 *   Required parameters:
 *     params.files_list           - Same input CSV file used in first run
 *     params.basedir              - Same basedir used in first run
 *     params.runID                - Same runID used in first run
 *
 *   Usage:
 *     nextflow run elden.nf -entry run_accelsearch_on_timeseries -params-file params.config
 *
 * MODE 2: Single-file mode (explicit paths)
 *   For processing a single filterbank with its time series.
 *
 *   Required parameters:
 *     params.timeseries_input_dir  - Directory containing pre-dumped .dat and .inf files
 *     params.filterbank_file       - Original filterbank file for folding
 *
 *   Usage:
 *     nextflow run elden.nf -entry run_accelsearch_on_timeseries \
 *         --timeseries_input_dir /path/to/TIMESERIES \
 *         --filterbank_file /path/to/file.fil
 *
 * Common parameters:
 *   params.presto.fold_backend   - 'pulsarx' or 'presto' (default: 'pulsarx')
 *   params.presto.fft_size       - FFT size used in original search (for period correction)
 */
workflow run_accelsearch_on_timeseries {
    main:
    // Determine which mode to use
    def single_file_mode = params.timeseries_input_dir && params.filterbank_file
    def multi_file_mode = params.files_list && params.basedir && params.runID

    if (!single_file_mode && !multi_file_mode) {
        error """ERROR: Missing required parameters.

For multi-file mode (recommended), provide:
  - params.files_list  (same input CSV as first run)
  - params.basedir     (same basedir as first run)
  - params.runID       (same runID as first run)

For single-file mode, provide:
  - params.timeseries_input_dir  (directory with .dat/.inf files)
  - params.filterbank_file       (original filterbank file)
"""
    }

    if (single_file_mode) {
        // ===== SINGLE FILE MODE =====
        def dat_channel = Channel.fromPath("${params.timeseries_input_dir}/*.dat")
        def inf_channel = Channel.fromPath("${params.timeseries_input_dir}/*.inf")
        def fil_channel = Channel.fromPath(params.filterbank_file)

        def meta_channel = fil_channel.map { fil ->
            def basename = fil.baseName
            tuple(
                basename,                    // pointing
                fil,                         // filterbank path
                params.cluster ?: 'unknown', // cluster
                basename,                    // beam_name
                '0',                         // beam_id
                'unknown',                   // utc
                '00:00:00.0',               // ra
                '+00:00:00.0',              // dec
                params.cdm ?: '0.0'         // cdm
            )
        }

        run_accelsearch_single(dat_channel, inf_channel, fil_channel, meta_channel)

    } else {
        // ===== MULTI-FILE MODE =====
        // Parse input CSV (same format as first run)
        def input_channel = Channel.fromPath("${params.files_list}")
            .splitCsv(header: true, sep: ',')
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
                tuple(pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, cdm)
            }

        // For each input file, find its time series directory
        // Structure from run 1: ${basedir}/${runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/TIMESERIES/
        // We search for all segment directories with TIMESERIES subfolder
        def timeseries_channel = input_channel.flatMap { pointing, fits_file, cluster, beam_name, beam_id, utc, ra, dec, cdm ->
            // Find all TIMESERIES directories for this beam
            def base_path = "${params.basedir}/${params.runID}/${beam_name}"
            def timeseries_dirs = []

            // Search for TIMESERIES directories in all segment subdirectories
            def base_dir = file(base_path)
            if (base_dir.exists()) {
                base_dir.eachDirRecurse { dir ->
                    if (dir.name == 'TIMESERIES' && dir.isDirectory()) {
                        def dat_files = file("${dir}/*.dat")
                        def inf_files = file("${dir}/*.inf")
                        if (dat_files || inf_files) {
                            timeseries_dirs << tuple(
                                pointing, fits_file, cluster, beam_name, beam_id, utc, ra, dec, cdm,
                                dir.toString()
                            )
                        }
                    }
                }
            }

            if (timeseries_dirs.isEmpty()) {
                log.warn "No TIMESERIES directories found for beam ${beam_name} in ${base_path}"
            }

            return timeseries_dirs
        }

        // Process each beam's time series
        timeseries_channel.map { pointing, fits_file, cluster, beam_name, beam_id, utc, ra, dec, cdm, ts_dir ->
            def dat_files = file("${ts_dir}/*.dat")
            def inf_files = file("${ts_dir}/*.inf")
            tuple(
                pointing, file(fits_file), cluster, beam_name, beam_id, utc, ra, dec, cdm,
                dat_files, inf_files, ts_dir
            )
        }.set { beam_timeseries }

        // Process each beam separately
        beam_timeseries.each { pointing, fil_file, cluster, beam_name, beam_id, utc, ra, dec, cdm, dat_files, inf_files, ts_dir ->
            log.info "Processing time series for beam ${beam_name} from ${ts_dir}"

            def dat_ch = Channel.fromList(dat_files instanceof List ? dat_files : [dat_files])
            def inf_ch = Channel.fromList(inf_files instanceof List ? inf_files : [inf_files])
            def fil_ch = Channel.of(fil_file)
            def meta_ch = Channel.of(tuple(pointing, fil_file, cluster, beam_name, beam_id, utc, ra, dec, cdm))

            run_accelsearch_single(dat_ch, inf_ch, fil_ch, meta_ch)
        }
    }
}

/*
 * Helper workflow to run accelsearch on a single set of time series files
 */
workflow run_accelsearch_single {
    take:
    dat_channel   // Channel of .dat files
    inf_channel   // Channel of .inf files
    fil_channel   // Filterbank file channel
    meta_channel  // Metadata channel for PulsarX

    main:
    // Run acceleration search on pre-dumped time series
    presto_accelsearch(
        dat_channel.flatten(),
        inf_channel.flatten(),
        file('NO_ZAPLIST')  // No zaplist for pre-dumped time series
    ).set { accel_out }

    // Sift candidates
    presto_sift_candidates(accel_out.accel_files.collect(), fil_channel)
        .set { sifted_out }

    // Choose folding backend
    def fold_backend = params.presto?.fold_backend ?: 'pulsarx'

    if (fold_backend == 'pulsarx') {
        // Generate meta file for PulsarX folding
        def fft_size = params.presto?.fft_size ?: params.peasoup?.fft_size ?: 134217728
        generate_fold_meta(meta_channel, fft_size, fold_backend)
            .set { meta_out }

        // Fold with PulsarX
        presto_fold_pulsarx(fil_channel, sifted_out.sifted_csv, meta_out.meta_file)
            .set { fold_out }

        log.info "Folding complete with PulsarX. Output files are in the publishDir."

    } else {
        // Fold with prepfold (no meta file needed)
        presto_prepfold_batch(fil_channel, sifted_out.sifted_csv, file('NO_MASK'))
            .set { fold_out }

        // Post-processing: PNG conversion and tarball
        presto_pfd_to_png(fold_out.pfd_files.collect())
            .set { png_out }

        // Create merged results
        presto_fold_merge(
            fold_out.pfd_files.collect(),
            png_out.png_files.collect(),
            fold_out.bestprof_files.collect().ifEmpty([]),
            sifted_out.sifted_csv
        ).set { merged_out }

        // Create tarball
        presto_create_tarball(
            merged_out.merged_csv,
            merged_out.all_pfd,
            merged_out.all_png,
            fil_channel
        ).set { tarball_out }
    }

    emit:
    sifted_csv = sifted_out.sifted_csv
}

// Default to `full` if no --entry is given
workflow {
    full()
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

    // Safe access to params (defaults defined at top of this file)
    def outputDir = params.basedir ?: '.'
    def runIDStr = params.runID ?: 'N/A'
    def notificationConfig = params.notification

    // ========================================================================
    // Cumulative runtime tracking across resume runs
    // ========================================================================
    def runtimeFile = file("${outputDir}/.cumulative_runtime_${workflow.sessionId}.txt")
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
    Output Dir:    ${outputDir}
    Run ID:        ${runIDStr}

    Command:       ${workflow.commandLine}
    Profile:       ${workflow.profile}
    Container:     ${workflow.container ?: 'N/A'}
    """.stripIndent()

    println summary

    // Write summary to file
    def summaryFile = file("${outputDir}/pipeline_summary_${workflow.runName}.txt")
    try {
        summaryFile.text = summary
        println "Pipeline summary written to: ${summaryFile}"
    } catch (Exception e) {
        println "WARNING: Could not write summary file: ${e.message}"
    }

    // Send email notification if enabled
    if (notificationConfig?.enabled && notificationConfig?.email) {
        def subject = "[ELDEN-RING] Pipeline ${success ? 'COMPLETED' : 'FAILED'}: ${runIDStr != 'N/A' ? runIDStr : workflow.runName}"
        def shouldSend = (success && notificationConfig?.on_complete) ||
                         (!success && notificationConfig?.on_fail)

        if (shouldSend) {
            try {
                sendMail(
                    to: notificationConfig.email,
                    subject: subject,
                    body: summary
                )
                println "Email notification sent to: ${notificationConfig.email}"
            } catch (Exception e) {
                println "WARNING: Could not send email notification: ${e.message}"
                println "Make sure your system has sendmail configured or SMTP settings are correct."
            }
        }
    }
}

workflow.onError {
    // Safe access to params (defaults defined at top of this file)
    def runIDStr = params.runID ?: 'N/A'
    def notificationConfig = params.notification

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
    if (notificationConfig?.enabled && notificationConfig?.email && notificationConfig?.on_error) {
        def subject = "[ELDEN-RING] Pipeline ERROR: ${runIDStr != 'N/A' ? runIDStr : workflow.runName}"

        try {
            sendMail(
                to: notificationConfig.email,
                subject: subject,
                body: errorMessage
            )
            println "Error notification sent to: ${notificationConfig.email}"
        } catch (Exception e) {
            println "WARNING: Could not send error notification: ${e.message}"
        }
    }
}
