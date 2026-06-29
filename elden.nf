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
include { palClean } from './processes'
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
include { generateRfiFilter as generateRfiFilterCleaned } from './processes'
include { readfile as readfileCleaned } from './processes'
include { alpha_beta_gamma_test } from './processes'
include { syncFiles } from './processes'
include { create_candyjar_tarball } from './processes'
include { parfold } from './processes'
include { candypolice_pulsarx} from './processes'
include { extract_candidates } from './processes'
include { dada_to_fits } from './processes'
include { merge_filterbanks } from './processes'
include { split_filterbank } from './processes.nf'

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
params.entry = 'full'
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
            def fits_pattern = row.fits_files.trim()
            def fits_files = file(fits_pattern)
            if (fits_files instanceof List) {
                if (fits_files.isEmpty()) {
                    error "No files matched fits_files pattern '${fits_pattern}' for ${row.cluster}/${row.beam_name}"
                }
                fits_files = fits_files.sort()
            }
            def cluster = row.cluster.trim()
            def beam_name = row.beam_name.trim()
            def beam_id = row.beam_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-")
            def ra = row.ra.trim()
            def dec = row.dec.trim()
            def cdm = row.cdm.trim()
            def filename = (fits_files instanceof List ? fits_files[0] : fits_files).getName()
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
    // readfile now outputs path files for numeric values — read contents here
    def rdout = readfile(orig_fits_channel).map { p,f,c,bn,bi,u,ra,dec,cdm,tpf_f,ts_f,ns_f,si_f ->
        tuple(p,f,c,bn,bi,u,ra,dec,cdm,tpf_f.text.trim(),ts_f.text.trim(),ns_f.text.trim(),si_f.text.trim())
    }

    def suffix = params.generateRfiFilter.suffix ?: ""
    if (params.generateRfiFilter.run_rfi_filter) {
        fil_input = generateRfiFilter(rdout, suffix).map { p,f,c,bn,bi,u,ra,dec,cdm,rfi_f,ts,ns,si,png ->
            return tuple(p,f,c,bn,bi,u,ra,dec,cdm,rfi_f.text.trim(),ts,ns,si)
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

// rfi_clean: read metadata → (opt) generateRfiFilter → filtool → (opt) QC RFI filter
workflow rfi_clean {
    take:
    fil_input

    main:
    if (params.palClean.use_pal_clean) {
        new_fil = palClean(fil_input, params.threads, params.telescope)
    } else if (params.filtool.run_filtool) {
        new_fil = filtool(fil_input, params.threads, params.telescope)
    } else {
        new_fil = fil_input.map{ p,f,c,bn,bi,u,ra,dec,cdm,rfi,ts,ns,si ->
            tuple(p,f,c,bn,bi,u,ra,dec,cdm)
        }
    }

    // Run RFI filter on cleaned data for QC verification
    if (params.generateRfiFilter.run_rfi_filter && params.generateRfiFilter.repeat_rfi_filter) {
        cleaned_with_meta = new_fil.map { p,f,c,bn,bi,u,ra,dec,cdm ->
            tuple(p,f,c,bn,bi,u,ra,dec,cdm,f.getName())
        }
        cleaned_rdout = readfileCleaned(cleaned_with_meta).map { p,f,c,bn,bi,u,ra,dec,cdm,tpf_f,ts_f,ns_f,si_f ->
            tuple(p,f,c,bn,bi,u,ra,dec,cdm,tpf_f.text.trim(),ts_f.text.trim(),ns_f.text.trim(),si_f.text.trim())
        }
        generateRfiFilterCleaned(cleaned_rdout, "_cleaned")
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
    // segmented_params outputs tsamp/nsamples as path files — read contents here
    def seg_ch = segmented_params(split_params).map { p,fp,c,bn,bi,u,ra,dec,cdm,ts_f,ns_f,seg,segf ->
        tuple(p,fp,c,bn,bi,u,ra,dec,cdm,ts_f.text.trim(),ns_f.text.trim(),seg,segf)
    }

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
        .filter { p,c,bn,bi,u,ra,dec,cdm,fft_size,seg,seg_id,fil_base,fil_file,start_sample,filtered_candidate_csv,candfile,metafile ->
            // Skip candfiles with no candidates (header-only lines start with #)
            candfile.readLines().any { !it.startsWith('#') && it.trim() }
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

    def joined = abg.join(pics, by: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11])
        .map { row ->
            def pointing = row[1]
            def cluster  = row[2]
            tuple(pointing, cluster, row.join(','))
        }

    // Collect lines into one CSV per (pointing, cluster) pair.
    // collectFile drops the upstream tuple, so we track keys separately via an
    // MD5 index and re-join — no parsing of filenames needed.
    def keys = joined
        .map { pointing, cluster, line -> tuple(pointing, cluster) }
        .unique()
        .map { pointing, cluster ->
            def idx = "${pointing}_${cluster}".md5()
            tuple(idx, pointing, cluster)
        }

    def csv_files = joined
        .map { pointing, cluster, line ->
            def idx = "${pointing}_${cluster}".md5()
            tuple(idx, line)
        }
        .collectFile { idx, line -> [ "${idx}_combined.csv", line + '\n' ] }
        .map { f -> tuple(f.getName().replace('_combined.csv', ''), f) }

    keys
        .join(csv_files)
        .map { idx, pointing, cluster, f -> tuple(pointing, cluster, f) }
        .set{ abg_pics_combined_csv }

    emit:
    abg_pics_combined_csv
}

workflow candyjar_tarball {
    take:
    abg_pics_combined_csv  // tuple(pointing, cluster, csv_file)

    main:
    if (params.alpha_beta_gamma.create_candyjar_tarball) {
        def tar_input = abg_pics_combined_csv.map { pointing, cluster, f ->
            def safe_pointing = pointing.replaceAll(':', '-')
            def tarball_name = "${cluster}_${safe_pointing}_${params.runID}.tar.gz"
            tuple(cluster, f, tarball_name)
        }

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
    def sf = dada_to_fits(dada_intake())
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
    dada_to_fits(dada_intake()).set{ digifits_out }
}

// -------------DADA->SF->FIL->STACK--------------
workflow run_dada_clean_stack {
    main:
    def sf = dada_to_fits(dada_intake())
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
    main:
    def intake_ch = intake()
    rfi_filter(intake_ch)
}

//-------------Filtool the files ----------------
workflow run_rfi_clean {
    main:
    def intake_ch = intake()
    def rfi_ch    = rfi_filter(intake_ch)
    rfi_clean(rfi_ch)
}

// ---------- Run search and fold on filtooled files -----
// run_search assumes rfi_cleaned files inside the files_lists
workflow run_search_fold {
    main:
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
    main:
    def parfile_ch    = Channel.fromPath("${params.parfold.parfile_path}")
    def intake_ch     = intake()
    def rfi_ch        = rfi_filter(intake_ch)
    def cleaned_ch    = rfi_clean(rfi_ch)
    def parfold_input = cleaned_ch.map { p, fil, c, bn, bi, u, ra, dec, dm ->
        tuple(p, fil, c, bn, bi, u, ra, dec)
    }
    def parfold_out   = parfold(parfold_input, parfile_ch)

    emit:
    parfold_out
}



workflow candypolice {
    main:
    def intake_ch    = intake()
    def rdout        = readfile(intake_ch).map { p,f,c,bn,bi,u,ra,dec,tpf,ts,ns,si ->
        tuple(p,f,c,bn,bi,u,ra,dec,ts,ns,si)
    }
    def candyjar_csv = Channel.fromPath("${params.candypolice.input_csv}")
    def candfiles    = extract_candidates(candyjar_csv).flatMap { cd, cf ->
        def cList = cf instanceof List ? cf : [cf]
        return tuple(cList)
    }
    candypolice_pulsarx(rdout, candfiles)
}

// Select workflow via params.entry (replaces the deprecated -entry CLI flag in NF26+)
// Default: full
workflow {
    def entry = params.entry ?: 'full'
    if      (entry == 'full')                 full()
    else if (entry == 'run_search_fold')      run_search_fold()
    else if (entry == 'run_dada_search')      run_dada_search()
    else if (entry == 'run_dada_clean_stack') run_dada_clean_stack()
    else if (entry == 'run_digifits')         run_digifits()
    else if (entry == 'generate_rfi_filter')  generate_rfi_filter()
    else if (entry == 'run_rfi_clean')        run_rfi_clean()
    else if (entry == 'fold_par')             fold_par()
    else if (entry == 'candypolice')          candypolice()
    else error "Unknown entry workflow: '${entry}'. Valid options: full, run_search_fold, run_dada_search, run_dada_clean_stack, run_digifits, generate_rfi_filter, run_rfi_clean, fold_par, candypolice"

    workflow.onComplete { onComplete() }
    workflow.onError    { onError()    }
}

// ============================================================================
// WORKFLOW COMPLETION HANDLERS - Monitoring and Notifications
// ============================================================================

def formatDuration(long millis) {
    def totalSeconds = (millis / 1000L) as long
    def seconds = totalSeconds % 60
    def minutes = Math.floorDiv(totalSeconds, 60L) % 60
    def hours   = Math.floorDiv(totalSeconds, 3600L) % 24
    def days    = Math.floorDiv(totalSeconds, 86400L)
    if      (days > 0)    return String.format("%dd %dh %dm %ds", days, hours, minutes, seconds)
    else if (hours > 0)   return String.format("%dh %dm %ds", hours, minutes, seconds)
    else if (minutes > 0) return String.format("%dm %ds", minutes, seconds)
    else                  return String.format("%ds", seconds)
}

def onComplete() {
    def duration   = workflow.duration
    def success    = workflow.success
    def exitStatus = workflow.exitStatus
    def outputDir  = params.basedir ?: '.'
    def runIDStr   = params.runID ?: 'N/A'
    def notificationConfig = params.notification

    // Cumulative runtime tracking across resume runs
    def runtimeFile      = file("${outputDir}/.cumulative_runtime_${workflow.sessionId}.txt")
    def cumulativeMillis = 0L
    def runCount         = 1
    try {
        if (runtimeFile.exists()) {
            def lines = runtimeFile.readLines()
            if (lines.size() >= 2) {
                cumulativeMillis = lines[0].toLong()
                runCount         = lines[1].toInteger() + 1
            }
        }
    } catch (Exception e) {
        println "WARNING: Could not read cumulative runtime file: ${e.message}"
    }
    def currentMillis = duration.toMillis()
    def totalMillis   = cumulativeMillis + currentMillis
    try {
        runtimeFile.text = "${totalMillis}\n${runCount}\n"
    } catch (Exception e) {
        println "WARNING: Could not write cumulative runtime file: ${e.message}"
    }

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
    This run:               ${formatDuration(currentMillis)}
    Total (across resumes): ${formatDuration(totalMillis)}
    Number of runs:         ${runCount}

    Completed:     ${workflow.complete}

    Work Dir:      ${workflow.workDir}
    Output Dir:    ${outputDir}
    Run ID:        ${runIDStr}

    Command:       ${workflow.commandLine}
    Profile:       ${workflow.profile}
    Container:     ${workflow.container ?: 'N/A'}
    """.stripIndent()

    println summary

    def summaryFile = file("${outputDir}/pipeline_summary_${workflow.runName}.txt")
    try {
        summaryFile.text = summary
        println "Pipeline summary written to: ${summaryFile}"
    } catch (Exception e) {
        println "WARNING: Could not write summary file: ${e.message}"
    }

    if (notificationConfig?.enabled && notificationConfig?.email) {
        def subject    = "[ELDEN-RING] Pipeline ${success ? 'COMPLETED' : 'FAILED'}: ${runIDStr != 'N/A' ? runIDStr : workflow.runName}"
        def shouldSend = (success && notificationConfig?.on_complete) ||
                         (!success && notificationConfig?.on_fail)
        if (shouldSend) {
            try {
                sendMail(to: notificationConfig.email, subject: subject, body: summary)
                println "Email notification sent to: ${notificationConfig.email}"
            } catch (Exception e) {
                println "WARNING: Could not send email notification: ${e.message}"
            }
        }
    }
}

def onError() {
    def runIDStr           = params.runID ?: 'N/A'
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

    if (notificationConfig?.enabled && notificationConfig?.email && notificationConfig?.on_error) {
        def subject = "[ELDEN-RING] Pipeline ERROR: ${runIDStr != 'N/A' ? runIDStr : workflow.runName}"
        try {
            sendMail(to: notificationConfig.email, subject: subject, body: errorMessage)
            println "Error notification sent to: ${notificationConfig.email}"
        } catch (Exception e) {
            println "WARNING: Could not send error notification: ${e.message}"
        }
    }
}
