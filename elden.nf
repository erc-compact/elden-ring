#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
    
    orig_fits_channel.view()

    emit:
    orig_fits_channel
}

workflow dada_intake {
    main:
    // Parse the CSV file to get the list of dada files and parameters
    dada_file_channel_and_meta = Channel.fromPath("${params.dada.dada_csv}")
        .splitCsv(header : true, sep : ',')
        .map { row -> 
            def pointing = row.pointing.trim()
            def source = row.dada_files
            def cluster = row.cluster.trim()
            def beam_name = row.beam_name.trim()
            def beam_id = row.beam_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-")
            def ra = row.ra.trim()
            def dec = row.dec.trim()
            def cdms = new groovy.json.JsonSlurper().parseText(row.cdm_list)

            def dada_files = source.isDirectory()
                ? source.listFiles().findAll { it.name.startsWith(${params.dada.dada_prefix}) && it.name.endsWith('.dada') }
                : source.text.readLines().collect { new File(it) }

            cdms.collect { cdm ->
                tuple(pointing, dada_files, cluster, beam_name, beam_id, utc_start, ra, dec, cdm)
            }
        }

    dada_file_channel_and_meta.view()

    emit:
    dada_file_channel_and_meta
}

workflow dm {
    main:
    generateDMFiles()
        .flatMap {it}
        .set{ dm_file }
    emit:
    dm_file
}

// generate_rfi_filter 
workflow rfi_filter {
    take:
    orig_fits_channel
    
    main:
    readfile(orig_fits_channel).set{ rdout }

    if (params.generateRfiFilter.run_rfi_filter) {
        fil_input = generateRfiFilter(rdout).map { p,f,c,bn,bi,u,ra,dec,cdm,rfi,ts,ns,si -> 
            return tuple(p,f,c,bn,bi,u,ra,dec, cdm,rfi,ts,ns,si)   
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
        tuple(p,f,c,bn,bi,u,ra,dec,cdm,ts,ns,si)
    }

    emit:
    new_fil
}

// segmentation: break into segments → format for peasoup
workflow segmentation {
    take:
    new_fil

    main:
    split_params = new_fil
        .flatMap { p,fp,c,bn,bi,u,ra,dec,cdm,ts,ns,si -> 
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
            
// search: birdies → search → grouping
workflow search {
    take:
    segmentation
    dm_file

    main:
    birdies(segmentation)
        .set{ bird_out }

    peasoup(bird_out, dm_file)
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
        .groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        .map { p,c,bn,bi,u,ra,dec,fft_size,seg,seg_id,fil_base,fil,f_csv,candfile,metatext,pngs,ar,cands ->
            def f_csv_sort = f_csv instanceof List ? f_csv.sort() : [f_csv]
            def first_f_csv = (f_csv_sort instanceof List) ? f_csv_sort[0] : f_csv_sort
            ar = ar instanceof List ? ar : [ar]
            cands = cands instanceof List ? cands : [cands]
            ar = ar.flatten()
            cands = cands.flatten()
            tuple(p,c,bn,bi,u,ra,dec,fft_size,seg,seg_id,fil_base,first_f_csv,ar,cands)
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

    abg.flatMap{ it }
        .join(pics.flatMap{ it }, by: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
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
}
}

// ----------------Main workflow ------------------
workflow full {
    main:
    def dm_ch        = dm()
    def intake_ch    = intake()
    def rfi_ch       = rfi_filter(intake_ch)
    def cleaned_ch   = rfi_clean(rfi_ch)
    def seg_ch       = segmentation(cleaned_ch)
    def search_ch    = search(seg_ch, dm_ch)
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
    intake()
    dm()
    readfile(intake.out).map{ p,f,c,bn,bi,u,ra,dec,tpf,ts,ns,si ->
        tuple(p,f,c,bn,bi,u,ra,dec,ts,ns,si)
    }.set{rdout}
    segmentation(rdout)
    search(segmentation.out,dm.out)
    xml_parse(search.out)
    fold(xml_parse.out)
    fold_merge(fold.out)
    classify(fold_merge.out)
    candyjar_tarball(classify.out)
}

//------------- parfold workflow ---------------
workflow fold_par {
    parfile_ch = Channel.fromPath("${params.parfold.parfile_path}")
    intake()
    generate_rfi_filter(intake.out)
    rfi_clean(generate_rfi_filter.out)
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

// Default to `full` if no --entry is given
workflow {
    full()
}