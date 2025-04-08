#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filtool } from './processes'
// include { nearest_power_of_two_calculator } from './processes'
include { generateDMFiles } from './processes'
include { segmented_params } from './processes'
include { peasoup } from './processes'
include { birdies } from './processes'
include { parse_xml } from './processes'
// include { splitcands } from './processes'
include { search_fold_merge } from './processes'
include { psrfold } from './processes'
include { pics_classifier } from './processes'
include { readfile } from './processes'
include { generateRfiFilter } from './processes'
include { alpha_beta_gamma_test } from './processes'
include { syncFiles } from './processes'
include { create_candyjar_tarball } from './processes'


workflow {
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
            // extract file name from path
            def filename = fits_files.tokenize("/")[-1]
            return tuple(pointing,fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, filename)
        }

    // Copy from tape? 
    def orig_fits_channel
    if (params.copy_from_tape.run_copy) {
        orig_fits_channel = syncFiles(fits_file_channel_and_meta)
    } else {
        orig_fits_channel = fits_file_channel_and_meta
    }

    // Generate DM files
    dm_file = generateDMFiles().flatMap {it}

    // def new_fil_file_channel  // Declare the variable here so it can be used for each case

    if (params.filtool.run_filtool) {

        // Run readfile process to get metadata
        readfile_output = readfile(orig_fits_channel)

        if (params.generateRfiFilter.run_rfi_filter) {
            // Run generateRfiFilter process using the time_per_file
            generateRfiFilter_output = generateRfiFilter(readfile_output).map { pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, rfi_filter_string, tsamp, nsamples, subintlength, pngs, rfitxt -> 
                return tuple(pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, rfi_filter_string, tsamp, nsamples, subintlength)
            }

            // RFI removal with filtool using the generated rfi_filter_string   
            new_fil_file_channel = filtool(generateRfiFilter_output, params.threads, params.telescope)

        } else {
            // Skip RFI filtering processes
            // Use default RFI filters based on telescope
            filtool_input = readfile_output.map { pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, time_per_file, tsamp, nsamples, subintlength ->
                def default_rfi_filter = params.filtool.rfi_filter_list[params.telescope]
                return tuple(pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, default_rfi_filter, tsamp, nsamples, subintlength)
            }
            // RFI removal with filtool using default rfi_filter_string
            new_fil_file_channel = filtool(filtool_input, params.threads, params.telescope)
        } 

        // def cleanup_input = orig_fits_channel.combine(new_fil_file_channel).map { tuple -> tuple[0] }.view()

        // if (params.filtool.run_filtool_cleanup) {
        //     // Cleanup branch: Remove the original files after processing
        //     cleanupResult = dataCleanup(cleanup_input)
        // } else {
        //     // No cleanup branch: Do not remove the original files after processing
        //     cleanupResult = cleanup_input
        // }

    // no filtool

    } else {
        // Run readfile process to get metadata
        readfile_output = readfile(orig_fits_channel)

        // Create a new channel with the metadata and new file path
        new_fil_file_channel = readfile_output.map { pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, time_per_file, tsamp, nsamples, subintlength -> 
            return tuple(pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, tsamp, nsamples, subintlength)
        }
    }

    // segmented search
    // Split the data into segments
    split_params_input = new_fil_file_channel.flatMap { pointing, filepath, cluster, beam_name, beam_id, utc_start, ra, dec, tsamp, nsamples, subintlength -> params.peasoup.segments.collect { segments -> 
        return tuple(pointing, filepath, cluster, beam_name, beam_id, utc_start, ra, dec, segments)}}

    segmented_params = segmented_params(split_params_input)

    // Create a channel with the segmented parameters for peasoup
    peasoup_input = segmented_params.flatMap { pointing, filepath, cluster, beam_name, beam_id, utc_start, ra, dec, tsamp, nsamples, segments, segment_file -> 
        // parse the segment file
        segment_file.splitCsv(header : true, sep : ',').collect { row -> 
            def segment_id = row.i.trim()
            def fft_size = row.fft_size.trim()
            def start_sample = row.start_sample.trim()
            def nsamples_per_segment = row.nsamples_per_segment.trim()
            return tuple(pointing, filepath, cluster, beam_name, beam_id, utc_start, ra, dec, tsamp, nsamples_per_segment, segments, segment_id, fft_size, start_sample)
        }
    }

    // you wanted to work from here onwards by making the changes to this channel.

    // // Generate birdies file
    birdies_output = birdies(peasoup_input)

    // Peasoup processing
    peasoup_output = peasoup(birdies_output, dm_file)

    // Aggregate the output from peasoup
    peasoup_output_grouped = peasoup_output.map { pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, dm_file, fil_file, xml_path, birdies_file, start_sample, nsamples -> 
        def fil_base_name = fil_file.getBaseName()
        return tuple(pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, dm_file, fil_base_name, fil_file, xml_path, start_sample)
    }
    .groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 14])
    //Group by the pointing, cluster, beam_name, beam_id, utc_start,ra, dec, fft_size, segments, segment_id, fil_base_name, and start_sample
    .map { pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, dm_file, fil_base_name, fil_file, xml_paths, start_sample -> 
        def combined = []
         // select only one fil file name from the list of same fil_file names.
        def first_fil_file = (fil_file instanceof List) ? fil_file[0] : fil_file
        // Combine and sort to keep the resume functionality
        (0..<dm_file.size()).each { i -> 
            combined << [dm: dm_file[i], xml: xml_paths[i]]
        }
        def combined_sorted = combined.sort { a, b -> a.dm <=> b.dm }
        def sorted_dm_file = combined_sorted.collect { it.dm }
        def sorted_xml_file = combined_sorted.collect { it.xml}
        return tuple(pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, sorted_dm_file, fil_base_name, first_fil_file, sorted_xml_file, start_sample)
    }

    // Parse the XML files and adjust the tuple
    parse_xml_results = parse_xml(peasoup_output_grouped)

    //add sifting here
    // splitcands = splitcands(parse_xml_results)

    // Split the candidates and adjust the tuple
    splitcands_channel = parse_xml_results
        .flatMap { it -> 
            def (pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, dm_file, fil_base_name, fil_file, xml_file, start_sample, filtered_candidate_csv, unfiltered_candidates_csv, candfiles, metafile, allCands) = it
            def cList = candfiles instanceof List ? candfiles : [candfiles]
            cList.collect { candfile ->
            return tuple(pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, fil_base_name, fil_file, start_sample, filtered_candidate_csv, candfile, metafile)
            }
        }

    // Fold the candidates
    // pulsarx_output = psrfold(splitcands_channel)

    pulsarx_output_grouped = psrfold(splitcands_channel)
        .groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 7, 9])
        .map { pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, fil_base_name, fil_file, filtered_candidate_csv, candfile, metatext, pngs, archives, cands -> 
            pngs = pngs instanceof List ? pngs : [pngs]
            archives = archives instanceof List ? archives : [archives]
            cands = cands instanceof List ? cands : [cands]
            pngs = pngs.flatten()
            archives = archives.flatten()
            cands = cands.flatten()
            return tuple(pointing, cluster, beam_name, beam_id, utc_start, ra, dec, fft_size, segments, segment_id, fil_base_name, filtered_candidate_csv, candfile, metatext, pngs, archives, cands)
        }

    search_fold_merged = search_fold_merge(pulsarx_output_grouped)
    
    // collect all the output from pulsarx for a single filterbank file.

    alpha_beta_gamma = alpha_beta_gamma_test(search_fold_merged)

    // Classify the candidates (optional)
    classified_cands = pics_classifier(search_fold_merged)

    //Wait for all the alpha_beta_output and pics_output to be ready
    alpha_beta_output_combined = alpha_beta_gamma.collect(flat: false)
    pics_output_combined = classified_cands.collect(flat: false)
    
    alpha_beta_channel   = alpha_beta_output_combined.flatMap{ it }
    pics_channel = pics_output_combined.flatMap { it }
    alpha_beta_pics_combined = alpha_beta_channel.join(pics_channel, by: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    
    alpha_beta_pics_results_file = alpha_beta_pics_combined
        .map { row -> row.join(',') } // Convert each list to a CSV line
        .collectFile(name: 'alpha_beta_pics_combined.csv', newLine: true, storeDir: params.basedir)  // Append to the CSV file
     if (params.alpha_beta_gamma.create_candyjar_tarball) {

    def output_tarballname = "${params.alpha_beta_gamma.output_dir_alpha_pics_results}.tar.gz"

    candyjar_tarball_input = alpha_beta_pics_results_file.map { it -> 
        return tuple(it, output_tarballname)
    }
    create_candyjar_tarball(candyjar_tarball_input)
    }
    
    // Prepare candidates for candyjar and create the tarball
    // candyjar_output = candyjar(classified_cands)
}
