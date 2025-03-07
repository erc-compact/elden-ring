#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filtool } from './processes'
// include { nearest_power_of_two_calculator } from './processes'
include { generateDMFiles } from './processes'
include { segmented_params } from './processes'
include { peasoup } from './processes'
include { birdies } from './processes'
include { parse_xml } from './processes'
include { splitcands } from './processes'
include { psrfold } from './processes'
include { pics_classifier } from './processes'
include { readfile } from './processes'
include { generateRfiFilter } from './processes'
include { alpha_beta_gamma_test } from './processes'
include { syncFiles } from './processes'


workflow {
    // Parse the CSV file to get the list of FITS files and parameters
    fits_file_channel_and_meta = Channel.fromPath("${params.files_list}")
        .splitCsv(header : true, sep : ',')
        .map { row -> 
            def fits_files = row.fits_files.trim()
            def cluster = row.cluster.trim()
            def beam_name = row.beam_name.trim()
            def beam_id = row.beam_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-")
            // extract file name from path
            def filename = fits_files.tokenize("/")[-1]
            return tuple(fits_files, cluster, beam_name, beam_id, utc_start, filename)
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
            generateRfiFilter_output = generateRfiFilter(readfile_output).map { fits_files, cluster, beam_name, beam_id, utc_start, rfi_filter_string, tsamp, nsamples, subintlength, pngs, rfitxt -> 
                return tuple(fits_files, cluster, beam_name, beam_id, utc_start, rfi_filter_string, tsamp, nsamples, subintlength)
            }

            // RFI removal with filtool using the generated rfi_filter_string   
            new_fil_file_channel = filtool(generateRfiFilter_output, params.threads, params.telescope)

        } else {
            // Skip RFI filtering processes
            // Use default RFI filters based on telescope
            filtool_input = readfile_output.map { fits_files, cluster, beam_name, beam_id, utc_start, time_per_file, tsamp, nsamples, subintlength ->
                def default_rfi_filter = params.filtool.rfi_filter_list[params.telescope]
                return tuple(fits_files, cluster, beam_name, beam_id, utc_start, default_rfi_filter, tsamp, nsamples, subintlength)
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
        new_fil_file_channel = readfile_output.map { fits_files, cluster, beam_name, beam_id, utc_start, time_per_file, tsamp, nsamples, subintlength -> 
            return tuple(fits_files, cluster, beam_name, beam_id, utc_start, tsamp, nsamples, subintlength)
        }
    }

    // segmented search
    // Split the data into segments
    split_params_input = new_fil_file_channel.flatMap { filepath, cluster, beam_name, beam_id, utc_start, tsamp, nsamples, subintlength -> params.peasoup.segments.collect { segments -> 
        return tuple(filepath, cluster, beam_name, beam_id, utc_start, segments)}}

    segmented_params = segmented_params(split_params_input)

    // Create a channel with the segmented parameters for peasoup
    peasoup_input = segmented_params.flatMap { filepath, cluster, beam_name, beam_id, utc_start, tsamp, nsamples, segments, segment_file -> 
        // parse the segment file
        segment_file.splitCsv(header : true, sep : ',').collect { row -> 
            def segment_id = row.i.trim()
            def fft_size = row.fft_size.trim()
            def start_sample = row.start_sample.trim()
            def nsamples_per_segment = row.nsamples_per_segment.trim()
            return tuple(filepath, cluster, beam_name, beam_id, utc_start, tsamp, nsamples_per_segment, segments, segment_id, fft_size, start_sample)
        }
    }

    // you wanted to work from here onwards by making the changes to this channel.

    // // Generate birdies file
    birdies_file = birdies(peasoup_input)

    // Peasoup processing
    peasoup_output = peasoup(peasoup_input, dm_file, birdies_file)

    // Aggregate the output from peasoup
    basename_ch = peasoup_output.map { cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id, dm_file, fil_file, xml_path, birdies_file, start_sample, nsamples -> 
        def fil_base_name = fil_file.getBaseName()
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id, dm_file, fil_base_name, fil_file, xml_path, start_sample)
    }

    // Group peasoup outputs xml files of the same cluster, beam and segment chunk
    peasoup_output_grouped = basename_ch.groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 8, 11]).map { cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id, dm_file, fil_base_name, fil_files, xml_paths, start_sample -> 
        def first_fil_file = fil_files[0]
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id, dm_file, fil_base_name, first_fil_file, xml_paths, start_sample)
    }

    // Parse the XML files and adjust the tuple
    parse_xml_results = parse_xml(peasoup_output_grouped)

    // // add sifting here
    splitcands = splitcands(parse_xml_results)

    // Split the candidates and adjust the tuple
    splitcands_channel = splitcands
        .flatMap { it -> 
            def (cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id, dm_file, fil_base_name, fil_file, xml_file, start_sample, candidate_csv, metafile, allCands, candfiles, metatext) = it
            def cList = candfiles instanceof List ? candfiles : [candfiles]
            cList.collect { candfile ->
                tuple(cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id, fil_base_name, fil_file, start_sample, candfile, metatext)
            }
        }

    // Fold the candidates
    pulsarx_output = psrfold(splitcands_channel)

    alpha_beta_gamma = alpha_beta_gamma_test(pulsarx_output)

    // Classify the candidates (optional)
    classified_cands = pics_classifier(pulsarx_output)
}
