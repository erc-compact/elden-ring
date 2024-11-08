#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filtool } from './processes'
include { nearest_power_of_two_calculator } from './processes'
include { generateDMFiles } from './processes'
include { peasoup } from './processes'
include { parse_xml } from './processes'
include { splitcands } from './processes'
include { psrfold } from './processes'
include { pics_classifier } from './processes'
include { readfile } from './processes'
include { generateRfiFilter } from './processes'

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
            return tuple(fits_files, cluster, beam_name, beam_id, utc_start)
        }

    if (params.generateRfiFilter.run_rfi_filter) {
        // Run readfile process to get time_per_file
        readfile_output = readfile(fits_file_channel_and_meta)

        // Run generateRfiFilter process using the time_per_file
        generateRfiFilter_output = generateRfiFilter(readfile_output)

        // RFI removal with filtool using the generated rfi_filter_string
        processed_fil_file = filtool(generateRfiFilter_output, params.threads, params.telescope)
    } else {
        // Skip RFI filtering processes
        // Use default RFI filters based on telescope
        filtool_input = fits_file_channel_and_meta.map { fits_file_meta ->
            def default_rfi_filter = params.filtool.rfi_filter_list[params.telescope]
            return tuple(fits_file_meta, default_rfi_filter)
        }

        // RFI removal with filtool using default rfi_filter_string
        processed_fil_file = filtool(filtool_input, params.threads, params.telescope)
    }

    