#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filtool } from './processes'
include { generateRfiFilter } from './rfi-test-processes'
include { generateRfiFilterSecond } from './rfi-test-processes'
include { readfile } from './processes'

workflow {

    // Parse the CSV file to get the list of FITS files and parameters
    fits_file_channel_and_meta = Channel.fromPath("${params.files_list}")
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            def fits_files = row.fits_files.trim()
            def cluster = row.cluster.trim()
            def beam_name = row.beam_name.trim()
            def beam_id = row.beam_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-")
            return tuple(fits_files, cluster, beam_name, beam_id, utc_start)
        }

    // Step 1: First RFI Filter Generation
    readfile_output = readfile(fits_file_channel_and_meta)
    first_generateRfiFilter_output = generateRfiFilter(readfile_output)

    // // Step 2: Run filtool using the first generated RFI filter
    // processed_fil_file = filtool(first_generateRfiFilter_output, params.threads, params.telescope)

    // // Step 3: Second RFI Filter Generation on Cleaned File
    // second_generateRfiFilter_output = generateRfiFilterSecond(processed_fil_file)
}
