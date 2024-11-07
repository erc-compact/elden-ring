#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filtool as filtool } from './processes'
include { nearest_power_of_two_calculator as nearest_power_of_two_calculator } from './processes'
include { generateDMFiles as generateDMFiles } from './processes'
include { peasoup as peasoup } from './processes'
include {parse_xml as parse_xml} from './processes'
include {splitcands as splitcands} from './processes'
include {psrfold as psrfold} from './processes'
include {pics_classifier as pics_classifier} from './processes'


workflow {
    // parse the csv file to get the list of fits files and params
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

    // RFI removal with filtool
    processed_fil_file = filtool(fits_file_channel_and_meta, params.threads, params.telescope)

    // create a new channel with the metadata and new file path
    new_fil_file_channel = processed_fil_file.map { metadata, filepath -> 
        def (raw_file, cluster, beam_name, beam_id, utc_start) = metadata
                    // Trim the string values
        cluster = cluster.trim()
        beam_name = beam_name.trim()
        beam_id = beam_id.trim()
        utc_start = utc_start.trim()
        return tuple(filepath, cluster, beam_name, beam_id, utc_start)
        }

    //find the nearest power of two
    nearest_power_of_two = nearest_power_of_two_calculator( new_fil_file_channel )//.map { fil_file, cluster, beam_name, beam_id, utc_start, nearest_power_of_2 -> 
    //     return tuple(fil_file, cluster, beam_name, beam_id, utc_start, nearest_power_of_2)
    // }

    // Generate DM files
    dm_file = generateDMFiles()
    // dm_file.view()

    // Peasoup
    // peasoup_input_channel = nearest_power_of_two.combine(dm_files)
    // peasoup_input_channel.view()
    peasoup_output = peasoup(nearest_power_of_two, dm_file)

    // Aggregate the output from peasoup
    // gives me a list of xml files for each cluster, beam_name, utc and fil_file

    basename_ch = peasoup_output.map { cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_file, xml_path -> 
        def fil_base_name = fil_file.getBaseName()
        
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_base_name, fil_file, xml_path)
        }

    peasoup_output_grouped = basename_ch.groupTuple(by: [0, 1, 2, 3, 4,6]).map { cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_base_name, fil_file, xml_path -> 
        def first_fil_file = fil_file[0]
        
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_base_name, first_fil_file, xml_path)
        }
    // peasoup_output_grouped.view()

    // Parse the XML files
    parse_xml_results = parse_xml(peasoup_output_grouped).map { cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_base_name, fil_file, xml_path, candidate_csv, metafile -> 
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_base_name, fil_file, xml_path, candidate_csv, metafile)
        }
    // parse_xml_results.view()
    
    // Fold the candidates after splitting the csv file into candfiles
    splitcands_channel = splitcands(parse_xml_results).map { cluster, beam_name, beam_id, utc_start, fft_size, dm_file, fil_base_name, fil_file, xml_path, candidate_csv, metafile, allCands, candfile, metatext -> 
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, fil_base_name, fil_file, candfile, metatext)  
    }.transpose().view()
    // splitcands_channel.view()

    // Fold the candidates
    folded_cands = psrfold(splitcands_channel).map { cluster, beam_name, beam_id, utc_start, fft_size, fil_base_name, fil_file, candfile, metatext, pngs, ars, cands -> 
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, fil_base_name, fil_file, candfile, metatext, pngs, ars, cands)
        }.view() // just doing this to keep pics dependent on psrfold

    // Classify the candidates
    // classified_cands = pics_classifier(folded_cands)
}
    