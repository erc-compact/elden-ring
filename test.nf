#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow {
    // Create a sample tuple that mimics your splitcands output
    def sample = tuple(
        "CLUSTER1",            // cluster
        "BEAM_A",              // beam_name
        1,                     // beam_id
        "2024-04-05T03:23:38",  // utc_start
        33554432,              // fft_size
        2,                     // segments
        1,                     // segment_id
        "DM1",                 // dm_file
        "FIL_BASE",            // fil_base_name
        "fil_file.fil",        // fil_file (path)
        "xml_file.xml",        // xml_files (path)
        12345,                 // start_sample
        "candidate.csv",       // candidate_csv (path)
        "metafile.txt",        // metafile (path)
        "allCands.txt",        // allCands (path)
        [ "cand1.candfile", "cand2.candfile", "cand3.candfile" ], // candfiles array (paths)
        "meta.txt"             // meta_txt (path)
    )
    
    // Create a channel with the sample tuple
    splitcands_channel = Channel.of(sample).view()
    
    // Flatten the candidate files
    flattened = splitcands_channel
        .flatMap { t ->
            // Destructure the tuple into individual components
            def (cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id,
                 dm_file, fil_base_name, fil_file, xml_files, start_sample, candidate_csv,
                 metafile, allCands, candfiles, meta_txt) = t
                 
            // For each candidate file, emit a new tuple with all the other metadata intact
            candfiles.collect { candfile ->
                tuple(cluster, beam_name, beam_id, utc_start, fft_size, segments, segment_id,
                      dm_file, fil_base_name, fil_file, xml_files, start_sample, candidate_csv,candfile, meta_txt)
            }
         }
         
    // View the output tuples  
    flattened.view { it -> println it }
}
