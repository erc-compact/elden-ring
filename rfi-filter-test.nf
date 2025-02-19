#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process readfile {
    label 'readfile'
    container "${params.presto_image}"

    input:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start)

    output:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), env(tsamp), env(nsamples), env(subintlength)

    script:
    """
    #!/bin/bash
    output=\$(readfile ${fits_files})
    echo "\$output"
    time_per_file=\$(echo "\$output" | grep "Time per file (sec)" | awk '{print \$6}')
    tsamp=\$(echo "\$output" | grep "Sample time (us)" | awk '{print \$5}')
    nsamples=\$(echo "\$output" | grep "Spectra per file" | awk '{print \$5}')
    subintlength=\$(echo "scale=10; \$nsamples * \$tsamp * (1/1000000) / 64.0" | bc -l | awk '{print int(\$0)}')
    echo "\${time_per_file}" > time_per_file.txt
    """
}

process filtool {
    label 'filtool'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/CLEANEDFIL/", pattern: "*.fil", mode: 'copy'

    input:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples) , val(subintlength)
    val(rfi_filter_string_id)
    each val(rfi_filter_string)
    val threads
    val telescope

    output:
    tuple path("*.fil"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples), val(subintlength), val(rfi_filter_string_id), val(rfi_filter_string)
    
    script:
    def outputFile = "${cluster.trim()}_${utc_start.trim()}_${beam_name.trim()}_rfi_${rfi_filter_string_id}_clean"
    def source_name = "${cluster.trim()}"

    // Prepare the rfi_filter option
    def zaplist = ''
    if (rfi_filter_string) {
        zaplist = "-z ${rfi_filter_string}"
    }
    """
    #!/bin/bash
    # Get the first file from the inputFile string
    # This is used to determine the file extension
    first_file=\$(echo ${fits_files} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    if [[ ${telescope} == "effelsberg" ]]; then
        if [[ "\${file_extension}" == "fits" ]]; then
            filtool --psrfits --flip -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f ${fits_files} -s ${source_name}
        else 
            filtool -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f ${fits_files} -s ${source_name}
        fi

    elif [[ ${telescope} == "meerkat" ]]; then
        if [[ "\${file_extension}" == "sf" ]]; then
            filtool --psrfits -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f ${fits_files} -s ${source_name}
        else 
            filtool -t ${threads} --telescope ${telescope} ${zaplist} -o "${outputFile}" -f ${fits_files} -s ${source_name}
        fi
    fi
    """
}

process segmented_params {
    label 'segmented_params'
    container "${params.presto_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/SEGPARAMS/", pattern: "*.csv", mode: 'copy'

    input:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(segments), val(rfi_filter_string_id), val(rfi_filter_string)
    
    output:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), env(tsamp), env(nsamples), val(segments), val(rfi_filter_string_id), val(rfi_filter_string), path("*.csv")

    script:
    """
    #!/bin/bash
    output=\$(readfile ${fil_file})
    echo "\$output"
    tsamp=\$(echo "\$output" | grep "Sample time (us)" | awk '{print \$5}')
    nsamples=\$(echo "\$output" | grep "Spectra per file" | awk '{print \$5}')

    # Calculate nsample/segments
    nsamples_per_segment=\$((\${nsamples}/${segments}))

    log2=\$(echo "l(\${nsamples_per_segment})/l(2)" | bc -l)
    decimal_part=\$(echo "\$log2" | awk -F"." '{print "0."\$2}')
    rounded_log2=\$(echo "\$log2" | awk -F"." '{second_num = substr(\$2, 1, 2) + 0; if (second_num >= 35) print \$1+1; else print \$1}')

    fft_size=\$((2**\$rounded_log2))

    # Generate params.csv 
    output_file="${beam_name}_segments_${segments}_params.csv"
    echo "i,fft_size,start_sample,nsamples_per_segment" > \${output_file}
    start_sample=0
    for ((i=0; i<=${segments}-1; i++ )); do
        echo "\$i,\$fft_size,\$start_sample,\$nsamples_per_segment" >> \${output_file}
        start_sample=\$((start_sample + nsamples_per_segment))
    done
    """
}

process peasoup {
    label 'peasoup'
    container "${params.peasoup_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/SEARCH/", pattern: "*.xml", mode: 'copy'
    scratch true

    input:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples),val(rfi_filter_string_id), val(rfi_filter_string), val(segments), val(segment_id), val(fft_size), val(start_sample)
    each path(dm_file) 

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string_id), val(rfi_filter_string), val(segments), val(segment_id), path(dm_file), path(fil_file, followLinks: false), path("*.xml"), path(birdies_file), val(start_sample), val(nsamples)

    script:
    """
    #!/bin/bash

    # check if birdies files is empty
    if [ ! -s ${birdies_file} ]; then
        echo "Birdies file is empty"
        birdies_string=""
    else
        birdies_string="--zapfile ${birdies_file}"
    fi

    peasoup -p -v -i ${fil_file} --fft_size ${fft_size} --limit ${params.peasoup.total_cands_limit} -m ${params.peasoup.min_snr} -t ${params.peasoup.ngpus} -n ${params.peasoup.nharmonics} --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file} \${birdies_string} --start_sample ${start_sample} 

    #Rename the output file
    mv **/*.xml ${beam_name}_${dm_file.baseName}_ck${segments}${segment_id}_rfi__overview.xml
    """
}

process parse_xml {
    label 'parse_xml'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/", pattern: "*{csv,meta}", mode: 'copy'
    stageInMode 'symlink'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample), path("*candidates.csv"), path("*.meta")
    
    script:
    """
    python3 ${params.parse_xml.script} ${xml_files} --chunk_id ${segments}${segment_id} --outfile ${utc_start}_${beam_name}_ck${segments}${segment_id}_candidates.csv --metafile ${utc_start}_${beam_name}_ck${segments}${segment_id}_metafile.meta
    """
}

process rfi_filter_efficiency {
    label 'rfi_filter_efficiency'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/RFI_FILTER/", pattern: "*{fil,meta}", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), val(start_sample), path(candidate_csv), path(metafile)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), val(start_sample), path(candidate_csv), path(metafile), path('*allCands.txt'),path('*candfile'),path('*meta.txt')

    script:
    """
    
    """
}

process splitcands {
    label 'splitcands'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/SPLITS/", pattern: "*{allCands.txt,candfile,meta.txt}", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), val(start_sample), path(candidate_csv), path(metafile)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), val(start_sample), path(candidate_csv), path(metafile), path('*allCands.txt'),path('*candfile'),path('*meta.txt')

    script:
    """
    echo "Running splitcands"
    python3 ${params.splitcands.script} -i ${candidate_csv} -fil ${fil_file} -t pulsarx -n ${params.splitcands.nh} -b ${beam_name}  -threads ${params.splitcands.threads} -p ${params.template_dir} --snr_min ${params.splitcands.snr_min} -ncands ${params.splitcands.ncands} -clfd ${params.splitcands.clfd} -cpn ${params.splitcands.cands_per_node} --beam_id ${beam_id} --utc ${utc_start} --metafile ${metafile}
    """
}

process psrfold {
    label "psrfold"
    container "${params.pulsarx_image}"
    scratch true
    // maxForks 100
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "*.png", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "*.ar", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "*.cands", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "search_fold_cands.csv", mode: 'copy'
    
    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(fil_base_name), path(fil_file), val(start_sample), path(candfile), path(metatext)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(fil_base_name), path(fil_file, followLinks: false), path(candfile), path(metatext), path("*.png"), path("*.ar"), path("*.cands"), path("search_fold_cands.csv")

    script:
    """
    python3 ${baseDir}/scripts/pulsarx_fold.py -meta ${metatext} -cands ${candfile}

    # take candfile number and multiply that with the number of candidates 

    fold_cands=\$(ls -v *.ar)
    pulsarx_cands_file=\$(ls -v *.cands)

    python3 ${baseDir}/scripts/fold_cands_to_csv.py -f \${fold_cands} -c \${pulsarx_cands_file}
    """
}

process pics_classifier {
    label "pics_classifier"
    container "/hercules/scratch/fkareem/singularity_img/trapum_pulsarx_fold_docker_20220411.sif"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/CLASSIFICATION/", pattern: "*scored.csv", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands), path(search_fold_cands_csv)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(search_fold_cands_csv), path("*scored.csv") 

    script:
    output_csv = "${cluster}_${beam_name}_ck${segments}${segment_id}_scored.csv"
    """
    python2 ${baseDir}/scripts/pics_classifier_multiple_models.py -m ${params.pics_model_dir} -o ${output_csv}
    """
}

process alpha_beta_gamma_test {
label "alpha_beta_gamma_test"
container "${params.pulsarx_image}"
publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/ZERODM/", pattern: "DM0*.png", mode: 'copy'
publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/ABG", pattern: "*alpha_beta_gamma.csv", mode: 'copy'

input:
tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands), path(search_fold_cands_csv)

output:
tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(rfi_filter_string), val(segments), val(segment_id), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands), path(search_fold_cands_csv), path("*alpha_beta_gamma.csv")

script:
"""
#!/bin/bash
publish_dir="${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/ABG"
mkdir -p \${publish_dir}
python3 ${baseDir}/scripts/calculate_alpha_beta_gamma_dmffdot.py -i ${search_fold_cands_csv} -o ${cluster}_${beam_name}_ck${segments}${segment_id}_alpha_beta_gamma.csv -t ${params.alpha_beta_gamma.snr_min} -c -p \${publish_dir} 
"""
}


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
    // each row on rfi_filter file is a filter
    rfi_filter_list = Channel.fromPath("${params.rfi_filter_test.rfi_filter_file}")
        .splitCsv(header: false, sep: ',')
        .map { row -> 
            def id = row[0].trim()
            def filter = row[1].trim()
            return tuple(id, filter)
        }.flatMap{ it }.view()

    readfile_output = readfile(fits_file_channel_and_meta)

    // RFI removal with filtool using default rfi_filter_string
    new_fil_file_channel = filtool(readfile_output, rfi_filter_list[0], rfi_filter_list[1], params.threads, params.telescope).view()

    // Split the data into segments
    split_params_input = new_fil_file_channel.flatMap { filepath, cluster, beam_name, beam_id, utc_start, tsamp, nsamples, subintlength, rfi_filter_string_id, rfi_filter_string -> params.peasoup.segments.collect { segments -> 
        return tuple(filepath, cluster, beam_name, beam_id, utc_start, segments, rfi_filter_string_id, rfi_filter_string)}}

    segmented_params = segmented_params(split_params_input)

    // Create a channel with the segmented parameters for peasoup
    peasoup_input = segmented_params.flatMap { filepath, cluster, beam_name, beam_id, utc_start, tsamp, nsamples, segments, rfi_filter_string_id, rfi_filter_string, segment_file -> 
        // parse the segment file
        segment_file.splitCsv(header : true, sep : ',').collect { row -> 
            def segment_id = row.i.trim()
            def fft_size = row.fft_size.trim()
            def start_sample = row.start_sample.trim()
            def nsamples_per_segment = row.nsamples_per_segment.trim()
            return tuple(filepath, cluster, beam_name, beam_id, utc_start, tsamp, nsamples_per_segment, rfi_filter_string_id, rfi_filter_string, segments, segment_id, fft_size, start_sample)
        }
    }

    dm_file = generateDMFiles(fits_file_channel_and_meta).flatMap {it}

    // Peasoup processing
    peasoup_output = peasoup(peasoup_input, dm_file)

    // Aggregate the output from peasoup
    basename_ch = peasoup_output.map { cluster, beam_name, beam_id, utc_start, fft_size, rfi_filter_string, segments, segment_id, dm_file, fil_file, xml_path, birdies_file, start_sample, nsamples -> 
        def fil_base_name = fil_file.getBaseName()
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, rfi_filter_string, segments, segment_id, dm_file, fil_base_name, fil_file, xml_path, start_sample)
    }

    // Group peasoup outputs xml files of the same cluster, beam and segment chunk
    peasoup_output_grouped = basename_ch.groupTuple(by: [0, 1, 2, 3, 4, 5, 6, 7, 9, 12]).map { cluster, beam_name, beam_id, utc_start, fft_size, rfi_filter_string, segments, segment_id, dm_file, fil_base_name, fil_files, xml_paths, start_sample -> 
        def first_fil_file = fil_files[0]
        return tuple(cluster, beam_name, beam_id, utc_start, fft_size, rfi_filter_string, segments, segment_id, dm_file, fil_base_name, first_fil_file, xml_paths, start_sample)
    }

    // Parse the XML files and adjust the tuple
    parse_xml_results = parse_xml(peasoup_output_grouped)

    // // add sifting here
    splitcands = splitcands(parse_xml_results)

    // Split the candidates and adjust the tuple
    splitcands_channel = splitcands
        .flatMap { it -> 
            def (cluster, beam_name, beam_id, utc_start, fft_size, rfi_filter_string, segments, segment_id, dm_file, fil_base_name, fil_file, xml_file, start_sample, candidate_csv, metafile, allCands, candfiles, metatext) = it
            def cList = candfiles instanceof List ? candfiles : [candfiles]
            cList.collect { candfile ->
                tuple(cluster, beam_name, beam_id, utc_start, fft_size, rfi_filter_string, segments, segment_id, fil_base_name, fil_file, start_sample, candfile, metatext)
            }
        }.view()

    // Fold the candidates
    pulsarx_output = psrfold(splitcands_channel)

    alpha_beta_gamma = alpha_beta_gamma_test(pulsarx_output)

    // Classify the candidates (optional)
    classified_cands = pics_classifier(pulsarx_output)
}
