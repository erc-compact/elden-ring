process readfile {
    label 'readfile'
    container "${params.presto_image}"

    input:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start)

    output:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), env(time_per_file), env(tsamp), env(nsamples), env(subintlength)

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

process generateRfiFilter {
    label 'generate_rfi_filter'
    container "${params.rfi_mitigation_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/RFIFILTER/", pattern: "*.{png,txt}", mode: 'symlink'

    input:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(time_per_file), val(tsamp), val(nsamples), val(subintlength)

    output:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), env(rfi_filter_string) , val(tsamp), val(nsamples) , val(subintlength), path("*.png"), path("*.txt")

    script:
    def num_intervals = Math.floor(time_per_file.toFloat() / 400) as int
    // def num_intervals = 2
    """
    #!/bin/bash
    export MPLCONFIGDIR=/tmp
    export NUMBA_CACHE_DIR=/tmp
    python3 ${params.generateRfiFilter.script} ${fits_files} . --target_resolution_ms ${params.generateRfiFilter.target_resolution_ms} --num_intervals ${num_intervals}

    zap_commands=\$(grep -Eo '[0-9.]+ *- *[0-9.]+' combined_frequent_outliers.txt | \\
    awk -F '-' '{gsub(/ /,""); print "zap "\$1" "\$2}' | tr '\\n' ' ')

    rfi_filter_string="kadaneF 8 4 zdot \${zap_commands}"
    echo "\${rfi_filter_string}" > rfi_filter_string.txt

    mv combined_sk_heatmap_and_histogram.png ${beam_name}_rfi.png
    mv combined_frequent_outliers.txt combined_frequent_outliers_${beam_name}.txt
    mv block_bad_channel_percentages.txt block_bad_channel_percentages_${beam_name}.txt
    """
}

process filtool {
    label 'filtool'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/CLEANEDFIL/", pattern: "*.fil", mode: 'copy'

    input:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(rfi_filter_string), val(tsamp), val(nsamples) , val(subintlength)
    val threads
    val telescope

    output:
    tuple path("*.fil"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples), val(subintlength)
    
    script:
    def outputFile = "${cluster.trim()}_${utc_start.trim()}_${beam_name.trim()}_clean"
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

// process nearest_power_of_two_calculator {
//     label 'nearest_power_two'
//     container "${params.presto_image}"

//     input:
//     tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples), val(subintlength)
  
//     output:
//     tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples), env(nearest_power_of_2)

//     script:
//     """
//     #!/bin/bash

//     # Function to calculate nearest power of 2
//     nearest_power_of_2() {
//     python3 -c "import math; value = int(\${1}); log2 = math.log(value, 2); rounded_log2 = math.ceil(log2) if log2 - int(log2) >= 0.35 else int(log2); nearest_power_of_2 = 2 ** rounded_log2; print(nearest_power_of_2)"
//     }

//     # Calculate nearest power of 2 for nsamples, nsamples/2, and nsamples/4
//     nearest_power_of_2_value=\$(nearest_power_of_2 ${nsamples})

//     # Output the results to a file
//     echo "$nearest_power_of_2_value" > nearest_power_of_2.txt
//     """
// }

process generateDMFiles {
    label "generateDMFiles"
    container "${params.presto_image}"
    publishDir "${params.basedir}/${cluster}/DMFILES/", pattern: "*.dm", mode: 'copy'

    input:
    tuple path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start)

    output:
    path("*.dm")

    script:
    """
    #!/usr/bin/env python3
    import numpy as np

    # Generate the DM file
    dm_start = ${params.ddplan.dm_start}
    dm_end = ${params.ddplan.dm_end}
    dm_step = ${params.ddplan.dm_step}
    dm_sample = ${params.ddplan.dm_sample}

    # Create DM values with a step of dm_step
    dm_values = np.round(np.arange(dm_start, dm_end, dm_step), 3)

    # Split DM values into multiple files, each containing dm_sample number of lines
    for i in range(0, len(dm_values), dm_sample):
        chunk = dm_values[i:i + dm_sample]
        end_index = min(i + dm_sample, len(dm_values))
        filename = f'dm_{dm_values[i]}_{dm_values[end_index - 1]}.dm'
        np.savetxt(filename, chunk, fmt='%f')
    """
}

process birdies {
    label 'birdies'
    container "${params.peasoup_image}"
    publishDir "${params.basedir}/${cluster}/segment_${segments}/${segments}${segment_id}/${beam_name}/BIRDIES/", pattern: "*.{xml,txt}", mode: 'copy'

    input:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample)

    output:
    path("birdies.txt")

    script:
    """
    #!/bin/bash
    peasoup -p -i ${fil_file} --fft_size ${fft_size} -m 10.0 -t 1 -n 16 --acc_start 0.0 --acc_end 0.0 --ram_limit_gb 50.0 --dm_start 0.0 --dm_end 0.0  --start_sample ${start_sample} 

    #Rename the output file
    mv **/*.xml ${beam_name}_birdies.xml

    python3 ${projectDir}/scripts/birdies_parser.py --xml_file  *birdies.xml
    """
}

process segmented_params {
    label 'segmented_params'
    container "${params.presto_image}"
    publishDir "${params.basedir}/${cluster}/segment_${segments}/${beam_name}/SEGPARAMS/", pattern: "*.csv", mode: 'copy'

    input:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(segments)
    
    output:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), env(tsamp), env(nsamples), val(segments), path("*.csv")

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
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/SEARCH/", pattern: "*.xml", mode: 'copy'
    scratch true

    input:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample) 
    path(birdies_file)
    each path(dm_file)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), path(dm_file), path(fil_file, followLinks: false), path("*.xml"), path(birdies_file), val(start_sample), val(nsamples)

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

    peasoup -p -i ${fil_file} --fft_size ${fft_size} --limit ${params.peasoup.total_cands_limit} -m ${params.peasoup.min_snr} -t ${params.peasoup.ngpus} -n ${params.peasoup.nharmonics} --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file} \${birdies_string} --start_sample ${start_sample} 

    #Rename the output file
    mv **/*.xml ${beam_name}_${dm_file.baseName}_ck${segments}${segment_id}_overview.xml
    """
}

process parse_xml {
    label 'parse_xml'
    container "${params.presto_image}"
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/PARSEXML/", pattern: "*{csv,meta}", mode: 'copy'
    stageInMode 'symlink'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size),val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files),  val(start_sample), path("*candidates.csv"), path("*.meta")
    
    script:
    """
    python3 ${params.parse_xml.script} ${xml_files} --chunk_id ${segments}${segment_id} --outfile ${utc_start}_${beam_name}_ck${segments}${segment_id}_candidates.csv --metafile ${utc_start}_${beam_name}_ck${segments}${segment_id}_metafile.meta
    """
}

process splitcands {
    label 'splitcands'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/SPLITS/", pattern: "*{allCands.txt,candfile,meta.txt}", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), val(start_sample), path(candidate_csv), path(metafile)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), val(start_sample), path(candidate_csv), path(metafile), path('*allCands.txt'), path('*candfile'), path('*meta.txt')

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
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/FOLDING/", pattern: "*.png", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/FOLDING/", pattern: "*.ar", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/FOLDING/", pattern: "*.cands", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/FOLDING/", pattern: "search_fold_cands.csv", mode: 'copy'
    
    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), val(start_sample), path(candfile), path(metatext)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file, followLinks: false), path(candfile), path(metatext), path("*.png"), path("*.ar"), path("*.cands"), path("search_fold_cands.csv")


    
    script:
    """
    python3 ${baseDir}/scripts/pulsarx_fold.py -meta ${metatext} -cands ${candfile}

    # thought about adding the csv step here. 07/02/25

    fold_cands = \$(ls -v *.ar)
    pulsarx_cands_file = \$(ls -v *.cands)

    python3 ${baseDir}/scripts/fold_cands_to_csv.py -f \${fold_cands} -c \${pulsarx_cands_file}
    """
}

process pics_classifier {
    label "pics_classifier"
    container "/hercules/scratch/fkareem/singularity_img/trapum_pulsarx_fold_docker_20220411.sif"
    publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/CLASSIFICATION/", pattern: "*pics_scores.csv", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands), path(search_fold_cands_csv)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(search_fold_cands_csv), path("*pics_scores.csv") 

    script:
    output_csv = "${cluster}_${beam_name}_ck${segments}${segment_id}_scored.csv"
    """
    python2 ${baseDir}/scripts/pics_classifier_multiple_models.py -m ${params.pics_model_dir} -o ${output_csv}
    """
}

process alpha_beta_gamma_test {
label "alpha_beta_gamma_test"
container "${params.pulsarx_image}"
publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/ZERODM/", pattern: "DM0*.png", mode: 'copy'
publishDir "${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/ABG", pattern: "search_fold_alpha_beta_gamma_merged.csv", mode: 'copy'

input:
tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands), path(search_fold_cands_csv)

output:
tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands), path(search_fold_cands_csv), path("*alpha_beta_gamma.csv")

script:
"""
#!/bin/bash
publish_dir="${params.basedir}/${cluster}/segment_${segments}${segment_id}/${beam_name}/ABG"
mkdir -p \${publish_dir}
python3 ${baseDir}/scripts/calculate_alpha_beta_gamma_dmffdot.py -i ${search_fold_cands_csv} -o ${cluster}_${beam_name}_ck${segments}${segment_id}_alpha_beta_gamma.csv -t ${params.alpha_beta_gamma.snr_min} -c -p \${publish_dir} 
"""
}

