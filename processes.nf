process filtool {
    label 'filtool'
    container "${params.pulsarx_image}"

    publishDir "${params.filtool.output_path}/", pattern: "*.fil", mode: 'symlink'

    input:
    val(fits_file_channel_and_meta)
    val threads
    val telescope

    output:
    tuple val(fits_file_channel_and_meta) , path("*.fil")
    
    script:
    def outputFile = "${fits_file_channel_and_meta[1].trim()}_${fits_file_channel_and_meta[4].trim()}_${fits_file_channel_and_meta[2].trim()}_clean"
    def inputFile = "${fits_file_channel_and_meta[0].trim()}"
    def source_name = "${fits_file_channel_and_meta[1].trim()}"

    // Use appropriate rfi masks for the telescope
    def rfi_filter = "${params.filtool.rfi_filter_list[telescope]}"

    """
    # Get the first file from the input_data string
    first_file=\$(echo ${inputFile} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    if [[ ${telescope} == "effelsberg" ]]; then
        if [[ \${file_extension} == "fits" ]]; then
            filtool -psrfits --scloffs --flip -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
        else 
            filtool -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
        fi

    elif [[ ${telescope} == "meerkat" ]]; then
        if [[ \${file_extension} == "sf" ]]; then
            filtool --psrfits --scloffs -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
        else 
            filtool -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
        fi
    fi
    """
}



process nearest_power_of_two_calculator {
    label 'nearest_power_two'
    container "/hercules/scratch/vishnu/singularity_images/presto5_pddot.sif"

    input:
    tuple path(fil_file), val(cluster),val(beam_name), val(beam_id), val(utc_start)

    output:
    tuple path(fil_file), val(cluster),val(beam_name), val(beam_id), val(utc_start), env(nearest_power_of_2)

    script:
    """
    #!/bin/bash

    output=\$(readfile ${fil_file})
    echo "\$output"

    value=\$(echo "\$output" | grep "Spectra per file" | awk '{print \$5}')

    # Use Python to calculate the nearest power of 2
    nearest_power_of_2=\$(python3 -c "import math; value = int(\${value}); log2 = math.log(value, 2); rounded_log2 = math.ceil(log2) if log2 - int(log2) >= 0.35 else int(log2); nearest_power_of_2 = 2 ** rounded_log2; print(nearest_power_of_2)")

    echo \$nearest_power_of_2
    """
}



process generateDMFiles {
    label "generateDMFiles"
    container "${params.presto_image}"

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

    # Split DM values into multiple files, each containing dm_samp number of lines
    for i in range(0, len(dm_values), dm_sample):
        chunk = dm_values[i:i + dm_sample]
        end_index = min(i + dm_sample, len(dm_values))
        filename = f'dm_{dm_values[i]}_{dm_values[end_index - 1]}.dm'
        np.savetxt(filename, chunk, fmt='%f')
    """
}



process peasoup {
    label 'peasoup'
    container "${params.peasoup_image}"
    publishDir "${params.peasoup.output_path}/", pattern: "*.xml", mode: 'copy'

    input:
    tuple path(fil_file), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size)
    each path(dm_file)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), path(dm_file), path(fil_file, followLinks: false), path("*.xml")

    script:
    """
    #!/bin/bash

    peasoup -p -i ${fil_file} --fft_size ${fft_size} --limit ${params.peasoup.total_cands_limit} -m ${params.peasoup.min_snr} -t ${params.peasoup.ngpus} -n ${params.peasoup.nharmonics} --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file}

    # Rename the output file
    mv **/*.xml ${beam_name}_${dm_file.baseName}_overview.xml
    """
}

process parse_xml {
    label 'parse_xml'
    container "${params.presto_image}"
    publishDir "${params.parse_xml.output_path}/", pattern: "*csv", mode: 'copy'
    publishDir "${params.parse_xml.output_path}/", pattern: "*meta", mode: 'copy'
    stageInMode 'symlink'


    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), path("*.csv"), path("*.meta")
    
    script:
    """
    python3 ${params.parse_xml.script} ${xml_files} --outfile ${utc_start}_${beam_name}_candidates.csv --metafile ${utc_start}_${beam_name}_metafile.meta
    """
}

process splitcands {
    label 'splitcands'
    container "${params.pulsarx_image}"
    publishDir "${params.splitcands.output_path}/", pattern: "*{allCands.txt,candfile,meta.txt}", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), path(candidate_csv), path(metafile)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(dm_file), val(fil_base_name), path(fil_file), path(xml_files), path(candidate_csv), path(metafile), path('*allCands.txt'), path('*candfile'), path('*meta.txt')

    script:
    """
    echo "Running splitcands"
    python3 ${params.splitcands.script} -i ${candidate_csv} -fil ${fil_file} -t pulsarx -n ${params.splitcands.nh} -b ${beam_name}  -threads ${params.splitcands.threads} -p ${params.template_dir} --snr_min ${params.splitcands.snr_min} -ncands ${params.splitcands.ncands} -clfd ${params.splitcands.clfd} -cpn ${params.splitcands.cands_per_node} --beam_id ${beam_id} --utc ${utc_start} --metafile ${metafile}
    """
}

process psrfold {
    label "psrfold"
    container "${params.pulsarx_image}"
    maxForks 100
    publishDir "${params.psrfold.output_path}/${beam_name}/", pattern: "*.png", mode: 'copy'
    publishDir "${params.psrfold.output_path}/${beam_name}/", pattern: "*.ar", mode: 'copy'
    publishDir "${params.psrfold.output_path}/${beam_name}/", pattern: "*.cands", mode: 'copy'
    
    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(fil_base_name), path(fil_file), path(candfile), path(metatext)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path("*.png"), path("*.ar"), path("*.cands")
    
    script:
    """
    python3 ${baseDir}/scripts/pulsarx_fold.py -meta ${metatext} -cands ${candfile}
    """
}

process pics_classifier {
    label "pics_classifier"
    container "/hercules/scratch/fkareem/singularity_img/trapum_pulsarx_fold_docker_20220411.sif"
    publishDir "${params.pulsarx_fold.output_path}", pattern: "*.csv", mode: 'copy'

    input:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path(pngs), path(ars), path(cands)

    output:
    tuple val(cluster),val(beam_name), val(beam_id), val(utc_start), val(fft_size), val(fil_base_name), path(fil_file), path(candfile), path(metatext), path("*.csv")

    script:
    output_csv = "${cluster}_${beam_name}_scored.csv"
    """
    python /hercules/scratch/fkareem/New_Elden_Ring/scripts/pics_classifier_v2.0.py -m ${params.pics_model_dir} -o ${output_csv}
    """
}





    










