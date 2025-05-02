process syncFiles {
    executor 'local'
    maxForks 11
    
    input:
    tuple val(pointing), val(fitsfilepath), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(filename)

    output:
    tuple val(pointing), path("${filename}"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(filename)

    script:
    """
    #!/bin/bash
    mkdir -p ${params.basedir}/${cluster}/Data

    rsync -avz ${params.copy_from_tape.remoteUser}@${params.copy_from_tape.remoteHost}:${fitsfilepath} ${params.basedir}/${cluster}/Data/
     
    # Create symlink in the work directory
    ln -s ${params.basedir}/${cluster}/Data/${filename} ${filename}
    """
}

process readfile {
    label 'readfile'
    container "${params.presto_image}"

    input:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(filename)

    output:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), env(time_per_file), env(tsamp), env(nsamples), env(subintlength)

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

process generateDMFiles {
    label "generateDMFiles"
    container "${params.presto_image}"
    publishDir "${params.basedir}/DMFILES/", pattern: "*.dm", mode: 'copy'

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

process generateRfiFilter {
    label 'generate_rfi_filter'
    container "${params.rfi_mitigation_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/RFIFILTER/", pattern: "*.{png,txt}", mode: 'symlink'

    input:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(time_per_file), val(tsamp), val(nsamples), val(subintlength)

    output:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), env(rfi_filter_string) , val(tsamp), val(nsamples) , val(subintlength), path("*.png"), path("*.txt")

    script:
    def num_intervals = Math.floor(time_per_file.toFloat() / 200) as int
    // def num_intervals = 2
    """
    #!/bin/bash
    export MPLCONFIGDIR=/tmp
    export NUMBA_CACHE_DIR=/tmp
    python3 ${projectDir}/scripts/rfi_mitigation_modified.py ${fits_files} . --target_resolution_ms ${params.generateRfiFilter.target_resolution_ms} --num_intervals ${num_intervals}

    zap_commands=\$(grep -Eo '[0-9.]+ *- *[0-9.]+' combined_frequent_outliers.txt | \\
    awk -F '-' '{gsub(/ /,""); print "zap "\$1" "\$2}' | tr '\\n' ' ')

    rfi_filter_string="kadaneF 8 4 kadaneT 8 4 zdot \${zap_commands}"
    echo "\${rfi_filter_string}" > rfi_filter_string.txt

    mv combined_sk_heatmap_and_histogram.png ${beam_name}_rfi.png
    mv combined_frequent_outliers.txt combined_frequent_outliers_${beam_name}.txt
    mv block_bad_channel_percentages.txt block_bad_channel_percentages_${beam_name}.txt
    """
}

process filtool {
    label 'filtool'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/CLEANEDFIL/", pattern: "*.fil", mode: 'symlink'

    input:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(rfi_filter_string), val(tsamp), val(nsamples) , val(subintlength)
    val threads
    val telescope

    output:
    tuple val(pointing), path("*clean_01.fil"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(tsamp), val(nsamples), val(subintlength)
    
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

process segmented_params {
    label 'segmented_params'
    container "${params.presto_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/SEGPARAMS/", pattern: "*.csv", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(segments)
    
    output:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), env(tsamp), env(nsamples), val(segments), path("*.csv")

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

process birdies {
    label 'birdies'
    container "${params.peasoup_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/BIRDIES/", pattern: "*.{xml,txt}", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample)

    output:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path("birdies.txt")

    script:
    """
    #!/bin/bash
    echo 'Running birdies'
    echo 'What are the parameters?'

    
    peasoup -p -v -i ${fil_file} --fft_size ${fft_size} -m 8.5 -t 1 -n ${params.peasoup.nharmonics} --acc_start 0.0 --acc_end 0.0 --ram_limit_gb 200.0 --dm_start 0.0 --dm_end 0.0  --start_sample ${start_sample} 

    #Rename the output file
    mv **/*.xml ${beam_name}_birdies.xml

    python3 ${projectDir}/scripts/birdies_parser.py --xml_file  *birdies.xml
    """
}

process peasoup {
    label 'peasoup'
    container "${params.peasoup_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/SEARCH/", pattern: "*.xml", mode: 'copy'
    cache 'lenient'
    
    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path(birdies_file)
    each path(dm_file)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), path(dm_file), path(fil_file, followLinks: false), path("*.xml"), path(birdies_file), val(start_sample), val(nsamples)

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

    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${params.peasoup.total_cands_limit} -m ${params.peasoup.min_snr} -t ${params.peasoup.ngpus} -n ${params.peasoup.nharmonics} --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file} \${birdies_string} --start_sample ${start_sample} 

    #Rename the output file
    mv **/*.xml ${beam_name}_${dm_file.baseName}_ck${segments}${segment_id}_overview.xml
    """
}

process parse_xml {
    label 'parse_xml'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/", pattern: "*.{csv,meta,txt,candfile}", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size),val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample), path("filtered_candidates_file.csv"), path("unfiltered_for_folding.csv"), path("*.candfile"), path("*.meta"), path("*allCands.txt")
    
    script:
    def subintlengthstring = params.psrfold.subintlength && params.psrfold.subintlength != "None" ? "-sub ${params.psrfold.subintlength}" : ""
    """ 
    echo "Rerun this"
    #!/bin/bash
    echo "Running parse_xml"
    echo "What are the parameters?"
    python3 ${params.parse_xml.script} -i ${xml_files} --chunk_id ${segments}${segment_id} --fold_technique ${params.psrfold.fold_technique} --nbins_default ${params.psrfold.nbins} --binplan "${params.psrfold.binplan}" ${subintlengthstring} -nsub ${params.psrfold.nsub} -clfd ${params.psrfold.clfd} -b ${beam_name} -b_id ${beam_id} -utc ${utc_start} -threads ${params.psrfold.threads}  --template_dir ${params.psrfold.template_dir} --telescope ${params.telescope} --config_file ${params.parse_xml.config_file} --cdm ${params.psrfold.cdm} --cands_per_node ${params.psrfold.cands_per_node}
    """
}


process psrfold {
    label "psrfold"
    container "${params.pulsarx_image}"
    // maxForks 100
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "*.{png,ar,cands}", mode: 'copy'
    
    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), val(start_sample), path(filtered_candidate_csv), path(candfile), path(metafile)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file, followLinks: false), path(filtered_candidate_csv), path(candfile), path(metafile), path("*.png"), path("*.ar"), path("*.cands")

    script:
    """
    #!/bin/bash
    python3 ${baseDir}/scripts/pulsarx_fold.py -meta ${metafile} -cands ${candfile}

    fold_cands=\$(ls -v *.ar)

    #Run dmffdot if there are missing png files.
    for file in \$fold_cands; do
        png_file="\${file%.ar}.png"
        if [ ! -f "\$png_file" ]; then
            echo "Missing PNG file for \$png_file. Running dmffdot."
            dmffdot --telescope ${params.telescope} -f \$file
        fi
    done

    shopt -s nullglob # Enables automatic removal of non-matching patterns
    for file in *.ar *.png; do
        candfile_no=\$(basename \${file} | cut -d'_' -f2)
        # candidate number is the last number in the file name in the format 0000n.ar
        cand_no=\$(basename \${file} | cut -d'_' -f6 | cut -d'.' -f1)
        new_cand_no=\$(( (10#\${candfile_no} -1) * ${params.psrfold.cands_per_node} + 10#\${cand_no} ))
        new_cand_fmt=\$(printf "%05d" \${new_cand_no})
        name=\$(basename \${file})
        prefix=\${name%_*}
        suffix=\${name##*.}
        new_name="\${prefix}_\${new_cand_fmt}.\${suffix}"
        # Move the file only if the new name is different
        if [[ \${file} != \${new_name} ]]; then
            echo "Renaming \${file} to \${new_name}"
            mv \${file} \${new_name}
        else
            echo "No renaming needed for \${file}"
        fi
    done
    shopt -u nullglob  # Restore default behavior
    """
}

process search_fold_merge {
    label "search_fold_merge"
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "*{.csv,master.cands}", mode: 'copy'

    input:
    tuple val(pointing), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), path(ars), path(cands)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), env(publish_dir), path(ars), path("*master.cands"), path("search_fold_cands*.csv")
    
    script:
    """
    publish_dir="${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING"
    mkdir -p \${publish_dir}

    fold_cands=\$(ls -v *.ar)
    pulsarx_cands_file=\$(ls -v *.cands)
    
    python3 ${baseDir}/scripts/fold_cands_to_csv.py -f \${fold_cands} -c \${pulsarx_cands_file} -x ${filtered_candidate_csv} -o search_fold_cands_${beam_name}_ck${segments}${segment_id}.csv -p \${publish_dir} 
    """
}

process alpha_beta_gamma_test {
    label "alpha_beta_gamma_test"
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/ZERODM/", pattern: "DM0*.png", mode: 'copy'
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/ABG", pattern: "*alpha_beta_gamma.csv", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), val(png_source_dir), path(ars), path(master_cands), path(search_fold_cands_csv)

    output:
    tuple path("*alpha_beta_gamma.csv"), val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(segments), val(segment_id), val(png_source_dir)

    script:
    """
    #!/bin/bash
    publish_dir="${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/ABG"
    mkdir -p \${publish_dir}
    python3 ${baseDir}/scripts/calculate_alpha_beta_gamma_dmffdot.py -i ${search_fold_cands_csv} -o ${cluster}_${beam_name}_ck${segments}${segment_id}_alpha_beta_gamma.csv -t ${params.alpha_beta_gamma.snr_min} -p \${publish_dir} -s ${png_source_dir} -c
    """
}

process pics_classifier {
    label "pics_classifier"
    container "${params.pics_classifier_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/CLASSIFICATION/", pattern: "*scored.csv", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), val(png_source_dir), path(ars), path(master_cands), path(search_fold_cands_csv)

    output:
    tuple path("*scored.csv"), val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(segments), val(segment_id), val(png_source_dir)

    script:
    output_csv = "${cluster}_${beam_name}_ck${segments}${segment_id}_scored.csv"
    """
    python2 ${baseDir}/scripts/pics_classifier_multiple_models.py -m ${params.pics_model_dir} -o ${output_csv}
    """
}


process create_candyjar_tarball {
    executor 'local'
    container "${params.pulsarx_image}"
    //publishDir "${params.publish_dir_prefix}/${target}/CANDIDATE_TARBALLS", pattern: "*.tar.gz", mode: 'copy'

    input:
    tuple path(candidate_results_file), val(output_tarball_name)

    output:
    stdout

    script:
    """
    #!/bin/bash
    header="pointing,target,beam,beam_id,utc_start,ra,dec,segments,segment_id,fold_cands_filepath,alpha_beta_file,pics_file"

    # Extract basename without extension and append "_header.csv"
    candidate_results_file_with_header="\$(basename ${candidate_results_file} .csv)_header.csv"
    echo "\$header" > "\$candidate_results_file_with_header"
    cat "${candidate_results_file}" >> "\$candidate_results_file_with_header"

    python ${baseDir}/scripts/create_candyjar_tarball.py -i \$candidate_results_file_with_header -o ${output_tarball_name} --verbose --npointings 0 -m ${params.metafile_source_path} -d ${params.basedir} --threshold ${params.alpha_beta_gamma.threshold}
    """
}
