process syncFiles {
    executor 'local'
    maxForks 11
    
    input:
    tuple val(pointing), val(fitsfilepath), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(filename)

    output:
    tuple val(pointing), path("${filename}"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(filename)

    script:
    """
    #!/bin/bash
    mkdir -p ${params.basedir}/${cluster}/Data

    rsync -avz ${params.copy_from_tape.remoteUser}@${params.copy_from_tape.remoteHost}:${fitsfilepath} ${params.basedir}/${cluster}/Data/
     
    # Create symlink in the work directory
    ln -s ${params.basedir}/${cluster}/Data/${filename} ${filename}
    """
}

process dada_to_fits {
    label 'dada_to_fits'
    container "${params.edd_pulsar_image}"
    
    input:
    tuple val(pointing), val(dada_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm)

    output:
    tuple val(pointing), path(filename), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm)

    script:
    filename = "${cluster}_${utc_start}_${beam_name}_cdm_${cdm}.sf"
    publish_dir = "${params.basedir}/${cluster}/${beam_name}/FITS/"
    """
    #!/bin/bash
    set -euo pipefail
    work_dir=\$(pwd)
    mkdir -p ${publish_dir}

    # Base digifits command
    base_cmd=(
        digifits -v -cuda 0 -r -u
        -b "${params.dada.bits}"
        -p "${params.dada.npol}"
        -nsblk "${params.dada.nsblk}"
        -F "${params.dada.nchan}:D"
        -x "${params.dada.nfft}"
        -t "${params.dada.tsamp}"
        -do_dedisp -D "${cdm}"
        -o "${filename}"
    )
    
    # Run digifits with retry logic
    run_digifits() {
        local uval=\$1
        local extra_args=()
        [[ \$uval -ne 0 ]] && extra_args+=(-U "\$uval")
        
        echo "Running: \${base_cmd[@]} \${extra_args[@]} ${dada_files}"
        "\${base_cmd[@]}" "\${extra_args[@]}" ${dada_files} > digifits.log 2>&1
    }

    # First try with default settings
    if ! run_digifits 0; then
        if grep -q 'insufficient RAM: limit=' digifits.log && \
           grep -q 'a minimum of "-U' digifits.log
        then
            # FIXED: Properly escaped regex for Groovy
            recommended_u=\$(grep 'a minimum of' digifits.log | sed -E 's/.*-U ?"?([0-9]+).*/\\1/' | head -1)
            
            if [[ -n "\$recommended_u" && "\$recommended_u" =~ ^[0-9]+\$ ]]; then
                echo "Retrying with -U \$recommended_u"
                run_digifits "\$recommended_u" || {
                    echo "digifits failed after retry"
                    cat digifits.log
                    exit 1
                }
            else
                echo "ERROR: Failed to parse recommended U-value"
                echo "Tried to parse from:"
                grep 'a minimum of' digifits.log
                exit 1
            fi
        else
            echo "digifits failed with unexpected error"
            cat digifits.log
            exit 1
        fi
    fi

    # file check and renaming logic

    output_found=false
    for f in *.sf; do
        if [[ -f "\$f" && "\$f" != "${filename}" ]]; then
            echo "Renaming output file from \$f to ${filename}"
            mv "\$f" "${filename}"
            output_found=true
            break
        elif [[ -f "\$f" ]]; then
            output_found=true
            break
        fi
    done

    if ! \$output_found; then
        echo "ERROR: No .sf output file found in:"
        ls -la
        exit 1
    fi

    echo "digifits completed successfully."

    # Move the output file to the publish directory
    mv "${filename}" ${publish_dir}/${filename}

    ln -s ${publish_dir}/${filename} ${filename}
    """
}

process readfile {
    label 'readfile'
    container "${params.presto_image}"

    input:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(filename)

    output:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), env(time_per_file), env(tsamp), env(nsamples), env(subintlength)

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
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(time_per_file), val(tsamp), val(nsamples), val(subintlength)

    output:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), env(rfi_filter_string) , val(tsamp), val(nsamples) , val(subintlength), path("*.png"), path("*.txt")

    script:
    def num_intervals = Math.floor(time_per_file.toFloat()) as int
    """
    #!/bin/bash
    export MPLCONFIGDIR=/tmp
    export NUMBA_CACHE_DIR=/tmp
    python3 ${projectDir}/scripts/rfi_mitigation_modified.py ${fits_files} . --target_resolution_ms ${params.generateRfiFilter.target_resolution_ms} --num_intervals ${num_intervals}

    zap_commands=\$(grep -Eo '[0-9.]+ *- *[0-9.]+' combined_frequent_outliers.txt | \\
    awk -F '-' '{gsub(/ /,""); print "zap "\$1" "\$2}' | tr '\\n' ' ')

    rfi_filter_string="${params.generateRfiFilter.default_flag} \${zap_commands}"
    echo "\${rfi_filter_string}" > rfi_filter_string.txt

    mv combined_sk_heatmap_and_histogram.png ${beam_name}_rfi.png
    mv combined_frequent_outliers.txt combined_frequent_outliers_${beam_name}.txt
    mv block_bad_channel_percentages.txt block_bad_channel_percentages_${beam_name}.txt
    """
}

process filtool {
    label 'filtool'
    container "${params.pulsarx_image}"
    // publishDir "${params.basedir}/${cluster}/CLEANEDFIL/", pattern: "*.fil", mode: 'symlink'
    // removed publishDir to avoid symlinks, now using publishDir in the script

    input:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(rfi_filter_string), val(tsamp), val(nsamples) , val(subintlength)
    val threads
    val telescope

    output:
    tuple val(pointing), path("*clean_01.fil"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(subintlength)
    
    script:
    def outputFile = "${cluster.trim()}_${utc_start.trim()}_${beam_name.trim()}_cdm_${cdm}_clean"
    def source_name = "${cluster.trim()}"

    // Prepare the rfi_filter option
    def zaplist = ''
    if (rfi_filter_string) {
        zaplist = "-z ${rfi_filter_string}"
    }
    """
    #!/bin/bash
    workdir=\$(pwd)
    echo "Working directory: \${workdir}"
    publish_dir="${params.basedir}/${cluster}/${beam_name}/CLEANEDFIL"
    mkdir -p \${publish_dir}
    cd \${publish_dir}

    # Get the first file from the inputFile string
    # This is used to determine the file extension
    first_file=\$(echo ${fits_files} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    flip_flag=""
    if [[ ${params.filtool.flip} == true ]]; then
        flip_flag="--flip"
    fi

    if [[ "\${file_extension}" == "fits" || "\${file_extension}" == "sf" || "\${file_extension}" == "rf" ]]; then
        echo "Running: filtool --psrfits --scloffs \${flip_flag} --td ${params.filtool.td} --fd ${params.filtool.fd} -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f \${workdir}/${fits_files} -s ${source_name}"
        filtool --psrfits --scloffs \${flip_flag} --td ${params.filtool.td} --fd ${params.filtool.fd} -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f \${workdir}/${fits_files} -s ${source_name}
    else 
        echo "Running: filtool \${flip_flag} --td ${params.filtool.td} --fd ${params.filtool.fd} -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f \${workdir}/${fits_files} -s ${source_name}"
        filtool \${flip_flag} --td ${params.filtool.td} --fd ${params.filtool.fd} -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f \${workdir}/${fits_files} -s ${source_name}
    fi

    # create a symlink to the cleaned file in the work directory
    cd \${workdir}
    ln -s \${publish_dir}/${outputFile}_01.fil ${outputFile}_01.fil
    """
}

process merge_filterbanks {
    label 'merge_filterbanks'
    container "${params.filtools_sig_image}"
    // publishDir "${params.basedir}/${cluster}/${beam_name}/MERGED/", pattern: "*.fil", mode: 'copy'

    input:
    tuple val(pointing), val(cluster), val(utc), val(ra), val(dec), val(cdm), val(group_label), path(fil_files)

    output:
    tuple val(pointing), path("*stacked.fil"), val(cluster), env(beam_name), val(group_label), val(utc), val(ra), val(dec), val(cdm)

    script:
    def outputFile = "${cluster}.${utc}_cfbf${group_label}_stacked.fil"
    def filelist = fil_files.collect { it.getName() }.join(' ')
    def publishDir = "${params.basedir}/${cluster}/${beam_name}/MERGED"
    """
    #!/bin/bash
    workdir=\$(pwd)

    beam_name="cfbf${group_label}"

    mkdir -p ${publishDir}
    cd ${publishDir}
    echo "Merging files for cdm = ${cdm}, group_label = ${group_label}"
    python ${baseDir}/scripts/freq_stack.py -o ${outputFile} $filelist

    echo "Merged file created: ${outputFile}"
    cd \${workdir}
    ln -s ${publishDir}/${outputFile} ${outputFile}
    """

}
process segmented_params {
    label 'segmented_params'
    container "${params.presto_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/SEGPARAMS/", pattern: "*.csv", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec),val(cdm), val(segments)
    
    output:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), env(tsamp), env(nsamples), val(segments), path("*.csv")

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
    output_file="${beam_name}_cdm_${cdm}_segments_${segments}_params.csv"
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
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample)

    output:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path("birdies.txt")

    script:
    """
    #!/bin/bash
    echo 'Running birdies'
    echo 'What are the parameters?'

    
    peasoup -p -v -i ${fil_file} --cdm ${cdm} --fft_size ${fft_size} -m ${params.peasoup.birdies_min_snr} -t 1 -n ${params.peasoup.nharmonics} --acc_start 0.0 --acc_end 0.0 --ram_limit_gb 200.0 --dm_start 0.0 --dm_end 0.0  --start_sample ${start_sample} 

    mv **/*.xml ${beam_name}_cdm_${cdm}_birdies.xml

    python3 ${projectDir}/scripts/birdies_parser.py --xml_file  *birdies.xml
    """
}


process generateDMFiles {
    label "generateDMFiles"
    container "${params.presto_image}"
    publishDir "${params.basedir}/DMFILES/", pattern: "*.dm", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path(birdies_file)

    output:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path(birdies_file), path("*.dm")

    script:
    """
    #!/usr/bin/env python3
    import numpy as np

    # Generate the DM file
    dm_start = ${cdm} + ${params.ddplan.dm_start}
    dm_end = ${cdm} + ${params.ddplan.dm_end}
    dm_step = ${params.ddplan.dm_step}
    dm_sample = ${params.ddplan.dm_sample}

    # Create DM values with a step of dm_step
    dm_values = np.round(np.arange(dm_start, dm_end, dm_step), 3)

    # Split DM values into multiple files, each containing dm_sample number of lines
    for i in range(0, len(dm_values), dm_sample):
        chunk = dm_values[i:i + dm_sample]
        end_index = min(i + dm_sample, len(dm_values))
        filename = f'cdm_${cdm}_dm_{dm_values[i]}_{dm_values[end_index - 1]}.dm'
        np.savetxt(filename, chunk, fmt='%f')
    """
}

process peasoup {
    label 'peasoup'
    container "${params.peasoup_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/SEARCH/", pattern: "*.xml", mode: 'copy'
    cache 'lenient'
    
    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path(birdies_file), path(dm_file)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), path(dm_file), path(fil_file, followLinks: false), path("*.xml"), path(birdies_file), val(start_sample), val(nsamples)

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

    peasoup -i ${fil_file} --cdm ${cdm} --fft_size ${fft_size} --limit ${params.peasoup.total_cands_limit} -m ${params.peasoup.min_snr} -t ${params.peasoup.ngpus} -n ${params.peasoup.nharmonics} --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file} \${birdies_string} --start_sample ${start_sample} 

    #Rename the output file
    mv **/*.xml ${beam_name}_cdm_${cdm}_${dm_file.baseName}_ck${segments}${segment_id}_overview.xml
    """
}

process parse_xml {
    label 'parse_xml'
    container "${params.pulsarx_image}"
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/", pattern: "*.{csv,meta,txt,candfile}", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size),val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample), path("filtered_candidates_file.csv"), path("unfiltered_for_folding.csv"), path("*.candfile"), path("*.meta"), path("*allCands.txt")
    
    script:
    def subintlengthstring = params.psrfold.subintlength && params.psrfold.subintlength != "None" ? "-sub ${params.psrfold.subintlength}" : ""
    """ 
    #!/bin/bash
    python3 ${params.parse_xml.script} -i ${xml_files} --chunk_id ${segments}${segment_id} --fold_technique ${params.psrfold.fold_technique} --nbins_default ${params.psrfold.nbins} --binplan "${params.psrfold.binplan}" ${subintlengthstring} -nsub ${params.psrfold.nsub} -clfd ${params.psrfold.clfd} -b ${beam_name} -b_id ${beam_id} -utc ${utc_start} -threads ${params.psrfold.threads}  --template_dir ${params.psrfold.template_dir} --telescope ${params.telescope} --config_file ${params.parse_xml.config_file} --cdm ${cdm} --cands_per_node ${params.psrfold.cands_per_node}
    """
}


process psrfold {
    label "psrfold"
    container "${params.pulsarx_image}"
    // maxForks 100
    publishDir "${params.basedir}/${cluster}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/", pattern: "*.{png,ar,cands}", mode: 'copy'
    
    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), val(start_sample), path(filtered_candidate_csv), path(candfile), path(metafile)

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


// other pipelines 
process parfold {
    label "parfold"
    container "${params.pulsarx_image}"
    maxForks 100
    publishDir "${params.parfold.output_path}/", pattern: "*.png", mode: 'copy'
    publishDir "${params.parfold.output_path}/", pattern: "*.ar", mode: 'copy'
    publishDir "${params.parfold.output_path}/", pattern: "*.cands", mode: 'copy'
    
    input:
    tuple val(pointing), path(fil_file), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(tsamp), val(nsamples), val(subintlength)
    each path(parfile_channel)

    output:
    tuple path("*.png"), path("*.ar"), path("*.cands")
    
    script:
    def Outname = "${beam_name}_${parfile_channel.getName().replace(".par", "")}"
    """
    #!/bin/bash

    psrfold_fil --plotx --nosearch -v -t ${params.parfold.threads} --parfile ${parfile_channel} -n ${params.parfold.nsub} -b ${params.parfold.nbins} --nbinplan ${params.parfold.binplan} --template ${params.template_dir}/Effelsberg_${beam_id}.template --clfd ${params.parfold.clfd} -L ${subintlength} -f ${fil_file} -o ${Outname}
    """
}

// process for candypolice

process extract_candidates {
    container "${params.pulsarx_image}"
    publishDir "${params.candypolice.output_dir}", pattern: "*{txt,candfile}", mode: 'copy'
    input:
    path csv_file

    output:
    tuple path("candidate_details.txt"), path("*candfile")

    script:
    """
    #!/usr/bin/env python3
    import csv, os

    # columns to extract, in order
    fields = [
        'pointing_id',
        'beam_name',
        'beam_id',
        'source_name',
        'segment_id',
        'f0_opt',
        'f1_opt',
        'p0_new',
        'p1_new',
        'acc_opt',
        'dm_opt',
        'sn_fold',
        'pepoch',
        'fold_cands_filename',
        'classification'
    ]

    ## writing candfiles for each beam.
    groups = {}

    # Read + filter + write the big CSV *and* stash rows for the .candfile

    with open('${csv_file}', newline='') as csvf, \
         open('candidate_details.txt', 'w', newline='') as outf:

        reader = csv.DictReader(csvf)
        writer = csv.writer(outf, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # write header for candidate_details.txt
        writer.writerow(fields)

        for row in reader:
            if row['classification'] in ('T1_CAND', 'T2_CAND'):
                writer.writerow([row[f] for f in fields])

                # group it to later write .candfile
                key = (row['source_name'], row['beam_name'], row['pepoch'], row['segment_id'])
                groups.setdefault(key, []).append(row)

    # writing candfiles per group for pulsarx
    for (source, beam, pepoch, seg_id), rows in groups.items():
        fname = f"{source}_{beam}_{seg_id}_{pepoch}.candfile"
        with open(fname, 'w') as out:
            out.write('#id dm acc F0 F1 F2 S/N\\n')
            for r in rows:
                vals = [
                    r['id'],       # id
                    r['dm_opt'],   # dm
                    r['acc_opt'],  # acc
                    r['f0_opt'],   # F0
                    r['f1_opt'],   # F1
                    '0',           # F2 placeholder
                    r['sn_fold']   # S/N
                ]
                out.write(' '.join(vals) + '\\n')

    """
}


process candypolice_pulsarx {
    label "candypolice_pulsarx"
    container "${params.pulsarx_image}"
    publishDir "${params.candypolice.output_path}", pattern: "**/*.{ar,png,cands}", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file),val(cluster),val(beam_name),val(beam_id),val(utc_start),val(ra),val(dec),val(ts),val(ns),val(si)
    each path(candfile)

    output:

    script:
    """
    #!/bin/bash
    filename=\$(basename ${candfile})
    base=\${filename%.candfile}

    IFS='_' read -r source beam pepoch seg_id <<< "\$base"

    Outname="\${base}_${beam_id}"
    
    i="\${seg_id:0:1}"
    j="\${seg_id:1:1}"

    end_frac=\$(awk -v i=\$i -v j=\$j 'BEGIN { printf "%.2f", (1/i)*(j+1) }')
    start_frac=\$(awk -v i=\$i -v ef=\$end_frac 'BEGIN { printf "%.2f", ef - 1/i }')
    echo "seg_id = \$seg_id, i=\$i, j=\$j"
    echo "  end_frac = \$end_frac"
    echo "  start_frac = \$start_frac"

    psrfold_fil --plotx --nosearch -v -t ${params.candypolice.pulsarx_threads} --candfile ${candfile} -n ${params.candypolice.nsub} -b ${params.candypolice.nbins} --nbinplan ${params.candypolice.binplan} --template ${params.template_dir}/Effelsberg_${beam_id}.template --clfd ${params.candypolice.clfd} -L ${si} -f ${fil_file} -o \${Outname}
    """
}

process candypolice_presto {
    label "candypolice_presto"
    container "${params.presto_image}"
    scratch true
    publishDir "${params.candypolice.output_path}", pattern: "**/*.{ps,png,pfd}", mode: 'copy'
    
    input:
    each candidate_details_line
    path(input_file)

    output:
    path "**/*.{ps,png,pfd}"

    script:
    """
    #!/bin/env python3
    import os
    import subprocess

    # Parse candidate details from the input line
    line = "${candidate_details_line}".strip()
    parts = line.split('|')
    candidate = parts[0]
    candidate_details = {kv.split(':')[0]: kv.split(':')[1] for kv in parts[1:]}


    # Create a directory for the candidate
    
    # Expects input png_path in format: /fpra/timing/01/fazal/Scripts/Elden_Ring/work/9e/8421a89c0992bc364eec620d298e3a/TERZAN5_Band3_dm_file_48.0_48.9_1/60403.1458332005_cbf00000_00096.png

    #expects input file in format: /fpra/timing/01/fazal/Eff_Data_Proc/NGC6544/Filtool/NGC6544_Band3_01.fil

    
    name = os.path.basename(candidate_details['png_path']).replace('.png', '').replace('/', '_').split('_')[-1]
    bandname = os.path.basename("${input_file}").split('.')[0]
    detection_band = os.path.dirname(candidate_details['png_path']).split('/')[-1]
    dir_name = f"{detection_band}_{name}"

    os.makedirs(dir_name, exist_ok=True)
    
    # Define commands
    commands = [
        f"prepfold -topo -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name} ${input_file}",
        f"prepfold -topo -dm 0 -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_0dm ${input_file}",
        f"prepfold -topo -fine -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_fine ${input_file}",
        f"prepfold -topo -nosearch -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_nosearch ${input_file}",
        f"prepfold -topo -pfact 2 -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_pfact2 ${input_file}",
        f"prepfold -topo -ffact 2 -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_ffact2 ${input_file}"
    ]

    # Execute all commands in parallel
    processes = [subprocess.Popen(cmd, shell=True) for cmd in commands]
    for p in processes:
        p.wait()
    """
}
