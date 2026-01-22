process syncFiles {
    executor 'local'
    maxForks 11
    tag "${cluster}_${filename}"

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

    //TODO this should be universal
    container "${params.edd_pulsar_image}"
    tag "${cluster}_${beam_name}_cdm_${cdm}"
    cache 'lenient'

    input:
    tuple val(pointing), val(dada_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm)

    output:
    tuple val(pointing), path(filename), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(filename)

    script:
    filename = "${cluster}_${utc_start}_${beam_name}_cdm_${cdm}.sf"
    // Use a shared location independent of runID for caching
    publish_dir = "${params.basedir}/shared_cache/${cluster}/FITS/"
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
    mv "${filename}" "${publish_dir}/${filename}"

    # Create symlink in runID-specific directory for easy access
    if [[ -n "${params.runID}" ]]; then
        runid_dir="${params.basedir}/${params.runID}/${beam_name}/FITS"
        mkdir -p "\${runid_dir}"
        ln -sf "${publish_dir}/${filename}" "\${runid_dir}/${filename}"

        # Verify symlink was created
        if [[ ! -L "\${runid_dir}/${filename}" ]]; then
            echo "WARNING: Failed to create symlink at \${runid_dir}/${filename}"
        fi
    fi

    # Create symlink in work directory
    ln -s "${publish_dir}/${filename}" "${filename}"
    """
}

process readfile {
    label 'readfile'
    container "${params.presto_image}"
    tag "${cluster}_${beam_name}_cdm_${cdm}"
    cache 'lenient'

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
    tag "${cluster}_${beam_name}_cdm_${cdm}"
    cache 'lenient'
    // Don't use runID in publishDir for caching - use shared cache location
    publishDir "${params.basedir}/shared_cache/${cluster}/${beam_name}/RFIFILTER/", pattern: "*.{png,txt}", mode: 'copy'

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
    if awk 'BEGIN{exit !((${cdm} <= 60))}'; then
      echo "not using zdot for cdm = ${cdm}"
      default_flag="${params.generateRfiFilter.default_flag}"
    else
      echo "cdm = ${cdm}; using zdot"
      default_flag="${params.generateRfiFilter.default_flag} zdot"
    fi
    rfi_filter_string="\${default_flag} \${zap_commands}"
    echo "\${rfi_filter_string}" > rfi_filter_string_cdm_${cdm}.txt

    mv combined_sk_heatmap_and_histogram.png ${beam_name}_cdm_${cdm}_rfi.png
    mv combined_frequent_outliers.txt combined_frequent_outliers_${beam_name}_${cdm}.txt
    mv block_bad_channel_percentages.txt block_bad_channel_percentages_${beam_name}_${cdm}.txt

    # Also create symlinks in runID-specific directory for easy access
    if [[ -n "${params.runID}" ]]; then
        runid_dir="${params.basedir}/${params.runID}/${beam_name}/RFIFILTER"
        mkdir -p "\${runid_dir}"

        ln -sf "${params.basedir}/shared_cache/${cluster}/${beam_name}/RFIFILTER/${beam_name}_cdm_${cdm}_rfi.png" "\${runid_dir}/"
        ln -sf "${params.basedir}/shared_cache/${cluster}/${beam_name}/RFIFILTER/combined_frequent_outliers_${beam_name}_${cdm}.txt" "\${runid_dir}/"
        ln -sf "${params.basedir}/shared_cache/${cluster}/${beam_name}/RFIFILTER/block_bad_channel_percentages_${beam_name}_${cdm}.txt" "\${runid_dir}/"

        # Verify symlinks were created
        if [[ ! -L "\${runid_dir}/${beam_name}_cdm_${cdm}_rfi.png" ]]; then
            echo "WARNING: Failed to create symlink for RFI filter plots"
        fi
    fi
    """
}

process filtool {
    label 'filtool'
    container "${params.pulsarx_image}"
    tag "${cluster}_${beam_name}_cdm_${cdm}"
    cache 'lenient'
    // publishDir "${params.basedir}/${cluster}/CLEANEDFIL/", pattern: "*.fil", mode: 'symlink'
    // removed publishDir to avoid symlinks, now using publishDir in the script

    input:
    tuple val(pointing), path(fits_files), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(rfi_filter_string), val(tsamp), val(nsamples) , val(subintlength)
    val threads
    val telescope

    output:
    tuple val(pointing), path("*clean_01.fil"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm)

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
    # Use shared cache location independent of runID
    publish_dir="${params.basedir}/shared_cache/${cluster}/${beam_name}/CLEANEDFIL"
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

    # Also create symlink in runID-specific directory for easy access
    if [[ -n "${params.runID}" ]]; then
        runid_dir="${params.basedir}/${params.runID}/${beam_name}/CLEANEDFIL"
        mkdir -p "\${runid_dir}"
        ln -sf "\${publish_dir}/${outputFile}_01.fil" "\${runid_dir}/${outputFile}_01.fil"

        # Verify symlink was created
        if [[ ! -L "\${runid_dir}/${outputFile}_01.fil" ]]; then
            echo "WARNING: Failed to create symlink at \${runid_dir}/${outputFile}_01.fil"
        fi
    fi

    # create a symlink to the cleaned file in the work directory
    cd \${workdir}
    ln -s "\${publish_dir}/${outputFile}_01.fil" "${outputFile}_01.fil"
    """
}

process split_filterbank {
    label 'split_filterbank'
    container "${params.filtools_sig_image}"
    tag "${cluster}_${beam_name}_cdm_${cdm}"
    cache 'lenient'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm)

    output:
    tuple val(pointing), path("*cut.fil"), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm)

    script:
    """
    #!/bin/bash
    outputFile="${cluster.trim()}_${utc_start.trim()}_${beam_name.trim()}_cdm_${cdm}_clean"
    # Use shared cache location independent of runID
    publish_dir="${params.basedir}/shared_cache/${cluster}/${beam_name}/CLEANEDFIL"
    mkdir -p "\${publish_dir}"

    if [[ ${params.split_fil} == true ]]; then
      if [[ ${beam_id} == 1 ]]; then
        echo "Splitting band 1"
        python ${baseDir}/scripts/cut_filterbank.py -i ${fil_file} -c ${params.split_freq} -l \${outputFile}_low.fil -u \${outputFile}_cut.fil
        cp \${outputFile}_cut.fil "\${publish_dir}/"
      else
        echo "File good. Skipping"
        mv ${fil_file} \${outputFile}_cut.fil
      fi
    fi

    # Create symlink in runID-specific directory for easy access
    if [[ -n "${params.runID}" && -f "\${outputFile}_cut.fil" ]]; then
        runid_dir="${params.basedir}/${params.runID}/${beam_name}/CLEANEDFIL"
        mkdir -p "\${runid_dir}"
        ln -sf "\${publish_dir}/\${outputFile}_cut.fil" "\${runid_dir}/\${outputFile}_cut.fil"

        # Verify symlink was created
        if [[ ! -L "\${runid_dir}/\${outputFile}_cut.fil" ]]; then
            echo "WARNING: Failed to create symlink at \${runid_dir}/\${outputFile}_cut.fil"
        fi
    fi
    """
  }


process merge_filterbanks {
    label 'merge_filterbanks'
    container "${params.filtools_sig_image}"
    tag "${cluster}_cfbf${group_label}_cdm_${cdm}"
    cache 'lenient'
    // publishDir "${params.basedir}/${cluster}/${beam_name}/MERGED/", pattern: "*.fil", mode: 'copy'

    input:
    tuple val(pointing), val(cluster), val(utc), val(ra), val(dec), val(cdm), val(group_label), val(fil_files)

    output:
    tuple val(pointing), path("*stacked.fil"), val(cluster), env(beam_name), val(group_label), val(utc), val(ra), val(dec), val(cdm)

    script:
    def beam_name="cfbf${group_label}"
    def outputFile = "${cluster}.${utc}_cfbf${group_label}_cdm_${cdm}_stacked.fil"
    def filelist = fil_files.collect { it }.join(' ')
    // Use shared cache location independent of runID
    def publishDir = "${params.basedir}/shared_cache/${cluster}/${beam_name}/MERGED"
    """
    #!/bin/bash
    beam_name="cfbf${group_label}"
    workdir=\$(pwd)
    mkdir -p ${publishDir}
    cd ${publishDir}
    echo "Merging files for cdm = ${cdm}, group_label = ${group_label}"
    echo "python ${baseDir}/scripts/freq_stack.py -o ${outputFile} $filelist"
    python ${baseDir}/scripts/freq_stack.py -o ${outputFile} ${filelist}

    echo "Merged file created: ${outputFile}"

    # Also create symlink in runID-specific directory for easy access
    if [[ -n "${params.runID}" ]]; then
        runid_dir="${params.basedir}/${params.runID}/cfbf${group_label}/MERGED"
        mkdir -p "\${runid_dir}"
        ln -sf "${publishDir}/${outputFile}" "\${runid_dir}/${outputFile}"

        # Verify symlink was created
        if [[ ! -L "\${runid_dir}/${outputFile}" ]]; then
            echo "WARNING: Failed to create symlink at \${runid_dir}/${outputFile}"
        fi
    fi

    cd \${workdir}
    ln -s "${publishDir}/${outputFile}" "${outputFile}"
    """
}

process segmented_params {
    label 'segmented_params'
    container "${params.presto_image}"
    tag "${cluster}_${beam_name}_cdm_${cdm}_seg${segments}"
    cache 'lenient'
    // Don't use publishDir directive - handle publishing in script for shared cache approach

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

    # Use shared cache location based on segments (not runID)
    publish_dir="${params.basedir}/shared_cache/${cluster}/${beam_name}/segment_${segments}/SEGPARAMS"
    mkdir -p "\${publish_dir}"

    # Generate CSV in shared cache location
    cd "\${publish_dir}"
    echo "i,fft_size,start_sample,nsamples_per_segment" > "\${output_file}"
    start_sample=0
    for ((i=0; i<=${segments}-1; i++ )); do
        echo "\$i,\$fft_size,\$start_sample,\$nsamples_per_segment" >> "\${output_file}"
        start_sample=\$((start_sample + nsamples_per_segment))
    done

    # Create symlink in runID-specific directory for easy access
    if [[ -n "${params.runID}" ]]; then
        runid_dir="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/SEGPARAMS"
        mkdir -p "\${runid_dir}"
        ln -sf "\${publish_dir}/\${output_file}" "\${runid_dir}/\${output_file}"

        # Verify symlink was created
        if [[ ! -L "\${runid_dir}/\${output_file}" ]]; then
            echo "WARNING: Failed to create symlink at \${runid_dir}/\${output_file}"
        fi
    fi

    # Create symlink in work directory
    workdir="\$(pwd | sed 's|/shared_cache/.*||')/\$(basename \$(pwd))"
    cd "\${workdir}" 2>/dev/null || cd -
    ln -s "\${publish_dir}/\${output_file}" "\${output_file}"
    """
}

process birdies {
    label 'birdies'
    container "${params.peasoup_image}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}"
    cache 'lenient'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/BIRDIES/", pattern: "*.{xml,txt}", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample)

    output:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path("*birdies.txt"), path("*birdies.xml")

    script:
    """
    #!/bin/bash
    echo 'Running birdies'
    echo 'What are the parameters?'

    
    peasoup -p -v -i ${fil_file} --cdm ${cdm} --fft_size ${fft_size} -m ${params.peasoup.birdies_min_snr} -t 1 -n ${params.peasoup.birdies_harmonics} --acc_start 0.0 --acc_end 0.0 --ram_limit_gb ${params.peasoup.birdies_ram_limit_gb} --dm_start 0.0 --dm_end 0.0  --start_sample ${start_sample} --max_freq ${params.peasoup.birdies_max_freq}

    mv **/*.xml ${beam_name}_cdm_${cdm}_birdies.xml

    python3 ${projectDir}/scripts/birdies_parser.py --default_birdies_file ${params.peasoup.default_birdies} --xml_file  *birdies.xml

    mv birdies.txt ${beam_name}_cdm_${cdm}_birdies.txt
    """
}


process generateDMFiles {
    label "generateDMFiles"
    container "${params.presto_image}"
    tag "cdm_${cdm}_seg${segments}${segment_id}"
    cache 'lenient'
    publishDir "${params.basedir}/${params.runID}/DMFILES/", pattern: "*.dm", mode: 'copy'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path(birdies_file), path(birdies_xml)

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

        if ${params.ddplan.use_zero_dm ? 'True' : 'False'}:
            chunk = np.concatenate(([0.0], chunk))
            chunk = np.unique(chunk)   # remove duplicates (in case 0.0 already present)
            chunk.sort()

        filename = f'cdm_${cdm}_dm_{dm_values[i]}_{dm_values[end_index - 1]}.dm'
        np.savetxt(filename, chunk, fmt='%f')
    """
}

process peasoup {
    label 'peasoup'
    container "${params.peasoup_image}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}_${dm_file.baseName}"
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/SEARCH/", pattern: "*.xml", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/TIMESERIES/", pattern: "timeseries_dump/*.{dat,inf}", mode: 'copy', enabled: params.peasoup?.dump_timeseries ?: false, saveAs: { filename -> filename.split('/').last() }
    cache 'lenient'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(tsamp), val(nsamples), val(segments), val(segment_id), val(fft_size), val(start_sample), path(birdies_file), path(dm_file)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), path(dm_file), path(fil_file, followLinks: false), path("*.xml"), path(birdies_file), val(start_sample), val(nsamples), emit: search_results
    tuple path("timeseries_dump/*.dat"), path("timeseries_dump/*.inf"), emit: timeseries_data, optional: true

    script:
    def dump_timeseries = params.peasoup?.dump_timeseries ?: false
    """
    #!/bin/bash

    # check if birdies files is empty
    if [ ! -s ${birdies_file} ]; then
        echo "Birdies file is empty"
        birdies_string=""
    else
        birdies_string="--zapfile ${birdies_file}"
    fi

    # Run normal peasoup search
    peasoup -i ${fil_file} --cdm ${cdm} --fft_size ${fft_size} --limit ${params.peasoup.total_cands_limit} -m ${params.peasoup.min_snr} -t ${params.peasoup.ngpus} -n ${params.peasoup.nharmonics} --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file} \${birdies_string} --start_sample ${start_sample}

    # Rename the output file
    mv **/*.xml ${beam_name}_cdm_${cdm}_${dm_file.baseName}_ck${segments}${segment_id}_overview.xml

    # If dump_timeseries is enabled, run a second peasoup to dump .dat/.inf files
    if [ "${dump_timeseries}" == "true" ]; then
        echo "Dumping time series with peasoup -d timeseries_dump"
        mkdir -p timeseries_dump
        peasoup -i ${fil_file} --cdm ${cdm} --fft_size ${fft_size} -t ${params.peasoup.ngpus} --ram_limit_gb ${params.peasoup.ram_limit_gb} --dm_file ${dm_file} -d timeseries_dump --start_sample ${start_sample}
    fi
    """
}

process parse_xml {
    label 'parse_xml'
    container "${params.rusty_candypicker}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}"
    cache 'lenient'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/", pattern: "*.{csv,meta,txt,candfile}", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size),val(segments), val(segment_id), val(dm_file), val(fil_file_base), path(fil_file), path(xml_files), val(start_sample), path("filtered_candidates_file*.csv"), path("unfiltered_for_folding*.csv"), path("*.candfile"), path("*.meta"), path("*allCands.txt")
    
    script:
    def subintlengthstring = params.psrfold.subintlength && params.psrfold.subintlength != "None" ? "-sub ${params.psrfold.subintlength}" : ""
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob

    echo "running parse xml"

    # prefer overview XMLs; fallback to any *.xml
    xmls=( *overview.xml )
    if (( \${#xmls[@]} == 0 )); then
        xmls=( *xml )
    fi
    if (( \${#xmls[@]} == 0 )); then
        echo "[ERROR] No XML files found (tried *overview.xml then *xml)" >&2
        exit 1
    fi

    if [[ ${params.parse_xml.pick_candies} == true ]]; then
        echo "Picking candies"

        # Building optional flags as an array so word-splitting is correct
        birdie_flag=()
        if [[ ${params.parse_xml.candy_picker_remove_birdies} == true ]]; then
            birdie_flag+=( --birdies "${params.parse_xml.candy_picker_birdies_file}" )
            birdie_flag+=( --birdie-harmonics "${params.parse_xml.birdies_harmonics}" )
            if [[ ${params.parse_xml.scale_birdie_width} == true ]]; then
                birdie_flag+=( --scale-birdie-width )
            fi
        fi

        dm_flag=()
        if [[ ${params.parse_xml.candy_picker_dm_tolerance} -gt 0 ]]; then
            dm_flag+=( -d "${params.parse_xml.candy_picker_dm_tolerance}" )
        fi

        # Run the picker (period threshold is required; other flags optional)
        candy_picker_rs -p "${params.parse_xml.candy_picker_period_threshold}" "\${birdie_flag[@]}" "\${dm_flag[@]}" "\${xmls[@]}"

        # Collect outputs; if none produced, fail loudly so downstream doesn’t break silently
        picked_xml_files=( *overview_picked.xml )
        if (( \${#picked_xml_files[@]} == 0 )); then
            echo "[ERROR] No *_overview_picked.xml produced by candy_picker_rs" >&2
            exit 1
        fi

        # rename pivots if present
        if [[ -f pivots.csv ]]; then
            mv pivots.csv "pivots_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}.csv"
        fi

        PICKED_XML_DIR="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/XML"
        mkdir -p "\${PICKED_XML_DIR}"
        cp "\${picked_xml_files[@]}" "\${PICKED_XML_DIR}/"
        cp -f *pivots*.csv "\${PICKED_XML_DIR}/" 2>/dev/null || true

    else
        echo "Not picking candies"
        picked_xml_files=( *overview.xml )
        PICKED_XML_DIR="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/XML"
        mkdir -p "\${PICKED_XML_DIR}"
        cp "\${picked_xml_files[@]}" "\${PICKED_XML_DIR}/"
    fi

    # Optional candidate filtering flag
    if [[ "${params.parse_xml.filter_cands}" == true ]]; then
        echo "Filtering candidates using config file: ${params.parse_xml.config_file}"
        config_flag=( --config_file "${params.parse_xml.config_file}" )
    else
        echo "Not filtering candidates"
        config_flag=()
    fi
    
    python3 ${params.parse_xml.script} -i \${picked_xml_files[@]} --chunk_id ${segments}${segment_id} --fold_technique ${params.psrfold.fold_technique} --nbins_default ${params.psrfold.nbins} --binplan "${params.psrfold.binplan}" ${subintlengthstring} -nsub ${params.psrfold.nsub} -clfd ${params.psrfold.clfd} -b ${beam_name} -b_id ${beam_id} -utc ${utc_start} -threads ${params.psrfold.threads}  --template_dir ${params.psrfold.template_dir} --telescope ${params.telescope} \${config_flag[@]} --cdm ${cdm} --cands_per_node ${params.psrfold.cands_per_node}

    mv filtered_candidates_file* filtered_candidates_file_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}.csv
    mv unfiltered_for_folding* unfiltered_for_folding_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}.csv
    """
}


process psrfold {
    label "psrfold"
    container "${params.pulsarx_image}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}"
    cache 'lenient'
    // maxForks 100
    // Organized output: separate directories for PNG, AR, and CANDS files
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/PNG/", pattern: "*.png", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/AR/", pattern: "*.ar", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/CANDS/", pattern: "*.cands", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file), val(start_sample), path(filtered_candidate_csv), path(candfile), path(metafile)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(fil_file, followLinks: false), path(filtered_candidate_csv), path(candfile), path(metafile), path("*.png"), path("*.ar"), path("*.cands")

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
        cand_no=\$(basename \${file} | cut -d'_' -f8 | cut -d'.' -f1)
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
    container "${params.rusty_candypicker}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}"
    cache 'lenient'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/CSV/", pattern: "*{.csv,master.cands}", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/PROVENANCE/", pattern: "*_provenance.csv", mode: 'copy'

    input:
    tuple val(pointing), val(cluster), val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), path(ars), path(cands)

    output:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), env(publish_dir), path(ars), path("*master.cands"), path("search_fold_cands*picked.csv")

    script:
    """
    # Base publish directory for PNG files (needed by downstream processes)
    publish_dir="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/PNG"
    mkdir -p \${publish_dir}
    mkdir -p "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/CSV"
    mkdir -p "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/PROVENANCE"

    fold_cands=\$(ls -v *.ar)
    pulsarx_cands_file=\$(ls -v *.cands)

    python3 ${baseDir}/scripts/fold_cands_to_csv.py -f \${fold_cands} -c \${pulsarx_cands_file} -x ${filtered_candidate_csv} -o search_fold_cands_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}.csv --cands_per_node ${params.psrfold.cands_per_node} -p \${publish_dir}

    echo "Number of candidates before candy picking:"
    cat search_fold_cands_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}.csv | wc -l

    if [[ ${params.psrfold.cluster_folded} == true ]]; then
        dm_flag=()
        if [[ ${params.psrfold.dm_tolerance} -gt 0 ]]; then
            dm_flag+=( --dmtol "${params.psrfold.dm_tolerance}" )
        fi
        echo "Picking candies"
        csv_candypicker --ptol ${params.parse_xml.candy_picker_period_threshold} "\${dm_flag[@]}" --tobs 0 -i search_fold_cands*.csv -o search_fold_cands_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}_picked.csv
    else
        echo "Not picking candies"
        cp search_fold_cands_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}.csv search_fold_cands_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}_picked.csv
    fi

    # ========================================================================
    # Generate provenance tracking file
    # This links each candidate to its source candfile, .cands file, and #id
    # ========================================================================
    provenance_file="${beam_name}_cdm_${cdm}_ck${segments}${segment_id}_provenance.csv"

    # Get paths for XML and candfile directories
    xml_dir="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML/XML"
    parsexml_dir="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/PARSEXML"
    ar_dir="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/FOLDING/AR"

    # Find XML files
    xml_files=""
    if [[ -d "\${xml_dir}" ]]; then
        xml_files=\$(ls "\${xml_dir}"/*_picked.xml 2>/dev/null | tr '\\n' ';' | sed 's/;\$//')
        if [[ -z "\${xml_files}" ]]; then
            xml_files=\$(ls "\${xml_dir}"/*.xml 2>/dev/null | tr '\\n' ';' | sed 's/;\$//')
        fi
    fi
    xml_files=\${xml_files:-"N/A"}

    # Write provenance header with comprehensive tracking columns
    echo "candidate_id,png_file,ar_file,cands_file,cands_file_id,candfile,candfile_id,xml_source,dm,period,f0,f1,snr,beam,segment_id,cdm,cluster,utc_start,ra,dec" > "\${provenance_file}"

    # Parse the master.cands file to extract provenance for each candidate
    master_cands_file=\$(ls *_master.cands 2>/dev/null | head -1)

    if [[ -f "\${master_cands_file}" ]]; then
        # Skip header line and process each row
        tail -n +2 "\${master_cands_file}" | while IFS=',' read -r id f0_new f1_new dm_new sn_new f0_old f1_old dm_old sn_old f0_err f1_err candidate_name filename rest; do
            # candidate_name format: CLUSTER_CANDFILENUM_cdm_CDM_ckSEGMENT_PEPOCH_BEAM_CANDNUM.png
            # filename format: CLUSTER_CANDFILENUM_cdm_CDM_ckSEGMENT_PEPOCH_BEAM.cands

            # Extract info from candidate_name
            if [[ -n "\${candidate_name}" ]]; then
                png_file="\${publish_dir}/\${candidate_name}"
                ar_basename=\${candidate_name%.png}.ar
                ar_file="\${ar_dir}/\${ar_basename}"

                # Extract candfile number from filename
                # The candfile number is the second field in the filename
                candfile_num=\$(echo "\${filename}" | cut -d'_' -f2)

                # Find the matching candfile
                candfile_path=""
                if [[ -d "\${parsexml_dir}" ]]; then
                    # Look for candfile with matching pattern
                    candfile_path=\$(ls "\${parsexml_dir}"/*_\${candfile_num}_*.candfile 2>/dev/null | head -1)
                    if [[ -z "\${candfile_path}" ]]; then
                        candfile_path=\$(ls "\${parsexml_dir}"/*.candfile 2>/dev/null | head -1)
                    fi
                fi
                candfile_path=\${candfile_path:-"N/A"}

                # The #id in the original candfile is: final_id - (candfile_num - 1) * cands_per_node
                # But we need the original candidate ID within the candfile
                # This is computed from the candidate_name's last field (the 5-digit number)
                final_cand_num=\$(echo "\${candidate_name}" | rev | cut -d'_' -f1 | cut -d'.' -f2 | rev)
                # Convert to integer (remove leading zeros)
                final_cand_num=\$((10#\${final_cand_num}))
                # Calculate original ID in candfile
                candfile_id=\$(( (final_cand_num - 1) % ${params.psrfold.cands_per_node} + 1 ))

                # Calculate period from f0
                period="N/A"
                if [[ -n "\${f0_new}" && "\${f0_new}" != "0" ]]; then
                    period=\$(echo "scale=12; 1.0 / \${f0_new}" | bc -l 2>/dev/null || echo "N/A")
                fi

                # cands_file is the filename column
                cands_file="\${filename}"
                # The ID in the cands file is the same as the #id column (row number)
                cands_file_id="\${id}"

                echo "\${id},\${png_file},\${ar_file},\${cands_file},\${cands_file_id},\${candfile_path},\${candfile_id},\${xml_files},\${dm_new},\${period},\${f0_new},\${f1_new},\${sn_new},${beam_name},${segments}${segment_id},${cdm},${cluster},${utc_start},${ra},${dec}" >> "\${provenance_file}"
            fi
        done
    else
        echo "WARNING: No master.cands file found, provenance file will be incomplete"
    fi

    echo "Provenance file created: \${provenance_file}"
    echo "Total candidates tracked: \$(wc -l < \${provenance_file})"
    """
}

process alpha_beta_gamma_test {
    label "alpha_beta_gamma_test"
    container "${params.pulsarx_image}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}"
    cache 'lenient'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/ZERODM/", pattern: "DM0*.png", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/ABG", pattern: "*alpha_beta_gamma.csv", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), val(png_source_dir), path(ars), path(master_cands), path(search_fold_cands_csv)

    output:
    tuple path("*alpha_beta_gamma.csv"), val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(segments), val(segment_id), val(png_source_dir)

    script:
    """
    #!/bin/bash
    publish_dir="${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/ABG"
    mkdir -p \${publish_dir}
    python3 ${baseDir}/scripts/calculate_alpha_beta_gamma_dmffdot.py -i ${search_fold_cands_csv} -o ${cluster}_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}_alpha_beta_gamma.csv -t ${params.alpha_beta_gamma.snr_min} -p \${publish_dir} -s ${png_source_dir} -c
    """
}

process pics_classifier {
    label "pics_classifier"
    container "${params.pics_classifier_image}"
    tag "${cluster}_${beam_name}_seg${segments}${segment_id}_cdm_${cdm}"
    cache 'lenient'
    publishDir "${params.basedir}/${params.runID}/${beam_name}/segment_${segments}/${segments}${segment_id}/CLASSIFICATION/", pattern: "*scored.csv", mode: 'copy'

    input:
    tuple val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(fft_size), val(segments), val(segment_id), val(fil_base_name), path(filtered_candidate_csv), val(png_source_dir), path(ars), path(master_cands), path(search_fold_cands_csv)

    output:
    tuple path("*scored.csv"), val(pointing), val(cluster),val(beam_name), val(beam_id), val(utc_start), val(ra), val(dec), val(cdm), val(segments), val(segment_id), val(png_source_dir)

    script:
    output_csv = "${cluster}_${beam_name}_cdm_${cdm}_ck${segments}${segment_id}_scored.csv"
    """
    python2 ${baseDir}/scripts/pics_classifier_multiple_models.py -m ${params.pics_model_dir} -o ${output_csv}
    """
}


process create_candyjar_tarball {
    executor 'local'
    container "${params.pulsarx_image}"
    tag "${output_tarball_name}"
    // Publish CSVs to dedicated directory for easy access
    publishDir "${params.basedir}/${params.runID}/TARBALL_CSV/", pattern: "*_header.csv", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/TARBALL_CSV/", pattern: "candidates.csv", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/TARBALL_CSV/", pattern: "candidates_alpha_below_one.csv", mode: 'copy'
    publishDir "${params.basedir}/${params.runID}/TARBALL_CSV/", pattern: "candidates_pics_above_threshold.csv", mode: 'copy'

    input:
    tuple path(candidate_results_file), val(output_tarball_name)

    output:
    tuple path("*_header.csv"), path("candidates.csv"), path("candidates_alpha_below_one.csv"), path("candidates_pics_above_threshold.csv")

    script:
    """
    #!/bin/bash
    header="pointing,target,beam,beam_id,utc_start,ra,dec,cdm,segments,segment_id,fold_cands_filepath,alpha_beta_file,pics_file"

    # Extract basename without extension and append "_header.csv"
    candidate_results_file_with_header="\$(basename ${candidate_results_file} .csv)_header.csv"
    echo "\$header" > "\$candidate_results_file_with_header"
    cat "${candidate_results_file}" >> "\$candidate_results_file_with_header"

    publish_dir="${params.basedir}/${params.runID}/CANDIDATE_TARBALLS"
    csv_dir="${params.basedir}/${params.runID}/TARBALL_CSV"
    mkdir -p \${publish_dir}
    mkdir -p \${csv_dir}

    # Run the tarball creation script
    # This creates candidates.csv, candidates_alpha_below_one.csv, and candidates_pics_above_threshold.csv
    # in the current working directory
    python ${baseDir}/scripts/create_candyjar_tarball.py -i \$candidate_results_file_with_header -o ${output_tarball_name} --verbose --npointings 0 -m ${params.metafile_source_path} -d \${publish_dir} --threshold ${params.alpha_beta_gamma.threshold} --snr_threshold ${params.alpha_beta_gamma.snr_threshold}
    """
}


// other pipelines
process parfold {
    label "parfold"
    container "${params.pulsarx_image}"
    tag "${beam_name}_${parfile_channel.getName()}"
    cache 'lenient'
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
    tag "extract_${csv_file.baseName}"
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
    tag "${cluster}_${beam_name}_${candfile.baseName}"
    cache 'lenient'
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
