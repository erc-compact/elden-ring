process readfile {
    label 'readfile'
    container "${params.presto_image}"

    input:
    val(fits_file_channel_and_meta)

    output:
    tuple val(fits_file_channel_and_meta), env(time_per_file)

    script:
    def inputFile = "${fits_file_channel_and_meta[0].trim()}"
    """
    #!/bin/bash
    output=\$(readfile ${inputFile})
    echo "\$output"
    time_per_file=\$(echo "\$output" | grep "Time per file (sec)" | awk '{print \$6}')
    echo "\${time_per_file}" > time_per_file.txt
    """
}

process generateRfiFilter {
    label 'generate_rfi_filter'
    container "${params.rfi_mitigation_image}"
    publishDir "${params.generateRfiFilter.fits_rfi_output_dir}/", pattern: "*.{png,txt}", mode: 'symlink'

    input:
    tuple val(fits_file_channel_and_meta), val(time_per_file)

    output:
    tuple val(fits_file_channel_and_meta), env(rfi_filter_string)

    script:
    def inputFile = "${fits_file_channel_and_meta[0].trim()}"
    def beam_name = "${fits_file_channel_and_meta[2].trim()}"
    def num_intervals = Math.floor(time_per_file.toFloat() / 200) as int
    // def num_intervals = 2
    """
    #!/bin/bash
    export MPLCONFIGDIR=/tmp
    export NUMBA_CACHE_DIR=/tmp
    python3 ${params.generateRfiFilter.script} ${inputFile} . --target_resolution_ms ${params.generateRfiFilter.target_resolution_ms} --num_intervals ${num_intervals}

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

    publishDir "${params.filtool.output_path}/", pattern: "*.fil", mode: 'symlink'

    input:
    tuple val(fits_file_channel_and_meta), val(rfi_filter_string)
    val threads
    val telescope

    output:
    tuple val(fits_file_channel_and_meta), path("*.fil")
    
    script:
    def (fits_file, cluster, beam_name, beam_id, utc_start) = fits_file_channel_and_meta
    def outputFile = "${cluster.trim()}_${utc_start.trim()}_${beam_name.trim()}_clean"
    def inputFile = "${fits_file.trim()}"
    def source_name = "${cluster.trim()}"

    // Prepare the rfi_filter option
    def zaplist = ''
    if (rfi_filter_string) {
        zaplist = "-z ${rfi_filter_string}"
    }
    """
    #!/bin/bash
    # Get the first file from the inputFile string
    first_file=\$(echo ${inputFile} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    if [[ ${telescope} == "effelsberg" ]]; then
        if [[ "\${file_extension}" == "fits" ]]; then
            filtool --psrfits --scloffs --flip -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f ${inputFile} -s ${source_name}
        else 
            filtool -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f ${inputFile} -s ${source_name}
        fi

    elif [[ ${telescope} == "meerkat" ]]; then
        if [[ "\${file_extension}" == "sf" ]]; then
            filtool --psrfits --scloffs -t ${threads} --telescope ${telescope} ${zaplist} -o ${outputFile} -f ${inputFile} -s ${source_name}
        else 
            filtool -t ${threads} --telescope ${telescope} ${zaplist} -o "${outputFile}" -f ${inputFile} -s ${source_name}
        fi
    fi
    """
}

process psrfold {
    label "parfold"
    container "${params.pulsarx_image}"
    maxForks 100
    publishDir "${params.parfold.output_path}/", pattern: "*.png", mode: 'copy'
    publishDir "${params.parfold.output_path}/", pattern: "*.ar", mode: 'copy'
    publishDir "${params.parfold.output_path}/", pattern: "*.cands", mode: 'copy'
    
    input:
    tuple path(fil_file), val(cluster), val(beam_name), val(beam_id), val(utc_start)
    each path(parfile_channel)

    output:
    tuple path("*.png"), path("*.ar"), path("*.cands")
    
    script:
    def Outname = "${beam_name}_${parfile_channel.getName().replace(".par", "")}"

    """
    #!/bin/bash
    psrfold_fil2 --plotx -v -t ${params.splitcands.threads} --parfile ${parfile_channel} -n 64 -b 64 --nbinplan 0.1 128 --template ${params.template_dir}/Effelsberg_${beam_id}.template --clfd 2.0 -L 112 -f ${fil_file} -o ${Outname}
    """
}

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

    parfile_channel = Channel.fromPath("${params.parfold.parfile_path}/*.par").transpose()

    // Create a new channel with the metadata and new file path
    new_fil_file_channel = processed_fil_file.map { metadata, filepath -> 
        def (raw_file, cluster, beam_name, beam_id, utc_start) = metadata
        return tuple(filepath, cluster, beam_name, beam_id, utc_start)
    }   

    folded_cands = psrfold(new_fil_file_channel, parfile_channel)
}