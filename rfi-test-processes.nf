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
    // def num_intervals = Math.floor(time_per_file.toFloat() / 400) as int
    def num_intervals = 2
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

process generateRfiFilterSecond {
    label 'generate_rfi_filter'
    container "${params.rfi_mitigation_image}"

    publishDir "${params.generateRfiFilter.fil_rfi_output_dir}/", pattern: "*.png", mode: 'copy'
    publishDir "${params.generateRfiFilter.fil_rfi_output_dir}/", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(fits_file_channel_and_meta), path(cleaned_fil_file)

    output:
    tuple val(fits_file_channel_and_meta), env(rfi_filter_string)

    script:
    def inputFile = "${cleaned_fil_file}"
    def beam_name = "${fits_file_channel_and_meta[2].trim()}"
    """
    #!/bin/bash
    export MPLCONFIGDIR=/tmp
    export NUMBA_CACHE_DIR=/tmp
    python3 ${params.generateRfiFilter.script} ${inputFile} . --target_resolution_ms ${params.generateRfiFilter.target_resolution_ms}

    zap_commands=\$(grep -Eo '[0-9.]+ *- *[0-9.]+' combined_frequent_outliers.txt | \\
    awk -F '-' '{gsub(/ /,""); print "zap "\$1" "\$2}' | tr '\\n' ' ')

    rfi_filter_string="kadaneF 8 4 zdot \${zap_commands}"
    echo "\${rfi_filter_string}" > rfi_filter_string.txt

    mv combined_sk_heatmap_and_histogram.png ${beam_name}_second.png
    mv combined_frequent_outliers.txt combined_frequent_outliers_${beam_name}_second.txt
    mv block_bad_channel_percentages.txt block_bad_channel_percentages_${beam_name}_second.txt
    """
}