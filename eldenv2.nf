nextflow.enable.dsl=2

process find_files {
    label "find_files"

    output:
    path ("fits_files.txt")

    script:
    """
    #!/bin/bash
    baseDir=\$(dirname "${params.filtool.fits_files}")
    fileName=\$(basename "${params.filtool.fits_files}")
    find \${baseDir} -name "\${fileName}" 2>/dev/null | while read line; do 
    if [[ -f "\$line" ]] ; then
        echo "\$line" >> fits_files.txt
    fi
    done
    """
}

process filtool {
    label "filtool"
    container "${params.pulsarx_image}"
    publishDir "${params.filtool.output_path}", pattern: "*.fil", mode: 'copy' 
    
    input:
    path file_path

    output:
    path("*.fil")

    when:
    def file_path_str = file_path.toString()
    // check the file name and see if _Band_ is present, don't proceed for that file
    file_path_str.contains("_Band_")
    script:
    def file_path_str = file_path.toString()
    def file_name = new File(file_path_str).name.replaceAll(/\.fits$|\.sf$|\.rf$/, "")
    def GC = new File(file_path_str).name.replaceAll(/_Band_.*$/, "")
    def Band = new File(file_path_str).name.replaceAll(/.*_Band_(\d+).*/, '$1')
    def mask_param_name = "Band${Band}_Mask"
    // println "mask_param_name: $mask_param_name
    // Get the value of the mask using the dynamic mask_param_name
    def mask_value = params.filtool[mask_param_name]
    if (mask_value == null) {
        throw new IllegalStateException("Mask parameter ${mask_param_name} is not defined in params.filtool")
    }
    // println "mask_value: $mask_value"

    """
    #!/bin/bash
    # Run filtool with extracted values
    filtool -v --flip -o ${file_name} --psrfits ${file_path} -z ${mask_value}
    """
}

process nearest_power_of_two_calculator {
    label 'nearest_power_two'
    container "${params.presto_image}"

    input:
    path(fil_file)

    output:
    tuple path(fil_file), env(nearest_power_of_2)

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

process dmfile_gen {
    label "dmfile_gen"
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
    label "peasoup"
    container "${params.peasoup_image}"
    publishDir "${params.peasoup.output_path}", pattern: "**/*.xml", mode: 'copy'
    

    input:
    tuple val(fil_file), val(fft_size)
    each path(dm_file)

    output:
    tuple path("**/*.xml") , env(dm_name)

    script:
    def file_path_str = fil_file.toString()
    def GC = new File(file_path_str).name.replaceAll(/_Band_.*$/, "")
    def bandName = new File(file_path_str).name.replaceAll(/.*_Band_(\d+).*/, '$1')
    """
    #!/bin/bash
    dm_name=\$(basename "${dm_file}" .dm)

    echo "Running peasoup with DM file: \${dm_name}" and FFT size: ${fft_size} on ${fil_file}

    peasoup -p -i ${fil_file} -o ${GC}_${bandName}_{\${dm_name} --limit 100000 -m ${params.peasoup.min_snr} -t 1 --acc_start ${params.peasoup.acc_start} --acc_end ${params.peasoup.acc_end} --dm_file ${dm_file} --ram_limit_gb ${params.peasoup.ram_limit_gb} -n ${params.peasoup.nharmonics} --fft_size ${fft_size}
    """
    
}

process parse_xml{
    label "xml_parse"
    container "${params.pulsarx_image}"
    publishDir "${params.parse_xml.output_path}", pattern: "*.candfile", mode: 'copy'
    publishDir "${params.parse_xml.output_path}", pattern: "*_meta.txt", mode: 'copy'

    input:
    tuple path (xml_file), val(dm_name)


    output:
    tuple path("*.candfile"), path("*_meta.txt")

    script:
    """
    #!/bin/bash
    #need to add sub and chan mask later. 
    # Running the XML Parsing task
    python3 ${baseDir}/candidate_parser.py -i ${xml_file} -dm_name ${dm_name} -f ${params.parse_xml.fold_technique} -cpn ${params.parse_xml.pulsarx_cpn} -n ${params.parse_xml.nh} -nsub ${params.parse_xml.nsub} -clfd ${params.parse_xml.clfd} -fn ${params.parse_xml.fast_nbins} -sn ${params.parse_xml.slow_nbins}
    """
}

process pulsarx_fold{
    label "pulsarx_fold"
    container "${params.pulsarx_image}"
    publishDir "${params.pulsarx_fold.output_path}", pattern: "*.png", mode: 'copy'
    publishDir "${params.pulsarx_fold.output_path}", pattern: "*.ar", mode: 'copy'
    publishDir "${params.pulsarx_fold.output_path}", pattern: "*.cands", mode: 'copy'
    
    input:
    val xml_parser_results

    output:
    tuple path("*.png"), path("*.ar"), path("*.cands")

    script:
    """
    #!/bin/bash
    python3 ${baseDir}/pulsarx_fold.py -threads ${params.pulsarx_fold.ncpus} -p ${params.template_dir} -meta ${xml_parser_results[1]} -cands ${xml_parser_results[0]} -isitfits ${params.pulsarx_fold.isitfits} -t ${params.parse_xml.fold_technique} 
    """
}

process pics_classifier {
    label "pics_classifier"
    container "${params.psrchive_image}"
    publishDir "${params.pulsarx_fold.output_path}", pattern: "*.csv", mode: 'copy'

    input:
    path (fold_out)

    output:
    path("*.csv")
    script:
    """
    #!/bin/bash
    # Running the PICS Classifier
    python3 ${baseDir}/pics_classifier_v2.0.py -i ${params.pulsarx_fold.output_path} -m ${params.pics_model_dir}
    """
}


workflow {
    fits_files_ch = find_files()
    fits_queue = fits_files_ch..map splitText(){ it.trim() }
    fits_queue.view()
    filtool_channel = filtool(fits_queue)
    // filtool_channel = Channel.fromPath("/hercules/scratch/fkareem/NGC7099/Filtool/NGC7099_Band_*fil") //This is a quick fix to avoid the filtool process for fil files
    nearest_power_two_results = nearest_power_of_two_calculator(filtool_channel).transpose()
    filtool_output = nearest_power_two_results.map { item -> 
        def (fil_file, fft_size) = item
        return [fil_file, fft_size]
        }
    dmfiles = dmfile_gen()
    peasoup_results = peasoup(filtool_output, dmfiles)
    xml_parser_results= parse_xml(peasoup_results).transpose()
    xml_parser_results.view()
    fold_out = pulsarx_fold(xml_parser_results).collectFile(name: 'finished.txt', newLine : true)
    pics_classifier(fold_out)
}

