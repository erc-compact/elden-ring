nextflow.enable.dsl=2

process parse_xml{
    label "xml_parse"
    container "${params.pulsarx_image}"
    publishDir "${params.output_path}/${params.beam_name}/", pattern: "**/*.candfile", mode: 'copy'
    publishDir "${params.output_path}/${params.beam_name}/", pattern: "**/*_meta.txt", mode: 'copy'

    input:
    path xml_file
    val output_path


    output:
    tuple path("candidates/*.candfile"), path("candidates/*_meta.txt")

    script:
    """
    #!/bin/bash
    # Ensure the output directory exists
    mkdir -p ${output_path}/${params.beam_name}/candidates

    # Running the XML Parsing task
    python3 ${baseDir}/candidate_parser.py -i ${xml_file} -o candidates/  -f ${params.fold_technique} -cpn ${params.pulsarx_cpn} -b ${params.beam_name} -rfi ${params.rfi_filter} -n ${params.nh} -nsub ${params.nsub} -clfd ${params.clfd} -fn ${params.fast_nbins} -sn ${params.slow_nbins}
    """
}


process pulsarx_fold{
    label "pulsarx_fold"
    container "${params.pulsarx_image}"
    publishDir "${params.output_path}/Band_${params.band}", pattern: "**/*.png", mode: 'copy'
    publishDir "${params.output_path}/Band_${params.band}", pattern: "**/*.ar", mode: 'copy'
    
    input:
    path xml_parser_results
    val output_path

    output:
    tuple path("**/*.png"), path("**/*.ar")

    script:
    """
    #!/bin/bash
    python3 ${baseDir}/pulsarx_fold.py -threads ${params.threads} -p ${params.template_dir} -meta ${xml_parser_results[1]} -cands ${xml_parser_results[0]} -fits ${params.isitfits} -band ${params.band}
    """
}

workflow {
    xml_queue = Channel.fromPath(params.xml_files)
    xml_parser_results= parse_xml(xml_queue, params.output_path).transpose()
    
    pulsarx_fold(xml_parser_results, params.output_path)
}
