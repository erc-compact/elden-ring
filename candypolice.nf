nextflow.enable.dsl=2

process parse_csv {
    label "candcsv_parse"
    container "${params.presto_image}"

    input:
    val sorted_candy_csv_file

    output:
    path ("candidates.txt")

    script:
    """
    #!/bin/bash
    if [ -f ${sorted_candy_csv_file} ]; then
        awk -F, '\$NF == "T1_CAND" || \$NF == "T2_CAND" {print \$3}' ${sorted_candy_csv_file} > candidates.txt
    else
        echo "${sorted_candy_csv_file}" > candidates.txt
    fi
    """
}

process candypolice {
    label "candypolice"
    container "${params.presto_image}"
    scratch true
    // publishDir "${params.candypolice.output_path}", pattern: "**/*.{ps,png,pfd}", mode: 'copy', saveAs: { path -> "${path.parent}/${path.name}" }
    publishDir "${params.candypolice.output_path}", pattern: "**/*.{ps,png,pfd}", mode: 'copy'
    input:
    val candidates
    path output_path

    output:
    tuple path("**/*.ps"), path("**/*.png"), path("**/*.pfd"), val(candidates)

    script:
    """
    #!/bin/bash
    python3 ${baseDir}/CandyPolice.py -csv ${params.candypolice.candidates_csv} -cand ${candidates} -fil ${params.candypolice.fold_fil_path} -num_cores ${params.candypolice.num_cores}
    """
}

workflow {
    candidates = parse_csv(params.candypolice.sorted_candy_csv_file)
    candidates_channel = candidates.splitText().map { it.trim() }
    candypolice(candidates_channel, params.candypolice.output_path)
}



