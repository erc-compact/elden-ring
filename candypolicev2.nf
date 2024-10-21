nextflow.enable.dsl=2

process extract_candidates {
    input:
    path csv_file // Path to the candidates CSV file
    path candidate_input // Path to the candidate input file

    output:
    path "candidate_details.txt" // Output the details to a text file

    script:
    """
    #!/usr/bin/env python3
    import csv
    import os

    csv_path = '${csv_file}'
    candidate_input_path = '${candidate_input}'

    # Read candidate input file and filter based on classification
    candidates = []
    with open(candidate_input_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=['beamid', 'utc', 'png', 'classification'])
        next(reader)  # Skip header
        for row in reader:
            if row['classification'] in ['T1_CAND', 'T2_CAND']:
                candidates.append(row['png'])

    # Read the CSV file and find matching candidates
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        with open('candidate_details.txt', 'w') as out_file:
            for row in reader:
                for candidate in candidates:
                    if candidate in row['png_path']:
                        details = {
                            'pointing_id': row['pointing_id'],
                            'beam_id': row['beam_id'],
                            'beam_name': row['beam_name'],
                            'source_name': row['source_name'],
                            'ra': row['ra'],
                            'dec': row['dec'],
                            'gl': row['gl'],
                            'gb': row['gb'],
                            'mjd_start': row['mjd_start'],
                            'utc_start': row['utc_start'],
                            'f0_user': row['f0_user'],
                            'f0_opt': row['f0_opt'],
                            'f0_opt_err': row['f0_opt_err'],
                            'f1_user': row['f1_user'],
                            'f1_opt': row['f1_opt'],
                            'f1_opt_err': row['f1_opt_err'],
                            'acc_user': row['acc_user'],
                            'acc_opt': row['acc_opt'],
                            'acc_opt_err': row['acc_opt_err'],
                            'dm_user': row['dm_user'],
                            'dm_opt': row['dm_opt'],
                            'dm_opt_err': row['dm_opt_err'],
                            'sn_fft': row['sn_fft'],
                            'sn_fold': row['sn_fold'],
                            'pepoch': row['pepoch'],
                            'maxdm_ymw16': row['maxdm_ymw16'],
                            'dist_ymw16': row['dist_ymw16'],
                            'pics_trapum_ter5': row['pics_trapum_ter5'],
                            'pics_palfa': row['pics_palfa'],
                            'pics_meerkat_l_sband_combined_best_recall': row['pics_meerkat_l_sband_combined_best_recall'],
                            'pics_palfa_meerkat_l_sband_best_fscore': row['pics_palfa_meerkat_l_sband_best_fscore'],
                            'png_path': row['png_path'],
                            'metafile_path': row['metafile_path'],
                            'filterbank_path': row['filterbank_path'],
                            'candidate_tarball_path': row['candidate_tarball_path']
                        }
                        # Write the candidate and details to the output file in one line
                        out_file.write(f"{candidate}|" + "|".join(f"{key}:{value}" for key, value in details.items()) + "\\n")
                        break  # Exit after finding the first match
    """
}

process candypolice {
    label "candypolice"
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

workflow {
    // Path to the candidates CSV file
    csv_file = params.candypolice.candidates_csv

    // Specify either a candidate string or path to cand_csv
    candidate_input = params.candypolice.sorted_candy_csv_file

    // Use the extract_candidates process to get candidate details
    details_file = extract_candidates(csv_file, candidate_input)

    // Read the candidate details file and split into lines
    candidate_details_lines = details_file.splitText().map { it.trim() }.collect()

    // // Input files for the fil files
    input_files = Channel.fromPath(params.candypolice.fold_fil_path).view() // Channel of input files
    
    candypolice(candidate_details_lines, input_files)
}