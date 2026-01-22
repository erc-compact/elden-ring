/*
 * ============================================================================
 * PRESTO Pipeline - Processes and Workflows
 * ============================================================================
 *
 * This file contains all PRESTO-related processes and workflows for pulsar
 * searching. The pipeline consists of 6 stages that can be run together or
 * independently using the unified presto_pipeline entry point.
 *
 * Pipeline Stages:
 *   1. RFI Detection      (rfifind)
 *   2. Birdie Detection   (zero-DM accelsearch)
 *   3. Dedispersion       (prepsubband)
 *   4. Acceleration Search (accelsearch)
 *   5. Sift and Fold      (candidate sifting + prepfold/pulsarx)
 *   6. Post-processing    (PNG conversion + tarball creation)
 *
 * Usage Examples:
 *
 *   # Run full pipeline (all 6 stages):
 *   nextflow run elden.nf -entry presto_pipeline --input_fil /path/to/file.fil
 *
 *   # Run only stages 1-3 (RFI through Dedispersion):
 *   nextflow run elden.nf -entry presto_pipeline \
 *       --input_fil file.fil --presto.end_stage 3
 *
 *   # Resume from stage 4 using state file:
 *   nextflow run elden.nf -entry presto_pipeline \
 *       --presto.start_stage 4 --state_file dedisperse_state.json
 *
 *   # Run with CSV input (multiple files):
 *   nextflow run elden.nf -entry presto_full_entry --files_list input.csv
 *
 * State files are saved to: ${params.output_dir}/${params.target_name}/presto/state/
 * ============================================================================
 */

// ============================================================================
// PRESTO PROCESSES
// ============================================================================

/*
 * PRESTO RFI Detection - rfifind
 * Creates RFI mask for the observation
 */
process presto_rfifind {
    label 'presto'
    label 'short'

    publishDir "${params.output_dir}/${params.target_name}/presto/rfi", mode: 'copy', pattern: "*.mask"

    input:
    path input_file

    output:
    path "*.mask", emit: rfi_mask
    path "*_rfifind.stats", emit: rfi_stats
    path "*_rfifind.inf", emit: rfi_inf
    path "*_rfifind.out", emit: rfi_outfile
    path "*rfifind*", emit: rfi_files

    script:
    def basename = input_file.baseName
    def time_int = params.presto?.rfifind_time ?: 2.0
    def chanfrac = params.presto?.rfifind_chanfrac ?: 0.5
    def intfrac = params.presto?.rfifind_intfrac ?: 0.3
    def timesig = params.presto?.rfifind_timesig ?: 10
    def freqsig = params.presto?.rfifind_freqsig ?: 4
    """
    rfifind -time ${time_int} -chanfrac ${chanfrac} -intfrac ${intfrac} \
        -timesig ${timesig} -freqsig ${freqsig} \
        -o ${basename} ${input_file} 2>&1 | tee ${basename}_rfifind.out
    """
}

/*
 * PRESTO prepdata for zero-DM time series (birdie detection)
 */
process presto_prepdata_zerodm {
    label 'presto'
    label 'short'

    input:
    path input_file
    path rfi_mask
    path rfi_stats
    path rfi_inf

    output:
    path "*.dat", emit: dat_file
    path "*.inf", emit: inf_file

    script:
    def basename = input_file.baseName
    """
    # Ensure rfifind sidecars are available for prepdata
    prepdata -dm 0 -mask ${rfi_mask} -o ${basename}_DM0 ${input_file}
    """
}

/*
 * PRESTO accelsearch at zero DM for birdie list
 */
process presto_accelsearch_zerodm {
    label 'presto'
    label 'medium'

    input:
    path dat_file
    path inf_file

    output:
    path "*.zaplist", emit: zaplist, optional: true
    path "*_ACCEL_*", emit: accel_files

    script:
    def basename = dat_file.baseName
    def zmax = params.presto?.birdie_zmax ?: 0
    def numharm = params.presto?.numharm ?: 8
    """
    realfft ${dat_file}
    accelsearch -zmax ${zmax} -numharm ${numharm} ${basename}.fft

    # Create zaplist from significant zero-DM candidates (potential birdies)
    if [ -f ${basename}_ACCEL_${zmax} ]; then
        awk 'NR>3 && \$1!="#" && \$4>5 {print \$2}' ${basename}_ACCEL_${zmax} > birdies.zaplist || true
    fi
    touch birdies.zaplist
    """
}

/*
 * PRESTO prepsubband - dedispersion across DM range
 */
process presto_prepsubband {
    label 'presto'
    label 'long'

    publishDir "${params.output_dir}/${params.target_name}/presto/subbands/${dm_low}_${dm_high}", mode: 'symlink'

    input:
    path input_file
    path input_inf
    path rfi_mask
    path rfi_stats
    tuple val(dm_low), val(dm_high), val(dm_step), val(downsamp)

    output:
    path "*.dat", emit: dat_files
    path "*.inf", emit: inf_files
    tuple val(dm_low), val(dm_high), val(dm_step), val(downsamp), path("*.dat"), path("*.inf"), emit: subband_data

    script:
    def basename = input_file.baseName
    def nsub = params.presto?.nsub ?: 128
    def ndm = Math.max(1, ((dm_high.toFloat() - dm_low.toFloat()) / dm_step.toFloat()).toInteger() + 1)
    """
    # Verify input files exist
    if [ ! -f "${input_file}" ]; then
        echo "ERROR: Input file ${input_file} not found"
        exit 1
    fi

    if [ ! -f "${rfi_mask}" ]; then
        echo "ERROR: RFI mask ${rfi_mask} not found"
        exit 1
    fi

    if [ ! -f "${rfi_stats}" ]; then
        echo "ERROR: RFI stats ${rfi_stats} not found"
        exit 1
    fi

    # Optional zerodm .inf file (prepsubband can operate without it for .fil)
    if [ ! -f "${input_inf}" ]; then
        echo "WARNING: zerodm .inf file ${input_inf} not found; continuing without it"
    fi

    # Run prepsubband with correct flags (-lodm/-dmstep/-numdms)
    prepsubband -lodm ${dm_low} -dmstep ${dm_step} -numdms ${ndm} \
        -nsub ${nsub} -downsamp ${downsamp} \
        -mask ${rfi_mask} -o ${basename} ${input_file}
    """
}

/*
 * PRESTO accelsearch - acceleration/jerk search
 * Can use GPU acceleration if available
 */
process presto_accelsearch {
    label 'presto_gpu'
    label 'long'

    publishDir "${params.output_dir}/${params.target_name}/presto/accelsearch/${(params.presto?.segment_chunks ?: 1) > 1 ? 'segmented' : 'full'}", mode: 'symlink'

    input:
    path dat_file
    path inf_file
    path zaplist

    output:
    path "*_ACCEL_*", emit: accel_files
    path "*_ACCEL_*.cand", emit: cand_files, optional: true
    path "*_ACCEL_*.inf", emit: accel_inf, optional: true

    script:
    def basename = dat_file.baseName
    def zmax = params.presto?.zmax ?: 200
    def wmax = params.presto?.wmax ?: 0
    def numharm = params.presto?.numharm ?: 8
    def sigma = params.presto?.sigma_threshold ?: 2.0
    def use_gpu = params.presto?.use_gpu ?: false
    def accel_bin = (use_gpu ? (params.presto?.accelsearch_cuda_bin ?: 'accelsearch_cu') : 'accelsearch')
    def zap_flag = zaplist.name != 'NO_ZAPLIST' ? "-zaplist ${zaplist}" : ""
    def chunks = params.presto?.segment_chunks ?: 1
    """
    set -euo pipefail

    run_accel () {
        local dat="\$1"
        local inf="\$2"
        local outbase="\$3"

        realfft "\${dat}"
        rednoise "\${outbase}.fft"
        cp "\${inf}" "\${outbase}_red.inf"
        rm -f "\${outbase}.fft"

        if [ -f "${zaplist}" ] && [ -s "${zaplist}" ]; then
            zapbirds -zap -zapfile ${zaplist} "\${outbase}_red.fft"
            mv "\${outbase}_red.fft" "\${outbase}.fft"
            mv "\${outbase}_red.inf" "\${outbase}.inf"
        else
            mv "\${outbase}_red.fft" "\${outbase}.fft"
            mv "\${outbase}_red.inf" "\${outbase}.inf"
        fi

        ${accel_bin} -zmax ${zmax} -wmax ${wmax} -numharm ${numharm} -sigma ${sigma} ${zap_flag} "\${outbase}.fft"

        # Keep matching inf with ACCEL output for downstream
        cp "\${outbase}.inf" "\${outbase}_ACCEL_${zmax}${wmax > 0 ? "_JERK_${wmax}" : ""}.inf"
        rm -f "\${outbase}.fft" "\${outbase}.inf"
    }

    # Handle possibly multiple dat/inf pairs (one per DM)
    dat_files=( ${dat_file} )
    inf_files=( ${inf_file} )

    if [ \${#dat_files[@]} -ne \${#inf_files[@]} ]; then
        echo "ERROR: dat/inf count mismatch: \${#dat_files[@]} vs \${#inf_files[@]}"
        exit 1
    fi

    for idx in "\${!dat_files[@]}"; do
        dat=\${dat_files[\$idx]}
        inf=\${inf_files[\$idx]}
        base=\$(basename "\${dat}" .dat)

        # Ensure the directory referenced by the .inf "Data file name" exists
        red_dat=\$(awk -F'[:=]' '/Data file name/ {print \$2; exit}' "\${inf}" | xargs)
        if [ -n "\${red_dat}" ]; then
            red_dir=\$(dirname "\${red_dat}")
            if [ "\${red_dir}" != "." ]; then
                mkdir -p "\${red_dir}"
            fi
        fi

        if [ ${chunks} -le 1 ]; then
            run_accel "\${dat}" "\${inf}" "\${base}"
        else
            nsamples=\$(grep "Number of bins in the time series" "\${inf}" | awk -F'=' '{print \$2}' | tr -d ' ')
            samples_per_chunk=\$(( nsamples / ${chunks} ))
            if [ \$((samples_per_chunk % 2)) -ne 0 ]; then
                samples_per_chunk=\$((samples_per_chunk - 1))
            fi
            for i in \$(seq 1 ${chunks}); do
                start_frac=\$(python - << 'PY'
import sys
chunks=int(sys.argv[1]); idx=int(sys.argv[2])
print(f"{(idx-1)/chunks:.6f}")
PY
 ${chunks} \$i)
                seg_base="\${base}_seg\$i"
                prepdata -nobary -dm 0 -start \$start_frac -numout \$samples_per_chunk -o "\${seg_base}" "\${dat}"
                run_accel "\${seg_base}.dat" "\${seg_base}.inf" "\${seg_base}"
            done
        fi
    done
    """
}

/*
 * PRESTO candidate sifting
 */
process presto_sift_candidates {
    label 'presto'
    label 'short'

    publishDir "${params.output_dir}/${params.target_name}/presto/sifted", mode: 'copy'

    input:
    path cand_files
    path input_file

    output:
    path "sifted_candidates.csv", emit: sifted_csv
    path "sifted_candidates.candfile", emit: candfile
    path "sifted_candidates.provenance.csv", emit: provenance_csv

    script:
    def short_period = params.presto?.sift_short_period ?: 0.0005
    def long_period = params.presto?.sift_long_period ?: 15.0
    def sigma_threshold = params.presto?.sift_sigma_threshold ?: 4.0
    def max_cands = params.presto?.sift_max_cands ?: 500
    def min_num_dms = params.presto?.sift_min_num_dms ?: 1
    def low_dm_cutoff = params.presto?.sift_low_dm_cutoff ?: 2.0
    def zmax = params.presto?.zmax ?: 200
    def wmax = params.presto?.wmax ?: 0
    """
    shopt -s nullglob
    # Collect only main ACCEL files (without extensions like .cand, .txtcand, .inf)
    cand_list=()
    for file in *_ACCEL_*; do
        # Only include files that don't have common extensions
        if [[ ! "\$file" =~ \\.(cand|txtcand|inf)\$ ]]; then
            if [ ${wmax} -gt 0 ]; then
                if [[ "\$file" == *_ACCEL_${zmax}_JERK_${wmax} ]]; then
                    cand_list+=("\$file")
                fi
            else
                if [[ "\$file" == *_ACCEL_${zmax} ]]; then
                    cand_list+=("\$file")
                fi
            fi
        fi
    done
    shopt -u nullglob

    if [ \${#cand_list[@]} -eq 0 ]; then
        echo "No ACCEL files found; creating empty outputs"
        echo "id,dm,f0,f1,f2,snr,sigma,sn_fft,acc,accel,z,w,period_ms,file,cand_num,accel_file" > sifted_candidates.csv
        echo "#id dm acc F0 F1 F2 S/N" > sifted_candidates.candfile
        echo "id,accel_file,cand_num,dm,f0,f1,f2,sigma,snr" > sifted_candidates.provenance.csv
        exit 0
    fi

    echo "Found \${#cand_list[@]} ACCEL files to sift"

    python3 ${projectDir}/scripts/presto_accel_sift.py \\
        --sigma-threshold ${sigma_threshold} \\
        --min-period ${short_period} \\
        --max-period ${long_period} \\
        --min-num-dms ${min_num_dms} \\
        --low-dm-cutoff ${low_dm_cutoff} \\
        --max-cands-to-fold ${max_cands} \\
        --remove-duplicates \\
        --remove-harmonics \\
        --output sifted_candidates \\
        \${cand_list[@]}
    """
}

/*
 * PRESTO prepfold - fold single candidate
 */
process presto_prepfold {
    label 'presto'
    label 'medium'

    publishDir "${params.output_dir}/${params.target_name}/presto/folded", mode: 'copy'

    input:
    path input_file
    path candidate_csv
    path rfi_mask
    path rfi_stats
    val candidate_line

    output:
    path "*.pfd", emit: pfd_files
    path "*.pfd.bestprof", emit: bestprof_files
    path "*.pfd.ps", emit: ps_files

    script:
    def basename = input_file.baseName
    def nsub = params.presto?.fold_nsub ?: 128
    def npart = params.presto?.fold_npart ?: 64
    def pstep = params.presto?.fold_pstep ?: 1
    def pdstep = params.presto?.fold_pdstep ?: 2
    def dmstep = params.presto?.fold_dmstep ?: 2
    """
    #!/usr/bin/env python3
    import subprocess
    import os

    line = "${candidate_line}".strip()
    parts = line.split(',')

    # Expected CSV format: filename,cand_num,dm,period,f0,f1,accel,sigma
    if len(parts) >= 5:
        dm = parts[2]
        period = parts[3]
        f0 = parts[4]
        f1 = parts[5] if len(parts) > 5 else "0"
        accel = parts[6] if len(parts) > 6 else "0"
        cand_id = parts[1]

        outname = f"${basename}_cand_{cand_id}"

        cmd = f"prepfold -dm {dm} -p {period} -pd {f1} " \
              f"-mask ${rfi_mask} -nsub ${nsub} -npart ${npart} " \
              f"-pstep ${pstep} -pdstep ${pdstep} -dmstep ${dmstep} " \
              f"-o {outname} ${input_file}"

        subprocess.run(cmd, shell=True, check=True)
    """
}

/*
 * PRESTO prepfold batch - fold multiple candidates
 *
 * This process applies period correction using the formulas from foldx.py:
 *   pdot = period * accel / LIGHT_SPEED
 *   fold_period = period - pdot * fft_size * tsamp / 2
 *
 * The corrected period is used with prepfold's -topo flag for topocentric folding.
 */
process presto_prepfold_batch {
    label 'presto'
    label 'long'

    publishDir "${params.presto?.output_dir ?: params.basedir}/presto/folded", mode: 'copy'

    input:
    path input_file
    path candidate_csv
    path rfi_mask
    path rfi_stats

    output:
    path "*.pfd", emit: pfd_files, optional: true
    path "*.pfd.bestprof", emit: bestprof_files, optional: true
    path "*.pfd.ps", emit: ps_files, optional: true

    script:
    def basename = input_file.baseName
    def npart = params.presto?.fold_npart ?: 64
    def max_cands = params.presto?.max_fold_cands ?: 100
    def mask_opt = rfi_mask.name != 'NO_MASK' ? "-mask ${rfi_mask}" : ""
    def start_frac = params.presto?.fold_start_frac ?: 0.0
    def end_frac = params.presto?.fold_end_frac ?: 1.0
    def extra_flags = params.presto?.fold_extra_flags ?: ""
    """
    #!/usr/bin/env python3
    import subprocess
    import csv
    from concurrent.futures import ProcessPoolExecutor, as_completed

    def fold_candidate(args):
        f0, f1, f2, dm, cand_id, basename, input_file, mask_opt = args
        outname = f"{basename}_cand_{cand_id}"
        cmd = (
            f"prepfold -noxwin "
            f"-f {f0:.15f} -fd {f1:.15f} -fdd {f2:.15f} "
            f"-dm {dm:.6f} -npart ${npart} {mask_opt} "
            f"-start ${start_frac} -end ${end_frac} "
            f"${extra_flags} "
            f"-o {outname} {input_file}"
        )

        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            return True, cand_id
        except Exception as e:
            return False, str(e)

    # Read candidates from CSV
    # Expected columns: f0,f1,f2,dm (from sifter output)
    candidates = []
    with open("${candidate_csv}", 'r') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= ${max_cands}:
                break

            try:
                f0 = float(row.get('f0', 0.0))
                f1 = float(row.get('f1', 0.0))
                f2 = float(row.get('f2', 0.0))
                dm = float(row.get('dm', row.get('DM', 0)))
                cand_id = row.get('cand_id', row.get('id', i))

                candidates.append((
                    f0,
                    f1,
                    f2,
                    dm,
                    cand_id,
                    "${basename}",
                    "${input_file}",
                    "${mask_opt}"
                ))
            except (KeyError, ValueError) as e:
                print(f"Warning: Could not parse candidate row: {e}")
                continue

    print(f"Folding {len(candidates)} candidates with prepfold")
    if not candidates:
        raise SystemExit(0)

    # Fold in parallel
    with ProcessPoolExecutor(max_workers=min(8, len(candidates))) as executor:
        futures = [executor.submit(fold_candidate, c) for c in candidates]
        for future in as_completed(futures):
            success, result = future.result()
            if not success:
                print(f"Warning: Folding failed: {result}")
    """
}

/*
 * PRESTO PFD to PNG conversion
 * Uses show_pfd to render plots from .pfd files, then converts to PNG.
 */
process presto_pfd_to_png {
    label 'presto'
    label 'short'

    publishDir "${params.output_dir}/${params.target_name}/presto/plots", mode: 'copy'

    input:
    path pfd_files

    output:
    path "*.png", emit: png_files, optional: true

    script:
    """
    set -euo pipefail

    for pfd in ${pfd_files}; do
        if [ -f "\$pfd" ]; then
            show_pfd -noxwin "\$pfd"
        else
            echo "WARNING: PFD file not found: \$pfd" >&2
        fi
    done

    shopt -s nullglob
    for ps in *.ps; do
        base="\${ps%.ps}"
        if command -v gs >/dev/null 2>&1; then
            gs -dSAFER -dBATCH -dNOPAUSE -dAutoRotatePages=/None \
                -sDEVICE=png16m -r150 \
                -sOutputFile="\${base}.png" \
                -c "<</Orientation 3>> setpagedevice" \
                -f "\$ps"
        elif command -v convert >/dev/null 2>&1; then
            convert -density 150 -rotate 90 "\$ps" "\${base}.png"
        else
            echo "WARNING: gs/convert not found; leaving \$ps" >&2
        fi
    done
    """
}

/*
 * PRESTO PS to PNG conversion (legacy)
 * Converts PostScript files (.pfd.ps) directly to PNG.
 */
process presto_ps_to_png {
    label 'presto'
    label 'short'

    publishDir "${params.output_dir}/${params.target_name}/presto/plots", mode: 'copy'

    input:
    path ps_files

    output:
    path "*.png", emit: png_files

    script:
    """
    python ${projectDir}/scripts/presto_ps_to_png.py \
        --input_dir . \
        --output_dir . \
        --workers 4
    """
}

/*
 * PRESTO merge folded results
 */
process presto_fold_merge {
    label 'presto'
    label 'short'

    publishDir "${params.output_dir}/${params.target_name}/presto/merged", mode: 'copy'

    input:
    path pfd_files
    path png_files
    path bestprof_files
    path sifted_csv
    path provenance_csv
    val meta_info

    output:
    path "merged_results.csv", emit: merged_csv
    path "pfd_files/*", emit: all_pfd, optional: true
    path "png_files/*", emit: all_png, optional: true

    script:
    def pointing = meta_info[0]
    def cluster = meta_info[1]
    def beam_name = meta_info[2]
    def beam_id = meta_info[3]
    def utc_start = meta_info[4]
    def ra = meta_info[5]
    def dec = meta_info[6]
    def cdm = meta_info[7]
    def filterbank_path = meta_info[8]
    """
    #!/usr/bin/env python3
    import os
    import csv
    import shutil
    import glob

    os.makedirs("pfd_files", exist_ok=True)
    os.makedirs("png_files", exist_ok=True)

    # Copy pfd and png files
    for pfd in glob.glob("*.pfd"):
        shutil.copy(pfd, "pfd_files/")
    for png in glob.glob("*.png"):
        shutil.copy(png, "png_files/")

    meta_defaults = {
        "pointing_id": "${pointing}",
        "beam_id": "${beam_id}",
        "beam_name": "${beam_name}",
        "utc_start": "${utc_start}",
        "ra": "${ra}",
        "dec": "${dec}",
        "cdm": "${cdm}",
        "filterbank_path": "${filterbank_path}",
        "cluster": "${cluster}",
        "segment_id": "0",
        "source_name": "",
    }

    if meta_defaults["utc_start"]:
        meta_defaults["metafile_path"] = "meta/%s.meta" % meta_defaults["utc_start"]
    else:
        meta_defaults["metafile_path"] = ""

    # Read sifted candidates
    sifted_data = {}
    with open("${sifted_csv}", 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cand_key = row.get('cand_id') or row.get('id') or row.get('cand_num', '')
            if cand_key:
                sifted_data[str(cand_key)] = row
            cand_num = row.get('cand_num', '')
            if cand_num:
                sifted_data[str(cand_num)] = row

    # Read provenance data
    prov_data = {}
    with open("${provenance_csv}", 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            prov_id = row.get('id', '')
            if prov_id:
                prov_data[str(prov_id)] = row

    def parse_bestprof(path):
        info = {}
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith("# Epoch_bary (MJD)"):
                        val = line.split("=", 1)[1].strip()
                        try:
                            info["mjd_start"] = float(val)
                        except ValueError:
                            pass
                    elif line.startswith("# Best DM"):
                        val = line.split("=", 1)[1].strip()
                        try:
                            info["dm_opt"] = float(val)
                        except ValueError:
                            pass
        except IOError:
            pass
        return info

    results = []

    for bp in glob.glob("*.bestprof"):
        raw_cand_id = bp.split('_cand_')[1].replace('.pfd.bestprof', '') if '_cand_' in bp else ''
        cand_id = raw_cand_id
        if raw_cand_id:
            cand_id = raw_cand_id.split('_')[0]
        basename = os.path.basename(bp).replace('.pfd.bestprof', '')

        result = {
            'basename': basename,
            'cand_id': cand_id,
            'cand_label': raw_cand_id,
            'pfd_file': basename + '.pfd',
            'png_file': basename + '.png' if os.path.exists(basename + '.png') else ''
        }

        result.update(meta_defaults)

        if cand_id in sifted_data:
            result.update(sifted_data[cand_id])

        if cand_id in prov_data:
            result["accel_file"] = prov_data[cand_id].get("accel_file", "")
            result["cand_num"] = prov_data[cand_id].get("cand_num", "")
            result["sn_fft"] = prov_data[cand_id].get("sn_fft", prov_data[cand_id].get("sigma", ""))
            if result.get("sn_fft") not in (None, ""):
                try:
                    result["sn_fft"] = float(result["sn_fft"])
                except Exception:
                    pass

        bp_info = parse_bestprof(bp)
        result.update(bp_info)
        if (not result.get("utc_start")) and ("mjd_start" in result):
            result["utc_start"] = str(result.get("mjd_start"))

        sigma_val = result.get("sigma") or result.get("snr") or result.get("sn_fft")
        if sigma_val not in (None, ""):
            try:
                result["sn_fold"] = float(sigma_val)
            except Exception:
                result["sn_fold"] = sigma_val

        if "dm" in result and "dm_user" not in result:
            result["dm_user"] = result.get("dm")
        if "dm_opt" not in result and "dm" in result:
            result["dm_opt"] = result.get("dm")

        if "f0" in result and "f0_user" not in result:
            result["f0_user"] = result.get("f0")
        if "f1" in result and "f1_user" not in result:
            result["f1_user"] = result.get("f1")

        results.append(result)

    if results:
        fieldnames = list(results[0].keys())
        with open("merged_results.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
    else:
        with open("merged_results.csv", 'w') as f:
            f.write("basename,cand_id,pfd_file,png_file\\n")
    """
}

/*
 * PRESTO PICS classifier
 * Runs PICS classification on PRESTO .pfd files
 */
process presto_pics_classifier {
    label 'pics_classifier'
    container "${params.pics_classifier_image}"

    publishDir "${params.output_dir}/${params.target_name}/presto/classification", mode: 'copy', pattern: "*.csv"

    input:
    path pfd_files
    path merged_csv

    output:
    path "pics_scored.csv", emit: pics_csv
    path "presto_classified.csv", emit: classified_csv

    script:
    """
    #!/bin/bash
    shopt -s nullglob
    pfd_list=( *.pfd )
    shopt -u nullglob
    if [ \${#pfd_list[@]} -eq 0 ]; then
        echo "No PFD files found; creating empty classification outputs"
        cp "${merged_csv}" presto_classified.csv
        echo "filename" > pics_scored.csv
        exit 0
    fi

    # Run PICS classifier on PFD files
    python2 ${projectDir}/scripts/pics_classifier_multiple_models.py \
        -i . \
        -m ${params.pics_model_dir} \
        -o pics_scored.csv

    # Merge PICS scores with search results (python2-only container)
    python2 << 'MERGE_SCRIPT'
import pandas as pd
import os

# Read merged search results
search_df = pd.read_csv("${merged_csv}")

# Read PICS scores
pics_df = pd.read_csv("pics_scored.csv")

# Extract basename from pfd filename for matching
pics_df['basename'] = pics_df['filename'].apply(lambda x: os.path.splitext(x)[0])

# Merge on basename or pfd_file
if 'pfd_file' in search_df.columns:
    search_df['match_key'] = search_df['pfd_file'].apply(lambda x: os.path.splitext(str(x))[0])
elif 'basename' in search_df.columns:
    search_df['match_key'] = search_df['basename']
else:
    search_df['match_key'] = search_df.index.astype(str)

# Perform merge
merged = search_df.merge(
    pics_df,
    left_on='match_key',
    right_on='basename',
    how='left',
    suffixes=('', '_pics')
)

# Drop temporary columns
if 'match_key' in merged.columns:
    merged.drop(columns=['match_key'], inplace=True)
if 'basename_pics' in merged.columns:
    merged.drop(columns=['basename_pics'], inplace=True)
if 'filename' in merged.columns:
    merged.drop(columns=['filename'], inplace=True)

# Save classified results
merged.to_csv("presto_classified.csv", index=False)
print("Merged %d candidates with PICS scores" % len(merged))
MERGE_SCRIPT
    """
}

/*
 * PRESTO create tarball for CandyJar
 * Creates CandyJar-compatible tarball with the same structure and CSV fields
 * as create_candyjar_tarball.py (PNG + CSV + metafiles, no PFD/archive files)
 *
 * Tarball Structure:
 *   {tarball_prefix}_presto/
 *   ├── candidates.csv
 *   ├── candidates_pics_above_threshold.csv
 *   ├── metafiles/
 *   │   └── {utc_start}.meta
 *   └── plots/
 *       └── *.png
 */
process presto_create_tarball {
    label 'pulsarx'
    label 'short'
    container "${params.pulsarx_container}"

    publishDir "${params.output_dir}/${params.target_name}/presto/tarballs", mode: 'copy'

    input:
    path classified_csv
    path png_files
    val tarball_prefix

    output:
    path "*.tar.gz", emit: tarball
    path "candidates.csv", emit: final_csv
    path "candidates_pics_above_threshold.csv", emit: pics_filtered_csv

    script:
    def tarball_name = "${tarball_prefix}_presto"
    def metafile_dir = params.metafile_source_path ?: ""
    def pics_threshold = params.presto?.pics_threshold ?: 0.1
    def snr_threshold = params.presto?.snr_threshold ?: 4.0
    """
    # Create tarball using the CandyJar-compatible script
    python3 ${projectDir}/scripts/create_presto_tarball.py \\
        -i ${classified_csv} \\
        -o ${tarball_name}.tar.gz \\
        --png-dir . \\
        ${metafile_dir ? "--metafile-dir ${metafile_dir}" : ""} \\
        --tarball-prefix ${tarball_prefix} \\
        --pics-threshold ${pics_threshold} \\
        --snr-threshold ${snr_threshold} \\
        --verbose

    # The script emits candidates.csv and candidates_pics_above_threshold.csv
    # in the work directory as well as in the tarball
    """
}

/*
 * Fold PRESTO/peasoup candidates using PulsarX (psrfold_fil2)
 */
process presto_fold_pulsarx {
    label 'pulsarx'
    label 'long'

    publishDir "${params.presto?.output_dir ?: params.basedir}/pulsarx_fold", mode: 'copy'

    input:
    path input_file
    path sifted_csv
    path meta_file

    output:
    path "*.png", emit: png_files, optional: true
    path "*.ar", emit: ar_files, optional: true

    script:
    def fold_backend = params.presto?.fold_backend ?: 'pulsarx'
    def threads = params.presto?.fold_threads ?: 16
    """
    python ${projectDir}/scripts/pulsarx_fold.py \\
        -meta ${meta_file} \\
        -cands ${sifted_csv} \\
        -o . \\
        --csv \\
        --fold_backend ${fold_backend}
    """
}

/*
 * Generate meta file for CSV-based folding (PRESTO/peasoup candidates)
 */
process generate_fold_meta {
    label 'short'

    input:
    tuple val(pointing), path(fil_file), val(cluster), val(beam_name), val(beam_id),
          val(utc_start), val(ra), val(dec), val(cdm)
    val(fft_size)
    val(fold_backend)

    output:
    path "fold_meta.txt", emit: meta_file

    script:
    def telescope = params.telescope ?: 'meerkat'
    def template_dir = params.template_dir ?: params.psrfold?.template_dir ?: "${projectDir}/templates"
    def template = "${template_dir}/${telescope}_fold.template"
    def subint_length = params.presto?.fold_subint_length ?: params.psrfold?.subintlength ?: 10.0
    def nsubband = params.presto?.fold_nsub ?: params.psrfold?.nsub ?: 64
    def clfd = params.presto?.fold_clfd ?: params.psrfold?.clfd ?: 2.0
    def nbins = params.presto?.fold_nbins ?: params.psrfold?.nbins ?: 64
    def binplan = params.presto?.fold_binplan ?: params.psrfold?.binplan ?: "0.005 32 0.01 64 0.1 128"
    def threads = params.presto?.fold_threads ?: params.psrfold?.threads ?: 16
    def coherent_dm = cdm ?: 0.0
    """
    # Read filterbank header to get tstart, tsamp, nsamples
    python3 << 'PYEOF'
import struct
import os

def read_fil_header(filepath):
    header = {}
    with open(filepath, 'rb') as f:
        # Skip to find header keywords
        f.seek(0)
        content = f.read(4096)  # Read first 4KB

        # Parse common keywords
        if b'tstart' in content:
            idx = content.find(b'tstart')
            f.seek(idx + 6)
            header['tstart'] = struct.unpack('d', f.read(8))[0]

        if b'tsamp' in content:
            idx = content.find(b'tsamp')
            f.seek(idx + 5)
            header['tsamp'] = struct.unpack('d', f.read(8))[0]

        if b'nsamples' in content:
            idx = content.find(b'nsamples')
            f.seek(idx + 8)
            header['nsamples'] = struct.unpack('i', f.read(4))[0]

    return header

# Try to read header, use defaults if fails
try:
    header = read_fil_header('${fil_file}')
    tstart = header.get('tstart', 0.0)
    tsamp = header.get('tsamp', 64e-6)
    nsamples = header.get('nsamples', 0)
except:
    tstart = 0.0
    tsamp = 64e-6
    nsamples = 0

# Write meta file
with open('fold_meta.txt', 'w') as f:
    f.write(f"filterbank_file: ${fil_file}\\n")
    f.write(f"fft_size: ${fft_size}\\n")
    f.write(f"tsamp: {tsamp}\\n")
    f.write(f"tstart: {tstart}\\n")
    f.write(f"nsamples: {nsamples}\\n")
    f.write(f"source_name: ${cluster}\\n")
    f.write(f"telescope: ${telescope}\\n")
    f.write(f"beam_name: ${beam_name}\\n")
    f.write(f"beam_id: ${beam_id}\\n")
    f.write(f"utc_beam: ${utc_start}\\n")
    f.write(f"ra: ${ra}\\n")
    f.write(f"dec: ${dec}\\n")
    f.write(f"cdm: ${coherent_dm}\\n")
    f.write(f"subint_length: ${subint_length}\\n")
    f.write(f"nsubband: ${nsubband}\\n")
    f.write(f"clfd_q_value: ${clfd}\\n")
    f.write(f"nbins: ${nbins}\\n")
    f.write(f"binplan: ${binplan}\\n")
    f.write(f"threads: ${threads}\\n")
    f.write(f"template: ${template}\\n")
    f.write("start_fraction: 0.0\\n")
    f.write("end_fraction: 1.0\\n")
    f.write("chunk_id: 0\\n")
    f.write("cmask: \\n")
    f.write("rfi_filter: \\n")

print(f"Meta file created with tstart={tstart}, tsamp={tsamp}")
PYEOF
    """
}

/*
 * Peasoup dedispersion-only mode to dump PRESTO-compatible time series
 */
process peasoup_dump_timeseries {
    label 'peasoup'
    label 'long'

    publishDir "${params.output_dir}/${params.target_name}/timeseries/${dm_low}_${dm_high}", mode: 'symlink', saveAs: { filename -> filename.split('/').last() }

    input:
    path input_file
    path dm_file
    path birdies_file

    output:
    path "timeseries_out/*.dat", emit: dat_files
    path "timeseries_out/*.inf", emit: inf_files
    tuple path("timeseries_out/*.dat"), path("timeseries_out/*.inf"), emit: timeseries_data

    script:
    def basename = input_file.baseName
    def ram_limit = params.peasoup?.ram_limit_gb ?: 90
    def ngpus = params.peasoup?.ngpus ?: 1
    def fft_size = params.peasoup?.fft_size ?: 134217728
    """
    mkdir -p timeseries_out

    # Run peasoup with --nosearch and -d to dump time series
    # Peasoup generates both .dat and .inf files automatically
    peasoup -d timeseries_out \\
        -i ${input_file} \\
        --dm_file ${dm_file} \\
        -z ${birdies_file} \\
        --fft_size ${fft_size} \\
        --ram_limit_gb ${ram_limit} \\
        -t ${ngpus}
    """
}

// ============================================================================
// STATE FILE PROCESSES - For Workflow Chaining
// ============================================================================

/*
 * Save RFI state - output from presto_rfi, input to presto_birdies
 */
process save_presto_rfi_state {
    label 'short'
    publishDir "${params.output_dir}/${params.target_name}/presto/state", mode: 'copy'

    input:
    path input_file
    path rfi_mask
    path rfi_stats
    path rfi_inf

    output:
    path "rfi_state.json", emit: state_file

    script:
    """
    python3 << 'PYEOF'
import json
import os

state = {
    "stage": "rfi",
    "input_file": os.path.abspath("${input_file}"),
    "rfi_mask": os.path.abspath("${rfi_mask}"),
    "rfi_stats": os.path.abspath("${rfi_stats}"),
    "rfi_inf": os.path.abspath("${rfi_inf}"),
    "next_workflow": "presto_birdies"
}

with open("rfi_state.json", "w") as f:
    json.dump(state, f, indent=2)
print(f"RFI state saved: {state}")
PYEOF
    """
}

/*
 * Save birdies state - output from presto_birdies, input to presto_dedisperse
 */
process save_presto_birdies_state {
    label 'short'
    publishDir "${params.output_dir}/${params.target_name}/presto/state", mode: 'copy'

    input:
    path input_file
    path rfi_mask
    path zaplist
    path rfi_stats
    path rfi_inf
    path zerodm_inf

    output:
    path "birdies_state.json", emit: state_file

    script:
    """
    python3 << 'PYEOF'
import json
import os

zaplist_path = "${zaplist}"
# Handle NO_ZAPLIST case
if zaplist_path == "NO_ZAPLIST" or not os.path.exists(zaplist_path):
    zaplist_path = None
else:
    zaplist_path = os.path.abspath(zaplist_path)

state = {
    "stage": "birdies",
    "input_file": os.path.abspath("${input_file}"),
    "rfi_mask": os.path.abspath("${rfi_mask}"),
    "rfi_stats": os.path.abspath("${rfi_stats}"),
    "rfi_inf": os.path.abspath("${rfi_inf}"),
    "zerodm_inf": os.path.abspath("${zerodm_inf}"),
    "zaplist": zaplist_path,
    "next_workflow": "presto_dedisperse"
}

with open("birdies_state.json", "w") as f:
    json.dump(state, f, indent=2)
print(f"Birdies state saved: {state}")
PYEOF
    """
}

/*
 * Save dedisperse state - output from presto_dedisperse, input to presto_search
 */
process save_presto_dedisperse_state {
    label 'short'
    publishDir "${params.output_dir}/${params.target_name}/presto/state", mode: 'copy'

    input:
    path input_file
    path rfi_mask
    path rfi_stats
    path zaplist
    path dat_files
    path inf_files

    output:
    path "dedisperse_state.json", emit: state_file

    script:
    """
    python3 << 'PYEOF'
import json
import os
import glob

zaplist_path = "${zaplist}"
if zaplist_path == "NO_ZAPLIST" or not os.path.exists(zaplist_path):
    zaplist_path = None
else:
    zaplist_path = os.path.abspath(zaplist_path)

dat_files = sorted([os.path.abspath(f) for f in glob.glob("*.dat")])
inf_files = sorted([os.path.abspath(f) for f in glob.glob("*.inf")])

state = {
    "stage": "dedisperse",
    "input_file": os.path.abspath("${input_file}"),
    "rfi_mask": os.path.abspath("${rfi_mask}"),
    "rfi_stats": os.path.abspath("${rfi_stats}"),
    "zaplist": zaplist_path,
    "dat_files": dat_files,
    "inf_files": inf_files,
    "next_workflow": "presto_search"
}

with open("dedisperse_state.json", "w") as f:
    json.dump(state, f, indent=2)
print(f"Dedisperse state saved with {len(dat_files)} dat files")
PYEOF
    """
}

/*
 * Save search state - output from presto_search, input to presto_sift_fold
 */
process save_presto_search_state {
    label 'short'
    publishDir "${params.output_dir}/${params.target_name}/presto/state", mode: 'copy'

    input:
    path input_file
    path rfi_mask
    path rfi_stats
    path accel_files

    output:
    path "search_state.json", emit: state_file

    script:
    """
    python3 << 'PYEOF'
import json
import os
import glob

accel_all = sorted([os.path.abspath(f) for f in glob.glob("*_ACCEL_*")])
# Main ACCEL outputs (exclude sidecars like .cand/.txtcand/.inf)
accel_files = [f for f in accel_all if not any(f.endswith(ext) for ext in [".cand", ".txtcand", ".inf"])]
# Retain cand files if needed for compatibility
cand_files = [f for f in accel_all if f.endswith(".cand")]

state = {
    "stage": "search",
    "input_file": os.path.abspath("${input_file}"),
    "rfi_mask": os.path.abspath("${rfi_mask}"),
    "rfi_stats": os.path.abspath("${rfi_stats}"),
    "accel_files": accel_files,
    "cand_files": cand_files,
    "next_workflow": "presto_sift_fold"
}

with open("search_state.json", "w") as f:
    json.dump(state, f, indent=2)
print(f"Search state saved with {len(accel_files)} ACCEL files and {len(cand_files)} cand files")
PYEOF
    """
}

/*
 * Save sift_fold state - output from presto_sift_fold, input to presto_postprocess
 */
process save_presto_sift_fold_state {
    label 'short'
    publishDir "${params.output_dir}/${params.target_name}/presto/state", mode: 'copy'

    input:
    path input_file
    path sifted_csv
    path pfd_files
    path provenance_csv

    output:
    path "sift_fold_state.json", emit: state_file

    script:
    """
    python3 << 'PYEOF'
import json
import os
import glob

pfd_files = sorted([os.path.abspath(f) for f in glob.glob("*.pfd")])

state = {
    "stage": "sift_fold",
    "input_file": os.path.abspath("${input_file}"),
    "sifted_csv": os.path.abspath("${sifted_csv}"),
    "provenance_csv": os.path.abspath("${provenance_csv}"),
    "pfd_files": pfd_files,
    "next_workflow": "presto_postprocess"
}

with open("sift_fold_state.json", "w") as f:
    json.dump(state, f, indent=2)
print(f"Sift/fold state saved with {len(pfd_files)} PFD files")
PYEOF
    """
}

// ============================================================================
// PRESTO WORKFLOWS - Chainable with State Files
// ============================================================================

/*
 * PRESTO RFI Detection workflow
 * Input:  input_files (filterbank) OR params.input_fil
 * Output: rfi_state.json -> presto_birdies
 */
workflow presto_rfi {
    take:
    input_files  // Channel of filterbank files

    main:
    // Run RFI detection
    presto_rfifind(input_files)
        .set { rfi_out }

    // Save state file for chaining
    save_presto_rfi_state(input_files, rfi_out.rfi_mask, rfi_out.rfi_stats, rfi_out.rfi_inf)
        .set { state_out }

    emit:
    state_file = state_out.state_file  // rfi_state.json -> presto_birdies
    rfi_mask = rfi_out.rfi_mask
    rfi_files = rfi_out.rfi_files
    rfi_stats = rfi_out.rfi_stats
    rfi_inf = rfi_out.rfi_inf
}

/*
 * PRESTO Birdie Detection workflow (zero-DM search for RFI identification)
 * Input:  (input_file, rfi_mask) OR params.state_file (rfi_state.json)
 * Output: birdies_state.json -> presto_dedisperse
 */
workflow presto_birdies {
    take:
    input_file  // Filterbank file
    rfi_mask    // RFI mask from presto_rfi
    rfi_stats   // rfifind stats from presto_rfi
    rfi_inf     // rfifind inf from presto_rfi

    main:
    presto_prepdata_zerodm(input_file, rfi_mask, rfi_stats, rfi_inf)
        .set { zerodm_out }

    presto_accelsearch_zerodm(zerodm_out.dat_file, zerodm_out.inf_file)
        .set { birdie_out }

    // Handle empty zaplist
    zaplist_ch = birdie_out.zaplist.ifEmpty(file('NO_ZAPLIST'))

    // Save state file for chaining
    save_presto_birdies_state(input_file, rfi_mask, zaplist_ch, rfi_stats, rfi_inf, zerodm_out.inf_file)
        .set { state_out }

    emit:
    state_file = state_out.state_file  // birdies_state.json -> presto_dedisperse
    zaplist = birdie_out.zaplist
    rfi_mask_passthrough = rfi_mask
    rfi_stats_passthrough = rfi_stats
    rfi_inf_passthrough = rfi_inf
    zerodm_inf = zerodm_out.inf_file  // .inf file from prepdata zero-DM for prepsubband
}

/*
 * PRESTO Dedispersion workflow
 * Input:  (input_file, input_inf, rfi_mask, zaplist, dm_ranges) OR params.state_file (birdies_state.json)
 * Output: dedisperse_state.json -> presto_search
 */
workflow presto_dedisperse {
    take:
    input_file  // Filterbank file
    input_inf   // .inf file from prepdata zero-DM
    rfi_mask    // RFI mask
    rfi_stats   // RFI stats
    zaplist     // Birdie zaplist (or NO_ZAPLIST)
    dm_ranges   // Channel of tuples: (dm_low, dm_high, dm_step, downsamp)

    main:
    presto_prepsubband(input_file, input_inf, rfi_mask, rfi_stats, dm_ranges)
        .set { subband_out }

    // Collect all dat/inf files for state
    dat_collected = subband_out.dat_files.collect()
    inf_collected = subband_out.inf_files.collect()

    // Save state file for chaining
    save_presto_dedisperse_state(
        input_file, rfi_mask, rfi_stats, zaplist.ifEmpty(file('NO_ZAPLIST')),
        dat_collected, inf_collected
    ).set { state_out }

    emit:
    state_file = state_out.state_file  // dedisperse_state.json -> presto_search
    dat_files = subband_out.dat_files
    inf_files = subband_out.inf_files
    subband_data = subband_out.subband_data
}

/*
 * PRESTO Acceleration Search workflow
 * Input:  (input_file, rfi_mask, dat_files, inf_files, zaplist) OR params.state_file
 * Output: search_state.json -> presto_sift_fold
 */
workflow presto_search {
    take:
    input_file  // Original filterbank (for state tracking)
    rfi_mask    // RFI mask (for state tracking)
    rfi_stats   // RFI stats (for state tracking)
    dat_files   // Dedispersed time series
    inf_files   // Info files
    zaplist     // Birdie zaplist

    main:
    // Handle case where zaplist may not exist
    zaplist_ch = zaplist.ifEmpty(file('NO_ZAPLIST'))

    presto_accelsearch(dat_files, inf_files, zaplist_ch)
        .set { accel_out }

    // Collect accel files for state
    accel_collected = accel_out.accel_files.collect()

    // Save state file for chaining
    save_presto_search_state(input_file, rfi_mask, rfi_stats, accel_collected)
        .set { state_out }

    emit:
    state_file = state_out.state_file  // search_state.json -> presto_sift_fold
    accel_files = accel_out.accel_files
    cand_files = accel_out.cand_files
}

/*
 * PRESTO Sift and Fold workflow
 * Input:  (input_file, rfi_mask, accel_files) OR params.state_file (search_state.json)
 * Output: sift_fold_state.json -> presto_postprocess
 */
workflow presto_sift_fold {
    take:
    input_file   // Original filterbank
    rfi_mask     // RFI mask
    rfi_stats    // RFI stats
    accel_files  // ACCEL search results
    cand_files   // ACCEL cand files

    main:
    cand_files.collect().set { cand_files_all }

    presto_sift_candidates(cand_files_all, input_file)
        .set { sifted_out }

    presto_prepfold_batch(input_file, sifted_out.sifted_csv, rfi_mask, rfi_stats)
        .set { fold_out }

    // Collect pfd files for state
    pfd_collected = fold_out.pfd_files.flatten().collect()

    // Save state file for chaining
    save_presto_sift_fold_state(input_file, sifted_out.sifted_csv, pfd_collected, sifted_out.provenance_csv)
        .set { state_out }

    emit:
    state_file = state_out.state_file  // sift_fold_state.json -> presto_postprocess
    pfd_files = fold_out.pfd_files
    ps_files = fold_out.ps_files       // PostScript files for PNG conversion
    sifted_csv = sifted_out.sifted_csv
    provenance_csv = sifted_out.provenance_csv
}

/*
 * PRESTO Post-processing workflow (PNG conversion, PICS classification, and tarball creation)
 * Input:  (input_file, sifted_csv, pfd_files, ps_files) OR params.state_file (sift_fold_state.json)
 * Output: Final tarball and CSV (tarball contains only PNG + CSV, no PFD files)
 *
 * Pipeline: PFD -> PNG -> Merge -> PICS Classification (on PFD) -> Tarball
 */
workflow presto_postprocess {
    take:
    input_file   // Original filterbank
    sifted_csv   // Sifted candidates CSV
    pfd_files    // Folded profile files (.pfd) - needed for PICS
    ps_files     // PostScript files (.pfd.ps) - unused when using show_pfd
    provenance_csv // Candidate -> ACCEL provenance
    meta_info   // Metadata tuple for candidates.csv

    main:
    // Normalize PFD channel (process outputs are lists)
    def pfd_flat = pfd_files.flatten()

    // Step 1: Convert PFD to PNG using show_pfd
    presto_pfd_to_png(pfd_flat.collect())
        .set { png_out }

    // Collect bestprof files (may be empty)
    bestprof_ch = pfd_flat.map { pfd ->
        file("${pfd}.bestprof")
    }.filter { it.exists() }.collect().ifEmpty([])

    // Step 2: Merge fold results with search candidates
    presto_fold_merge(
        pfd_flat.collect(),
        png_out.png_files.collect(),
        bestprof_ch,
        sifted_csv,
        provenance_csv,
        meta_info
    ).set { merged_out }

    // Step 3: Run PICS classification on PFD files
    presto_pics_classifier(
        merged_out.all_pfd,
        merged_out.merged_csv
    ).set { pics_out }

    // Step 4: Create tarball (PNG + CSV only, no PFD files)
    // Use tarball prefix from config
    def tarball_prefix = params.tarball_prefix ?: params.target_name ?: input_file.baseName

    presto_create_tarball(
        pics_out.classified_csv,
        merged_out.all_png,
        tarball_prefix
    ).set { tarball_out }

    emit:
    tarball = tarball_out.tarball
    final_csv = tarball_out.final_csv
    pics_filtered_csv = tarball_out.pics_filtered_csv
    png_files = png_out.png_files
}

// ============================================================================
// UNIFIED PRESTO PIPELINE ENTRY POINT
// ============================================================================

/*
 * Unified PRESTO Pipeline with Stage Control
 *
 * This workflow replaces the 6 standalone entry points with a single
 * entry point controlled by parameters.
 *
 * Stages:
 *   1 = RFI detection (rfifind)
 *   2 = Birdie detection (zero-DM accelsearch)
 *   3 = Dedispersion (prepsubband)
 *   4 = Acceleration search (accelsearch)
 *   5 = Sift and fold candidates
 *   6 = Post-processing (PNG + tarball)
 *
 * Usage Examples:
 *   # Full pipeline (all 6 stages)
 *   nextflow run elden.nf -entry presto_pipeline --input_fil /path/to/file.fil
 *
 *   # Run only stages 1-3 (RFI through Dedispersion)
 *   nextflow run elden.nf -entry presto_pipeline --input_fil file.fil --presto.end_stage 3
 *
 *   # Resume from stage 4 using state file
 *   nextflow run elden.nf -entry presto_pipeline --presto.start_stage 4 --state_file dedisperse_state.json
 *
 * Parameters:
 *   --input_fil           Input filterbank file (required for start_stage 1)
 *   --state_file          State file to resume from (required for start_stage > 1)
 *   --presto.start_stage  First stage to run (default: 1)
 *   --presto.end_stage    Last stage to run (default: 6)
 */
workflow presto_pipeline {
    main:
    // Get stage parameters
    start_stage = params.presto?.start_stage ?: 1
    end_stage = params.presto?.end_stage ?: 6

    // Validate parameters
    if (start_stage < 1 || start_stage > 6) {
        error "Invalid start_stage: ${start_stage}. Must be between 1 and 6."
    }
    if (end_stage < 1 || end_stage > 6) {
        error "Invalid end_stage: ${end_stage}. Must be between 1 and 6."
    }
    if (start_stage > end_stage) {
        error "start_stage (${start_stage}) cannot be greater than end_stage (${end_stage})"
    }

    // Create DM ranges channel at workflow level (outside conditionals)
    // Parse dm_ranges if it's a string (from command line), otherwise use as-is
    def dm_ranges_list = params.presto?.dm_ranges ?: []
    if (dm_ranges_list instanceof String) {
        dm_ranges_list = new groovy.json.JsonSlurper().parseText(dm_ranges_list)
    }

    dm_ranges_ch = channel.from(dm_ranges_list).map { range ->
        tuple(range.dm_low, range.dm_high, range.dm_step, range.downsamp)
    }

    // Load state file if resuming from a previous stage
    if (start_stage > 1) {
        if (!params.state_file) {
            error "State file required when starting from stage ${start_stage}. Use --state_file parameter."
        }
        log.info "Loading state from: ${params.state_file}"
        state = new groovy.json.JsonSlurper().parse(file(params.state_file))
    } else {
        state = null
    }

    // Build metadata tuple for candidates.csv
    def input_fil_path = params.input_fil ?: (state ? state.input_file : null)
    def default_meta = tuple(
        params.get('pointing_id', "0"),
        params.runID ?: "",
        params.target_name ?: (input_fil_path ? new File(input_fil_path).name.replaceFirst(/\\.fil$/, "") : ""),
        params.get('beam_id', "0"),
        params.get('utc_start', ""),
        params.get('ra', ""),
        params.get('dec', ""),
        params.get('cdm', "0.0"),
        input_fil_path ?: ""
    )

    meta_ch = Channel.value(default_meta)
    if (params.files_list && input_fil_path) {
        def input_base = new File(input_fil_path).name
        meta_ch = channel.fromPath(params.files_list)
            .splitCsv(header: true, sep: ',')
            .filter { row ->
                def row_path = row.fits_files?.trim()
                row_path == input_fil_path || (row_path && new File(row_path).name == input_base)
            }
            .map { row ->
                tuple(
                    row.pointing.trim(),
                    row.cluster.trim(),
                    row.beam_name.trim(),
                    row.beam_id.trim(),
                    row.utc_start.trim().replace(" ", "-"),
                    row.ra.trim(),
                    row.dec.trim(),
                    row.cdm.trim(),
                    row.fits_files.trim()
                )
            }
            .ifEmpty(Channel.value(default_meta))
    }

    // Stage 1: RFI Detection
    if (start_stage <= 1 && end_stage >= 1) {
        log.info "Running Stage 1: RFI Detection"

        if (!params.input_fil) {
            error "Input filterbank file required for stage 1. Use --input_fil parameter."
        }

        input_ch = channel.fromPath(params.input_fil)
        rfi_out = presto_rfi(input_ch)

        rfi_mask_ch = rfi_out.rfi_mask
        rfi_stats_ch = rfi_out.rfi_stats
        rfi_inf_ch = rfi_out.rfi_inf

        if (end_stage == 1) {
            log.info "Pipeline complete. RFI detection finished."
        }
    }

    // Stage 2: Birdie Detection
    if (start_stage <= 2 && end_stage >= 2) {
        log.info "Running Stage 2: Birdie Detection"

        // Load from state if starting at stage 2
        if (start_stage == 2) {
            input_ch = channel.fromPath(state.input_file)
            rfi_mask_ch = channel.fromPath(state.rfi_mask)
            if (!state.rfi_stats || !state.rfi_inf) {
                error "State file missing rfi_stats or rfi_inf; rerun stage 1 to regenerate rfifind outputs."
            }
        rfi_stats_ch = channel.fromPath(state.rfi_stats)
        rfi_inf_ch = channel.fromPath(state.rfi_inf)
        if (!state.zerodm_inf) {
            error "State file missing zerodm_inf; rerun stage 2 to regenerate prepdata outputs."
        }
        zerodm_inf_ch = channel.fromPath(state.zerodm_inf)
    }

    birdie_out = presto_birdies(input_ch, rfi_mask_ch, rfi_stats_ch, rfi_inf_ch)

    zaplist_ch = birdie_out.zaplist
    zerodm_inf_ch = birdie_out.zerodm_inf

        if (end_stage == 2) {
            log.info "Pipeline complete. Birdie detection finished."
        }
    }

    // Stage 3: Dedispersion
    if (start_stage <= 3 && end_stage >= 3) {
        log.info "Running Stage 3: Dedispersion"

        // Load from state if starting at stage 3
        if (start_stage == 3) {
            input_ch = channel.fromPath(state.input_file)
            rfi_mask_ch = channel.fromPath(state.rfi_mask)
            zaplist_ch = state.zaplist ? channel.fromPath(state.zaplist) : channel.value(file('NO_ZAPLIST'))
            if (!state.zerodm_inf) {
                error "State file missing zerodm_inf; rerun stage 2 to regenerate prepdata outputs."
            }
            zerodm_inf_ch = channel.fromPath(state.zerodm_inf)
            if (!state.rfi_stats) {
                error "State file missing rfi_stats; rerun stage 1 to regenerate rfifind outputs."
            }
            rfi_stats_ch = channel.fromPath(state.rfi_stats)
        }

        dedisperse_out = presto_dedisperse(input_ch, zerodm_inf_ch, rfi_mask_ch, rfi_stats_ch, zaplist_ch, dm_ranges_ch)

        dat_ch = dedisperse_out.dat_files
        inf_ch = dedisperse_out.inf_files

        if (end_stage == 3) {
            log.info "Pipeline complete. Dedispersion finished."
        }
    }

    // Stage 4: Acceleration Search
    if (start_stage <= 4 && end_stage >= 4) {
        log.info "Running Stage 4: Acceleration Search"

        // Load from state if starting at stage 4
        if (start_stage == 4) {
            input_ch = channel.fromPath(state.input_file)
            rfi_mask_ch = channel.fromPath(state.rfi_mask)
            rfi_stats_ch = channel.fromPath(state.rfi_stats)
            zaplist_ch = state.zaplist ? channel.fromPath(state.zaplist) : channel.value(file('NO_ZAPLIST'))
            dat_ch = channel.fromPath(state.dat_files).flatten()
            inf_ch = channel.fromPath(state.inf_files).flatten()
        }

        search_out = presto_search(input_ch, rfi_mask_ch, rfi_stats_ch, dat_ch, inf_ch, zaplist_ch)

        accel_ch = search_out.accel_files
        cand_ch = search_out.accel_files  // use main ACCEL files for sifting

        if (end_stage == 4) {
            log.info "Pipeline complete. Acceleration search finished."
        }
    }

    // Stage 5: Sift and Fold
    if (start_stage <= 5 && end_stage >= 5) {
        log.info "Running Stage 5: Sift and Fold"

        // Load from state if starting at stage 5
        if (start_stage == 5) {
            input_ch = channel.fromPath(state.input_file)
            rfi_mask_ch = channel.fromPath(state.rfi_mask)
            rfi_stats_ch = channel.fromPath(state.rfi_stats)
            accel_ch = channel.fromPath(state.accel_files)
            cand_ch = accel_ch
        }

        sift_fold_out = presto_sift_fold(input_ch, rfi_mask_ch, rfi_stats_ch, accel_ch, cand_ch)

        sifted_csv_ch = sift_fold_out.sifted_csv
        provenance_ch = sift_fold_out.provenance_csv
        pfd_ch = sift_fold_out.pfd_files
        ps_ch = sift_fold_out.ps_files

        if (end_stage == 5) {
            log.info "Pipeline complete. Sifting and folding finished."
        }
    }

    // Stage 6: Post-processing
    if (start_stage <= 6 && end_stage >= 6) {
        log.info "Running Stage 6: Post-processing"

        // Load from state if starting at stage 6
        if (start_stage == 6) {
            input_ch = channel.fromPath(state.input_file)
            sifted_csv_ch = channel.fromPath(state.sifted_csv)
            provenance_ch = channel.fromPath(state.provenance_csv)
            pfd_ch = Channel.from(state.pfd_files).map { file(it) }
            // PS files are in same directory as PFD files, just with .ps extension
            ps_ch = Channel.from(state.pfd_files).map { file("${it}.ps") }
        }

        presto_postprocess(input_ch, sifted_csv_ch, pfd_ch, ps_ch, provenance_ch, meta_ch)

        log.info "Pipeline complete. All stages finished."
    }
}

// ============================================================================
// FULL PRESTO PIPELINES
// ============================================================================

/*
 * Full PRESTO search pipeline
 * Runs: RFI detection -> Birdie detection -> Dedispersion -> Acceleration search -> Sifting -> Folding -> Post-processing
 *
 * Set params.presto.fold_backend = 'presto' (default) to fold with prepfold
 * Set params.presto.fold_backend = 'pulsarx' to fold with psrfold_fil2
 */
workflow presto_full {
    take:
    intake_ch  // Channel from intake()

    main:
    // Extract filterbank files and metadata
    def fil_channel = intake_ch.map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        file(f)
    }

    // Also keep the full metadata for PulsarX folding
    def meta_channel = intake_ch.map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        tuple(p, f, c, bn, bi, u, ra, dec, cdm)
    }
    def meta_info_ch = meta_channel.map { p, f, c, bn, bi, u, ra, dec, cdm ->
        tuple(p, c, bn, bi, u, ra, dec, cdm, f)
    }

    // Run RFI detection
    presto_rfi(fil_channel)
        .set { rfi_out }

    // Run birdie detection
    presto_birdies(fil_channel, rfi_out.rfi_mask, rfi_out.rfi_stats, rfi_out.rfi_inf)
        .set { birdie_out }

    // Generate DM ranges from params
    def dm_ranges = Channel.from(params.presto.dm_ranges).map { range ->
        tuple(range.dm_low, range.dm_high, range.dm_step, range.downsamp)
    }

    // Combine input file, zerodm inf, and RFI artifacts for dedispersion
    def dedisperse_input = fil_channel
        .zip(birdie_out.zerodm_inf)
        .zip(rfi_out.rfi_mask)
        .zip(rfi_out.rfi_stats)
        .map { tuple4 -> tuple4[0][0][0], tuple4[0][0][1], tuple4[0][1], tuple4[1] }

    // Run dedispersion for each DM range
    presto_prepsubband(
        dedisperse_input.map { it[0] },   // input file
        dedisperse_input.map { it[1] },   // zerodm inf
        dedisperse_input.map { it[2] },   // rfi mask
        dedisperse_input.map { it[3] },   // rfi stats
        dm_ranges
    ).set { subband_out }

    // Run acceleration search on all dedispersed data
    presto_accelsearch(
        subband_out.dat_files.flatten(),
        subband_out.inf_files.flatten(),
        birdie_out.zaplist.ifEmpty(file('NO_ZAPLIST'))
    ).set { accel_out }

    // Collect all ACCEL files and sift
    presto_sift_candidates(
        accel_out.accel_files.collect(),
        fil_channel
    ).set { sifted_out }

    // Choose folding backend based on params.presto.fold_backend
    def fold_backend = params.presto?.fold_backend ?: 'presto'

    if (fold_backend == 'pulsarx') {
        // Fold with PulsarX (psrfold_fil2) using pulsarx_fold.py --csv
        // Generate meta file with observation parameters
        def fft_size = params.presto?.fft_size ?: 134217728
        generate_fold_meta(meta_channel, fft_size, fold_backend)
            .set { meta_out }

        // Fold with PulsarX
        presto_fold_pulsarx(
            fil_channel,
            sifted_out.sifted_csv,
            meta_out.meta_file
        ).set { fold_out }

        // Note: PulsarX folding produces .ar and .png files directly
        // No additional tarball creation needed for PulsarX output
    } else {
        // Default: Fold with PRESTO prepfold
        // prepfold doesn't need a meta file - it uses command-line parameters
        // The period correction is handled by the presto_prepfold_batch process
        presto_prepfold_batch(
            fil_channel,
            sifted_out.sifted_csv,
            rfi_out.rfi_mask,
            rfi_out.rfi_stats
        ).set { fold_out }

        // Post-processing: PNG conversion from PFD files using show_pfd
        presto_pfd_to_png(fold_out.pfd_files.flatten().collect())
            .set { png_out }

        // Create merged results
        presto_fold_merge(
            fold_out.pfd_files.flatten().collect(),
            png_out.png_files.collect(),
            fold_out.bestprof_files.flatten().collect().ifEmpty([]),
            sifted_out.sifted_csv,
            sifted_out.provenance_csv,
            meta_info_ch
        ).set { merged_out }

        // Run PICS classification on PFD files (PICS requires .pfd files)
        presto_pics_classifier(
            merged_out.all_pfd,
            merged_out.merged_csv
        ).set { pics_out }

        // Create tarball (PNG + CSV only, no PFD files)
        def tarball_prefix = params.tarball_prefix ?: params.target_name ?: "presto_full"
        presto_create_tarball(
            pics_out.classified_csv,
            merged_out.all_png,
            tarball_prefix
        ).set { tarball_out }
    }
}

/*
 * PRESTO search on pre-cleaned filterbanks
 * Similar to run_search_fold but using PRESTO instead of peasoup+pulsarx
 */
workflow run_presto_search {
    take:
    cleaned_ch  // Channel from intake()

    main:
    def fil_ch = cleaned_ch.map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        file(f)
    }
    def meta_info_ch = cleaned_ch.map { p, f, c, bn, bi, u, ra, dec, cdm, fname ->
        tuple(p, c, bn, bi, u, ra, dec, cdm, f)
    }

    // Run RFI detection
    presto_rfi(fil_ch)
        .set { rfi_out }

    // Run birdie detection
    presto_birdies(fil_ch, rfi_out.rfi_mask, rfi_out.rfi_stats, rfi_out.rfi_inf)
        .set { birdie_out }

    // Generate DM ranges
    def dm_ranges = Channel.from(params.presto.dm_ranges).map { range ->
        tuple(range.dm_low, range.dm_high, range.dm_step, range.downsamp)
    }

    // Run dedispersion
    presto_prepsubband(fil_ch, birdie_out.zerodm_inf, rfi_out.rfi_mask, rfi_out.rfi_stats, dm_ranges)
        .set { subband_out }

    // Run acceleration search
    presto_accelsearch(
        subband_out.dat_files.flatten(),
        subband_out.inf_files.flatten(),
        birdie_out.zaplist.ifEmpty(file('NO_ZAPLIST'))
    ).set { accel_out }

    // Sift candidates
    presto_sift_candidates(accel_out.accel_files.collect(), fil_ch)
        .set { sifted_out }

    // Fold candidates
    presto_prepfold_batch(fil_ch, sifted_out.sifted_csv, rfi_out.rfi_mask, rfi_out.rfi_stats)
        .set { fold_out }

    // Convert PFD files to PNG using show_pfd
    presto_pfd_to_png(fold_out.pfd_files.flatten().collect())
        .set { png_out }

    // Merge fold results
    presto_fold_merge(
        fold_out.pfd_files.flatten().collect(),
        png_out.png_files.collect(),
        fold_out.bestprof_files.flatten().collect().ifEmpty([]),
        sifted_out.sifted_csv,
        sifted_out.provenance_csv,
        meta_info_ch
    ).set { merged_out }

    // Run PICS classification on PFD files (PICS requires .pfd files)
    presto_pics_classifier(
        merged_out.all_pfd,
        merged_out.merged_csv
    ).set { pics_out }

    // Create tarball (PNG + CSV only, no PFD files)
    def tarball_prefix = params.tarball_prefix ?: params.target_name ?: "presto_search"
    presto_create_tarball(
        pics_out.classified_csv,
        merged_out.all_png,
        tarball_prefix
    ).set { tarball_out }
}

/*
 * Peasoup time series dumping workflow
 * Dumps PRESTO-compatible .dat/.inf files from peasoup for subsequent PRESTO processing
 */
workflow peasoup_timeseries_dump {
    take:
    input_channel  // tuple(pointing, fil_file, cluster, beam_name, beam_id, utc, ra, dec, cdm)
    dm_file
    birdies_file

    main:
    def fil_channel = input_channel.map { p, f, c, bn, bi, u, ra, dec, cdm ->
        file(f)
    }

    peasoup_dump_timeseries(fil_channel, dm_file, birdies_file)
        .set { timeseries_out }

    emit:
    dat_files = timeseries_out.dat_files
    inf_files = timeseries_out.inf_files
    timeseries_data = timeseries_out.timeseries_data
}

/*
 * PRESTO search on peasoup-dumped time series
 * This workflow creates a separate PRESTO tarball with PNG + CSV (no PFD files)
 */
workflow presto_on_peasoup_timeseries {
    take:
    dat_files   // Channel of .dat files from peasoup
    inf_files   // Channel of .inf files from peasoup
    fil_channel // Original filterbank files for folding
    meta_channel // Metadata channel for PulsarX folding

    main:
    // Run acceleration search on peasoup time series
    presto_accelsearch(
        dat_files.flatten(),
        inf_files.flatten(),
        file('NO_ZAPLIST')  // No zaplist for peasoup time series
    ).set { accel_out }

    // Sift candidates
    presto_sift_candidates(accel_out.accel_files.collect(), fil_channel)
        .set { sifted_out }

    // Choose folding backend
    def fold_backend = params.presto?.fold_backend ?: 'pulsarx'
    def tarball_prefix = params.tarball_prefix ?: params.target_name ?: "hybrid"

    def meta_info_ch = meta_channel.map { p, f, c, bn, bi, u, ra, dec, cdm ->
        tuple(p, c, bn, bi, u, ra, dec, cdm, f)
    }

    if (fold_backend == 'pulsarx') {
        // Generate meta file for PulsarX folding
        def fft_size = params.peasoup?.fft_size ?: 134217728
        generate_fold_meta(meta_channel, fft_size, fold_backend)
            .set { meta_out }

        // Fold with PulsarX
        presto_fold_pulsarx(fil_channel, sifted_out.sifted_csv, meta_out.meta_file)
            .set { fold_out }

        // For PulsarX fold output, we need to run PICS on .ar files
        // and create a separate PRESTO tarball
        // Note: PulsarX outputs .ar files, not .pfd files
        log.info "PRESTO candidates folded with PulsarX. Output: .ar files in publishDir"

    } else {
        // Fold with prepfold (creates .pfd, .pfd.ps, .pfd.bestprof files)
        presto_prepfold_batch(fil_channel, sifted_out.sifted_csv, file('NO_MASK'), file('NO_STATS'))
            .set { fold_out }

        // Post-processing: PNG conversion from PFD files using show_pfd
        presto_pfd_to_png(fold_out.pfd_files.flatten().collect())
            .set { png_out }

        // Merge fold results
        presto_fold_merge(
            fold_out.pfd_files.flatten().collect(),
            png_out.png_files.collect(),
            fold_out.bestprof_files.flatten().collect().ifEmpty([]),
            sifted_out.sifted_csv,
            sifted_out.provenance_csv,
            meta_info_ch
        ).set { merged_out }

        // Run PICS classification on PFD files (PICS requires .pfd files)
        presto_pics_classifier(
            merged_out.all_pfd,
            merged_out.merged_csv
        ).set { pics_out }

        // Create PRESTO tarball (PNG + CSV only, no PFD files)
        presto_create_tarball(
            pics_out.classified_csv,
            merged_out.all_png,
            tarball_prefix
        ).set { tarball_out }
    }

    emit:
    sifted_csv = sifted_out.sifted_csv
}

/*
 * Helper workflow to run accelsearch on a single set of time series files
 */
workflow run_accelsearch_single {
    take:
    dat_channel   // Channel of .dat files
    inf_channel   // Channel of .inf files
    fil_channel   // Filterbank file channel
    meta_channel  // Metadata channel for PulsarX

    main:
    // Run acceleration search on pre-dumped time series
    presto_accelsearch(
        dat_channel.flatten(),
        inf_channel.flatten(),
        file('NO_ZAPLIST')  // No zaplist for pre-dumped time series
    ).set { accel_out }

    // Sift candidates
    presto_sift_candidates(accel_out.accel_files.collect(), fil_channel)
        .set { sifted_out }

    // Choose folding backend
    def fold_backend = params.presto?.fold_backend ?: 'pulsarx'

    def meta_info_ch = meta_channel.map { p, f, c, bn, bi, u, ra, dec, cdm ->
        tuple(p, c, bn, bi, u, ra, dec, cdm, f)
    }

    if (fold_backend == 'pulsarx') {
        // Generate meta file for PulsarX folding
        def fft_size = params.presto?.fft_size ?: params.peasoup?.fft_size ?: 134217728
        generate_fold_meta(meta_channel, fft_size, fold_backend)
            .set { meta_out }

        // Fold with PulsarX
        presto_fold_pulsarx(fil_channel, sifted_out.sifted_csv, meta_out.meta_file)
            .set { fold_out }

        log.info "Folding complete with PulsarX. Output files are in the publishDir."

    } else {
        // Fold with prepfold (creates .pfd, .pfd.ps, .pfd.bestprof files)
        presto_prepfold_batch(fil_channel, sifted_out.sifted_csv, file('NO_MASK'), file('NO_STATS'))
            .set { fold_out }

        // Post-processing: PNG conversion from PFD files using show_pfd
        presto_pfd_to_png(fold_out.pfd_files.flatten().collect())
            .set { png_out }

        // Create merged results
        presto_fold_merge(
            fold_out.pfd_files.flatten().collect(),
            png_out.png_files.collect(),
            fold_out.bestprof_files.flatten().collect().ifEmpty([]),
            sifted_out.sifted_csv,
            sifted_out.provenance_csv,
            meta_info_ch
        ).set { merged_out }

        // Run PICS classification on PFD files (PICS requires .pfd files)
        presto_pics_classifier(
            merged_out.all_pfd,
            merged_out.merged_csv
        ).set { pics_out }

        // Create tarball (PNG + CSV only, no PFD files)
        def tarball_prefix = params.tarball_prefix ?: params.target_name ?: "accelsearch"
        presto_create_tarball(
            pics_out.classified_csv,
            merged_out.all_png,
            tarball_prefix
        ).set { tarball_out }
    }

    emit:
    sifted_csv = sifted_out.sifted_csv
}

/*
 * Entry point for running accelsearch on pre-dumped time series from peasoup
 */
workflow run_accelsearch_on_timeseries {
    main:
    // Determine which mode to use
    def single_file_mode = params.timeseries_input_dir && params.filterbank_file
    def multi_file_mode = params.files_list && params.basedir && params.runID

    if (!single_file_mode && !multi_file_mode) {
        error """ERROR: Missing required parameters.

For multi-file mode (recommended), provide:
  - params.files_list  (same input CSV as first run)
  - params.basedir     (same basedir as first run)
  - params.runID       (same runID as first run)

For single-file mode, provide:
  - params.timeseries_input_dir  (directory with .dat/.inf files)
  - params.filterbank_file       (original filterbank file)
"""
    }

    if (single_file_mode) {
        // ===== SINGLE FILE MODE =====
        def dat_channel = Channel.fromPath("${params.timeseries_input_dir}/*.dat")
        def inf_channel = Channel.fromPath("${params.timeseries_input_dir}/*.inf")
        def fil_channel = Channel.fromPath(params.filterbank_file)

        def meta_channel = fil_channel.map { fil ->
            def basename = fil.baseName
            tuple(
                basename,                    // pointing
                fil,                         // filterbank path
                params.cluster ?: 'unknown', // cluster
                basename,                    // beam_name
                '0',                         // beam_id
                'unknown',                   // utc
                '00:00:00.0',               // ra
                '+00:00:00.0',              // dec
                params.cdm ?: '0.0'         // cdm
            )
        }

        run_accelsearch_single(dat_channel, inf_channel, fil_channel, meta_channel)

    } else {
        // ===== MULTI-FILE MODE =====
        // Parse input CSV (same format as first run)
        def input_channel = Channel.fromPath("${params.files_list}")
            .splitCsv(header: true, sep: ',')
            .map { row ->
                def pointing = row.pointing.trim()
                def fits_files = row.fits_files.trim()
                def cluster = row.cluster.trim()
                def beam_name = row.beam_name.trim()
                def beam_id = row.beam_id.trim()
                def utc_start = row.utc_start.trim().replace(" ", "-")
                def ra = row.ra.trim()
                def dec = row.dec.trim()
                def cdm = row.cdm.trim()
                tuple(pointing, fits_files, cluster, beam_name, beam_id, utc_start, ra, dec, cdm)
            }

        // For each input file, find its time series directory
        def timeseries_channel = input_channel.flatMap { pointing, fits_file, cluster, beam_name, beam_id, utc, ra, dec, cdm ->
            def base_path = "${params.basedir}/${params.runID}/${beam_name}"
            def timeseries_dirs = []

            def base_dir = file(base_path)
            if (base_dir.exists()) {
                base_dir.eachDirRecurse { dir ->
                    if (dir.name == 'TIMESERIES' && dir.isDirectory()) {
                        def dat_files = file("${dir}/*.dat")
                        def inf_files = file("${dir}/*.inf")
                        if (dat_files || inf_files) {
                            timeseries_dirs << tuple(
                                pointing, fits_file, cluster, beam_name, beam_id, utc, ra, dec, cdm,
                                dir.toString()
                            )
                        }
                    }
                }
            }

            if (timeseries_dirs.isEmpty()) {
                log.warn "No TIMESERIES directories found for beam ${beam_name} in ${base_path}"
            }

            return timeseries_dirs
        }

        // Process each beam's time series
        timeseries_channel.map { pointing, fits_file, cluster, beam_name, beam_id, utc, ra, dec, cdm, ts_dir ->
            def dat_files = file("${ts_dir}/*.dat")
            def inf_files = file("${ts_dir}/*.inf")
            tuple(
                pointing, file(fits_file), cluster, beam_name, beam_id, utc, ra, dec, cdm,
                dat_files, inf_files, ts_dir
            )
        }.set { beam_timeseries }

        // Process each beam separately
        beam_timeseries.each { pointing, fil_file, cluster, beam_name, beam_id, utc, ra, dec, cdm, dat_files, inf_files, ts_dir ->
            log.info "Processing time series for beam ${beam_name} from ${ts_dir}"

            def dat_ch = Channel.fromList(dat_files instanceof List ? dat_files : [dat_files])
            def inf_ch = Channel.fromList(inf_files instanceof List ? inf_files : [inf_files])
            def fil_ch = Channel.of(fil_file)
            def meta_ch = Channel.of(tuple(pointing, fil_file, cluster, beam_name, beam_id, utc, ra, dec, cdm))

            run_accelsearch_single(dat_ch, inf_ch, fil_ch, meta_ch)
        }
    }
}
