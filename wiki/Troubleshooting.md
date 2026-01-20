# Troubleshooting

This page covers common issues and their solutions when running ELDEN-RING.

## General Issues

### Pipeline Fails to Start

**Symptom**: Nextflow throws an error immediately without processing any data.

**Common causes and solutions**:

1. **Missing configuration file**
   ```
   Error: Cannot find configuration file: params.config
   ```
   **Solution**: Ensure `params.config` exists in the working directory or provide the full path:
   ```bash
   nextflow run elden.nf -entry full -c /full/path/to/params.config
   ```

2. **Invalid Nextflow version**
   ```
   Error: Nextflow version 21.10.0 or later is required
   ```
   **Solution**: Update Nextflow:
   ```bash
   nextflow self-update
   ```

3. **Invalid entry point**
   ```
   Error: Unknown workflow entry point: fulll
   ```
   **Solution**: Check spelling. Valid entries: `full`, `run_dada_search`, `run_rfi_clean`, etc.

---

### Container Issues

**Symptom**: Process fails with container-related errors.

1. **Container not found**
   ```
   Error: Container image not found: peasoup.sif
   ```
   **Solution**:
   - Verify container paths in configuration
   - Check `params.containers.*` settings
   - Ensure Singularity cache is accessible

2. **Singularity permission denied**
   ```
   FATAL: could not open image /path/to/container.sif: permission denied
   ```
   **Solution**:
   ```bash
   # Check file permissions
   ls -la /path/to/container.sif

   # Fix permissions if needed
   chmod 644 /path/to/container.sif
   ```

3. **Singularity cache issues**
   ```
   FATAL: Unable to create build: while extracting...
   ```
   **Solution**: Clear and recreate cache:
   ```bash
   rm -rf ~/.singularity/cache
   export NXF_SINGULARITY_CACHEDIR=/scratch/singularity_cache
   ```

---

### GPU Issues

**Symptom**: Peasoup search fails with GPU errors.

1. **No GPU available**
   ```
   Error: CUDA error: no CUDA-capable device is detected
   ```
   **Solution**:
   - Verify GPU allocation in cluster submission
   - Check SLURM options: `--gres=gpu:1`
   - Ensure `nvidia-smi` works on the node

2. **GPU memory exhausted**
   ```
   Error: CUDA error: out of memory
   ```
   **Solution**: Reduce search parameters:
   ```groovy
   params.peasoup {
       nharmonics = 2          // Reduce from 4
       segments = [1, 2]       // Remove 4
   }
   ```

3. **Wrong CUDA version**
   ```
   Error: CUDA driver version is insufficient
   ```
   **Solution**: Check driver compatibility and load correct module:
   ```bash
   module load cuda/11.0
   ```

---

## Input File Issues

### CSV Parsing Errors

**Symptom**: Pipeline fails during intake stage.

1. **Missing columns**
   ```
   Error: Required column 'cdm' not found
   ```
   **Solution**: Verify all required columns are present:
   ```
   pointing,cluster,beam_name,beam_id,utc_start,ra,dec,fits_files,cdm
   ```

2. **Invalid CSV format**
   ```
   Error: CSV parsing failed at line 3
   ```
   **Solution**: Check for:
   - Trailing commas
   - Missing values
   - Special characters
   - Inconsistent column counts

3. **File not found**
   ```
   Error: File does not exist: /data/beam3.fil
   ```
   **Solution**:
   - Use absolute paths
   - Verify file exists: `ls -la /data/beam3.fil`
   - Check for typos in paths

---

### Filterbank Issues

**Symptom**: readfile or filtool fails.

1. **Corrupt filterbank**
   ```
   Error: Cannot read header from filterbank file
   ```
   **Solution**: Verify file integrity:
   ```bash
   header /path/to/file.fil
   readfile /path/to/file.fil
   ```

2. **Unsupported format**
   ```
   Error: Unknown file format
   ```
   **Solution**: Ensure file is Sigproc filterbank format, not PSRFITS or other.

---

## Processing Issues

### RFI Filtering Fails

**Symptom**: `generateRfiFilter` or `filtool` process fails.

1. **Memory error**
   ```
   MemoryError: Unable to allocate array
   ```
   **Solution**: Increase memory allocation in profile:
   ```groovy
   withName: 'generateRfiFilter' {
       memory = '64 GB'
   }
   ```

2. **All channels masked**
   ```
   Warning: 100% of channels masked by RFI filter
   ```
   **Solution**: Adjust RFI thresholds:
   ```groovy
   params.generateRfiFilter {
       bandpass_sigma = 5.0    // Increase from 3.0
       spectral_kurtosis_sigma = 5.0
   }
   ```

---

### Search Fails

**Symptom**: Peasoup search terminates unexpectedly.

1. **Timeout**
   ```
   Error: Process exceeded time limit
   ```
   **Solution**: Increase time limit in profile:
   ```groovy
   withName: 'peasoup' {
       time = '24h'
   }
   ```

2. **No candidates found**
   ```
   Warning: Zero candidates in XML output
   ```
   **Solution**:
   - Lower S/N threshold: `params.peasoup.min_snr = 6.0`
   - Check DM range covers expected values
   - Verify data quality with diagnostic plots

---

### Folding Fails

**Symptom**: PulsarX folding errors.

1. **Invalid candidate parameters**
   ```
   Error: Period out of range
   ```
   **Solution**: Check candidate filtering:
   ```groovy
   params.parsexml {
       min_period = 0.0001    // Adjust limits
       max_period = 100.0
   }
   ```

2. **Template not found**
   ```
   Error: Cannot find template file
   ```
   **Solution**: Verify template path:
   ```bash
   ls -la templates/Effelsberg_0.template
   ```

---

## Cluster-Specific Issues

### SLURM Issues

1. **Job not starting**
   ```
   squeue shows job in pending state
   ```
   **Solution**: Check resource requests match available resources:
   ```bash
   sinfo -p <partition>
   sacctmgr show qos
   ```

2. **Wrong partition**
   ```
   Error: Invalid partition specified
   ```
   **Solution**: Update profile:
   ```groovy
   process {
       queue = 'normal'  // Change to valid partition
   }
   ```

### HTCondor Issues

1. **Held jobs**
   ```
   condor_q shows jobs on hold
   ```
   **Solution**: Check hold reason:
   ```bash
   condor_q -hold -af HoldReason
   ```

---

## Resume Issues

### Cache Invalidation

**Symptom**: `-resume` re-runs completed processes.

**Common causes**:
1. Input files modified
2. Configuration changed
3. Work directory cleaned

**Solution**: Check what changed:
```bash
nextflow log <run_name> -f hash,name,status
```

### Corrupt Work Directory

**Symptom**: Strange errors after resume.

**Solution**: Clean and restart:
```bash
rm -rf work/
nextflow run elden.nf -entry full -c params.config
```

---

## Memory Management

### Out of Memory Errors

**Symptom**: Process killed with OOM.

**Solution**: Configure process-specific memory:

```groovy
process {
    withName: 'generateRfiFilter' {
        memory = '64 GB'
    }
    withName: 'peasoup' {
        memory = '32 GB'
    }
    withName: 'psrfold' {
        memory = '16 GB'
    }
}
```

---

## Logging and Debugging

### Enable Verbose Logging

```bash
nextflow run elden.nf -entry full -c params.config \
    -with-trace trace.txt \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-dag flowchart.png
```

### Check Process Logs

```bash
# Find failed process work directory
nextflow log <run_name> -f workdir,name,status | grep FAILED

# View process logs
cat work/xx/yyyyyy/.command.log
cat work/xx/yyyyyy/.command.err
```

### Debug Mode

Run a single process manually:

```bash
cd work/xx/yyyyyy/
bash .command.sh
```

---

## Getting Help

If you cannot resolve an issue:

1. **Check existing issues** on GitHub
2. **Gather information**:
   - Nextflow version: `nextflow -version`
   - Full error message
   - Relevant configuration
   - Log files from failed process
3. **Open a GitHub issue** with the above information
