#!/bin/bash
#
# ELDEN-RING Pipeline Test Suite
# ==============================
#
# Test file: J0514-4002A_3HM_000_52.0_01_noise_t4_15min_inj_FAZAL_PRESTO.fil
# Injected pulsar at DM = 30, CDM = 30
#
# Tests cover:
#   1. PRESTO pipeline with prepfold
#   2. PRESTO pipeline with PulsarX fold
#   3. PRESTO pipeline with Riptide FFA search
#   4. Standalone Riptide FFA search
#   5. Hybrid Peasoup + PRESTO accelsearch + prepfold
#   6. Hybrid Peasoup + PRESTO accelsearch + PulsarX fold
#   7. Main Peasoup pipeline + Riptide FFA search
#   8. PRESTO full entry (CSV input)
#   9. PRESTO search+fold from state file (resumption)
#  10. Accelsearch on Peasoup timeseries
#
# Usage:
#   ./run_tests.sh           # Run all tests
#   ./run_tests.sh 1         # Run test 1 only
#   ./run_tests.sh 1 3 5     # Run tests 1, 3, and 5
#   ./run_tests.sh 9         # Run state file test (requires test 1 first)
#   ./run_tests.sh 10        # Run timeseries test (requires test 7 first)
#

set -e

# Change to pipeline directory
cd /hercules/results/fkareem/ELDEN/elden-ring

# Test configuration
TEST_FIL="/hercules/results/rsenzel/signal_inject/J0514-4002A_3HM_000_52.0_01_noise_t4_15min_inj_FAZAL_TEST_2.fil"
BASE_OUTPUT="/hercules/results/fkareem/ELDEN/elden-ring/test_presto"
CONFIG_FILE="test_presto/presto_test.config"
INPUT_CSV="${BASE_OUTPUT}/test_presto_inputfile.csv"

# Parse command line arguments for selective test execution
TESTS_TO_RUN=("$@")
if [ ${#TESTS_TO_RUN[@]} -eq 0 ]; then
    TESTS_TO_RUN=(1 2 3 4 5 6 7 8 9 10)
fi

should_run_test() {
    local test_num=$1
    for t in "${TESTS_TO_RUN[@]}"; do
        if [ "$t" -eq "$test_num" ]; then
            return 0
        fi
    done
    return 1
}

# ==============================================================================
# Test 1: PRESTO Pipeline with prepfold
# ==============================================================================
if should_run_test 1; then
    echo "=============================================="
    echo "Test 1: PRESTO search + PRESTO prepfold"
    echo "=============================================="

    nextflow run elden.nf \
        -entry presto_pipeline \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --input_fil ${TEST_FIL} \
        --basedir ${BASE_OUTPUT}/output_presto_prepfold \
        --target_name J0514_presto_prepfold \
        --tarball_prefix presto_prepfold_test \
        --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
        --presto.fold_backend presto \
        --presto.max_fold_cands 50 \
        --presto.zmax 100 \
        --presto.wmax 0 \
        -resume
fi

# ==============================================================================
# Test 2: PRESTO Pipeline with PulsarX fold
# ==============================================================================
if should_run_test 2; then
    echo "=============================================="
    echo "Test 2: PRESTO search + PulsarX fold"
    echo "=============================================="

    nextflow run elden.nf \
        -entry presto_pipeline \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --input_fil ${TEST_FIL} \
        --basedir ${BASE_OUTPUT}/output_presto_pulsarx \
        --target_name J0514_presto_pulsarx \
        --tarball_prefix presto_pulsarx_test \
        --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
        --presto.fold_backend pulsarx \
        --presto.max_fold_cands 50 \
        --presto.zmax 100 \
        --presto.wmax 0 \
        --presto.fold_threads 8 \
        --presto.fold_nbins 64 \
        --presto.fold_subint_length 10.0 \
        -resume
fi

# ==============================================================================
# Test 3: PRESTO Pipeline with Riptide FFA Search
# ==============================================================================
if should_run_test 3; then
    echo "=============================================="
    echo "Test 3: PRESTO search + Riptide FFA"
    echo "=============================================="

    nextflow run elden.nf \
        -entry presto_pipeline \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --input_fil ${TEST_FIL} \
        --basedir ${BASE_OUTPUT}/output_presto_ffa \
        --target_name J0514_presto_ffa \
        --tarball_prefix presto_ffa_test \
        --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
        --presto.fold_backend pulsarx \
        --presto.max_fold_cands 50 \
        --presto.zmax 100 \
        --riptide.run_ffa_search true \
        -resume
fi

# ==============================================================================
# Test 4: Standalone Riptide FFA Search
# ==============================================================================
if should_run_test 4; then
    echo "=============================================="
    echo "Test 4: Standalone Riptide FFA (PRESTO backend)"
    echo "=============================================="

    nextflow run elden.nf \
        -entry run_riptide \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --input_fil ${TEST_FIL} \
        --basedir ${BASE_OUTPUT}/output_riptide_standalone \
        --target_name J0514_riptide \
        --riptide.backend presto \
        --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
        -resume
fi

# ==============================================================================
# Test 5: Hybrid Peasoup + PRESTO Accelsearch + prepfold
# ==============================================================================
if should_run_test 5; then
    echo "=============================================="
    echo "Test 5: Peasoup + PRESTO accelsearch + prepfold"
    echo "=============================================="

    nextflow run elden.nf \
        -entry peasoup_with_presto_search \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --files_list ${INPUT_CSV} \
        --basedir ${BASE_OUTPUT}/output_hybrid_prepfold \
        --target_name J0514_hybrid_prepfold \
        --tarball_prefix hybrid_prepfold_test \
        --search_backend peasoup \
        --peasoup.dump_timeseries true \
        --peasoup.min_snr 6.0 \
        --peasoup.acc_start -50 \
        --peasoup.acc_end 50 \
        --presto.fold_backend presto \
        --presto.max_fold_cands 50 \
        --presto.zmax 100 \
        --presto.wmax 0 \
        --ddplan.dm_start -2 \
        --ddplan.dm_end 2 \
        --ddplan.dm_step 0.2 \
        -resume
fi

# ==============================================================================
# Test 6: Hybrid Peasoup + PRESTO Accelsearch + PulsarX fold
# ==============================================================================
if should_run_test 6; then
    echo "=============================================="
    echo "Test 6: Peasoup + PRESTO accelsearch + PulsarX fold"
    echo "=============================================="

    nextflow run elden.nf \
        -entry peasoup_with_presto_search \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --files_list ${INPUT_CSV} \
        --basedir ${BASE_OUTPUT}/output_hybrid_pulsarx \
        --target_name J0514_hybrid_pulsarx \
        --tarball_prefix hybrid_pulsarx_test \
        --search_backend peasoup \
        --peasoup.dump_timeseries true \
        --peasoup.min_snr 6.0 \
        --peasoup.acc_start -50 \
        --peasoup.acc_end 50 \
        --presto.fold_backend pulsarx \
        --presto.max_fold_cands 50 \
        --presto.fold_threads 8 \
        --presto.zmax 100 \
        --presto.wmax 0 \
        --ddplan.dm_start -2 \
        --ddplan.dm_end 2 \
        --ddplan.dm_step 0.2 \
        -resume
fi

# ==============================================================================
# Test 7: Main Peasoup Pipeline + Riptide FFA Search
# ==============================================================================
if should_run_test 7; then
    echo "=============================================="
    echo "Test 7: Peasoup full pipeline + Riptide FFA"
    echo "=============================================="

    nextflow run elden.nf \
        -entry full \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --files_list ${INPUT_CSV} \
        --basedir ${BASE_OUTPUT}/output_peasoup_ffa \
        --target_name J0514_peasoup_ffa \
        --tarball_prefix peasoup_ffa_test \
        --peasoup.dump_timeseries true \
        --peasoup.min_snr 6.0 \
        --peasoup.acc_start -50 \
        --peasoup.acc_end 50 \
        --riptide.run_ffa_search true \
        --ddplan.dm_start -2 \
        --ddplan.dm_end 2 \
        --ddplan.dm_step 0.2 \
        -resume
fi

# ==============================================================================
# Test 8: PRESTO Full Entry (CSV input)
# ==============================================================================
if should_run_test 8; then
    echo "=============================================="
    echo "Test 8: PRESTO full pipeline via CSV input"
    echo "=============================================="

    nextflow run elden.nf \
        -entry presto_full_entry \
        -profile hercules \
        -c ${CONFIG_FILE} \
        --files_list ${INPUT_CSV} \
        --basedir ${BASE_OUTPUT}/output_presto_csv \
        --target_name J0514_presto_csv \
        --tarball_prefix presto_csv_test \
        --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
        --presto.fold_backend pulsarx \
        --presto.max_fold_cands 50 \
        --presto.zmax 100 \
        -resume
fi

# ==============================================================================
# Test 9: PRESTO Search and Fold from State File (state file resumption)
# ==============================================================================
if should_run_test 9; then
    echo "=============================================="
    echo "Test 9: PRESTO search+fold from state file"
    echo "=============================================="

    # This test requires running Test 1 first to generate the state file
    STATE_FILE="${BASE_OUTPUT}/output_presto_prepfold/PRESTO_STATE/dedisperse_state.json"

    if [ ! -f "${STATE_FILE}" ]; then
        echo "WARNING: State file not found at ${STATE_FILE}"
        echo "Run Test 1 first to generate the dedisperse state file"
        echo "Skipping Test 9..."
    else
        nextflow run elden.nf \
            -entry presto_search_fold \
            -profile hercules \
            -c ${CONFIG_FILE} \
            --state_file ${STATE_FILE} \
            --basedir ${BASE_OUTPUT}/output_presto_statefile \
            --target_name J0514_statefile \
            --tarball_prefix statefile_test \
            --presto.fold_backend pulsarx \
            --presto.max_fold_cands 50 \
            --presto.zmax 100 \
            -resume
    fi
fi

# ==============================================================================
# Test 10: Accelsearch on Peasoup Timeseries
# ==============================================================================
if should_run_test 10; then
    echo "=============================================="
    echo "Test 10: PRESTO accelsearch on Peasoup timeseries"
    echo "=============================================="

    # This test requires running Test 7 first to generate time series files
    TIMESERIES_DIR="${BASE_OUTPUT}/output_peasoup_ffa"

    if [ ! -d "${TIMESERIES_DIR}" ]; then
        echo "WARNING: Time series directory not found at ${TIMESERIES_DIR}"
        echo "Run Test 7 first (with peasoup.dump_timeseries true) to generate time series"
        echo "Skipping Test 10..."
    else
        nextflow run elden.nf \
            -entry run_accelsearch_on_timeseries \
            -profile hercules \
            -c ${CONFIG_FILE} \
            --files_list ${INPUT_CSV} \
            --basedir ${TIMESERIES_DIR} \
            --runID J0514_peasoup_ffa \
            --target_name J0514_accelsearch \
            --tarball_prefix accelsearch_test \
            --presto.fold_backend pulsarx \
            --presto.max_fold_cands 50 \
            --presto.zmax 100 \
            -resume
    fi
fi

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "=============================================="
echo "Test Suite Complete"
echo "=============================================="
echo ""
echo "Output directories:"
echo "  Test 1  - PRESTO + prepfold:            ${BASE_OUTPUT}/output_presto_prepfold"
echo "  Test 2  - PRESTO + PulsarX:             ${BASE_OUTPUT}/output_presto_pulsarx"
echo "  Test 3  - PRESTO + FFA:                 ${BASE_OUTPUT}/output_presto_ffa"
echo "  Test 4  - Standalone Riptide:           ${BASE_OUTPUT}/output_riptide_standalone"
echo "  Test 5  - Hybrid Peasoup + prepfold:    ${BASE_OUTPUT}/output_hybrid_prepfold"
echo "  Test 6  - Hybrid Peasoup + PulsarX:     ${BASE_OUTPUT}/output_hybrid_pulsarx"
echo "  Test 7  - Peasoup full + FFA:           ${BASE_OUTPUT}/output_peasoup_ffa"
echo "  Test 8  - PRESTO CSV input:             ${BASE_OUTPUT}/output_presto_csv"
echo "  Test 9  - State file resumption:        ${BASE_OUTPUT}/output_presto_statefile"
echo "  Test 10 - Accelsearch on timeseries:    ${BASE_OUTPUT}/output_peasoup_ffa (reuses)"
echo ""
echo "Expected outputs:"
echo "  - *_presto.tar.gz: PRESTO pipeline tarball with PNG + CSV"
echo "  - *_peasoup.tar.gz: Peasoup pipeline tarball"
echo "  - RIPTIDE_SEARCH/: FFA candidates (tests 3, 4, 7)"
echo "  - PRESTO_STATE/*.json: State files for resumption (tests 1-3, 8)"
echo ""
echo "Test dependencies:"
echo "  - Test 9 requires Test 1 to generate state file"
echo "  - Test 10 requires Test 7 to generate timeseries files"
echo ""
