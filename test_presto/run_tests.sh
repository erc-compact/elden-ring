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
#   5. Hybrid Peasoup + PRESTO accelsearch
#
# Usage:
#   ./run_tests.sh           # Run all tests
#   ./run_tests.sh 1         # Run test 1 only
#   ./run_tests.sh 1 3 5     # Run tests 1, 3, and 5
#

set -e

# Change to pipeline directory
cd /hercules/results/fkareem/ELDEN/elden-ring

# Test configuration
TEST_FIL="/hercules/results/rsenzel/signal_inject/J0514-4002A_3HM_000_52.0_01_noise_t4_15min_inj_FAZAL_PRESTO.fil"
BASE_OUTPUT="/hercules/results/fkareem/ELDEN/elden-ring/test_presto"
CONFIG_FILE="test_presto/presto_test.config"
INPUT_CSV="${BASE_OUTPUT}/test_presto_inputfile.csv"

# Parse command line arguments for selective test execution
TESTS_TO_RUN=("$@")
if [ ${#TESTS_TO_RUN[@]} -eq 0 ]; then
    TESTS_TO_RUN=(1 2 3 4 5)
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
        --presto.dm_ranges '[{"dm_low": 28.0, "dm_high": 32.0, "dm_step": 0.2, "downsamp": 1}]' \
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
        --presto.dm_ranges '[{"dm_low": 28.0, "dm_high": 32.0, "dm_step": 0.2, "downsamp": 1}]' \
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
        --presto.dm_ranges '[{"dm_low": 28.0, "dm_high": 32.0, "dm_step": 0.2, "downsamp": 1}]' \
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
        --presto.dm_ranges '[{"dm_low": 28.0, "dm_high": 32.0, "dm_step": 0.2, "downsamp": 1}]' \
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
# Summary
# ==============================================================================
echo ""
echo "=============================================="
echo "Test Suite Complete"
echo "=============================================="
echo ""
echo "Output directories:"
echo "  Test 1 - PRESTO + prepfold:        ${BASE_OUTPUT}/output_presto_prepfold"
echo "  Test 2 - PRESTO + PulsarX:         ${BASE_OUTPUT}/output_presto_pulsarx"
echo "  Test 3 - PRESTO + FFA:             ${BASE_OUTPUT}/output_presto_ffa"
echo "  Test 4 - Standalone Riptide:       ${BASE_OUTPUT}/output_riptide_standalone"
echo "  Test 5 - Hybrid Peasoup + PRESTO:  ${BASE_OUTPUT}/output_hybrid_prepfold"
echo ""
echo "Expected outputs:"
echo "  - *_presto.tar.gz: PRESTO pipeline tarball with PNG + CSV"
echo "  - RIPTIDE_SEARCH/: FFA candidates (tests 3, 4)"
echo ""
