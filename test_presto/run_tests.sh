#!/bin/bash
#
# PRESTO Pipeline Test Commands
# =============================
#
# Test file: /hercules/results/rsenzel/signal_inject/J0514-4002A_3HM_000_52.0_01_noise_t4_15min_inj_FAZAL_PRESTO.fil
# Injected pulsar at DM = 30
#
# Pipeline combinations tested:
#   1. Pure PRESTO search + PRESTO prepfold
#   2. Pure PRESTO search + PulsarX fold
#   3. Peasoup search + PRESTO accelsearch (hybrid) + prepfold
#   4. Peasoup search + PRESTO accelsearch (hybrid) + PulsarX fold
#

# Change to pipeline directory
cd /hercules/results/fkareem/ELDEN/elden-ring

# Common test file
TEST_FIL="/hercules/results/rsenzel/signal_inject/J0514-4002A_3HM_000_52.0_01_noise_t4_15min_inj_FAZAL_PRESTO.fil"
BASE_OUTPUT="/hercules/results/fkareem/ELDEN/elden-ring/test_presto"
WMAX_JERK=50

# ==============================================================================
# SECTION A: PURE PRESTO PIPELINE TESTS
# ==============================================================================

# ==============================================================================
# # Test 1: Full PRESTO Pipeline with prepfold (PRESTO search + PRESTO fold)
# # ==============================================================================
# echo "=============================================="
# echo "Test 1: PRESTO search + PRESTO prepfold (wmax=0)"
# echo "=============================================="
#
# nextflow run elden.nf \
#   -entry presto_pipeline \
#   -profile hercules \
#   -c test_presto/presto_test.config \
#   --input_fil $TEST_FIL \
#   --basedir ${BASE_OUTPUT}/output_presto_prepfold \
#   --target_name J0514_presto_prepfold \
#   --tarball_prefix presto_prepfold_test \
#   --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
#   --presto.fold_backend presto \
#   --presto.max_fold_cands 50 \
#   --presto.zmax 100 \
#   --presto.wmax 0 \
#   -resume
#
# # ==============================================================================
# # Test 2: PRESTO search + PRESTO prepfold (jerk enabled)
# # ==============================================================================
# echo "=============================================="
# echo "Test 2: PRESTO search + PRESTO prepfold (wmax=${WMAX_JERK})"
# echo "=============================================="
#
# nextflow run elden.nf \
#   -entry presto_pipeline \
#   -profile hercules \
#   -c test_presto/presto_test.config \
#   --input_fil $TEST_FIL \
#   --basedir ${BASE_OUTPUT}/output_presto_prepfold_jerk \
#   --target_name J0514_presto_prepfold_jerk \
#   --tarball_prefix presto_prepfold_jerk_test \
#   --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.2, "downsamp": 1}]' \
#   --presto.fold_backend presto \
#   --presto.max_fold_cands 50 \
#   --presto.zmax 100 \
#   --presto.wmax ${WMAX_JERK} \
#   -resume
#
#
# # ==============================================================================
# echo "=============================================="
# echo "Test 3: PRESTO search + PulsarX fold"
# echo "=============================================="
#
# nextflow run elden.nf \
#     -entry presto_pipeline \
#     -profile hercules \
#     -c test_presto/presto_test.config \
#     --input_fil $TEST_FIL \
#     --basedir ${BASE_OUTPUT}/output_presto_pulsarx \
#     --target_name J0514_presto_pulsarx \
#     --tarball_prefix presto_pulsarx_test \
#     --presto.dm_ranges '[{"dm_low": 29.0, "dm_high": 31.0, "dm_step": 0.3, "downsamp": 1}]' \
#     --presto.fold_backend pulsarx \
#     --presto.max_fold_cands 50 \
#     --presto.zmax 100 \
#     --presto.wmax 0 \
#     --presto.fold_threads 8 \
#     --presto.fold_nbins 64 \
#     --presto.fold_subint_length 10.0 \
#     -resume
#
# ==============================================================================
# SECTION B: HYBRID PEASOUP + PRESTO PIPELINE TESTS
# ==============================================================================

# ==============================================================================
# Test 4: Peasoup search + PRESTO accelsearch + prepfold (Hybrid pipeline)
# ==============================================================================
# This runs peasoup for initial search, dumps time series, then runs PRESTO
# accelsearch on the dumped time series, and folds with prepfold
echo "=============================================="
echo "Test 4: Peasoup + PRESTO accelsearch + prepfold (Hybrid)"
echo "=============================================="

nextflow run elden.nf \
  -entry peasoup_with_presto_search \
  -profile hercules \
  -c test_presto/presto_test.config \
  --files_list ${BASE_OUTPUT}/test_presto_inputfile.csv \
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
  --ddplan.dm_start -1 \
  --ddplan.dm_end 1 \
  --ddplan.dm_step 0.2 \
  -resume

# ==============================================================================
# Test 5: Peasoup search + PRESTO accelsearch + PulsarX fold (Hybrid pipeline)
# ==============================================================================
echo "=============================================="
echo "Test 5: Peasoup + PRESTO accelsearch + PulsarX fold (Hybrid)"
echo "=============================================="

nextflow run elden.nf \
  -entry peasoup_with_presto_search \
  -profile hercules \
  -c test_presto/presto_test.config \
  --files_list ${BASE_OUTPUT}/test_presto_inputfile.csv \
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
  -resume

# # ==============================================================================
# # SECTION C: STAGED EXECUTION TESTS
# # ==============================================================================
#
# # ==============================================================================
# # Test 5: Run PRESTO stages 1-3 only (RFI -> Birdies -> Dedispersion)
# # ==============================================================================
# echo "=============================================="
# echo "Test 5: PRESTO stages 1-3 only (stop after dedispersion)"
# echo "=============================================="
#
# nextflow run elden.nf \
#     -entry presto_pipeline \
#     -profile hercules \
#     --input_fil $TEST_FIL \
#     --basedir ${BASE_OUTPUT}/output_stages_1_3 \
#     --target_name J0514_stages_1_3 \
#     --presto.start_stage 1 \
#     --presto.end_stage 3 \
#     --presto.dm_ranges '[{"dm_low": 25.0, "dm_high": 35.0, "dm_step": 0.3, "downsamp": 1}]' \
#     -resume
#
# # ==============================================================================
# # Test 6: Resume from stage 4 using state file
# # ==============================================================================
# # (Run this after Test 5 completes - uses the state file from stage 3)
# echo "=============================================="
# echo "Test 6: Resume from stage 4 using state file"
# echo "=============================================="
#
# # NOTE: Uncomment and update the state_file path after Test 5 completes
# # STATE_FILE="${BASE_OUTPUT}/output_stages_1_3/J0514_stages_1_3/presto/state/dedisperse_state.json"
# #
# # nextflow run elden.nf \
# #     -entry presto_pipeline \
# #     -profile hercules \
# #     --state_file $STATE_FILE \
# #     --basedir ${BASE_OUTPUT}/output_resumed \
# #     --target_name J0514_resumed \
# #     --presto.start_stage 4 \
# #     --presto.end_stage 6 \
# #     --presto.fold_backend presto \
# #     --tarball_prefix resumed_test \
# #     -resume
#
# # ==============================================================================
# # SECTION D: QUICK VALIDATION TESTS
# # ==============================================================================
#
# # ==============================================================================
# # Test 7: Quick PRESTO + prepfold (narrow DM range)
# # ==============================================================================
# echo "=============================================="
# echo "Test 7: Quick PRESTO + prepfold (narrow DM)"
# echo "=============================================="
#
# nextflow run elden.nf \
#     -entry presto_pipeline \
#     -profile hercules \
#     --input_fil $TEST_FIL \
#     --basedir ${BASE_OUTPUT}/output_quick_prepfold \
#     --target_name J0514_quick_prepfold \
#     --tarball_prefix quick_prepfold \
#     --presto.dm_ranges '[{"dm_low": 28.0, "dm_high": 32.0, "dm_step": 0.2, "downsamp": 1}]' \
#     --presto.fold_backend presto \
#     --presto.max_fold_cands 20 \
#     --presto.zmax 50 \
#     -resume
#
# # ==============================================================================
# # Test 8: Quick PRESTO + PulsarX fold (narrow DM range)
# # ==============================================================================
# echo "=============================================="
# echo "Test 8: Quick PRESTO + PulsarX fold (narrow DM)"
# echo "=============================================="
#
# nextflow run elden.nf \
#     -entry presto_pipeline \
#     -profile hercules \
#     --input_fil $TEST_FIL \
#     --basedir ${BASE_OUTPUT}/output_quick_pulsarx \
#     --target_name J0514_quick_pulsarx \
#     --tarball_prefix quick_pulsarx \
#     --presto.dm_ranges '[{"dm_low": 28.0, "dm_high": 32.0, "dm_step": 0.2, "downsamp": 1}]' \
#     --presto.fold_backend pulsarx \
#     --presto.max_fold_cands 20 \
#     --presto.zmax 50 \
#     -resume
#
# # ==============================================================================
# # SECTION E: CSV INPUT TESTS (Multiple files)
# # ==============================================================================
#
# # ==============================================================================
# # Test 9: Full PRESTO pipeline with CSV input + prepfold
# # ==============================================================================
# echo "=============================================="
# echo "Test 9: PRESTO with CSV input + prepfold"
# echo "=============================================="
#
# nextflow run elden.nf \
#     -entry presto_full_entry \
#     -profile hercules \
#     -c test_presto/presto_test.config \
#     --files_list ${BASE_OUTPUT}/test_presto_inputfile.csv \
#     --basedir ${BASE_OUTPUT}/output_csv_prepfold \
#     --target_name J0514_csv_prepfold \
#     --tarball_prefix csv_prepfold \
#     --presto.fold_backend presto \
#     -resume
#
# # ==============================================================================
# # Test 10: Full PRESTO pipeline with CSV input + PulsarX fold
# # ==============================================================================
# echo "=============================================="
# echo "Test 10: PRESTO with CSV input + PulsarX fold"
# echo "=============================================="
#
# nextflow run elden.nf \
#     -entry presto_full_entry \
#     -profile hercules \
#     -c test_presto/presto_test.config \
#     --files_list ${BASE_OUTPUT}/test_presto_inputfile.csv \
#     --basedir ${BASE_OUTPUT}/output_csv_pulsarx \
#     --target_name J0514_csv_pulsarx \
#     --tarball_prefix csv_pulsarx \
#     --presto.fold_backend pulsarx \
#     -resume
#
# echo "=============================================="
# echo "All tests complete!"
# echo "=============================================="
# echo ""
# echo "Output directories:"
# echo "  - Pure PRESTO + prepfold:    ${BASE_OUTPUT}/output_presto_prepfold"
# echo "  - Pure PRESTO + PulsarX:     ${BASE_OUTPUT}/output_presto_pulsarx"
# echo "  - Hybrid + prepfold:         ${BASE_OUTPUT}/output_hybrid_prepfold"
# echo "  - Hybrid + PulsarX:          ${BASE_OUTPUT}/output_hybrid_pulsarx"
# echo "  - Staged (1-3):              ${BASE_OUTPUT}/output_stages_1_3"
# echo "  - Quick prepfold:            ${BASE_OUTPUT}/output_quick_prepfold"
# echo "  - Quick PulsarX:             ${BASE_OUTPUT}/output_quick_pulsarx"
# echo "  - CSV + prepfold:            ${BASE_OUTPUT}/output_csv_prepfold"
# echo "  - CSV + PulsarX:             ${BASE_OUTPUT}/output_csv_pulsarx"
# echo ""
# echo "Expected tarballs in each output directory:"
# echo "  - *_presto.tar.gz (PRESTO pipeline tarball with PNG + CSV)"
# echo "  - *_peasoup.tar.gz (Peasoup pipeline tarball, for hybrid tests)"
