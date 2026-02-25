#!/usr/bin/env bash
set -u

# Run all 16 PhyloTracer example modules with low CPU pressure and per-module outputs.
# Usage:
#   cd /Users/apple/Documents/GitHub/PhyloTracer/example_data
#   bash run_all_examples.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TS="$(date +%Y%m%d_%H%M%S)"
BASE_OUT="${SCRIPT_DIR}/_example_runs/${TS}"
mkdir -p "${BASE_OUT}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR="${BASE_OUT}/mplcache"
mkdir -p "${MPLCONFIGDIR}"

PASS_COUNT=0
FAIL_COUNT=0
FAIL_LIST=()

run_cmd() {
  local tag="$1"
  shift
  local out="${BASE_OUT}/${tag}"
  mkdir -p "${out}"

  {
    echo "===== ${tag} ====="
    echo "CMD: python -m phylotracer.Phylo_Tracer $* --output_dir ${out}"
  } | tee -a "${BASE_OUT}/run.log"

  (
    cd "${REPO_DIR}" || exit 1
    python -m phylotracer.Phylo_Tracer "$@" --output_dir "${out}" >"${out}/stdout.log" 2>"${out}/stderr.log"
  )
  local code=$?

  if [[ ${code} -eq 0 ]]; then
    PASS_COUNT=$((PASS_COUNT + 1))
    echo "EXIT: ${code} (PASS)" | tee -a "${BASE_OUT}/run.log"
  else
    FAIL_COUNT=$((FAIL_COUNT + 1))
    FAIL_LIST+=("${tag}")
    echo "EXIT: ${code} (FAIL)" | tee -a "${BASE_OUT}/run.log"
  fi

  echo "TOP FILES:" | tee -a "${BASE_OUT}/run.log"
  find "${out}" -maxdepth 2 -type f | sed "s|${out}/||" | sort | sed -n '1,15p' | tee -a "${BASE_OUT}/run.log"
  echo | tee -a "${BASE_OUT}/run.log"
}

run_cmd m01 Phylo_Rooter \
  --input_GF_list "${SCRIPT_DIR}/1_Phylo_Rooter/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/1_Phylo_Rooter/gene2sps.imap" \
  --input_gene_length "${SCRIPT_DIR}/1_Phylo_Rooter/gene2length.imap" \
  --input_sps_tree "${SCRIPT_DIR}/1_Phylo_Rooter/sptree.nwk"

run_cmd m02 PhyloTree_CollapseExpand \
  --input_GF_list "${SCRIPT_DIR}/2_PhyloTree_CollapseExpand/GF_ID2path.imap" \
  --support_value 50

run_cmd m03 PhyloSupport_Scaler \
  --input_GF_list "${SCRIPT_DIR}/3_PhyloSupport_Scaler/GF_ID2path.imap" \
  --scale_to 100

run_cmd m04 BranchLength_NumericConverter \
  --input_GF_list "${SCRIPT_DIR}/4_BranchLength_NumericConverter/GF_ID2path.imap"

run_cmd m05 OrthoFilter_LB \
  --input_GF_list "${SCRIPT_DIR}/5_OrthoFilter_LB/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/5_OrthoFilter_LB/gene2sps.imap" \
  --absolute_branch_length 5 \
  --relative_branch_length 2.5

run_cmd m06 OrthoFilter_Mono \
  --input_GF_list "${SCRIPT_DIR}/6_OrthoFilter_Mono/GF_ID2path.imap" \
  --input_taxa "${SCRIPT_DIR}/6_OrthoFilter_Mono/gene2clade.imap" \
  --input_imap "${SCRIPT_DIR}/6_OrthoFilter_Mono/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/6_OrthoFilter_Mono/sptree.nwk"

run_cmd m07 TreeTopology_Summarizer \
  --input_GF_list "${SCRIPT_DIR}/7_TreeTopology_Summarizer/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/7_TreeTopology_Summarizer/gene2sps.imap"

run_cmd m08 Tree_Visualizer \
  --input_GF_list "${SCRIPT_DIR}/8_Tree_Visualizer/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/8_Tree_Visualizer/gene2sps.imap"

run_cmd m09 GD_Detector \
  --input_GF_list "${SCRIPT_DIR}/9_GD_Detector/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/9_GD_Detector/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/9_GD_Detector/sptree.nwk" \
  --gd_support 50 \
  --subclade_support 50 \
  --dup_species_proportion 0 \
  --dup_species_num 2 \
  --deepvar 1

run_cmd m10 GD_Visualizer \
  --input_sps_tree "${SCRIPT_DIR}/10_GD_Visualizer/numed_sptree.nwk" \
  --gd_result "${SCRIPT_DIR}/10_GD_Visualizer/gd_result_relaxed.txt" \
  --input_imap "${SCRIPT_DIR}/10_GD_Visualizer/gene2sps.imap"

run_cmd m11 GD_Loss_Tracker \
  --input_GF_list "${SCRIPT_DIR}/11_GD_Loss_Tracker/GF_ID2path.imap" \
  --input_sps_tree "${SCRIPT_DIR}/11_GD_Loss_Tracker/sptree.nwk" \
  --input_imap "${SCRIPT_DIR}/11_GD_Loss_Tracker/gene2sps.imap"

run_cmd m12 GD_Loss_Visualizer \
  --gd_loss_result "${SCRIPT_DIR}/12_GD_Loss_Visualizer/gd_loss_summary.txt" \
  --input_sps_tree "${SCRIPT_DIR}/12_GD_Loss_Visualizer/numed_sptree.nwk"

run_cmd m13 Ortho_Retriever \
  --input_GF_list "${SCRIPT_DIR}/13_Ortho_Retriever/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/13_Ortho_Retriever/gene2sps.imap" \
  --input_gene_length "${SCRIPT_DIR}/13_Ortho_Retriever/gene2length.imap"

run_cmd m14 Hybrid_Tracer \
  --input_GF_list "${SCRIPT_DIR}/14_Hybrid_Tracer/gf.txt" \
  --input_Seq_GF_list "${SCRIPT_DIR}/14_Hybrid_Tracer/gf_aln.txt" \
  --input_sps_tree "${SCRIPT_DIR}/14_Hybrid_Tracer/sptree.nwk" \
  --input_imap "${SCRIPT_DIR}/14_Hybrid_Tracer/gene2sps.imap"

run_cmd m15 Hybrid_Visualizer \
  --hyde_out "${SCRIPT_DIR}/15_Hybrid_Visualizer/hyde_out.txt" \
  --input_sps_tree "${SCRIPT_DIR}/15_Hybrid_Visualizer/sptree.nwk"

# Important: species labels in this example are lowercase (arh/ard).
run_cmd m16 HaploFinder \
  --mode haplofinder \
  --input_GF_list "${SCRIPT_DIR}/16_Haplofinder/gf.txt" \
  --input_imap "${SCRIPT_DIR}/16_Haplofinder/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/16_Haplofinder/sptree.nwk" \
  --species_a arh \
  --species_b ard \
  --species_a_gff "${SCRIPT_DIR}/16_Haplofinder/arh.gff" \
  --species_b_gff "${SCRIPT_DIR}/16_Haplofinder/ard.gff" \
  --species_a_lens "${SCRIPT_DIR}/16_Haplofinder/arh.lens" \
  --species_b_lens "${SCRIPT_DIR}/16_Haplofinder/ard.lens" \
  --gd_support 50 \
  --pair_support 50

echo "===== SUMMARY =====" | tee -a "${BASE_OUT}/run.log"
echo "PASS: ${PASS_COUNT}" | tee -a "${BASE_OUT}/run.log"
echo "FAIL: ${FAIL_COUNT}" | tee -a "${BASE_OUT}/run.log"
if [[ ${FAIL_COUNT} -gt 0 ]]; then
  echo "FAILED MODULES: ${FAIL_LIST[*]}" | tee -a "${BASE_OUT}/run.log"
fi
echo "OUTPUT_DIR: ${BASE_OUT}" | tee -a "${BASE_OUT}/run.log"

if [[ ${FAIL_COUNT} -gt 0 ]]; then
  exit 1
fi

