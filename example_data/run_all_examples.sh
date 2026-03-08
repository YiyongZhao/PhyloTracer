#!/usr/bin/env bash
set -u

# Run all 17 PhyloTracer example modules with low CPU pressure and per-module outputs.
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
  local use_output_dir="$2"
  shift 2
  local out="${BASE_OUT}/${tag}"
  mkdir -p "${out}"

  {
    echo "===== ${tag} ====="
    if [[ "${use_output_dir}" == "1" ]]; then
      echo "CMD: python -m phylotracer.Phylo_Tracer $* --output_dir ${out}"
    else
      echo "CMD: python -m phylotracer.Phylo_Tracer $*"
    fi
  } | tee -a "${BASE_OUT}/run.log"

  if [[ "${use_output_dir}" == "1" ]]; then
    (
      cd "${REPO_DIR}" || exit 1
      python -m phylotracer.Phylo_Tracer "$@" --output_dir "${out}" >"${out}/stdout.log" 2>"${out}/stderr.log"
    )
  else
    (
      cd "${REPO_DIR}" || exit 1
      python -m phylotracer.Phylo_Tracer "$@" >"${out}/stdout.log" 2>"${out}/stderr.log"
    )
  fi
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
  find "${out}" -maxdepth 2 -type f | sed "s|${out}/||" | sort | sed -n '1,20p' | tee -a "${BASE_OUT}/run.log"
  echo | tee -a "${BASE_OUT}/run.log"
}

run_cmd m01 1 Phylo_Rooter \
  --input_GF_list "${SCRIPT_DIR}/01_Phylo_Rooter/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/01_Phylo_Rooter/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/01_Phylo_Rooter/sptree.nwk"

# MulRF_Distance uses --output (not --output_dir), so run without extra output_dir injection.
run_cmd m02_mode1 0 MulRF_Distance \
  --mode 1 \
  --input_GF_list "${SCRIPT_DIR}/02_MulRF_Distance/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/02_MulRF_Distance/gene2sps.imap" \
  --output "${BASE_OUT}/m02_mode1/mulrf_mode1.tsv"

run_cmd m02_mode2 0 MulRF_Distance \
  --mode 2 \
  --input_GF_list "${SCRIPT_DIR}/02_MulRF_Distance/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/02_MulRF_Distance/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/02_MulRF_Distance/sptree.nwk" \
  --output "${BASE_OUT}/m02_mode2/mulrf_mode2.tsv"

run_cmd m03 1 PhyloTree_CollapseExpand \
  --input_GF_list "${SCRIPT_DIR}/03_PhyloTree_CollapseExpand/GF_ID2path.imap" \
  --support_value 50

run_cmd m04 1 PhyloSupport_Scaler \
  --input_GF_list "${SCRIPT_DIR}/04_PhyloSupport_Scaler/GF_ID2path.imap" \
  --scale_to 100

run_cmd m05 1 BranchLength_NumericConverter \
  --input_GF_list "${SCRIPT_DIR}/05_BranchLength_NumericConverter/GF_ID2path.imap"

run_cmd m06 1 OrthoFilter_LB \
  --input_GF_list "${SCRIPT_DIR}/06_OrthoFilter_LB/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/06_OrthoFilter_LB/gene2sps.imap" \
  --rrbr_cutoff 5 \
  --srbr_cutoff 2.5

run_cmd m07 1 OrthoFilter_Mono \
  --input_GF_list "${SCRIPT_DIR}/07_OrthoFilter_Mono/GF_ID2path.imap" \
  --input_taxa "${SCRIPT_DIR}/07_OrthoFilter_Mono/Clade.imap" \
  --input_imap "${SCRIPT_DIR}/07_OrthoFilter_Mono/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/07_OrthoFilter_Mono/sptree.nwk"

run_cmd m08 1 TreeTopology_Summarizer \
  --input_GF_list "${SCRIPT_DIR}/08_TreeTopology_Summarizer/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/08_TreeTopology_Summarizer/gene2sps.imap"

run_cmd m09 1 Tree_Visualizer \
  --input_GF_list "${SCRIPT_DIR}/09_Tree_Visualizer/GF_ID2path.visual20.imap" \
  --input_imap "${SCRIPT_DIR}/09_Tree_Visualizer/gene2sps.imap" \
  --gene_categories "${SCRIPT_DIR}/09_Tree_Visualizer/Family.imap" "${SCRIPT_DIR}/09_Tree_Visualizer/Order.imap" "${SCRIPT_DIR}/09_Tree_Visualizer/Clade.imap" \
  --input_sps_tree "${SCRIPT_DIR}/09_Tree_Visualizer/sptree.nwk" \
  --heatmap_matrix "${SCRIPT_DIR}/09_Tree_Visualizer/heatmap_matrix.txt" \
  --keep_branch 1 \
  --tree_style r \
  --visual_gd \
  --gd_support 50 \
  --subclade_support 0 \
  --dup_species_proportion 0.2 \
  --dup_species_num 2 \
  --deepvar 1

run_cmd m10 1 GD_Detector \
  --input_GF_list "${SCRIPT_DIR}/10_GD_Detector/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/10_GD_Detector/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/10_GD_Detector/sptree.nwk" \
  --gd_support 50 \
  --subclade_support 50 \
  --dup_species_proportion 0 \
  --dup_species_num 2 \
  --deepvar 1

run_cmd m11 1 GD_Visualizer \
  --input_sps_tree "${SCRIPT_DIR}/11_GD_Visualizer/numed_sptree.nwk" \
  --gd_result "${SCRIPT_DIR}/11_GD_Visualizer/gd_result_relaxed.txt" \
  --input_imap "${SCRIPT_DIR}/11_GD_Visualizer/gene2sps.imap"

run_cmd m12 1 GD_Loss_Tracker \
  --input_GF_list "${SCRIPT_DIR}/12_GD_Loss_Tracker/GF_ID2path.imap" \
  --input_sps_tree "${SCRIPT_DIR}/12_GD_Loss_Tracker/sptree.nwk" \
  --input_imap "${SCRIPT_DIR}/12_GD_Loss_Tracker/gene2sps.imap"

run_cmd m13 1 GD_Loss_Visualizer \
  --gd_loss_result "${SCRIPT_DIR}/13_GD_Loss_Visualizer/gd_loss_summary.txt" \
  --input_sps_tree "${SCRIPT_DIR}/13_GD_Loss_Visualizer/numed_sptree.nwk"

run_cmd m14 1 Ortho_Retriever \
  --input_GF_list "${SCRIPT_DIR}/14_Ortho_Retriever/GF_ID2path.imap" \
  --input_imap "${SCRIPT_DIR}/14_Ortho_Retriever/gene2sps.imap" \
  --input_gene_length "${SCRIPT_DIR}/14_Ortho_Retriever/gene2length.imap"

run_cmd m15 1 Hybrid_Tracer \
  --input_GF_list "${SCRIPT_DIR}/15_Hybrid_Tracer/gf.txt" \
  --input_Seq_GF_list "${SCRIPT_DIR}/15_Hybrid_Tracer/gf_aln.txt" \
  --input_sps_tree "${SCRIPT_DIR}/15_Hybrid_Tracer/sptree.nwk" \
  --input_imap "${SCRIPT_DIR}/15_Hybrid_Tracer/gene2sps.imap"

run_cmd m16 1 Hybrid_Visualizer \
  --hyde_out "${SCRIPT_DIR}/16_Hybrid_Visualizer/hyde_out.txt" \
  --input_sps_tree "${SCRIPT_DIR}/16_Hybrid_Visualizer/sptree.nwk"

run_cmd m17 1 HaploFinder \
  --mode haplofinder \
  --input_GF_list "${SCRIPT_DIR}/17_Haplofinder/gf.txt" \
  --input_imap "${SCRIPT_DIR}/17_Haplofinder/gene2sps.imap" \
  --input_sps_tree "${SCRIPT_DIR}/17_Haplofinder/sptree.nwk" \
  --species_a arh \
  --species_b ard \
  --species_a_gff "${SCRIPT_DIR}/17_Haplofinder/arh.gff" \
  --species_b_gff "${SCRIPT_DIR}/17_Haplofinder/ard.gff" \
  --species_a_lens "${SCRIPT_DIR}/17_Haplofinder/arh.lens" \
  --species_b_lens "${SCRIPT_DIR}/17_Haplofinder/ard.lens" \
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
