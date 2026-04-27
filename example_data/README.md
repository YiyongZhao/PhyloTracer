# Example Data

This directory contains runnable example inputs for the PhyloTracer modules.
Each numbered subdirectory corresponds to one command-line workflow and includes
the minimal files needed to reproduce the example.

## Quick Start

Run the full example suite from the repository root:

```bash
bash example_data/run_all_examples.sh
```

The script writes timestamped outputs under `example_data/_example_runs/`. Those
runtime outputs are local artifacts and are excluded from version control.

## Running Individual Examples

Use paths relative to the repository root so the commands work on any machine.
`12_GD_Loss_Tracker/gd_loss_tracker/` and `13_GD_Loss_Visualizer/` separate
`parsimony/` and `accumulate/` outputs so the two node-counting modes can be
compared directly.

### 01. `Phylo_Rooter`

```bash
python -m phylotracer.Phylo_Tracer Phylo_Rooter \
  --input_GF_list example_data/01_Phylo_Rooter/GF_ID2path.imap \
  --input_imap example_data/01_Phylo_Rooter/gene2sps.imap \
  --input_sps_tree example_data/01_Phylo_Rooter/sptree.nwk
```

### 02. `MulRF_Distance`

```bash
python -m phylotracer.Phylo_Tracer MulRF_Distance \
  --mode 1 \
  --input_GF_list example_data/02_MulRF_Distance/GF_ID2path.imap \
  --input_imap example_data/02_MulRF_Distance/gene2sps.imap \
  --output example_data/02_MulRF_Distance/mulrf_mode1.tsv
```

```bash
python -m phylotracer.Phylo_Tracer MulRF_Distance \
  --mode 2 \
  --input_GF_list example_data/02_MulRF_Distance/GF_ID2path.imap \
  --input_imap example_data/02_MulRF_Distance/gene2sps.imap \
  --input_sps_tree example_data/02_MulRF_Distance/sptree.nwk \
  --output example_data/02_MulRF_Distance/mulrf_mode2.tsv
```

### 03-17. Remaining Modules

For the remaining examples, use `example_data/run_all_examples.sh` as the
canonical reference. It provides portable commands for:

- `PhyloTree_CollapseExpand`
- `PhyloSupport_Scaler`
- `BranchLength_NumericConverter`
- `OrthoFilter_LB`
- `OrthoFilter_Mono`
- `TreeTopology_Summarizer`
- `Tree_Visualizer`
- `GD_Detector`
- `GD_Visualizer`
- `GD_Loss_Tracker`
- `GD_Loss_Visualizer`
- `Ortho_Retriever`
  Supports default ortholog splitting, optional synteny refinement, and optional strict outgroup attachment via `--add_outgroup`. When enabled, the run also produces `ortholog_outgroup_report.tsv`.
- `Hybrid_Tracer`
- `Hybrid_Visualizer`
- `HaploFinder`

## Notes

- Some example directories include generated outputs for demonstration.
- Large temporary outputs should be written to `example_data/_example_runs/`.
- If you want a clean rerun, delete or archive previous output directories first.
