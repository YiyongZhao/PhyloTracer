# PhyloTracer

PhyloTracer is a modular phylogenomics toolkit for gene-tree rooting, orthology refinement, gene-duplication/loss analysis, topology summarization, and hybridization signal screening.

It is designed for users who want practical command-line workflows for large gene-family datasets, while keeping each module independent and reusable.

## What PhyloTracer Does

PhyloTracer provides end-to-end utilities for common comparative-genomics tasks:

- Root gene trees with species-tree guidance (`Phylo_Rooter`)
- Clean noisy trees by long-branch or monophyly filters (`OrthoFilter_LB`, `OrthoFilter_Mono`)
- Detect and classify gene-duplication events (`GD_Detector`)
- Track and visualize duplication-associated loss patterns (`GD_Loss_Tracker`, `GD_Loss_Visualizer`)
- Summarize topology frequencies (`TreeTopology_Summarizer`)
- Build HyDe-ready inputs from GD-aware gene sets (`Hybrid_Tracer`)
- Visualize trees, GD patterns, and hybridization outputs (`Tree_Visualizer`, `GD_Visualizer`, `Hybrid_Visualizer`)
- Retrieve putative ortholog sets (`Ortho_Retriever`)

## Key Features

- **Modular command design**: each tool can run independently.
- **Reproducible CLI workflows**: clear inputs/outputs for pipeline integration.
- **Tree-centric analysis**: consistent support for Newick gene/species trees.
- **GD-focused analytics**: duplication type statistics plus downstream loss profiling.
- **Practical visualization outputs**: publication-ready tree and summary plots.

## Installation

### Option 1: Conda (recommended)

```bash
git clone https://github.com/YiyongZhao/PhyloTracer.git
cd PhyloTracer
conda env create -f environment.yml
conda activate phylotracer
```

### Option 2: PyPI

```bash
pip install PhyloTracer
```

### Optional Qt fix for headless environments

```bash
export QT_QPA_PLATFORM=offscreen
```

## Quick Start

### 1. Show global help

```bash
phylotracer -h
```

### 2. Detect gene duplications

```bash
phylotracer GD_Detector \
  --input_GF_list example_data/GD_Detector/GF_ID2path.imap \
  --input_imap example_data/GD_Detector/gene2sps.imap \
  --input_sps_tree example_data/GD_Detector/sptree.nwk \
  --gd_support 50 \
  --subclade_support 50 \
  --dup_species_proportion 0 \
  --dup_species_num 2 \
  --deepvar 1
```

### 3. Root gene trees

```bash
phylotracer Phylo_Rooter \
  --input_GF_list example_data/Phylo_Rooter/GF_ID2path.imap \
  --input_imap example_data/Phylo_Rooter/gene2sps.imap \
  --input_gene_length example_data/Phylo_Rooter/gene2length.imap \
  --input_sps_tree example_data/Phylo_Rooter/sptree.nwk
```

## Common Workflow

A typical analysis can be run in this order:

1. **Preprocess trees**
   - `PhyloSupport_Scaler`, `BranchLength_NumericConverter`, `PhyloTree_CollapseExpand`
2. **Root and clean trees**
   - `Phylo_Rooter`, then `OrthoFilter_LB` and/or `OrthoFilter_Mono`
3. **Infer duplication patterns**
   - `GD_Detector` → `GD_Visualizer`
4. **Track duplication-associated loss**
   - `GD_Loss_Tracker` → `GD_Loss_Visualizer`
5. **Optional downstream analyses**
   - `TreeTopology_Summarizer`, `Hybrid_Tracer`, `Hybrid_Visualizer`, `Ortho_Retriever`, `HaploFinder`

## Core Outputs

The exact outputs depend on the module. Frequently used files include:

- `gd_result_relaxed.txt` / `gd_result_strict.txt`: detailed GD events (`GD_Detector`)
- `gd_type_relaxed.tsv` / `gd_type_strict.tsv`: GD type counts by species-tree node
- `gd_loss_summary.txt`: event-level GD-loss records
- `gd_loss_count_summary.txt`: aggregated loss-path counts
- `numed_sptree.nwk`: numbered species tree used by visualization modules

## Input File Formats

Most modules use two-column TSV mapping files:

- `GF_ID2path.imap`: `<gene_family_id><TAB><tree_path>`
- `gene2sps.imap`: `<gene_id><TAB><species_name>`
- `gene2length.imap`: `<gene_id><TAB><gene_length>`
- Additional annotation maps (optional): family/order/clade/expression, depending on module

## Command Overview

PhyloTracer currently includes these subcommands:

- `BranchLength_NumericConverter`
- `GD_Detector`
- `GD_Loss_Tracker`
- `GD_Loss_Visualizer`
- `GD_Visualizer`
- `HaploFinder`
- `Hybrid_Tracer`
- `Hybrid_Visualizer`
- `OrthoFilter_LB`
- `OrthoFilter_Mono`
- `Ortho_Retriever`
- `PhyloSupport_Scaler`
- `PhyloTree_CollapseExpand`
- `Phylo_Rooter`
- `TreeTopology_Summarizer`
- `Tree_Visualizer`

For module-specific arguments:

```bash
phylotracer <SUBCOMMAND> -h
```

## Example Data

Example inputs are provided under `example_data/` and can be used directly to test each module.

## Documentation and Support

- GitHub: https://github.com/YiyongZhao/PhyloTracer
- PyPI: https://pypi.org/project/PhyloTracer
- Issues: https://github.com/YiyongZhao/PhyloTracer/issues

## Citation

If you use PhyloTracer in your work, please cite the corresponding paper and software release.

## License

PhyloTracer is released under the MIT License. See `LICENSE` for details.
