<div align="center">
  
# <img src="logo/PhyloTracer_logo.png" width="80" height="80" align="center"> PhyloTracer </div> 

```
###############################################################################################
                                                                                             
 ██████╗ ██████╗ ██╗   ██╗██╗      ██████╗ ████████╗██████╗  █████╗  ██████╗███████╗██████╗  
 ██╔══██╗██╔══██╗╚██╗ ██╔╝██║     ██╔═══██╗╚══██╔══╝██╔══██╗██╔══██╗██╔════╝██╔════╝██╔══██╗ 
 ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║   ██║   ██████╔╝███████║██║     █████╗  ██████╔╝ 
 ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║   ██║   ██╔══██╗██╔══██║██║     ██╔══╝  ██╔══██╗ 
 ██║     ██║  ██║   ██║   ███████╗╚██████╔╝   ██║   ██║  ██║██║  ██║╚██████╗███████╗██║  ██║ 
 ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝    ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚══════╝╚═╝  ╚═╝                             
                                                                                             
   PhyloTracer: A Versatile Toolkit for Comparative Genomics and Phylogenomics Analysis.
                                                                                             
    Pypi: https://pypi.org/project/PhyloTracer                                               
    Github: https://github.com/YiyongZhao/PhyloTracer                                        
    License: MIT license                                                                     
    Release Date: 2023-7                                                                     
    Contacts: Tao Li(l948777439@gmail.com); Yiyong Zhao(yiyong.zhao@yale.edu)
                                                                         
###############################################################################################
```
![Version](https://img.shields.io/badge/Version-1.0.3-blue)
[![CI](https://github.com/YiyongZhao/PhyloTracer/actions/workflows/ci.yml/badge.svg)](https://github.com/YiyongZhao/PhyloTracer/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/phylotracer/badge/?version=latest)](https://phylotracer.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/PhyloTracer.svg)](https://pypi.python.org/pypi/PhyloTracer)
![Python](https://img.shields.io/badge/Python-3.8--3.12-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

---
# PhyloTracer

An integrated phylogenomics toolkit designed for comprehensive post–gene-tree analysis, addressing systematic detection, quantification, and visualization of gene duplication, loss, hybridization, and subgenome fractionation signals from large-scale gene tree collections.

---
## What does PhyloTracer do?

`PhyloTracer` provides a reproducible workflow that bridges gene tree inference and biological interpretation. It employs a six-metric composite scoring framework for gene tree rooting, identifies gene duplication (GD) events and classifies them into evolutionary models (AABB/AXBB/AABX/Complex), traces per-species copy-state trajectories (2→2, 2→1, 2→0) for lineage-specific loss profiling, screens hybridization signals via HyDe D-statistics on GD-derived gene sets, and detects ancient recombination events on chromosomal synteny dotplots with subgenome-aware ortholog partitioning.

All 17 modules share a consistent data model (gene-to-species mapping, Newick trees, tab-delimited annotations) and can be used independently or combined in larger phylogenomic pipelines.

---
## Table of Contents

- [What does PhyloTracer do?](#what-does-phylotracer-do)
- [Module Features](#module-features)
- [Features](#features)
- [Getting started with PhyloTracer](#getting-started-with-phylotracer)
- [Installation](#installation)
- [Glossary](#glossary)
- [Shared Parameters](#shared-parameters)
- [Example input files](#example-input-files)
- [PhyloTracer Results Files](#phylotracer-results-files)
- [Command line options](#command-line-options)
- [Bug Reports](#bug-reports)
- [Contributing](#contributing)
- [Version History](#version-history)
- [License](#license)

---
## Module Features

PhyloTracer integrates 17 modular tools covering phylogenetic preprocessing, rooting, orthology refinement, duplication/loss detection, and hybridization analysis. Each module can run independently or be incorporated into larger evolutionary pipelines.

| # | Module | One-line summary | Key required inputs |
|---|--------|-----------------|---------------------|
| 1 | **Phylo_Rooter** | Root gene trees via six-metric composite scoring | `--input_GF_list`, `--input_imap`, `--input_sps_tree` |
| 2 | **MulRF_Distance** | Quantify gene-tree / species-tree topological conflict | `--mode`, `--input_GF_list`, `--input_imap` |
| 3 | **PhyloTree_CollapseExpand** | Collapse low-support branches into polytomies, or re-expand | `--input_GF_list`, `--support_value` |
| 4 | **PhyloSupport_Scaler** | Rescale branch support values between [0,1] and [0,100] | `--input_GF_list`, `--scale_to` |
| 5 | **BranchLength_NumericConverter** | Standardize branch-length precision | `--input_GF_list` |
| 6 | **OrthoFilter_LB** | Remove long-branch outliers (LBA artifacts) | `--input_GF_list`, `--input_imap`, `--rrbr_cutoff`, `--srbr_cutoff` |
| 7 | **OrthoFilter_Mono** | Remove monophyly-breaking tips iteratively | `--input_GF_list`, `--input_taxa`, `--input_imap`, `--input_sps_tree` |
| 8 | **TreeTopology_Summarizer** | Count and visualize topology frequencies | `--input_GF_list`, `--input_imap` |
| 9 | **Tree_Visualizer** | Render annotated gene/species trees as vector PDFs | `--input_GF_list`, `--input_imap` |
| 10 | **GD_Detector** | Detect and classify gene duplication events | `--input_GF_list`, `--input_imap`, `--input_sps_tree`, `--gd_support`, `--subclade_support`, `--dup_species_proportion`, `--dup_species_num`, `--deepvar` |
| 11 | **GD_Visualizer** | Map GD counts onto the species tree | `--input_sps_tree`, `--gd_result`, `--input_imap` |
| 12 | **GD_Loss_Tracker** | Trace per-species copy-state trajectories after GD | `--input_GF_list`, `--input_sps_tree`, `--input_imap` |
| 13 | **GD_Loss_Visualizer** | Pie charts of 2→2/2→1/2→0 distributions on species tree | `--gd_loss_result`, `--input_sps_tree` |
| 14 | **Ortho_Retriever** | Extract single-copy orthologs from paralogous gene families | `--input_GF_list`, `--input_imap`, `--input_gene_length` |
| 15 | **Hybrid_Tracer** | HyDe D-statistic hybridization screening on GD-derived gene sets | `--input_GF_list`, `--input_Seq_GF_list`, `--input_imap`, `--input_sps_tree` |
| 16 | **Hybrid_Visualizer** | Heatmaps of γ and D-statistic on species tree | `--hyde_out`, `--input_sps_tree` |
| 17 | **HaploFinder** | Subgenome dotplot analysis and FASTA partitioning | `--mode`, `--input_GF_list`/`--input_fasta`, `--input_imap` |

*Together, these modules provide a comprehensive workflow for constructing, refining, and interpreting large-scale phylogenomic data.*

---
## Features

- **Multi-criterion rooting** (via `Phylo_Rooter`):  
  Evaluates candidate roots from four strategies (outgroup, GD-node, MAD, MinVar) using a six-metric composite score (OD, BLV, GD, SO, GDC, MulRF) with direction-aware normalization and configurable weights, producing more robust roots than single-metric approaches.

- **Topology statistics** (via `TreeTopology_Summarizer`):  
  Computes absolute and relative topology frequencies for single-copy gene trees, with optional clade-level grouping and merged vector PDF visualization.

- **GD detection with type classification** (via `GD_Detector`):  
  Identifies gene duplication events via species-overlap reconciliation and classifies each into evolutionary models (AABB, AXBB, AABX, Complex) in both relaxed and strict modes.

- **Per-species copy-state tracking** (via `GD_Loss_Tracker`):  
  Traces post-duplication gene retention along the species tree with per-species 2→2 / 2→1 / 2→0 trajectories, target-species filtering, and MRCA-restricted analyses.

- **GD-informed hybridization screening** (via `Hybrid_Tracer`):  
  Extracts duplication-derived gene sets, constructs node-partitioned concatenated alignments, and runs HyDe (D-statistic) to detect introgression and hybrid speciation signals.

- **Haplotype analysis & subgenome splitting** (via `HaploFinder`):  
  Overlays GD pairs onto chromosomal synteny dotplots to detect gene conversion and crossover events; split mode reads `color_label.txt`, converts pair-level red/blue evidence into gene-level assignments, and writes subgenome-specific ortholog FASTA files.

---
## Getting started with PhyloTracer

### Option A (recommended): clone + conda + editable install

```bash
git clone https://github.com/YiyongZhao/PhyloTracer.git
cd PhyloTracer
conda env create -f environment.yml
conda activate PhyloTracer
python -m pip install -e .
```

Why this is recommended:
- `conda env create` resolves and installs all dependencies declared in `environment.yml`.
- `python -m pip install -e .` installs the current source in editable mode and registers the `PhyloTracer` command on your `PATH`.

Verify installation:
```bash
PhyloTracer -h
```

If `PhyloTracer` is not found in `PATH`, run via module entry:
```bash
PYTHONPATH=$(pwd) python -m phylotracer.Phylo_Tracer -h
```

### Option B: install released package from PyPI

```bash
python -m pip install PhyloTracer
PhyloTracer -h
```

### Quick start from GitHub ZIP (download + extract)

```bash
# 1) Download and unzip PhyloTracer-main.zip from GitHub
cd PhyloTracer-main

# 2) Create environment and install editable source
conda env create -f environment.yml
conda activate PhyloTracer
python -m pip install -e .

# 3) Run help
PhyloTracer -h

# 4) Example run (parameters intentionally set to permissive values for demonstration)
PhyloTracer GD_Detector \
  --input_GF_list example_data/10_GD_Detector/GF_ID2path.imap \
  --input_imap example_data/10_GD_Detector/gene2sps.imap \
  --input_sps_tree example_data/10_GD_Detector/sptree.nwk \
  --gd_support 50 \
  --subclade_support 50 \
  --dup_species_proportion 0 \
  --dup_species_num 2 \
  --deepvar 1
```

Optional (Linux headless / server environments — required when no display is available):
```bash
export QT_QPA_PLATFORM=offscreen
```

---
## Installation

> **Important: Python 3.13 is NOT supported.** PhyloTracer depends on `ete3`, which uses the `cgi` module that was removed in Python 3.13. Please use **Python 3.8–3.12**. We recommend creating a dedicated conda environment:
> ```bash
> conda create -n phylotracer python=3.12 -y
> conda activate phylotracer
> conda install -c conda-forge pyqt=5 -y
> pip install PhyloTracer
> export QT_QPA_PLATFORM=offscreen  # Required for headless/server environments
> ```

### Required dependencies

* Python 3.8–3.12 (Python 3.13+ is not supported due to `ete3` dependency)
* Core packages used by PhyloTracer:
  * `ete3`
  * `numpy`
  * `pandas`
  * `scipy`
  * `matplotlib`
  * `seaborn`
  * `tqdm`
  * `biopython`
  * `pypdf >= 3.0.0`
  * `pillow`
  * `pyqt5` (required by visualization modules)
  * `phyde >= 1.0.2` (HyDe Python interface, used by `Hybrid_Tracer`)

> Note: PhyloTracer uses `ete3` for tree analysis and rendering (https://etetoolkit.org/) and `HyDe` for ABBA-BABA hybridization testing (https://github.com/pblischak/HyDe).

---
## Glossary

| Term | Definition |
|------|-----------|
| **GD** | Gene Duplication — a duplication event inferred by reconciling a gene tree against the species tree |
| **MRCA** | Most Recent Common Ancestor — the deepest shared node connecting a set of taxa in a tree |
| **MulRF** | Multi-copy Robinson-Foulds distance — a topological distance metric that works at the species level for multi-copy gene families |
| **bipartition** | An internal branch that divides all leaf nodes into exactly two groups; the basic unit used when computing RF-based distances |
| **HyDe / D-statistic** | A statistical test (ABBA-BABA) that distinguishes incomplete lineage sorting (ILS) from true introgression/hybridization |
| **ILS** | Incomplete Lineage Sorting — ancestral polymorphism that is randomly resolved into discordant gene tree topologies, mimicking hybridization |
| **RRBR** | Root Relative Branch Ratio — deviation of a tip's root-to-tip distance from the tree-wide mean; used in `OrthoFilter_LB` |
| **SRBR** | Sister Relative Branch Ratio — asymmetry of a tip's root-to-tip distance relative to its sister subtree; used in `OrthoFilter_LB` |
| **LBA** | Long-Branch Attraction — a phylogenetic artifact in which two unrelated long-branch taxa are incorrectly grouped together |
| **monophyly** | A clade containing all descendants of a single common ancestor and no others |
| **polytomy** | A node with more than two immediate descendants (also called "multifurcation"); produced by collapsing low-support branches |
| **AABB / AXBB / AABX / Complex** | GD retention models: AABB = both daughter lineages retain genes from both sides of the species tree; AXBB or AABX = asymmetric loss on one side; Complex = other patterns |
| **γ (gamma)** | Admixture proportion in HyDe output — the estimated fraction of the hybrid genome derived from one parent lineage |

---
## Shared Parameters

The following parameter applies to **all 17 modules** and is not repeated in individual module sections below:

| Parameter | Description |
|-----------|-------------|
| `--output_dir DIR` | Write all output files directly into `DIR` (no extra module-named subfolder is created). If omitted, results are written to a module-specific subfolder in the current working directory. |

---
## Example input files

All `.imap` files are two-column TSV (tab-separated values). Column headers are optional; if present, the second column header is used as the display label in visualizations.

```
# GF_ID2path.imap — maps gene family IDs to gene tree file paths
# Format: <gene_family_ID> <TAB> <path_to_tree_file>
OG_104001  example_data/Phylo_Rooter/OG_104001.treefile   
OG_104002  example_data/Phylo_Rooter/OG_104002.treefile    
OG_104003  example_data/Phylo_Rooter/OG_104003.treefile

# gene2length.imap — maps gene IDs to sequence lengths (used by Ortho_Retriever)
# Format: <gene_id> <TAB> <sequence_length_in_bp>
AMTR_s00796p00010580  201
ATCG00500.1           1467
Glyma.07G273800.2     3417

# gene2sps.imap — maps gene IDs to species names (required by most modules)
# Format: <gene_id> <TAB> <species_name>
AMTR_s00796p00010580  Amborella_trichopoda
ATCG00500.1           Arabidopsis_thaliana
Glyma.07G273800.2     Glycine_max

# Family.imap — maps gene IDs to plant family labels (annotation layer for Tree_Visualizer)
# Format: <gene_id> <TAB> <family_name>
AMTR_s00796p00010580  Amborellaceae
ATCG00500.1           Brassicaceae
Glyma.07G273800.2     Fabaceae

# Order.imap — maps gene IDs to plant order labels
AMTR_s00796p00010580  Amborellales
ATCG00500.1           Brassicales
Glyma.07G273800.2     Fabales

# gene2taxa.imap — maps gene IDs to higher-level taxa
AMTR_s00796p00010580  Angiosperm
ATCG00500.1           Malvids
Glyma.07G273800.2     Fabids

# Clade.imap — maps gene IDs to functional clade labels
AMTR_s00796p00010580  Nitrogen-fixing
ATCG00500.1           Nitrogen-fixing
Glyma.07G273800.2     non-Nitrogen-fixing

# expression.csv — gene expression matrix (used by Tree_Visualizer heatmap overlay)
# Genes as row index; values should be normalized to [0, 100] before input
gene_id,SampleA,SampleB
AMTR_s00796p00010580,5.0,3.2
ATCG00500.1,12.0,9.4
Glyma.07G273800.2,0.0,0.8
```

> **Note:** You can supply any number of `.imap` annotation files to modules that accept `--gene_categories`. They are applied as sequential annotation layers from left to right alongside gene tip labels.

---
## PhyloTracer Results Files

Most modules write outputs to module-specific subdirectories (or to `--output_dir` if specified). Common outputs include:

- Rooted tree outputs: `rooted_trees/`
- MulRF distance outputs: `mulrf_distance.tsv`
- Long-branch filter outputs: `orthofilter_lb/pruned_tree/`, `orthofilter_lb/delete_gene/`
- Monophyly filter outputs: `orthofilter_mono/pruned_tree/`, `orthofilter_mono/delete_gene/`
- GD detection tables: `gd_result_*.txt`, `gd_type_*.tsv`
- GD-loss table: `gd_loss.csv`
- Hybridization outputs: `hyde_out.txt`, `hyde_filtered_out.txt`
- Ortholog retrieval outputs: `ortho_retriever_summary.txt`, `ortholog_trees.tsv`
- Topology summaries: `absolute_*.txt`, `relative_*.txt`, merged PDF summaries
- HaploFinder (haplofinder mode): `gd_pairs_dotplot.pdf`, `gd_pairs_dotplot.png`, `color_label.txt`
- HaploFinder (split mode): `haplofinder_split/split_assignment.tsv`, `haplofinder_split/split_subgenome_A.fasta`, `haplofinder_split/split_subgenome_B.fasta`, `haplofinder_split/split_subgenome_unknown.fasta`, `haplofinder_split/split_summary.txt`

---
## Command line options

> **Tip:** All modules accept `--output_dir DIR` to control where results are written. See [Shared Parameters](#shared-parameters).

---

### Rooting principle (for `Phylo_Rooter`)

`Phylo_Rooter` selects gene-tree roots by scoring multiple candidate roots and choosing the best-ranked result.

**1. Candidate root generation**

Candidates are generated from four complementary strategies:
- **Outgroup-based** — roots suggested by known outgroup taxa in the species tree
- **GD-node-based** — roots inferred from gene duplication node positions
- **MAD** (Minimal Ancestor Deviation) — a branch-length-based method
- **MinVar** (Minimum Variance) — minimizes root-to-tip distance variance

**2. Root quality metrics**

All six metrics are normalized to [0, 1] (higher = better) via direction-aware min-max normalization before weighting:

| Metric | Full name | Direction | Biological rationale |
|--------|-----------|-----------|----------------------|
| `OD` | Outgroup Depth | lower = better | Smaller topological distance from the species-tree root indicates a more basal outgroup position |
| `BLV` | Branch Length Variance | lower = better | More balanced root-to-tip distances after rooting |
| `GD` | Gene Duplication count | lower = better | Parsimony: fewer inferred duplication events is preferred |
| `SO` | Species Overlap (largest GD node) | higher = better | High species overlap at the GD node indicates a true duplication rather than an artifact |
| `GDC` | GD Consistency | higher = better | Mean of (size_symmetry × Jaccard similarity) across all GD nodes; higher = more internally consistent duplications |
| `MulRF` | Multi-copy RF distance | lower = better | Topological concordance with the species tree (see [Glossary](#glossary)) |

**3. Composite scoring**

$$
\text{score} = \sum_{m} w_m \cdot \widetilde{m}
$$

where $\widetilde{m} \in [0,1]$ is the direction-corrected normalized value of metric $m$ (higher is always better after normalization), and $w_m$ is the user-specified weight for that metric. The candidate with the **highest** composite score is selected as the optimal root.

---

### Phylo_Rooter

```
Description:
    Roots gene trees using a six-metric composite scoring framework (OD, BLV, GD, SO, GDC, MulRF)
    guided by the species tree. See "Rooting principle" above for metric definitions.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path); one tree per line
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_sps_tree        Species tree file in Newick format

Optional parameters:
    --weights               Weights for the six scoring metrics, in fixed order:
                              OD  BLV  GD  SO  GDC  MulRF
                            Enter exactly six floats that sum to 1.0.
                            Default: 0.30 0.10 0.30 0.10 0.10 0.10
                            The high default weights for OD and GD reflect a preference
                            for roots with a basal outgroup position and minimal inferred
                            duplications. If outgroup taxa are unavailable, consider
                            reducing OD and increasing GDC or MulRF accordingly.
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk

    PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk \
      --weights 0.30 0.10 0.30 0.10 0.10 0.10
```

---

### MulRF_Distance

```
Description:
    Computes species-level MulRF topological conflict distances between gene trees,
    or between gene trees and a species tree.
    Note: this is a topology-distance metric, not a sequence or genetic distance.
    See "MulRF" and "bipartition" in the Glossary for definitions.

    Supported modes:
      Mode 1 — Gene Tree vs Gene Tree: pairwise within one GF list (self-comparisons excluded)
      Mode 2 — Gene Tree vs Species Tree

    Practical uses:
      1) Quantify gene-tree/species-tree conflict intensity
         (higher distance = stronger discordance or more complex evolutionary history)
      2) Compare evolutionary pattern similarity among gene families for clustering/filtering
      3) Used internally by Phylo_Rooter to prefer roots with lower MulRF conflict

Required parameters:
    --mode                  Comparison mode: 1 = Gene Tree vs Gene Tree,
                                             2 = Gene Tree vs Species Tree
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)

Optional parameters:
    --input_sps_tree        Species tree file in Newick format (required when --mode 2)
    --output                Output TSV filename; default = mulrf_distance.tsv
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer MulRF_Distance --mode 1 --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap \
      [--output mulrf_mode1.tsv]
    PhyloTracer MulRF_Distance --mode 2 --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap \
      --input_sps_tree sptree.nwk [--output mulrf_mode2.tsv]

Output columns (mulrf_distance.tsv):
    tre_id_1                           Gene family ID of tree 1
    tre_id_2                           Gene family ID of tree 2
    gene_tree_1_leaf_count             Number of tip nodes in tree 1
    gene_tree_2_leaf_count             Number of tip nodes in tree 2
    gene_tree_1_species_count          Number of distinct species represented in tree 1
    gene_tree_2_species_count          Number of distinct species represented in tree 2
    shared_species_count               Number of species present in both trees
    mulrf_distance                     Raw MulRF topological conflict distance
    maximum_possible_mulrf_distance    Theoretical maximum MulRF distance for this pair
    normalized_mulrf_distance          mulrf_distance / maximum_possible_mulrf_distance; range [0, 1]
    shared_species_bipartition_count   Species-level bipartitions present in both trees
                                       (bipartition = an internal branch dividing leaves into two groups)
    gene_tree_1_only_bipartition_count Bipartitions unique to tree 1
    gene_tree_2_only_bipartition_count Bipartitions unique to tree 2
```

---

### PhyloTree_CollapseExpand

```
Description:
    Collapses internal branches whose support value falls below a user-defined threshold
    into polytomies (nodes with more than two children). Can also re-expand previously
    collapsed polytomies back to binary form for sensitivity analyses.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --support_value         Support value threshold for collapsing branches.
                            Branches with support BELOW this value are collapsed.
                            Use the same scale as your tree files:
                              bootstrap percentage → e.g., 50 (for 0-100 scale)
                              posterior probability → e.g., 0.5 (for 0-1 scale)
                            Default: 50

Optional parameters:
    --revert                If set, re-expand previously collapsed polytomies back to
                            binary form. Default: False
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50 \
      [--revert] [--output_dir DIR]
```

---

### PhyloSupport_Scaler

```
Description:
    Rescales branch support values (bootstrap percentages or posterior probabilities)
    to a standardized range for downstream computational compatibility.
    Useful when combining trees from different inference tools that use different scales.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --scale_to              Target support value range:
                              "1"   → rescale to [0, 1]   (posterior probability format)
                              "100" → rescale to [0, 100] (bootstrap percentage format)

Optional parameters:
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 100 [--output_dir DIR]
```

---

### BranchLength_NumericConverter

```
Description:
    Standardizes branch-length representation to a user-defined decimal precision.
    Some phylogenetic software outputs branch lengths in scientific notation (e.g., 1.23e-05)
    or with inconsistent decimal places, which can cause parsing errors in downstream tools.
    This module normalizes branch lengths while preserving the original notation style
    (decimal or scientific).

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)

Optional parameters:
    --decimal_place         Number of decimal places to retain for branch lengths.
                            Default: 10
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap \
      [--decimal_place 10] [--output_dir DIR]
```

---

### OrthoFilter_LB

**Description:** Removes tips with abnormally long branches from both single-copy and multi-copy gene family trees. Abnormally long branches can cause Long-Branch Attraction (LBA) artifacts — a well-known problem in phylogenetics where two unrelated long-branch taxa are incorrectly grouped together. See [Glossary](#glossary) for definitions of RRBR, SRBR, and LBA.

**Required parameters (CLI):** `--input_GF_list`, `--input_imap`, `--rrbr_cutoff`, `--srbr_cutoff`

**1. Root Relative Branch Ratio (RRBR)**

**Concept:** Measures how much a specific tip's root-to-tip path length deviates from the tree-wide average. A high RRBR indicates that a sequence is evolving significantly faster (or slower) than the rest of the gene family.

* **Formula:**

$$
\text{RRBR} = \frac{\text{Root-to-tip distance} - \text{Mean root-to-tip distance}}{\text{Mean root-to-tip distance}}
$$

* **Example:** RRBR = 5 means the tip's branch length is 6× the tree average (500% deviation).

**2. Sister Relative Branch Ratio (SRBR)**

**Concept:** Measures local branch-length asymmetry by comparing a tip's root-to-tip distance against its sister subtree. Even in fast-evolving gene families, a tip that is dramatically longer than its nearest neighbor is a strong candidate for removal.

* **Formula:**

$$
\text{SRBR} = \frac{\text{Root-to-tip distance} - \text{Sister root-to-tip distance}}{\text{Sister root-to-tip distance}}
$$

* **Where "Sister root-to-tip distance"** is: the root-to-tip distance of the sister tip (if sister is a leaf), or the mean root-to-tip distance across all descendant leaves of the sister subtree (if sister is an internal node).
* **Example:** SRBR = 2.5 means the tip's root-to-tip distance is 3.5× that of its sister.

**3. Long-branch decision mode**

Controls how RRBR and SRBR are combined:

* `or` (lenient): remove a tip when `RRBR >= rrbr_cutoff` **OR** `SRBR >= srbr_cutoff`. Recommended for most datasets.
* `and` (strict): remove a tip only when `RRBR >= rrbr_cutoff` **AND** `SRBR >= srbr_cutoff`. Use when false-positive removal is a concern.

```
Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --rrbr_cutoff           RRBR threshold for flagging outlier tips.
                            Default: 5 (tip must deviate >500% from the tree-wide mean)
                            For datasets with high divergence between lineages,
                            consider relaxing to 8-10. For conserved gene families, tighten to 3.
    --srbr_cutoff           SRBR threshold for flagging locally asymmetric tips.
                            Default: 2.5 (tip must be >250% longer than its sister subtree)

Optional parameters:
    --lb_mode               Long-branch decision mode: or | and. Default: or
    --visual                If set, export before/after tree visualization PDFs. Default: False
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap \
      --rrbr_cutoff 5 --srbr_cutoff 2.5 [--lb_mode or] [--visual] [--output_dir DIR]
```

---

### OrthoFilter_Mono

**Description:** Iteratively removes tips that break the expected monophyly of predefined taxonomic groups (e.g., a Brassicaceae gene nested inside a Fabaceae clade). This pattern can indicate horizontal gene transfer, contamination, or misannotation. Candidates for removal are ranked by a composite score and pruning stops once the target purity or removal cap is reached.

**Required parameters (CLI):** `--input_GF_list`, `--input_taxa`, `--input_imap`, `--input_sps_tree`

**Scoring logic:**

**1. Dominant Lineage Purity** — Measures how strongly the dominant lineage is composed of target-taxon tips.

$$\text{Purity}=\frac{N_{\text{target}}}{N_{\text{dominant tips}}}$$

**2. Phylogenetic Distance Score** — Alien lineages that are mapped farther from the target MRCA in the species tree receive a higher score and are more likely to be removed.

$$\text{PhyloDist}=\text{Depth}(\text{MRCA}(\text{target}))-\text{Depth}(\text{MRCA}(\text{target}\cup\text{alien}))$$

**3. Alien Coverage Score** — Alien lineages that occupy fewer tips within the dominant lineage (i.e., are more isolated) are more likely to be noise.

$$\text{AlienCov}=\frac{N_{\text{alien}}}{N_{\text{dominant tips}}}$$

**4. Alien Depth Score** — Alien lineages inserted more deeply relative to the dominant lineage root score higher for removal.

$$\text{AlienDepth}=\text{Depth}(\text{alien})-\text{Depth}(\text{MRCA}(\text{dom}))$$

**5. Combined Ranking Score** — Candidates are ranked by a multiplicative score using normalized components:

$$\text{Combined}=\text{Norm}(\text{PhyloDist})\times\text{Norm}(\text{AlienDepth})\times\left(-\log_{10}\left(\text{AlienCov} + 10^{-4}\right)\right)$$

**6. Removal Stopping Rules** — Pruning stops when the dominant-lineage purity reaches `purity_cutoff`, or when the removal cap `max_remove` is reached:

$$\text{max\_remove}=\max\left(\text{max\_remove\_fraction}\times N_{\text{dominant tips}},\ 1\right)$$

> **Caution:** The taxonomic groupings in `--input_taxa` should reflect well-established evolutionary relationships. Using poorly supported groupings may cause over-pruning of genuine diversity.

```
Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_taxa            Two-column mapping file (gene_id <TAB> clade_or_lineage_label)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_sps_tree        Species tree file in Newick format

Optional parameters:
    --purity_cutoff         Stop pruning when the dominant lineage reaches this purity fraction.
                            Default: 0.95  (i.e., stop when ≥95% of tips belong to the target taxon)
    --max_remove_fraction   Maximum fraction of tips that may be removed per gene tree.
                            Default: 0.5  (at most 50% of dominant-lineage tips removed)
    --visual                If set, export before/after pruning visualization PDFs. Default: False
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_taxa Clade.imap \
      --input_imap gene2sps.imap --input_sps_tree sptree.nwk \
      [--purity_cutoff 0.95] [--max_remove_fraction 0.5] [--visual] [--output_dir DIR]
```

---

### TreeTopology_Summarizer

```
Description:
    Enumerates and visualizes the frequency of both absolute and relative topologies
    across single-copy gene trees or user-defined predefined clades.
    Visualization outputs are vector PDFs: each topology panel is rendered at A4 width
    and merged in a single-column, top-to-bottom layout (no PNG rasterization).

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)

Optional parameters:
    --visual_top            Number of top-ranked (most frequent) topologies to include
                            in the visualization output. Default: 10
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap \
      --input_imap gene2sps.imap [--visual_top 10] [--output_dir DIR]
```

---

### Tree_Visualizer

```
Description:
    Renders annotated multi-layer gene trees with tip color annotations, optional
    expression heatmap overlays, and GD-node markers. Also maps gene duplication
    counts onto species trees. All outputs are publication-ready vector PDFs;
    category colors are assigned consistently by label across all output files.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)

Optional parameters:
    --keep_branch           Whether to preserve branch lengths in the plot: 1 = yes, 0 = no
    --tree_style            Tree layout style: r = rectangular, c = circular
    --gene_categories       One or more two-column annotation files (gene_id <TAB> category_label).
                            Each file adds one annotation layer displayed to the right of tip labels.
                            First line can optionally be a header ("gene_id <TAB> HeaderName");
                            HeaderName is used as the column title (e.g., Family, Order, Clade;
                            avoid spaces or special characters in header names).
                            Note: the FIRST file supplied is used as the family map for
                            species-tree duplication-count mapping.
                            Example files: Family.imap, Order.imap, Clade.imap
    --input_sps_tree        Species tree file in Newick format
    --heatmap_matrix        Gene expression (or other numeric) matrix file.
                            Recommended format: tab-delimited .txt or .tsv; also accepts .csv/.xls/.xlsx.
                            Genes as row index. Values should be in [0, 100]; normalize before input if needed.
    --visual_gd             If set, overlay predicted GD nodes on gene-tree figures. Default: False
    --gd_support            Minimum bootstrap/posterior support for a node to be considered a GD
                            candidate when --visual_gd is active. Range: 0-100. Default: 50
    --subclade_support      Minimum support required in both child subclades of a GD node
                            when --visual_gd is active. Range: 0-100. Default: 0
    --dup_species_proportion  Minimum overlap ratio of duplicated species between the two GD
                            child clades when --visual_gd is active. Range: 0-1. Default: 0.2
    --dup_species_num       Minimum number of species that must appear in both child clades of
                            a GD node when --visual_gd is active. Default: 2
    --deepvar               Maximum tolerated depth-variance score for GD node detection
                            when --visual_gd is active (see GD_Detector for full explanation).
                            Default: 1
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer Tree_Visualizer \
      --input_GF_list GF_ID2path.visual20.imap --input_imap gene2sps.imap \
      --gene_categories Family.imap Order.imap Clade.imap \
      --input_sps_tree sptree.nwk --heatmap_matrix heatmap_matrix.txt \
      --keep_branch 1 --tree_style r \
      --visual_gd --gd_support 50 --subclade_support 0 \
      --dup_species_proportion 0.2 --dup_species_num 2 --deepvar 1 \
      [--output_dir DIR]
```

---

### GD_Detector

```
Description:
    Identifies gene duplication (GD) events by reconciling gene trees against a species tree
    via species-overlap analysis. Each GD event is classified into one of four retention models:
      AABB    — both child lineages retain genes from both sides of the species tree (symmetric)
      AXBB    — asymmetric loss: one child lineage loses genes on the "A" side
      AABX    — asymmetric loss: one child lineage loses genes on the "B" side
      Complex — does not fit the above symmetric or simple asymmetric patterns
    Supports relaxed and strict detection modes (see --gdtype_mode).

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_sps_tree        Species tree file in Newick format
    --gd_support            Minimum support for a node to be considered a GD candidate.
                            Accepted range: 0-100. Typical values: 50-100. Default: 50
    --subclade_support      Minimum support required in both child subclades of a candidate
                            GD node. Accepted range: 0-100. Default: 0
                            Increase this to reduce false-positive GD calls in weakly
                            supported parts of the tree.
    --dup_species_proportion  Minimum overlap ratio of duplicated species between the two GD
                            child clades. Range: 0-1. Default: 0.2
                            A value of 0 allows GD nodes with no species overlap (permissive);
                            higher values (e.g., 0.5) require more symmetric species coverage.
    --dup_species_num       Minimum number of species that must appear in both child clades
                            of a GD node. Default: 2
    --deepvar               Maximum tolerated depth-variance score.
                            Depth variance measures how symmetrically the species under a GD
                            node are distributed in the species tree. A value of 0 requires
                            perfect depth symmetry; higher values allow more asymmetric GD
                            events to pass. Default: 1
                            For polyploid taxa or rapidly radiating lineages, consider
                            increasing to 2-3 to avoid discarding genuine GD events.

Optional parameters:
    --gdtype_mode           GD type assignment mode. Default: relaxed
                              relaxed — assigns type using species-overlap mapping only;
                                        higher recall, may include some false positives
                              strict  — additionally filters nodes using --deepvar depth
                                        consistency; higher precision, may miss some
                                        asymmetric GD events
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer GD_Detector \
      --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap \
      --input_sps_tree sptree.nwk \
      --gd_support 50 --subclade_support 50 \
      --dup_species_proportion 0 --dup_species_num 2 --deepvar 1 \
      [--gdtype_mode relaxed] [--output_dir DIR]
```

---

### GD_Visualizer

```
Description:
    Visualizes gene duplication detection results by mapping GD counts and types
    onto a species tree context.

Required parameters:
    --input_sps_tree        Numbered species tree in Newick format. Use the numbered tree
                            output from GD_Detector (e.g., numed_sptree.nwk), which contains
                            internal node labels required for GD mapping.
    --gd_result             GD result table produced by GD_Detector
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)

Optional parameters:
    --output                Output PDF file path. Default: <gd_result_basename>.pdf
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer GD_Visualizer \
      --input_sps_tree numed_sptree.nwk \
      --gd_result gd_result_relaxed.txt \
      --input_imap gene2sps.imap \
      [--output gd_result_relaxed.pdf]
```

---

### GD_Loss_Tracker

```
Description:
    Traces the fate of duplicated gene copies across species after each GD event.
    For every GD event and every species, assigns a copy-state label:
      2→2  — species retains both gene copies
      2→1  — species lost one of the two copies
      2→0  — species lost both copies
    Results are exported as detailed Excel reports and can be filtered by target species
    or restricted to GD events under a specified ancestral node.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_sps_tree        Species tree file in Newick format
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)

Optional parameters:
    --target_species        Only count loss paths ending at this species
                            (e.g., --target_species Arabidopsis_thaliana).
                            Can be specified multiple times for multiple species.
    --mrca_node             Only count loss paths passing through the MRCA of two named species.
                            Format: SpeciesA,SpeciesB (comma-separated, no spaces).
                            Can be specified multiple times.
    --include_unobserved_species
                            Policy for species absent from the current gene family tree.
                            If set: absent species are still assigned 2→2/2→1/2→0 labels
                              based on the presence/absence pattern in the left/right child clades.
                            If not set (default): absent species are labeled "missing_data".
                            Note: this flag only affects classification labels; it does not
                            change aggregate count statistics.
    --node_count_mode       How to count loss paths for path_count_* transition statistics.
                            Default: parsimony. Choices: parsimony | accumulate
                              parsimony  — if multiple species lost a gene in the same ancestor,
                                           that is counted as one shared loss event (reduces
                                           redundancy and is appropriate for most analyses)
                              accumulate — each independent loss path is counted separately
                                           (use when per-lineage resolution is needed)
    --parsimony_min_support_ratio
                            Minimum fraction of descendant species required to collapse
                            multiple losses into one shared parsimony node. Default: 0.5
                            Evaluated relative to all species under the current GD event node.
    --parsimony_min_support_species
                            Minimum number of descendant species required to collapse
                            multiple losses into one shared parsimony node. Default: 2
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer GD_Loss_Tracker \
      --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap \
      [--target_species Arabidopsis_thaliana] [--mrca_node SpeciesA,SpeciesB] \
      [--include_unobserved_species] [--node_count_mode parsimony] [--output_dir DIR]
```

---

### GD_Loss_Visualizer

```
Description:
    Renders pie charts of copy-state distributions (2→2 / 2→1 / 2→0) at each node
    of the species tree, enabling visual identification of lineage-specific gene loss patterns.

Required parameters:
    --gd_loss_result        CSV table generated by GD_Loss_Tracker (gd_loss.csv)
    --input_sps_tree        Numbered species tree in Newick format (numed_sptree.nwk from GD_Detector)

Optional parameters:
    --output                Output PDF file path. Default: gd_loss_pie_visualizer.pdf
    --output_dir            See Shared Parameters.

Usage:
    PhyloTracer GD_Loss_Visualizer \
      --input_sps_tree numed_sptree.nwk \
      --gd_loss_result gd_loss.csv \
      [--output gd_loss_pie_visualizer.pdf]
```

Example data layout:
`example_data/13_GD_Loss_Visualizer/parsimony/` and
`example_data/13_GD_Loss_Visualizer/accumulate/` each contain the rendered PDF
plus the exact `gd_loss.csv`, `gd_loss_node_summary.tsv`, and `numed_sptree.nwk`
files used for that mode.

---

### Ortho_Retriever

```
Description:
    Infers phylogenetically supported single-copy putative orthologs from large-scale
    gene family trees by recursively splitting paralogous clades. Uses gene length to
    resolve within-species paralog conflicts (longer genes are preferred as the ortholog).
    Optionally refines candidate sets using synteny block support, and can attach
    strict sister-clade outgroup genes to each written tree.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_gene_length     Two-column mapping file (gene_id <TAB> gene_length_in_bp)

Optional parameters:
    --input_synteny_blocks  Raw synteny block file for additional ortholog refinement.
                            Each block starts with "#"; each non-comment line contains one gene pair.
    --add_outgroup          Attach outgroup genes from the strict sister clade to each output tree.
                            Outgroup species must not overlap with ingroup species and must be
                            single-copy within the candidate sister clade.
    --output_dir            See Shared Parameters.

Outputs:
    ortho_retriever_summary.txt    Summary of inferred orthologs and written trees.
    ortholog_trees.tsv             Per-gene-family ortholog group assignments.
    ortholog_outgroup_report.tsv   (only when --add_outgroup) Selected outgroup genes per tree,
                                   with per-tree status (ok | skip_no_valid_outgroup).

Usage:
    PhyloTracer Ortho_Retriever \
      --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap \
      --input_gene_length gene2length.imap \
      [--input_synteny_blocks collinearity] [--add_outgroup] [--output_dir DIR]
```

---

### Hybrid_Tracer

```
Description:
    Screens for hybridization (introgression) signals using gene sets derived from GD events.
    For each GD node, relevant sequences are extracted, aligned, and tested with HyDe's
    D-statistic (ABBA-BABA test) to distinguish incomplete lineage sorting (ILS) from
    true hybridization. See Glossary for definitions of HyDe, D-statistic, and ILS.

    Before running HyDe, Hybrid_Tracer filters candidate species-tree nodes using two
    configurable heuristics (--min_gd_count and --min_asymmetric_ratio) to exclude nodes
    unlikely to carry hybridization signal.

    For each passing node, one outgroup species is selected per HyDe run by searching
    the sister branch in topological-distance order. GD events are grouped by the first
    valid outgroup candidate found; each group is tested in a separate HyDe run.
    Alignment columns with excessive missing data (fewer than 3 taxa with non-gap entries
    by default) are removed before testing.

Required parameters:
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_Seq_GF_list     Tab-delimited mapping file (GF_ID <TAB> alignment_path).
                            GF IDs must match those in --input_GF_list exactly.
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_sps_tree        Species tree file in Newick format

Optional parameters:
    --mrca_node             Restrict analysis to the subtree defined by the MRCA of two species.
                            Format: SpeciesA,SpeciesB (comma-separated, no spaces).
                            If multiple are provided, only the first valid pair is used.
    --outgroup_species      Explicitly specify the outgroup species. When used with --mrca_node,
                            only genes from this species are used as outgroup.
    --split_groups          Number of partitions for HyDe batch processing. Default: 1
                            Increase for large datasets to parallelize HyDe runs.
    --min_gd_count          Minimum number of classified GD events required at a species-tree
                            node to include it in HyDe testing. Default: 10
                            Nodes with too few GD events lack statistical power for HyDe.
    --min_asymmetric_ratio  Minimum fraction of asymmetric GD models (AXBB + AABX) required
                            at a species-tree node to include it in HyDe testing.
                            Accepted range: 0.0-1.0. Default: 0.1
                            Background: asymmetric gene loss after duplication is a hallmark
                            of hybridization/introgression-driven genome fractionation.
                            Nodes with very low asymmetric ratios are more consistent with
                            simple speciation and are excluded from hybridization testing.
    --output_dir            See Shared Parameters.

Outputs:
    hyde_out.txt              Full HyDe test results
    hyde_filtered_out.txt     Filtered HyDe results (significant tests only)
    hyde_summary.txt          Run summary
    hyde_event_assignments.tsv  GD event to HyDe test assignments
    (All HyDe tables include an "Outgroup" column indicating which outgroup species was used.)

Usage:
    PhyloTracer Hybrid_Tracer \
      --input_GF_list gf.txt --input_Seq_GF_list gf_aln.txt \
      --input_sps_tree sptree.nwk --input_imap gene2sps.imap \
      [--mrca_node SpeciesA,SpeciesB] [--outgroup_species SpeciesX] \
      [--split_groups 2] [--min_gd_count 10] [--min_asymmetric_ratio 0.1] \
      [--output_dir DIR]
```

---

### Hybrid_Visualizer

```
Description:
    Visualizes hybridization signals on the species tree by displaying admixture
    proportions (γ) and D-statistic support values as heatmaps.

Required parameters:
    --hyde_out              HyDe output table file (hyde_out.txt from Hybrid_Tracer)
    --input_sps_tree        Species tree file in Newick format

Optional parameters:
    --node                  Use node-mode heatmaps (monophyletic clade stacking) instead
                            of the default leaf-mode output.
    --output_dir            See Shared Parameters.

Legend:
    Red    = focal hybrid species or clade
    Blue   = target internal node label (--node mode only)
    Yellow = admixture proportion (γ) values
    White  = labels for tested hybridization combinations

Usage:
    PhyloTracer Hybrid_Visualizer \
      --hyde_out hyde_out.txt --input_sps_tree sptree.nwk \
      [--node] [--output_dir DIR]
```

---

### HaploFinder

```
Description:
    haplofinder mode — maps duplicated gene pairs onto chromosome-level synteny dotplots
      to identify regions of ancient gene conversion and crossover between subgenomes.
    split mode — reads color_label.txt produced by haplofinder mode, aggregates pair-level
      red/blue subgenome evidence into gene-level subgenome assignments by majority vote,
      and writes separate FASTA files for downstream polyploid genome evolution analysis.

    Red/blue colors in color_label.txt represent the two subgenomes inferred from the
    dotplot analysis. In split mode, each gene is assigned to subgenome A (if red_count >
    blue_count) or subgenome B (if blue_count > red_count).

Required parameters:
    --mode                  Run mode: "haplofinder" (default) or "split"

--- Mode = haplofinder required parameters ---
    --input_GF_list         Tab-delimited mapping file (GF_ID <TAB> gene_tree_path)
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_sps_tree        Species tree file in Newick format
    --species_a             Name of species A (must match names in --input_imap)
    --species_b             Name of species B (must match names in --input_imap)
    --species_a_gff         Genome annotation file for species A in GFF/GTF-compatible format
    --species_b_gff         Genome annotation file for species B in GFF/GTF-compatible format
    --species_a_lens        Chromosome-length file for species A; format: chr <TAB> length
    --species_b_lens        Chromosome-length file for species B; format: chr <TAB> length

--- Mode = haplofinder optional parameters ---
    --gd_support            Minimum support of GD nodes used for gene pair extraction.
                            Accepted range: 0-100. Default: 50
    --pair_support          Minimum support of ortholog/speciation pair nodes.
                            Accepted range: 0-100. Default: 50
    --visual_chr_a          File listing a subset of species A chromosomes to visualize
                            (one chromosome name per line). If omitted, all chromosomes are shown.
    --visual_chr_b          File listing a subset of species B chromosomes to visualize.
    --size                  Point size for dotplot rendering (positive float). Default: 0.0005
                            This is a relative unit: larger values produce bigger dots.
                            For short chromosomes or sparse gene pairs, try 0.001-0.005.
    --min_shared_pairs      Minimum number of gene pairs required between two chromosomes
                            to infer a collinear homolog relationship. Default: 5
                            Increase this threshold for noisy or fragmented assemblies.
    --min_conv_pairs        Minimum number of gene pairs on a chromosome pair required before
                            attempting gene conversion detection. Default: 10
    --n_permutations        Number of permutations for the gene conversion permutation test.
                            Default: 1000
    --p_threshold           Significance threshold for both the binomial test and the
                            permutation test in gene conversion detection. Default: 0.05
    --output_dir            See Shared Parameters.

--- Mode = split required parameters ---
    --input_imap            Two-column mapping file (gene_id <TAB> species_name)
    --input_fasta           Input FASTA file (.fa or .fasta)
    --color_label_file      Pair-level color label file from haplofinder mode (color_label.txt)
    --hyb_sps               Name of the hybrid species whose genes will be split by majority color
    --species_b_gff         Genome annotation file of the hybrid species, used to report
                            chromosome coordinates in the split output

--- Mode = split optional parameters ---
    --red_subgenome         Subgenome label assigned when red_count > blue_count. Default: A
    --blue_subgenome        Subgenome label assigned when blue_count > red_count. Default: B
    --cluster_file          Legacy compatibility field; optional and used only for logging.
    --chrs_per_subgenome    Legacy compatibility parameter; primary split assignment is now
                            determined by color-label majority vote, not this field.
    --output_dir            See Shared Parameters.

Split mode output (split_assignment.tsv) reports for each hybrid gene:
    red_count   — number of red (subgenome A) pair votes
    blue_count  — number of blue (subgenome B) pair votes
    status      — ok | tie | no_color

Usage (haplofinder mode):
    PhyloTracer HaploFinder \
      --mode haplofinder \
      --input_GF_list gf.txt --input_imap gene2sps.imap --input_sps_tree sptree.nwk \
      --species_a arh --species_b ard \
      --species_a_gff arh.gff --species_b_gff ard.gff \
      --species_a_lens arh.lens --species_b_lens ard.lens \
      [--gd_support 50] [--pair_support 50] \
      [--visual_chr_a chr_a.txt --visual_chr_b chr_b.txt --size 0.0001] \
      [--output_dir DIR]

Usage (split mode):
    PhyloTracer HaploFinder \
      --mode split \
      --input_imap gene2sps.imap --input_fasta proteins.fa \
      --color_label_file color_label.txt --hyb_sps Hybrid \
      --species_b_gff hybrid.gff \
      [--red_subgenome A --blue_subgenome B] [--output_dir DIR]
```

Bundled example data are available under `example_data/17_Haplofinder/`. Optional chromosome list files `chr_a.txt` and `chr_b.txt` are included for the haplofinder-mode example.

---
## Bug Reports

You can report bugs or request features through our [GitHub Issues page](https://github.com/YiyongZhao/PhyloTracer/issues). If you have any questions, suggestions, or issues, please do not hesitate to contact us.

## Contributing

If you're interested in contributing code or reporting bugs, we welcome your ideas and contributions to improve PhyloTracer! Please check out [Contribution Guidelines](https://docs.github.com/en/issues).

## Version History

Check the [Changelog](https://github.com/YiyongZhao/PhyloTracer/commits/PhyloTracer_v1.0.0) for details on different versions and updates.

## License

PhyloTracer is licensed under the [MIT LICENSE](LICENSE).
