<div align="center">
  
# <img src="logo/PhyloTracer_logo.png" width="80" height="80" align="center"> PhyloTracer </div> 

```
###############################################################################################
                                                                                             
 ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ████████╗██████╗  █████╗  ██████╗███████╗██████╗  
 ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗╚══██╔══╝██╔══██╗██╔══██╗██╔════╝██╔════╝██╔══██╗ 
 ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║   ██║   ██████╔╝███████║██║     █████╗  ██████╔╝ 
 ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║   ██║   ██╔══██╗██╔══██║██║     ██╔══╝  ██╔══██╗ 
 ██║     ██║  ██║   ██║   ███████╗╚██████╔╝   ██║   ██║  ██║██║  ██║╚██████╗███████╗██║  ██║ 
 ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝    ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚══════╝╚═╝  ╚═╝                             
                                                                                             
   PhyloTracer: A Versatile Toolkit for Comparative Genomics and Phylogenomics Analysis.
                                                                                             
    Pypi: https://pypi.org/project/PhyloTracer                                               
    Github: https://github.com/YiyongZhao/PhyloTracer                                        
    License: MIT license                                                                     
    Release Date: 2023-7                                                                     
    Contacts: Taoli(l948777439@163.com); Yiyong Zhao(yiyong.zhao@yale.edu)
                                                                         
###############################################################################################
```
![Version](https://img.shields.io/badge/Version-1.0.0-blue)
[![CI](https://github.com/YiyongZhao/PhyloTracer/actions/workflows/ci.yml/badge.svg)](https://github.com/YiyongZhao/PhyloTracer/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/phylotracer/badge/?version=latest)](https://phylotracer.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/PhyloTracer.svg)](https://pypi.python.org/pypi/PhyloTracer)
![Python](https://img.shields.io/badge/Python-3.8--3.12-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

---
# PhyloTracer

A user-friendly toolkit for gene tree rooting, topology summarization, species hybridization signal screening, gene-duplication (GD) and loss profiling, and subgenome-aware ortholog splitting, with practical utilities for tree format manipulation and visualization.

---
## What does PhyloTracer do?

`PhyloTracer` provides a reproducible workflow centered on accurate gene tree rooting and topology statistics. It further offers utilities for screening hybridization-like signals (via topology patterns such as ABAB/ABBA variants), summarizing GD and loss patterns, and subgenome-informed splitting of multi-copy families.
All modules are designed to be used independently or combined in larger phylogenomic pipelines. Where applicable, methods are documented with input assumptions and recommended validation steps to ensure rigorous interpretation.

---
## Table of Contents

- [What does PhyloTracer do?](#what-does-phylotracer-do)
- [Module Features](#module-features)
- [Features](#features)
- [Getting started with PhyloTracer](#getting-started-with-phylotracer)
- [Advanced installation notes](#advanced-installation-notes)
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

1. **Phylo_Rooter:** Implements an accurate, automated rooting algorithm for gene trees to enhance evolutionary inference.
2. **MulRF_Distance:** computes species-level topological conflict distance between gene trees (single-copy or multi-copy) and a species tree; useful for quantifying gene-tree/species-tree discordance and for comparing duplication/divergence patterns across gene families.
3. **PhyloTree_CollapseExpand:** Transforms a phylogenetic tree into a “comb-like” structure based on a predefined support threshold, and re-expands it back to binary form when needed.
4. **PhyloSupport_Scaler:** Recalibrates branch support values to standardized scales ([0–1] or [1–100]) for consistent computation.
5. **BranchLength_NumericConverter:** Converts branch-length strings to numeric format for downstream quantitative analyses.
6. **OrthoFilter_LB:** Removes tips with excessively long branches to eliminate phylogenomic noise.
7. **OrthoFilter_Mono:** Prunes non-monophyletic outliers and paralogs under user-defined taxonomic constraints.
8. **TreeTopology_Summarizer:** Summarizes frequencies of absolute and relative topologies across gene trees or predefined clades.
9. **Tree_Visualizer:** Visualizes duplication events, node labels, and gene-associated numeric profiles on gene and species trees.
10. **GD_Detector:** Identifies gene duplication events via reconciliation between gene and species trees.
11. **GD_Visualizer:** Displays detected duplication nodes in a species tree context.
12. **GD_Loss_Tracker:** Tracks duplication-loss patterns following major GD bursts across species tree nodes.
13. **GD_Loss_Visualizer:** Visualizes node/tip-specific gene loss summaries.
14. **Ortho_Retriever:** Extracts putative single-copy orthologs by recursively splitting paralogous clades.
15. **Hybrid_Tracer:** Detects hybridization signals from duplicated genes using coalescent-based phylogenetic invariants.
16. **Hybrid_Visualizer:** Highlights hybridization proportions (γ) and support values across the species tree.
17. **HaploFinder:** Identifies ancient recombination (conversion/crossover) events by tracing subgenome haplotypes.

*Together, these modules provide a comprehensive workflow for constructing, refining, and interpreting large-scale phylogenomic data.*

---
## Features

- Rooting (via `Phylo_Rooter`):  
  Automatically evaluates multiple rooting candidates and selects the most plausible root using a six-metric composite score (OD, BLV, GD, SO, GDC, MulRF) with direction-aware min-max normalization. Supports empirical or entropy (EWM) weight strategies.

- Topology statistics (via `TreeTopology_Summarizer`):  
  Computes absolute and relative topology frequencies for single-copy gene trees, with optional grouped summarization by user-defined clade labels.

- Hybridization screening (via `Hybrid_Tracer`):  
  Extracts GD-informed gene sets and performs node-focused hybridization signal screening in grouped or ungrouped modes.

- GD & loss profiling (via `GD_Detector` and `GD_Loss_Tracker`):  
  Reconciles gene–species trees to summarize duplication events and downstream lineage-specific loss patterns.

- HaploFinder:  
  Detects ancient recombination signals (including gene conversion/crossover) by tracing subgenome haplotype patterns.

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
- `conda env create` installs dependencies.
- `python -m pip install -e .` installs the current source in editable mode and registers the `PhyloTracer` command.

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

# 4) Example run
PhyloTracer GD_Detector \
  --input_GF_list example_data/GD_Detector/GF_ID2path.imap \
  --input_imap example_data/GD_Detector/gene2sps.imap \
  --input_sps_tree example_data/GD_Detector/sptree.nwk \
  --gd_support 50 \
  --subclade_support 50 \
  --dup_species_proportion 0 \
  --dup_species_num 2 \
  --deepvar 1
```

Optional (Linux headless visualization compatibility):
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

### Required dependencies:

* Python 3.8–3.12 (Python 3.13+ is not supported due to ete3 dependency)
* Core modules used by PhyloTracer:
  * ete3
  * numpy
  * pandas
  * scipy
  * matplotlib
  * seaborn
  * tqdm
  * biopython
  * pypdf>=3.0.0
  * pillow
  * pyqt5 (for visualization modules)
  * phyde>=1.0.2 (HyDe Python interface)
    
Note: PhyloTracer uses basic functions of analysis and visualization of trees from Python framework [ete3](https://etetoolkit.org/) and detects species hybridization signals using ABAB-BABA test by [HyDe](https://github.com/pblischak/HyDe).

---
## Example input files
The following input file should have two columns and be separated by a tab key.
```
Provide a two-column file in TSV format: each line contains <gene_family_IDs><TAB><file_paths>
------------GF_ID2path.imap-------------------------------------------------------------------------------------------------------
OG_104001  example_data/Phylo_Rooter/OG_104001.treefile   
OG_104002  example_data/Phylo_Rooter/OG_104002.treefile    
OG_104003  example_data/Phylo_Rooter/OG_104003.treefile

Provide a two-column file in TSV format: each line contains <gene_id><TAB><sequence_length>
------------gene2length.imap------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  201
ATCG00500.1           1467
Glyma.07G273800.2     3417

Provide a two-column file in TSV format: each line contains <gene_id><TAB><species>
------------gene2sps.imap---------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborella_trichopoda
ATCG00500.1           Arabidopsis_thaliana
Glyma.07G273800.2     Glycine_max

Provide a two-column file in TSV format: each line contains <gene_id><TAB><plant_family>
------------Family.imap------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborellaceae
ATCG00500.1           Brassicaceae
Glyma.07G273800.2     Fabaceae

Provide a two-column file in TSV format: each line contains <gene_id><TAB><plant_order>
------------Order.imap-------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborellales
ATCG00500.1           Brassicales
Glyma.07G273800.2     Fabales

Provide a two-column file in TSV format: each line contains <gene_id><TAB><higher-level_taxa>
------------gene2taxa.imap--------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Angiosperm
ATCG00500.1           Malvids
Glyma.07G273800.2     Fabids

Provide a two-column file in TSV format: each line contains <gene_id><TAB><functional_clade>
------------Clade.imap-------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Nitrogen-fixing
ATCG00500.1           Nitrogen-fixing
Glyma.07G273800.2     non-Nitrogen-fixing

Provide an expression matrix file (CSV/XLS/XLSX) with gene IDs as row index.
------------expression.csv---------------------------------------------------------------------------------------------------------
gene_id,SampleA,SampleB
AMTR_s00796p00010580,5.0,3.2
ATCG00500.1,12.0,9.4
Glyma.07G273800.2,0.0,0.8

#Note: You can add any number of imap files. They will sequentially provide annotations to the right of the gene tips according to the order of input.
```
---
## PhyloTracer Results Files

Most modules generate task-specific outputs in either the current working directory or module-specific subdirectories. Common outputs include:

- Rooted tree outputs: `rooted_trees/`
- MulRF distance outputs: `mulrf_distance/mulrf_distance.tsv`
- Long-branch filter outputs: `orthofilter_lb/pruned_tree/`, `orthofilter_lb/delete_gene/`
- Monophyly filter outputs: `orthofilter_mono/pruned_tree/`, `orthofilter_mono/delete_gene/`
- GD detection tables: `gd_result_*.txt`, `gd_type_*.tsv`
- GD-loss tables: `gd_loss_summary.txt`, `gd_loss_count_summary.txt`, `gd_loss.xlsx`
- Hybridization outputs: `hyde_out.txt`, `hyde_filtered_out.txt`
- Ortholog retrieval outputs: `ortho_retriever_summary.txt`, `ortholog_trees.tsv`
- Topology summaries: `absolute_*.txt`, `relative_*.txt`, merged PNG summaries
- HaploFinder (haplofinder mode): `gd_pairs_dotplot.pdf`, `gd_pairs_dotplot.png`, `color_label.txt`
- HaploFinder (split mode): `haplofinder_split/split_assignment.tsv`, `haplofinder_split/split_subgenome_A.fasta`, `haplofinder_split/split_subgenome_B.fasta`, `haplofinder_split/split_subgenome_unknown.fasta`, `haplofinder_split/split_summary.txt`

---
## Command line options

This section follows an OrthoFinder-like CLI reference style with compact layout and explicit parameter meanings.

### Rooting principle (for `Phylo_Rooter`)

`Phylo_Rooter` selects gene-tree roots by scoring multiple candidate roots and choosing the best-ranked result.

**1. Candidate root generation**

- Outgroup-based candidates
- GD-node-based candidates
- MAD candidates
- MinVar candidates

**2. Root quality metrics**

All six metrics are normalized to [0, 1] (higher = better) via direction-aware min-max normalization before weighting:

| Metric | Full name | Direction | Biological rationale |
|--------|-----------|-----------|----------------------|
| `OD` | Outgroup Depth | lower raw = better | Smaller topological distance from species-tree root = more basal outgroup |
| `BLV` | Branch Length Variance | lower raw = better | More balanced branch lengths after rooting |
| `GD` | Gene Duplication count | lower raw = better | Parsimony: fewer inferred GD events |
| `SO` | Species Overlap (largest GD node) | higher raw = better | High overlap indicates true duplications, not artefacts |
| `GDC` | GD Consistency | higher raw = better | Mean (size_symmetry × Jaccard) across GD nodes |
| `MulRF` | Multi-copy RF distance | lower raw = better | Topological concordance with species tree |

**3. Composite scoring**

$$
\text{score} = \sum_{m} w_m \cdot \widetilde{m}
$$

where $\widetilde{m} \in [0,1]$ is the direction-corrected normalized value of metric $m$ (higher is always better), and $w_m$ is the weight for that metric. The candidate with the **highest** composite score is selected as the optimal root.

**4. Weight strategies (`--weight-strategy`)**

- `empirical` (default): uses biologically-informed prior weights supplied via `--weights`; recommended when the number of candidate roots is small (< 10).
- `entropy`: Entropy Weight Method (EWM) — derives weights automatically from the candidate distribution; assigns higher weight to metrics with greater variance among candidates; recommended when many diverse candidates exist.

**5. Final root selection**

- The candidate with the **highest** composite score is selected as the optimal root.

### Phylo_Rooter
```
Description:
    Roots gene trees using a six-metric composite scoring framework (OD, BLV, GD, SO, GDC, MulRF)
    guided by the species tree. Supports empirical (prior) or entropy (data-driven) weight strategies.
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --input_sps_tree        Species tree file in Newick format
Optional parameter:
    --weights               Weights in fixed order: OD BLV GD SO GD_consistency RF; input exactly six
                            floats with sum = 1; used when --weight-strategy empirical
                            default = 0.30 0.10 0.30 0.10 0.10 0.10
    --weight-strategy       Scoring weight strategy: empirical (default) or entropy
                            empirical: biologically-informed prior weights from --weights;
                                       recommended for small candidate sets (< 10)
                            entropy:   Entropy Weight Method (EWM) -- data-driven; assigns higher
                                       weight to metrics with greater variance among candidates
    --output_dir            Output directory (default: current working directory)
Usage:
    # Default (empirical weights)
    PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk

    # Custom empirical weights (OD=0.20 BLV=0.05 GD=0.40 SO=0.15 GDC=0.10 RF=0.10)
    PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --weights 0.20 0.05 0.40 0.15 0.10 0.10

    # Entropy (data-driven) weighting
    PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --weight-strategy entropy
```
### MulRF_Distance
```
Description:
    To compute species-level MulRF topological conflict distances
    Supports two modes:
    1) mode 1: Gene Tree vs Gene Tree (pairwise within one GF list; self-comparisons are excluded)
    2) mode 2: Gene Tree vs Species Tree
    (this is a topology-distance metric, not a sequence/genetic distance)
    Practical uses:
    1) Quantify gene-tree/species-tree conflict intensity (higher distance = stronger discordance/complex history)
    2) Compare evolutionary pattern similarity among gene families for clustering/filtering/modeling
    3) Support rooting tie-break in Phylo_Rooter by preferring roots with lower normalized MulRF conflict
Required parameter:
    --mode                  Comparison mode: 1=Gene Tree vs Gene Tree, 2=Gene Tree vs Species Tree
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --input_sps_tree        Species tree file in Newick format (required when --mode 2)
    --output                Output TSV filename, default = mulrf_distance.tsv
Usage:
    PhyloTracer MulRF_Distance --mode 1 --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--output mulrf_mode1.tsv]
    PhyloTracer MulRF_Distance --mode 2 --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk [--output mulrf_mode2.tsv]
Output columns (mulrf_distance.tsv):
    tre_id_1                           Tree ID from GF1
    tre_id_2                           Tree ID from GF2
    gene_tree_1_leaf_count             Leaf count of tree 1
    gene_tree_2_leaf_count             Leaf count of tree 2
    gene_tree_1_species_count          Species count represented in tree 1
    gene_tree_2_species_count          Species count represented in tree 2
    shared_species_count               Number of shared species between the two trees
    mulrf_distance                     Raw MulRF topological conflict distance
    maximum_possible_mulrf_distance    Theoretical maximum MulRF distance for this pair
    normalized_mulrf_distance          mulrf_distance / maximum_possible_mulrf_distance
    shared_species_bipartition_count   Shared species-level bipartitions
    gene_tree_1_only_bipartition_count Bipartitions unique to tree 1
    gene_tree_2_only_bipartition_count Bipartitions unique to tree 2
```
### PhyloTree_CollapseExpand
```
Description:
    To transform a phylogenetic tree in Newick format into a 'comb' structure based on a predefined support value threshold. It can also revert this `comb` structure to a fully resolved binary tree, allowing dynamic topology adjustments
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --support_value         Node support cutoff used for collapsing internal branches, default = 50
Optional parameter:
    --revert                If set, expand previously collapsed comb structures back to binary form, default = False
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50 [--revert] [--output_dir DIR]
```
### PhyloSupport_Scaler
```
Description:
    To recalibrate support value from bootstrap or posterior probability in a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational compatibility, and vice versa to meet various analytical needs
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --scale_to              Target support scale: "1" for [0,1], "100" for [0,100]
Optional parameter:
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 100 [--output_dir DIR]
```
### BranchLength_NumericConverter
```
Description:
    To convert branch length values of a phylogenetic tree from string to numerical format
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
Optional parameter:
    --decimal_place         Number of decimal places to keep for branch lengths, default = 10
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap [--decimal_place 10] [--output_dir DIR]
```
### OrthoFilter_LB

**Description:** Prunes phylogenomic noise from both single-copy and multi-copy gene family trees by removing tips with abnormally long branches. This module helps eliminate potential artifacts such as Long Branch Attraction (LBA).

**Required parameters (CLI):** `--input_GF_list`, `--input_imap`, `--rrbr_cutoff`, `--srbr_cutoff`

**1. Root Relative Branch Ratio (RRBR)**

**Concept:** Measures the deviation of a specific gene's **root-to-tip path length** relative to the **global average** of root-to-tip distances in the same gene tree.

* **Purpose:** Detects outlier sequences evolving significantly faster or slower than the family norm.
* **Formula:**

$$
\text{RRBR} = \frac{\text{Root-to-tip distance} - \text{Average root-to-tip distance}}{\text{Average root-to-tip distance}}
$$

**2. Sister Relative Branch Ratio (SRBR)**

**Concept:** Measures the evolutionary distance of a gene relative to its **nearest neighbor** (sister branch).

* **Purpose:** Identifies local branch length asymmetry. A gene significantly longer than its "sister" is a high-risk candidate for phylogenetic noise, even in fast-evolving families.
* **Formula:**

$$
\text{SRBR} = \frac{\text{Branch Length} - \text{Sister Branch Length}}{\text{Sister Branch Length}}
$$

**Where:**

* **Root-to-tip distance:** Path length from the tip to the root.
* **Average root-to-tip distance:** Mean root-to-tip path length across all tips in the same gene tree.
* **Sister Branch Length:** The branch length of the nearest "neighbor" or "sister" gene.

**3. Long-branch decision mode**

**Concept:** Controls how RRBR and SRBR are combined during pruning.

* `or` (lenient): remove a tip when `RRBR >= rrbr_cutoff` **OR** `SRBR >= srbr_cutoff`.
* `and` (strict): remove a tip only when `RRBR >= rrbr_cutoff` **AND** `SRBR >= srbr_cutoff`.

```
Description:
    To prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --rrbr_cutoff           RRBR cutoff based on root-to-tip distance, default = 5
    --srbr_cutoff           SRBR cutoff based on sister-relative branch ratio, default = 2.5
Optional parameter:
    --lb_mode              Long-branch decision mode: or|and, default = or
    --visual                If set, export before/after tree visualization PDFs, default = False
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --rrbr_cutoff 5 --srbr_cutoff 2.5 [--lb_mode or] [--visual] [--output_dir DIR]
```
### OrthoFilter_Mono

**Required parameters (CLI):** `--input_GF_list`, `--input_taxa`, `--input_imap`, `--input_sps_tree`

**Scoring logic (current implementation concept):**

**1. Dominant Lineage Purity**

**Concept:** Measures how strongly a lineage is dominated by target taxa labels.

* **Formula:**

$$
\text{Purity}=\frac{\text{N}_{\text{target}}}{\text{N}_{\text{dominant tips}}}
$$

**Where:**

- $N_{\text{target}}$ = number of target taxa tips  
- $N_{\text{dominant tips}}$ = total tips in dominant lineage 

**2. Phylogenetic Distance Score**

**Concept:** Alien lineages mapped deeper and farther from the target lineage in the species tree are more likely to be removed.

* **Formula:**

$$
\text{PhyloDist}=\text{Depth}(\text{MRCA}(\text{alien}))-\text{Depth}(\text{MRCA}(\text{target}\cup\text{alien}))
$$

**3. Alien Coverage Score**

**Concept:** Alien lineages occupying fewer tips within a dominant lineage are more likely to be noise.

* **Formula:**

$$
\text{AlienCov}=\frac{\text{N}_{\text{alien}}}{\text{N}_{\text{dominant tips}}}
$$

**4. Alien Depth-Variation Score**

**Concept:** Alien lineages inserted more deeply relative to the dominant lineage root are more likely to be removed.

* **Formula:**

$$
\text{AlienDepth}=\text{Depth}(\text{alien})-\text{Depth}\left(\text{MRCA}(\text{dom})\right)
$$

**5. Combined Ranking Score**

**Concept:** Candidates are ranked by a multiplicative score using normalized components.

* **Formula:**
  
$$
\text{Combined}=\text{Norm}(\text{PhyloDist})\times\text{Norm}(\text{AlienDepth})\times\left(-\log_{10}\left(\text{AlienCov} + 10^{-4}\right)\right)
$$

**6. Removal Stopping Rules**

**Concept:** The iterative pruning stops as soon as the dominant-lineage purity reaches `purity_cutoff`, or when the removal cap `max_remove` is reached (including cases where removing the next candidate would exceed the cap).
 
* **Formula:**
  
$$
\text{max remove}=\max\left(\text{max remove fraction}\times\text{N}\_{\text{dominant tips}},1\right)
$$

```
Description:
    To prune phylogenomic noise from both single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing well-established relationships unless otherwise required for specific objectives
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_taxa            Two-column mapping file (gene_id<TAB>clade_or_lineage_label)
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --input_sps_tree        Species tree file in Newick format
Optional parameter:
    --purity_cutoff         Target purity for dominant lineage, default = 0.95
    --max_remove_fraction   Maximum fraction of tips allowed to be removed, default = 0.5
    --visual                If set, export before/after pruning visualization PDFs, default = False
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_taxa Clade.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk [--purity_cutoff 0.95] [--max_remove_fraction 0.5] [--visual] [--output_dir DIR]
```
### TreeTopology_Summarizer
```
Description:
    To enumerate and visualize the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades
    Visualization outputs are vector PDFs: each topology panel is rendered at A4 width and merged in a single-column, top-to-bottom layout (no PNG rasterization).
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --visual_top            Number of top-ranked topologies to visualize, default = 10
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--visual_top 10] [--output_dir DIR]
```
### Tree_Visualizer
```
Description:
    To mark tips of gene trees with provided tags, identify GD nodes, and integrate gene duplication results onto the species tree
    Tree figures are exported as vector PDFs. Category colors are assigned consistently by label across output files.
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --keep_branch           Whether to preserve branch lengths in plotting: 1=yes, 0=no
    --tree_style            Tree layout style: r=rectangular, c=circular
    --gene_categories       One or more two-column files (gene_id<TAB>category_label) for color annotations
                            Format: each line is "gene_id<TAB>label" (no header required)
                            Optional header naming rule: first line can be "gene_id<TAB>HeaderName";
                            HeaderName is used as the displayed column title (recommended: Family, Order, Clade; avoid spaces/special chars)
                            Meaning: each file is one categorical layer (e.g., family/order/clade)
                            Example files in 09_Tree_Visualizer: Family.imap, Order.imap, Clade.imap
                            Note: for species-tree family-duplication mapping, the first file in --gene_categories is used as the family map
    --input_sps_tree        Species tree file in Newick format
    --heatmap_matrix        Gene-associated numeric matrix file (recommended: .txt/.tsv tab-delimited; also supports .csv/.xls/.xlsx), genes as row index
    --visual_gd             If set, overlay predicted GD nodes on gene-tree figures, default = False
    --gd_support            Minimum support of a GD candidate node used by --visual_gd (range: 0-100), default = 50
    --subclade_support      Minimum support required in GD child subclades used by --visual_gd (range: 0-100), default = 0
    --dup_species_proportion  Minimum overlap ratio of duplicated species between GD child clades used by --visual_gd (range: 0-1), default = 0.2
    --dup_species_num       Minimum number of overlapping duplicated species under a GD node used by --visual_gd, default = 2
    --deepvar               Maximum tolerated depth-variance score used by --visual_gd, default = 1
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer Tree_Visualizer --input_GF_list GF_ID2path.visual20.imap --input_imap gene2sps.imap --gene_categories Family.imap Order.imap Clade.imap --input_sps_tree sptree.nwk --heatmap_matrix heatmap_matrix.txt --keep_branch 1 --tree_style r --visual_gd --gd_support 50 --subclade_support 0 --dup_species_proportion 0.2 --dup_species_num 2 --deepvar 1 [--output_dir DIR]
```
### GD_Detector
```
Description:
    To identify gene duplication events by reconciling gene trees with a species tree
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --gd_support            Minimum support of a GD candidate node (accepted range: 0-100; typical: 50-100), default = 50
    --subclade_support      Minimum support required in GD child subclades (accepted range: 0-100), default = 0
    --dup_species_proportion  Minimum overlap ratio of duplicated species between the two GD child clades (range: 0-1), default = 0.2
    --dup_species_num       Minimum number of overlapping duplicated species under a GD node, default = 2
    --input_sps_tree        Species tree file in Newick format
    --deepvar               Maximum tolerated depth-variance score for GD screening, default = 1
Optional parameter:
    --gdtype_mode           GD type assignment mode: relaxed uses species-overlap mapping only; strict additionally enforces depth-consistency filtering with --deepvar, default = relaxed
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0 --dup_species_num 2 --input_sps_tree sptree.nwk --deepvar 1 [--gdtype_mode relaxed] [--output_dir DIR]
```
### GD_Visualizer
```
Description:
    To visualize gene duplication detection results and integrate findings onto the species tree
Required parameter:
    --input_sps_tree        Numbered species tree file in Newick format (use the numbered tree output from GD_Detector, e.g., numed_sptree.nwk)
    --gd_result             GD result table produced by GD_Detector
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer GD_Visualizer --input_sps_tree numed_sptree.nwk --gd_result gd_result_relaxed.txt --input_imap gene2sps.imap [--output_dir DIR]
```
### GD_Loss_Tracker
```
Description:
    To analyze and summarize gene duplication loss events across nodes and tips in the species tree
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_sps_tree        Species tree file in Newick format
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --target_species        Only count loss paths ending in this species (e.g., Arabidopsis_thaliana). Can be used multiple times.
    --mrca_node             Only count loss paths passing through the MRCA of SP1 and SP2. Format: SpeciesA,SpeciesB (comma-separated, no space). Can be used multiple times.
    --include_unobserved_species  Classification policy for species absent from the current gene family, default = False. If set, unobserved species are still assigned 2-2/2-1/2-0 by left/right presence; if not set, they are labeled as missing_data. This affects classification labels only.
    --node_count_mode       Node counting mode for path_count_* transition statistics, default = nonaccumulate. Choices: nonaccumulate|accumulate.
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--target_species Arabidopsis_thaliana] [--mrca_node SpeciesA,SpeciesB] [--include_unobserved_species] [--node_count_mode nonaccumulate] [--output_dir DIR]
```
### GD_Loss_Visualizer
```
Description:
    To visualize the summary of gene duplication loss events on the context of species tree
Required parameter:
    --gd_loss_result        Detailed table generated by GD_Loss_Tracker (gd_loss_summary.txt)
    --input_sps_tree        Numbered species tree file in Newick format
Optional parameter:
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer GD_Loss_Visualizer --input_sps_tree numed_sptree.nwk --gd_loss_result gd_loss_summary.txt [--output_dir DIR]
```
### Ortho_Retriever
```
Description:
    To infer single-copy putative orthologs by splitting paralogs from large-scale gene family trees for multiple species
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --input_gene_length     Two-column mapping file (gene_id<TAB>gene_length)
Optional parameter:
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap [--output_dir DIR]
```
### Hybrid_Tracer
```
Description:
    To detect hybridization signals from GD-derived gene sets and run HyDe testing on species-tree mapped GD nodes
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_Seq_GF_list     Tab-delimited mapping file (GF_ID<TAB>alignment_path), aligned with --input_GF_list IDs
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --input_sps_tree        Species tree file in Newick format
Optional parameter:
    --mrca_node             Restrict Hybrid_Tracer to the MRCA of SP1 and SP2. Format: SpeciesA,SpeciesB (comma-separated, no space). If multiple are provided, only the first valid pair is used.
    --split_groups          Number of partitions for HYDE batch processing, default = 1
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer Hybrid_Tracer --input_GF_list gf.txt --input_Seq_GF_list gf_aln.txt --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--mrca_node SpeciesA,SpeciesB] [--split_groups 2] [--output_dir DIR]
```
### Hybrid_Visualizer
```
Description:
    To visualize hybridization signals, highlighting support from gene tree topologies and D-statistic signals
Required parameter:
    --hyde_out              HYDE output table file
    --input_sps_tree        Species tree file in Newick format
Optional parameter:
    --node                  Use node-mode heatmaps (monophyletic clade stacking) instead of leaf-mode output
    --output_dir            Output directory (default: current working directory)
Usage:
    PhyloTracer Hybrid_Visualizer --hyde_out hyde_out.txt --input_sps_tree sptree.nwk [--node] [--output_dir DIR]
```
### HaploFinder
```
Description:
    To detect haplotype-level GD signals and support split-mode FASTA partitioning
Required parameter:
    --mode                  Run mode: "haplofinder" for GD analysis, "split" for FASTA partitioning by color labels, default = haplofinder
Mode = haplofinder required:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); required in haplofinder mode
    --input_imap            Two-column mapping file (gene_id<TAB>species_name); required in both haplofinder and split modes
    --input_sps_tree        Species tree file in Newick format; required in haplofinder mode
    --species_a             Name of species A (required for haplofinder mode)
    --species_b             Name of species B (required for haplofinder mode)
    --species_a_gff         Genome annotation file for species A in GFF/GTF-compatible format
    --species_b_gff         Genome annotation file for species B in GFF/GTF-compatible format
    --species_a_lens        Chromosome-length file for species A (chr<TAB>length)
    --species_b_lens        Chromosome-length file for species B (chr<TAB>length)
Optional in haplofinder mode:
    --gd_support            Minimum support of GD nodes used for pair extraction (accepted range: 0-100, default = 50)
    --pair_support          Minimum support of ortholog/speciation pair nodes (accepted range: 0-100, default = 50)
    --visual_chr_a          Optional chromosome list file for species A visualization subset
    --visual_chr_b          Optional chromosome list file for species B visualization subset
    --size                  Point size in dotplot rendering (positive float, default = 0.0005)
Optional in both modes:
    --output_dir            Output directory (default: current working directory)
Mode = split required:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); required in haplofinder mode
    --input_imap            Two-column mapping file (gene_id<TAB>species_name); required in both haplofinder and split modes
    --input_fasta           Input FASTA file (.fa/.fasta), required in split mode
    --cluster_file          Split-mode cluster metadata file (legacy compatibility field; currently required by CLI checks)
    --hyb_sps               Hybrid species name used for subgenome assignment (required in split mode)
    --parental_sps          Parental species names used for split-mode assignment; provide as a single quoted, space-separated string
    --species_b_gff         Genome annotation file for species B in GFF/GTF-compatible format
Usage:
    PhyloTracer HaploFinder --mode haplofinder --input_GF_list gf.txt --input_imap gene2sps.imap --input_sps_tree sptree.nwk --species_a arh --species_b ard --species_a_gff arh.gff --species_b_gff ard.gff --species_a_lens arh.lens --species_b_lens ard.lens [--gd_support 50] [--pair_support 50] [--visual_chr_a chr_a.txt --visual_chr_b chr_b.txt --size 0.0001] [--output_dir DIR]
    PhyloTracer HaploFinder --mode split --input_GF_list gf.txt --input_imap gene2sps.imap --input_fasta proteins.fa --cluster_file cluster.tsv --hyb_sps Hybrid --parental_sps "P1 P2" --species_b_gff ard.gff [--output_dir DIR]
```

---
## Bug Reports

You can report bugs or request features through our [GitHub Issues page](https://github.com/YiyongZhao/PhyloTracer/issues). If you have any questions, suggestions, or issues, please do not hesitate to contact us.

## Contributing

If you're interested in contributing code or reporting bugs, we welcome your ideas and contributions to improve PhyloTracer! Please check out [Contribution Guidelines](https://docs.github.com/en/issues).

## Version History

Check the [Changelog](https://github.com/YiyongZhao/PhyloTracer/commits/PhyloTracer_v1.0.0) for details on different versions and updates.

## License

PhyloTracer is licensed under the [MIT LICENSE](LICENSE).
