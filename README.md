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
    Licence: MIT license                                                                     
    Release Date: 2023-7                                                                     
    Contacts: Taoli(l948777439@163.com); Yiyong Zhao(yiyongzhao1991@gmail.com)
                                                                         
###############################################################################################
```
![Version](https://img.shields.io/badge/Version-1.0.0-blue)
[![Documentation Status](http://readthedocs.org/projects/hybridization-detection/badge/?version=latest)](http://hybridization-detection.readthedocs.io)
[![PhyloTracer Issues](https://img.shields.io/badge/PhyloTracer-Issues-blue.svg)](https://github.com/YiyongZhao/PhyloTracer/issues)
![Build Status](https://travis-ci.org/YiyongZhao/PhyloTracer.svg?branch=master)
[![PyPI](https://img.shields.io/pypi/v/PhyloTracer.svg)](https://pypi.python.org/pypi/PhyloTracer)

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

PhyloTracer integrates 16 modular tools covering phylogenetic preprocessing, rooting, orthology refinement, duplication/loss detection, and hybridization analysis. Each module can run independently or be incorporated into larger evolutionary pipelines.

1. **PhyloTree_CollapseExpand:** Transforms a phylogenetic tree into a “comb-like” structure based on a predefined support threshold, and re-expands it back to binary form when needed.
2. **PhyloSupport_Scaler:** Recalibrates branch support values to standardized scales ([0–1] or [1–100]) for consistent computation.
3. **BranchLength_NumericConverter:** Converts branch-length strings to numeric format for downstream quantitative analyses.
4. **Phylo_Rooter:** Implements an accurate, automated rooting algorithm for gene trees to enhance evolutionary inference.
5. **OrthoFilter_LB:** Removes tips with excessively long branches to eliminate phylogenomic noise.
6. **OrthoFilter_Mono:** Prunes non-monophyletic outliers and paralogs under user-defined taxonomic constraints.
7. **TreeTopology_Summarizer:** Summarizes frequencies of absolute and relative topologies across gene trees or predefined clades.
8. **Tree_Visualizer:** Visualizes duplication events, node labels, and expression profiles on gene and species trees.
9. **GD_Detector:** Identifies gene duplication events via reconciliation between gene and species trees.
10. **GD_Visualizer:** Displays detected duplication nodes in a species tree context.
11. **GD_Loss_Tracker:** Tracks duplication-loss patterns following major GD bursts across species tree nodes.
12. **GD_Loss_Visualizer:** Visualizes node/tip-specific gene loss summaries.
13. **Ortho_Retriever:** Extracts putative single-copy orthologs by recursively splitting paralogous clades.
14. **Hybrid_Tracer:** Detects hybridization signals from duplicated genes using coalescent-based phylogenetic invariants.
15. **Hybrid_Visualizer:** Highlights hybridization proportions (γ) and support values across the species tree.
16. **HaploFinder:** Identifies ancient recombination (conversion/crossover) events by tracing subgenome haplotypes.

*Together, these modules provide a comprehensive workflow for constructing, refining, and interpreting large-scale phylogenomic data.*

---
## Getting started with PhyloTracer

### Clone and install environment

```bash
#A convenient one-click installation by using conda (https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) with the following commands:
git clone https://github.com/YiyongZhao/PhyloTracer.git
cd PhyloTracer
conda env create -f environment.yml
conda activate phylotracer

#Alternatively, a convenient one-click installation by using pip (the package installer for Python) with the following commands:
chmod +x install_packages.sh
bash install_packages.sh

#Reminder for potential visualization issues: qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found and this application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
#Alternative available platform plugins include: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-#egl, wayland-xcomposite-glx, webgl, xcb. before running PhyloTracer, please execute the following bash command:
export QT_QPA_PLATFORM=offscreen
```
### Install from PyPI with pip

```bash
pip install PhyloTracer
```

### Quick start from GitHub ZIP (download + extract)

```bash
# 1) Download and unzip PhyloTracer-main.zip from GitHub
# 2) Enter project directory
cd PhyloTracer-main

# 3) Create environment
conda env create -f environment.yml
conda activate phylotracer

# 4) Run help
phylotracer -h

# 5) Example run
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
---
## Features

- Rooting (core, via `Phylo_Rooter`):  
  Uses a two-stage strategy:
  1) **Stage 1 (fast screening, no RF)** computes `deep`, `var`, `GD`, and `species_overlap`;  
  2) **Stage 2 (fine screening)** computes `RF` only for top candidates, then selects the best root by sorting `RF` first, then Stage-1 score.

  Stage-1 scoring in code:
  $$
  score = w_{deep}\cdot norm(deep)+w_{var}\cdot norm(var)+w_{GD}\cdot norm(GD)-w_{SO}\cdot norm(species\_overlap)
  $$
  where `norm(x)` is min-max normalization in \([0,1]\).

  <details>
  <summary>📊 Stage-1 weights in <code>Phylo_Rooter</code> (from code)</summary>

  | Metric | Multi-copy GFs | Single-copy GFs (Stage-1) |
  |---|---:|---:|
  | `deep` | 0.30 | 0.70 |
  | `var` | 0.10 | 0.30 |
  | `GD` | 0.50 | 0.00 |
  | `species_overlap` (subtracted) | 0.10 | 0.00 |
  | `RF` | not used in Stage-1 | not used in Stage-1 |

  </details>

  Stage-2 selection in code:
  - keep top `max(20, ceil(0.8 * total_candidates))` candidates by Stage-1 score;
  - compute `RF` for them;
  - choose best by `sort_values(by=["RF", "score"], ascending=[True, True])`.

---

- Topology statistics (via `TreeTopology_Summarizer`):  
  Computes the absolute and relative topology frequencies for single-copy gene trees.  
  Supports grouped summarization by user-provided labels (e.g., family/order tags) when supplied.

---

- Hybridization screening (via `Hybrid_Tracer`):  
  `Hybrid_Tracer` extends conventional HyDe-style hybridization detection by leveraging gene duplication (GD)–based signal extraction, offering both grouped and ungrouped analytical modes for cleaner and more node-specific hybridization inference.

  Two complementary strategies are implemented:

  - Ungrouped mode (concatenation-based):  
    For a specific ancestral node, `Hybrid_Tracer` concatenates the alignment sequences from all duplicated genes descending from that node’s GD events, forming a targeted alignment matrix.  
    This enables direct HyDe-like inference of hybridization signals while minimizing unrelated noise, since only genes phylogenetically tied to that duplication origin are included.

  - Grouped mode (signal integration):  
    Alternatively, GD events can be partitioned into multiple groups according to their gene tree topology or taxonomic context.  
    Each group is analyzed independently to estimate its hybridization proportion (γ) and support.  
    These group-level results are then integrated (e.g., averaged or weighted) to infer the overall hybridization signal of that ancestral node.

  Compared with the traditional HyDe pipeline that concatenates all single-copy genes across the genome, `Hybrid_Tracer` focuses exclusively on GD-derived gene sets related to the evolutionary node of interest.  
  This design produces cleaner, more interpretable, and more localized hybridization signals, reducing interference from unrelated loci and improving the biological relevance of γ estimates.

---

- GD & loss profiling (via `GD_Detector` and `GD_Loss_Tracker`):  
  Reconciles gene–species trees to summarize gene duplication events and lineage-specific loss patterns.  
  Paired visualizers (`GD_Visualizer`, `GD_Loss_Visualizer`) assist in comparative interpretation.

---

- HaploFinder:  
  Detects ancient genome recombination signals, including gene conversions and crossover events, by tracing subgenome-specific haplotypes through phylogenomic profiling.  
  This module helps characterize the historical exchange of genetic material between subgenomes and provides insights into genome evolution following duplication or hybridization.


---
## Advanced installation notes

### Required dependencies:

* Python 3.0+ 
* Python modules:
  * ete3==3.1.3
  * HyDe (https://github.com/pblischak/HyDe)
  * pandas==2.2.1
  * numpy==1.26.4
  * tqdm==4.66.2
  * pypdf4==1.27.0
  * matplotlib==3.8.3
  * pyqt5==5.15.10
  * pyqt5-qt5==5.15.12
  * pyqt5-sip==12.13.0
  * pillow==10.2.0
  * fonttools==4.49.0
  * cycler==0.12.1
  * kiwisolver==1.4.5
  * packaging==23.2
  * pyparsing==3.1.1
  * python-dateutil==2.8.2
  * pytz==2024.1
  * six==1.16.0
  * tzdata==2024.1
    
Note: PhyloTracer uses basic functions of analysis and visualization of trees from Python framework [ete3](http://etetoolkit.org/) and detects species hybridization signals using ABAB-BABA test by [HyDe](https://github.com/pblischak/HyDe).

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
------------gene2family.imap------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborellaceae
ATCG00500.1           Brassicaceae
Glyma.07G273800.2     Fabaceae

Provide a two-column file in TSV format: each line contains <gene_id><TAB><plant_order>
------------gene2order.imap-------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborellales
ATCG00500.1           Brassicales
Glyma.07G273800.2     Fabales

Provide a two-column file in TSV format: each line contains <gene_id><TAB><higher-level_taxa>
------------gene2taxa.imap--------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Angiosperm
ATCG00500.1           Malvids
Glyma.07G273800.2     Fabids

Provide a two-column file in TSV format: each line contains <gene_id><TAB><functional_clade>
------------gene2clade.imap-------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Nitrogen-fixing
ATCG00500.1           Nitrogen-fixing
Glyma.07G273800.2     non-Nitrogen-fixing

Provide a two-column file in TSV format: each line contains <gene_id><TAB><expression_values>
------------gene2expression.imap--------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  5.0
ATCG00500.1           12.0
Glyma.07G273800.2     0.0

#Note: You can add any number of imap files. They will sequentially provide annotations to the right of the gene tips according to the order of input.
```
---
## PhyloTracer Results Files

Most modules generate task-specific outputs in either the current working directory or module-specific subdirectories. Common outputs include:

- Rooted tree outputs: `rooted_trees/`
- Long-branch filter outputs: `orthofilter_lb/`
- Monophyly filter outputs: `orthofilter_mono/`
- GD detection tables: `gd_result_*.txt`, `gd_type_*.tsv`
- GD-loss tables: `gd_loss_summary.txt`, `gd_loss_count_summary.txt`, `gd_loss.xlsx`
- Hybridization outputs: `hyde_out.txt`, `hyde_filtered_out.txt`
- Ortholog retrieval outputs: `ortho_retriever_summary.txt`, `ortholog_trees.tsv`
- Topology summaries: `absolute_*.txt`, `relative_*.txt`, merged PNG summaries

---
## Command line options

Compact reference format: each module uses one parameter table (`Parameter | Required | Description`) plus Input/Output/Usage.

### 1) PhyloTree_CollapseExpand
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path` (one tree per line). |
| `--support_value` | Yes | Support cutoff (`0-100`); nodes with support <= cutoff are collapsed. |
| `--revert` | No | Re-expand collapsed comb-like structures to binary form. |
Input: Gene trees from `--input_GF_list`.  
Output: `collapse_expand_tree/*.nwk`.  
Usage: `phylotracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50 [--revert]`

### 2) PhyloSupport_Scaler
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--scale_to` | Yes | Target support scale: `1` => `[0,1]`, `100` => `[0,100]`. |
Input: Gene trees with support values.  
Output: `support_scaler_tree/*.nwk`.  
Usage: `phylotracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 100`

### 3) BranchLength_NumericConverter
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--decimal_place` | No | Branch-length decimal precision (default `10`). |
Input: Gene trees.  
Output: `converter_tree/*.nwk`.  
Usage: `phylotracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap [--decimal_place 10]`

### 4) Phylo_Rooter
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--input_gene_length` | Yes | TSV `gene_id<TAB>gene_length`. |
| `--input_sps_tree` | Yes | Species tree in Newick format. |
Input: Gene trees + species map + gene length map + species tree.  
Output: `rooted_trees/*.nwk`, `stat_matrix.csv`.  
Usage: `phylotracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap --input_sps_tree sptree.nwk`

### 5) OrthoFilter_LB
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--absolute_branch_length` | Yes | Global long-branch multiplier threshold (integer). |
| `--relative_branch_length` | Yes | Sister-relative branch multiplier threshold (float). |
| `--visual` | No | Export before/after PDFs for each tree. |
Input: Gene trees + species map.  
Output: `orthofilter_lb/pruned_tree/*.nwk`, `orthofilter_lb/long_branch_gene/*_delete_gene.txt`, optional `orthofilter_lb/pruned_tree_pdf/*.pdf`.  
Usage: `phylotracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --absolute_branch_length 5 --relative_branch_length 2.5 [--visual]`

### 6) OrthoFilter_Mono
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_taxa` | Yes | TSV `gene_id<TAB>clade_label` for monophyly-guided pruning. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--input_sps_tree` | Yes | Species tree in Newick format. |
| `--purity_cutoff` | No | Target dominant-lineage purity (default `0.95`). |
| `--max_remove_fraction` | No | Maximum removable-tip fraction per dominant lineage (default `0.5`). |
| `--visual` | No | Export before/after visualization PDFs. |
Input: Gene trees + taxa labels + species map + species tree.  
Output: `orthofilter_mono/pruned_tree/*.nwk`, `orthofilter_mono/insert_gene/*_insert_gene.txt`, optional `orthofilter_mono/visual/*.pdf`.  
Usage: `phylotracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_taxa gene2clade.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk [--purity_cutoff 0.95 --max_remove_fraction 0.5 --visual]`

### 7) TreeTopology_Summarizer
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name` for label normalization. |
| `--visual_top` | No | Number of top topologies to visualize (default `10`). |
Input: Gene trees + species map.  
Output: `absolute_*.txt`, `relative_*.txt`, `merge_absolutely_top*.png`, `merge_relative_top*.png`.  
Usage: `phylotracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--visual_top 10]`

### 8) Tree_Visualizer
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--gene_categories` | No | One or more TSV files `gene_id<TAB>category`. |
| `--keep_branch` | No | `1` keep branch lengths, `0` normalize for display. |
| `--tree_style` | No | Tree style: `r` rectangular, `c` circular (default `r`). |
| `--gene_family` | No | TSV `gene_id<TAB>family_label` for family overlays. |
| `--input_sps_tree` | Cond. | Required when `--gene_family` is provided. |
| `--gene_expression` | No | Expression matrix (`.csv/.xls/.xlsx`), genes as row index. |
| `--visual_gd` | No | Overlay GD nodes on tree figures. |
Input: Gene trees + species map + optional metadata files.  
Output: `tree_visualizer/*.pdf`, optional `genefamily_map2_sptree.pdf`.  
Usage: `phylotracer Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--tree_style r] [--visual_gd]`

### 9) GD_Detector
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--gd_support` | Yes | Minimum support for GD candidate nodes (`0-100`). |
| `--subclade_support` | Yes | Minimum support in both GD child clades (`0-100`). |
| `--dup_species_proportion` | Yes | Minimum duplicated-species overlap ratio (`0-1`). |
| `--dup_species_num` | Yes | Minimum duplicated-species overlap count (`>=1`). |
| `--input_sps_tree` | Yes | Species tree Newick file. |
| `--deepvar` | Yes | Maximum tolerated species-depth variance (`>=0`). |
| `--gdtype_mode` | No | GD typing mode: `relaxed` (default) or `strict`. |
Input: Gene trees + species map + species tree.  
Output: `gd_result_relaxed.txt`/`gd_result_strict.txt`, `gd_type_relaxed.tsv`/`gd_type_strict.tsv`.  
Usage: `phylotracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0 --dup_species_num 2 --input_sps_tree sptree.nwk --deepvar 1 [--gdtype_mode relaxed]`

### 10) GD_Visualizer
| Parameter | Required | Description |
|---|---|---|
| `--input_sps_tree` | Yes | Numbered species tree in Newick format. |
| `--gd_result` | Yes | GD result table from `GD_Detector`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
Input: Numbered species tree + GD table + species map.  
Output: `phylotracer_gd_visualizer.pdf`.  
Usage: `phylotracer GD_Visualizer --input_sps_tree numed_sptree.nwk --gd_result gd_result_relaxed.txt --input_imap gene2sps.imap`

### 11) GD_Loss_Tracker
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_sps_tree` | Yes | Species tree in Newick format. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--target_species` | No | Restrict loss-path endpoints; repeatable argument. |
| `--mrca_node` | No | Restrict GD nodes by MRCA pair (`SP1,SP2`); repeatable. |
| `--include_unobserved_species` | No | Include family-unobserved species in classification output. |
Input: Gene trees + species tree + species map.  
Output: `gd_loss_summary.txt`, `gd_loss_count_summary.txt`, `gd_loss.xlsx`.  
Usage: `phylotracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap`

### 12) GD_Loss_Visualizer
| Parameter | Required | Description |
|---|---|---|
| `--gd_loss_result` | Yes | GD-loss summary table (recommended `gd_loss_summary.txt`). |
| `--input_sps_tree` | No | Numbered species tree in Newick format. |
Input: GD-loss table (+ species tree).  
Output: `gd_loss_pie_visualizer.PDF`.  
Usage: `phylotracer GD_Loss_Visualizer --gd_loss_result gd_loss_summary.txt --input_sps_tree numed_sptree.nwk`

### 13) Ortho_Retriever
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--input_gene_length` | Yes | TSV `gene_id<TAB>gene_length`. |
Input: Rooted gene trees + species map + gene lengths.  
Output: `ortho_retriever_summary.txt`, `ortholog_trees.tsv`.  
Usage: `phylotracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap`

### 14) Hybrid_Tracer
| Parameter | Required | Description |
|---|---|---|
| `--input_GF_list` | Yes | TSV file `GF_ID<TAB>gene_tree_path`. |
| `--input_Seq_GF_list` | Yes | TSV file `GF_ID<TAB>alignment_path` aligned to GF IDs. |
| `--input_sps_tree` | Yes | Species tree Newick file. |
| `--input_imap` | Yes | TSV `gene_id<TAB>species_name`. |
| `--mrca_node` | No | Restrict analysis to MRCA of `SpeciesA,SpeciesB`. |
| `--split_groups` | No | Number of GD groups for batched processing (default `1`). |
Input: Gene trees + alignments + species tree + species map.  
Output: `hyde_out.txt`, `hyde_filtered_out.txt`.  
Usage: `phylotracer Hybrid_Tracer --input_GF_list GF_ID2path.imap --input_Seq_GF_list Seq_GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--mrca_node SpeciesA,SpeciesB --split_groups 2]`

### 15) Hybrid_Visualizer
| Parameter | Required | Description |
|---|---|---|
| `--hyde_out` | Yes | HyDe result table from `Hybrid_Tracer`. |
| `--input_sps_tree` | Yes | Species tree in Newick format. |
| `--node` | No | Node-mode heatmap stacking; default is leaf-mode view. |
Input: HyDe output + species tree.  
Output: files matching `*_img_faces.png` and `*_hotmap.png`.  
Usage: `phylotracer Hybrid_Visualizer --hyde_out hyde_out.txt --input_sps_tree sptree.nwk [--node]`

### 16) HaploFinder
| Parameter | Required | Description |
|---|---|---|
| `--mode` | Yes | Run mode: `haplofinder` or `split`. |
| `--input_GF_list` | Cond. | Required in both modes; TSV `GF_ID<TAB>gene_tree_path`. |
| `--input_imap` | Cond. | Required in both modes; TSV `gene_id<TAB>species_name`. |
| `--input_sps_tree` | Cond. | Required in `haplofinder` mode. |
| `--species_a`, `--species_b` | Cond. | Required in `haplofinder` mode. |
| `--species_a_gff`, `--species_b_gff` | Cond. | Required in `haplofinder` mode (`species_b_gff` also required in `split` checks). |
| `--species_a_lens`, `--species_b_lens` | Cond. | Required in `haplofinder` mode. |
| `--gd_support` | No | GD support cutoff (`0-100`, default `50`). |
| `--pair_support` | No | Pair-support cutoff (`0-100`, default `50`). |
| `--visual_chr_a`, `--visual_chr_b` | No | Chromosome subset files for plotting. |
| `--size` | No | Dotplot point size (default `0.0005`). |
| `--input_fasta` | Cond. | Required in `split` mode. |
| `--cluster_file` | Cond. | Required in current `split` CLI checks. |
| `--hyb_sps` | Cond. | Required in `split` mode. |
| `--parental_sps` | Cond. | Required in `split` mode; quoted space-separated list. |
Input: Mode-specific trees/mappings/GFF/lens/FASTA files.  
Output: `*_dotplot.pdf`, `*_dotplot.png`, `gene_conversion_*.txt`, `color_label.txt`, and split FASTA outputs.  
Usage: `phylotracer HaploFinder --mode haplofinder ...` or `phylotracer HaploFinder --mode split ...`

---
## Bug Reports

You can report bugs or request features through our [GitHub Issues page](https://github.com/YiyongZhao/PhyloTracer/issues). If you have any questions, suggestions, or issues, please do not hesitate to contact us.

## Contributing

If you're interested in contributing code or reporting bugs, we welcome your ideas and contributions to improve PhyloTracer! Please check out [Contribution Guidelines](https://docs.github.com/en/issues).

## Version History

Check the [Changelog](https://github.com/YiyongZhao/PhyloTracer/commits/PhyloTracer_v1.0.0) for details on different versions and updates.

## License

PhyloTracer is licensed under the [MIT LICENSE](LICENSE).
