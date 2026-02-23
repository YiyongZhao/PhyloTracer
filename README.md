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
conda activate PhyloTracer

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
conda activate PhyloTracer

# 4) Run help
PhyloTracer -h

# 5) Example run
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
* Python modules:
  * ete3==3.1.3
  * HyDe (https://github.com/pblischak/HyDe)
  * pandas==2.2.1
  * numpy==1.26.4
  * tqdm==4.66.2
  * pypdf>=3.0.0
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
- HaploFinder (haplofinder mode): `gd_pairs_dotplot.pdf`, `gd_pairs_dotplot.png`, `color_label.txt`
- HaploFinder (split mode): `haplofinder_split/split_assignment.tsv`, `haplofinder_split/split_subgenome_A.fasta`, `haplofinder_split/split_subgenome_B.fasta`, `haplofinder_split/split_subgenome_unknown.fasta`, `haplofinder_split/split_summary.txt`

---
## Command line options

This section follows an OrthoFinder-like CLI reference style with compact layout and explicit parameter meanings.

### PhyloTree_CollapseExpand
```
Description:
    To transform a phylogenetic tree in Newick format into a 'comb' structure based on a predefined support value threshold. It can also revert this `comb` structure to a fully resolved binary tree, allowing dynamic topology adjustments
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --support_value         Node support cutoff used for collapsing internal branches, default=50
Optional parameter:
    --revert                If set, expand previously collapsed comb structures back to binary form, default=False
Usage:
    PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50 [--revert]
```
### PhyloSupport_Scaler
```
Description:
    To recalibrate support value from bootstrap or posterior probability in a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational compatibility, and vice versa to meet various analytical needs
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --scale_to              Target support scale: "1" for [0,1], "100" for [0,100]
Usage:
    PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 1
```
### BranchLength_NumericConverter
```
Description:
    To convert branch length values of a phylogenetic tree from string to numerical format
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
Optional parameter:
    --decimal_place         Number of decimal places to keep for branch lengths, default=10
Usage:
    PhyloTracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap [--decimal_place 10]
```
### Phylo_Rooter
```
Description:
    Enables an accurate method for gene tree rooting and enhancing the downstream evolutionary genomic analysis
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --input_gene_length     Two-column mapping file (gene_id<TAB>gene_length)
    --input_sps_tree        Species tree file in Newick format
Usage:
    PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap --input_sps_tree sptree.nwk
```
### OrthoFilter_LB

**Description:** Prunes phylogenomic noise from both single-copy and multi-copy gene family trees by removing tips with abnormally long branches. This module helps eliminate potential artifacts such as Long Branch Attraction (LBA).

**Required parameters (CLI):** `--input_GF_list`, `--input_imap`, `--absolute_branch_length`, `--relative_branch_length`

**1. Root Relative Branch Ratio (RRBR)**

**Concept:** Measures the deviation of a specific gene's branch length relative to the **global average** of the entire gene family tree.

* **Purpose:** Detects outlier sequences evolving significantly faster or slower than the family norm.
* **Formula:**

$$
\text{RRBR} = \frac{\text{Branch Length} - \text{Average Branch Length}}{\text{Average Branch Length}}
$$

**2. Sister Relative Branch Ratio (SRBR)**

**Concept:** Measures the evolutionary distance of a gene relative to its **nearest neighbor** (sister branch).

* **Purpose:** Identifies local branch length asymmetry. A gene significantly longer than its "sister" is a high-risk candidate for phylogenetic noise, even in fast-evolving families.
* **Formula:**

$$
\text{SRBR} = \frac{\text{Branch Length} - \text{Sister Branch Length}}{\text{Sister Branch Length}}
$$

**Where:**

* **Branch Length:** The branch length of the specified gene.
* **Average Branch Length:** The arithmetic mean of all branch lengths in the gene family tree.
* **Sister Branch Length:** The branch length of the nearest "neighbor" or "sister" gene.

```
Description:
    To prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --absolute_branch_length  Absolute branch-length multiplier threshold (integer), default=5
    --relative_branch_length  Relative branch-length multiplier threshold (float), default=2.5
Optional parameter:
    --visual                If set, export before/after tree visualization PDFs, default=False
Usage:
    PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --absolute_branch_length 5 --relative_branch_length 2.5 [--visual]
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
    --purity_cutoff         Target purity for dominant lineage, default=0.95
    --max_remove_fraction   Maximum fraction of tips allowed to be removed, default=0.5
    --visual                If set, export before/after pruning visualization PDFs, default=False
Usage:
    PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_taxa gene2clade.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk [--purity_cutoff 0.95 --max_remove_fraction 0.5 --visual]
```
### TreeTopology_Summarizer
```
Description:
    To enumerate and visualize the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --visual_top            Number of top-ranked topologies to visualize, default=10
Usage:
    PhyloTracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--visual_top 10]
```
### Tree_Visualizer
```
Description:
    To mark tips of gene trees with provided tags, identify GD nodes, and integrate gene duplication results onto the species tree
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
Optional parameter:
    --keep_branch           Whether to preserve branch lengths in plotting: 1=yes, 0=no
    --tree_style            Tree layout style: r=rectangular, c=circular
    --gene_categories       One or more two-column files (gene_id<TAB>category_label) for color annotations
    --gene_family           Two-column mapping file (gene_id<TAB>family_label)
    --input_sps_tree        Species tree file in Newick format (required with --gene_family)
    --gene_expression       Gene expression matrix file (.csv/.xls/.xlsx), genes as row index
    --visual_gd             If set, overlay predicted GD nodes on gene-tree figures, default=False
Usage:
    PhyloTracer Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--gene_categories gene2order.imap gene2taxa.imap gene2clade.imap --keep_branch {1,0} --tree_style {r,c} --gene_family gene2family.imap --input_sps_tree sptree.nwk --gene_expression gene2expression.csv --visual_gd]
```
### GD_Detector
```
Description:
    To identify gene duplication events by reconciling gene trees with a species tree
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --gd_support            Minimum support of a GD candidate node (accepted range: 0-100; typical: 50-100), default=50
    --subclade_support      Minimum support required in GD child subclades (accepted range: 0-100), default=0
    --dup_species_proportion  Minimum overlap ratio of duplicated species between the two GD child clades (range: 0-1), default=0.2
    --dup_species_num       Minimum number of overlapping duplicated species under a GD node, default=2
    --input_sps_tree        Species tree file in Newick format
    --deepvar               Maximum tolerated depth-variance score for GD screening, default=1
Optional parameter:
    --gdtype_mode           GD type mode: relaxed (species overlap only) or strict (overlap + depth constraint), default=relaxed
Usage:
    PhyloTracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0 --dup_species_num 2 --input_sps_tree sptree.nwk --deepvar 1 [--gdtype_mode relaxed]
```
### GD_Visualizer
```
Description:
    To visualize gene duplication detection results and integrate findings onto the species tree
Required parameter:
    --input_sps_tree        Numbered species tree file in Newick format
    --gd_result             GD result table produced by GD_Detector
    --input_imap            Two-column mapping file (gene_id<TAB>species_name), default=None
Usage:
    PhyloTracer GD_Visualizer --input_sps_tree sptree.nwk --gd_result gd_result.txt --input_imap gene2sps.imap
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
    --include_unobserved_species  If set, species unobserved in a gene family are still classified by left/right presence instead of labeled as missing_data.
Usage:
    PhyloTracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--target_species Arabidopsis_thaliana] [--mrca_node SpeciesA,SpeciesB] [--include_unobserved_species]
```
### GD_Loss_Visualizer
```
Description:
    To visualize the summary of gene duplication loss events on the context of species tree
Required parameter:
    --gd_loss_result        Detailed table generated by GD_Loss_Tracker (gd_loss_summary.txt)
    --input_sps_tree        Numbered species tree file in Newick format
Usage:
    PhyloTracer GD_Loss_Visualizer --input_sps_tree numbered_species_tree.nwk --gd_loss_result gd_loss_summary.txt
```
### Ortho_Retriever
```
Description:
    To infer single-copy putative orthologs by splitting paralogs from large-scale gene family trees for multiple species
Required parameter:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line
    --input_imap            Two-column mapping file (gene_id<TAB>species_name)
    --input_gene_length     Two-column mapping file (gene_id<TAB>gene_length)
Usage:
    PhyloTracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap
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
    --split_groups          Number of partitions for HYDE batch processing, default=1
Usage:
    PhyloTracer Hybrid_Tracer --input_GF_list GF_ID2path.imap --input_Seq_GF_list Seq_GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--mrca_node SpeciesA,SpeciesB --split_groups 2]
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
Usage:
    PhyloTracer Hybrid_Visualizer --hyde_out hyde.out --input_sps_tree sptree.nwk [--node]
```
### HaploFinder
```
Description:
    To detect haplotype-level GD signals and support split-mode FASTA partitioning
Required parameter:
    --mode                  Run mode: "haplofinder" for GD analysis, "split" for FASTA partitioning by color labels, default=haplofinder
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
    --gd_support            Minimum support of GD nodes used for pair extraction (accepted range: 0-100, default=50)
    --pair_support          Minimum support of ortholog/speciation pair nodes (accepted range: 0-100, default=50)
    --visual_chr_a          Optional chromosome list file for species A visualization subset
    --visual_chr_b          Optional chromosome list file for species B visualization subset
    --size                  Point size in dotplot rendering (positive float, default=0.0005)
Mode = split required:
    --input_GF_list         Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); required in haplofinder mode
    --input_imap            Two-column mapping file (gene_id<TAB>species_name); required in both haplofinder and split modes
    --input_fasta           Input FASTA file (.fa/.fasta), required in split mode
    --cluster_file          Split-mode cluster metadata file (legacy compatibility field; currently required by CLI checks)
    --hyb_sps               Hybrid species name used for subgenome assignment (required in split mode)
    --parental_sps          Parental species names used for split-mode assignment; provide as a single quoted, space-separated string
    --species_b_gff         Genome annotation file for species B in GFF/GTF-compatible format
Usage:
    PhyloTracer HaploFinder --mode haplofinder --input_GF_list GF.list --input_imap gene2sps.imap --input_sps_tree sptree.nwk --species_a A --species_b B --species_a_gff A.gff --species_b_gff B.gff --species_a_lens A.lens --species_b_lens B.lens --gd_support 50 --pair_support 50 [--visual_chr_a chr_a.txt --visual_chr_b chr_b.txt --size 0.0001]
    PhyloTracer HaploFinder --mode split --input_GF_list GF.list --input_imap gene2sps.imap --input_fasta proteins.fa --cluster_file cluster.txt --hyb_sps Hybrid --parental_sps "P1 P2" --species_b_gff B.gff
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
