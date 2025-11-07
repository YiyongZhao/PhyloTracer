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
## Introduction

`PhyloTracer` provides a reproducible workflow centered on accurate gene tree rooting and topology statistics. It further offers utilities for screening hybridization-like signals (via topology patterns such as ABAB/ABBA variants), summarizing GD and loss patterns, and subgenome-informed splitting of multi-copy families.
All modules are designed to be used independently or combined in larger phylogenomic pipelines. Where applicable, methods are documented with input assumptions and recommended validation steps to ensure rigorous interpretation.

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
### Clone and install environment:

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
---
### Install from PyPI with pip:

```bash
pip install PhyloTracer
```
---
## Features

- **Rooting (core, via `Phylo_Rooter`):**  
  Uses a **composite rooting score** combining five metrics, with **different weights** for multi-copy vs. single-copy gene families.  
  Metrics include: outgroup depth (OD), ingroup branch-length variance (BLV), Robinson–Foulds distance (RF), GD events count (GD), and GD-clade species overlap (SO).

  <details>
  <summary>📊 Show metric weights for <code>Phylo_Rooter</code></summary>

  | Metric | Multi-copy GFs | Single-copy GFs |
  |---|---:|---:|
  | Outgroup depth (OD) | 0.10 | 0.50 |
  | Branch length variance (BLV) | 0.05 | 0.10 |
  | Robinson–Foulds distance (RF) | 0.60 | 0.40 |
  | GD events count (GD) | 0.20 | 0.00 |
  | GD clade species overlap (SO) | 0.05 | 0.00 |

  </details>

---

- **Topology statistics (via `TreeTopology_Summarizer`):**  
  Computes the **absolute** and **relative** topology frequencies **for single-copy gene trees**.  
  Supports **grouped summarization** by user-provided **labels** (e.g., family/order tags) when supplied.

---

- **Hybridization screening (via `Hybrid_Tracer`):**  
  `Hybrid_Tracer` extends conventional HyDe-style hybridization detection by leveraging **gene duplication (GD)–based signal extraction**, offering both **grouped** and **ungrouped** analytical modes for cleaner and more node-specific hybridization inference.

  **Two complementary strategies are implemented:**

  - **Ungrouped mode (concatenation-based):**  
    For a specific ancestral node, `Hybrid_Tracer` concatenates the alignment sequences from **all duplicated genes descending from that node’s GD events**, forming a targeted alignment matrix.  
    This enables **direct HyDe-like inference** of hybridization signals while minimizing unrelated noise, since only genes phylogenetically tied to that duplication origin are included.

  - **Grouped mode (signal integration):**  
    Alternatively, GD events can be **partitioned into multiple groups** according to their gene tree topology or taxonomic context.  
    Each group is analyzed independently to estimate its **hybridization proportion (γ)** and support.  
    These group-level results are then **integrated (e.g., averaged or weighted)** to infer the overall hybridization signal of that ancestral node.

  Compared with the traditional HyDe pipeline that concatenates all single-copy genes across the genome, `Hybrid_Tracer` focuses exclusively on **GD-derived gene sets** related to the evolutionary node of interest.  
  This design produces **cleaner, more interpretable, and more localized hybridization signals**, reducing interference from unrelated loci and improving the biological relevance of γ estimates.

---

- **GD & loss profiling (via `GD_Detector` and `GD_Loss_Tracker`):**  
  Reconciles gene–species trees to summarize **gene duplication events** and **lineage-specific loss** patterns.  
  Paired visualizers (`GD_Visualizer`, `GD_Loss_Visualizer`) assist in comparative interpretation.

---

- **HaploFinder:**  
  Detects **ancient genome recombination signals**, including **gene conversions** and **crossover events**, by tracing **subgenome-specific haplotypes**    through phylogenomic profiling.  
  This module helps characterize the historical exchange of genetic material between subgenomes and provides insights into genome evolution following duplication or hybridization.



  
---
## Installation

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
Should sorted by second column with the recommended bash command: sort -k2,2 gene2sps.imap
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
## Usage
### PhyloTree_CollapseExpand
```
Description:
    To transform a phylogenetic tree in Newick format into a 'comb' structure based on a predefined support value threshold. It can also revert this `comb` structure to a fully resolved binary tree, allowing dynamic topology adjustments
Required parameter:
    --input_GF_list  File containing paths to gene tree files, one per line
    --support_value  Nodes whose support is less than or equal to 'support_value' will be converted and default=50
Optional parameter:
    --revert         Revert this 'comb' structure to a fully resolved binary tree
Usage:
    Phylo_Tracer.py PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50 [--revert]
```
### PhyloSupport_Scaler
```
Description:
    To recalibrate support value from bootstrap or posterior probability in a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational compatibility, and vice versa to meet various analytical needs
Required parameter:
    --input_GF_list  File containing paths to gene tree files, one per line
    --scale_to       Input '1' to scale support values from 1-100 to 0-1, or '100' to scale from 0-1 to 1-100
Usage:
    Phylo_Tracer.py PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 1
```
### BranchLength_NumericConverter
```
Description:
    To convert branch length values of a phylogenetic tree from string to numerical format
Required parameter:
    --input_GF_list  File containing paths to gene tree files, one per line
Optional parameter:
    --decimal_place  Return the branch length values to 10 decimal places and default = 10
Usage:
    Phylo_Tracer.py BranchLength_NumericConverter --input_GF_list GF_ID2path.imap [--decimal_place 10]
```
### Phylo_Rooter
```
Description:
    Enables an accurate method for gene tree rooting and enhancing the downstream evolutionary genomic analysis
Required parameter:
    --input_GF_list      File containing paths to gene tree files, one per line
    --input_imap         File with classification information of species corresponding to genes
    --input_gene_length  File with information corresponding to gene lengths
    --input_sps_tree     A species tree file with Newick format
Usage:
    Phylo_Tracer.py Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap --input_sps_tree sptree.nwk 
```
### OrthoFilter_LB
Root Relative Branch Ratio (RRBR): The relative branch length ratio based on the root. It represents the ratio between the branch length of a gene and the average branch length of all genes in the gene family tree.
$$
\text{RRBR} = \frac{\text{Branch Length} - \text{Average Branch Length}}{\text{Average Branch Length}}
$$

Sister Relative Branch Ratio (SRBR): The relative branch length ratio based on the nearest neighbor or "sister" gene. It represents the ratio between the branch length of a gene and the branch length of its nearest neighbor gene.

$$
\text{SRBR} = \frac{\text{Branch Length} - \text{Sister Branch Length}}{\text{Sister Branch Length}}
$$

Where:
- **Branch Length** : is the branch length of the specified gene.
- **Average Branch Length** : is the average branch length of all genes in the gene family tree.
- **Sister Branch Length** : is the branch length of the nearest "neighbor" or "sister" gene of the specified gene.

```
Description:
    To prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length
Required parameter:
    --input_GF_list             File containing paths to gene tree files, one per line
    --input_imap                File with classification information of species corresponding to genes
    --absolute_branch_length    Absolute branch length multiplier and default = 5
    --relative_branch_length    Relative branch length multiplier and default = 2.5
Optional parameter:
    --visual                    Visualize the results of gene family trees before and after removing long branches
Usage:
    Phylo_Tracer.py OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --absolute_branch_length 5 --relative_branch_length 2.5 [--visual]
```
### OrthoFilter_Mono
```
Description:
    To prune phylogenomic noise from both single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing well-established relationships unless otherwise required for specific objectives
Required parameter:
    --input_GF_list            File containing paths to gene tree files, one per line
    --input_taxa               File with taxonomic information for species
    --input_imap               File with classification information of species corresponding to genes
    --branch_length_multiples  Tips whose branch length is greater than or equal to the branch length multiples will be removed and default = 10
    --insert_branch_index      Nodes with insertion coverage ratio and insertion depth greater than or equal to the insert branch index will be removed and default = 10
Optional parameter:
    --visual                   Visualize the results of gene family trees before and after removing outliers and paralogs
Usage:
    Phylo_Tracer.py OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_taxa gene2clade.imap --branch_length_multiples 10 --insert_branch_index 10 [--visual]
```
### TreeTopology_Summarizer
```
Description:
    To enumerate and visualize the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades
Required parameter:
    --input_GF_list    File containing paths to gene tree files, one per line
    --input_imap       File with classification information of species corresponding to genes
Usage:
    Phylo_Tracer.py TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap 
```
### Tree_Visualizer
```
Description:
    To mark tips of gene trees with provided tags, identify GD nodes, and integrate gene duplication results onto the species tree
Required parameter:
    --input_GF_list       File containing paths to gene tree files, one per line
    --input_imap          File with classification information of species corresponding to genes
    --keep_branch         1 or 0 indicates whether or not to preserve branch length information
    --tree_style          The tree style, 'r' means rectangular, 'c' means circular
Optional parameter:
    --gene_categories     File with taxonomic information for species
    --gene_family         File with family classification information corresponding to genes
    --input_sps_tree      Species tree file in Newick format (required with --gene_family)
    --gene_expression     Gene expression level files
    --visual_gd           Visualize the gd node of gene family trees
Usage:
    Phylo_Tracer.py Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--gene_categories [ gene2order.imap  gene2taxa.imap gene2clade.imap] --keep_branch {1,0} --tree_style {r,c} --gene_family gene2family.imap --input_sps_tree sptree.nwk  --gene_expression gene2expression.imap --visual_gd]
```
### GD_Detector
```
Description:
    To identify gene duplication events by reconciling gene trees with a species tree
Required parameter:
    --input_GF_list            File containing paths to gene tree files, one per line
    --input_imap               File with classification information of species corresponding to genes
    --gd_support               GD node support and default = 50
    --subclade_support         The subclade support of GD node and default = 50
    --dup_species_proportion   The proportion of overlapped species from two subclades for a GD event with range [0-1] and default = 0.2
    --dup_species_num          The number of species with species duplications under the GD node
    --input_sps_tree           A species tree file with Newick format
    --deepvar                  Maximum variance of depth and default = 1
Usage:
    Phylo_Tracer.py GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0.5 --dup_species_num 1 --input_sps_tree sptree.nwk --deepvar 1
```
### GD_Visualizer
```
Description:
    To visualize gene duplication detection results and integrate findings onto the species tree
Required parameter:
    --input_sps_tree  A numbered species tree file with Newick format
    --gd_result       Result file of GD_Detector
    --input_imap         File with classification information of species corresponding to genes
Usage:
    Phylo_Tracer.py GD_Visualizer --input_sps_tree sptree.nwk --gd_result gd_result.txt --input_imap gene2sps.imap 
```
### GD_Loss_Tracker
```
Description:
    To analyze and summarize gene duplication loss events across nodes and tips in the species tree
Required parameter:
    --input_GF_list      File containing paths to gene tree files, one per line
    --input_sps_tree     A species tree file with Newick format
    --input_imap         File with classification information of species corresponding to genes
Usage:
    Phylo_Tracer.py GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap
```
### GD_Loss_Visualizer
```
Description:
    To visualize the summary of gene duplication loss events on the context of species tree
Required parameter:
    --gd_loss_result     Result file of gd loss count summary of GD_Loss_Tracker
    --input_sps_tree     A numbered species tree file with Newick format
Usage:
    Phylo_Tracer.py GD_Loss_Visualizer --input_sps_tree numbered_species_tree.nwk --gd_loss_result gd_loss_count_summary.txt 
```
### Ortho_Retriever
```
Description:
    To infer single-copy putative orthologs by splitting paralogs from large-scale gene family trees for multiple species
Required parameter:
    --input_GF_list     File containing paths to gene tree files, one per line
    --input_imap        File with classification information of species corresponding to genes
    --input_gene_length	File with information corresponding to gene lengths
Usage:
    Phylo_Tracer.py Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap
```
### Hybrid_Tracer
```
Description:
    To use the ABAB-BABA test to detect hybridization signals for each potential GD burst event across species tree detect species hybridization events for
Required parameter:
    --input_GF_list      File containing paths to gene tree files, one per line
    --input_Seq_GF_list  File containing paths to sequence alignment files corresponding to the gene trees
    --input_imap         File with classification information of species corresponding to genes
    --input_sps_tree     A species tree file with Newick format
Optional parameter:
    --target_node        File with species names mapped to specific GD nodes, one per line
    --split_groups       Split the GD data into a target number of groups for Hyde processing, and default = 1
Usage:
    Phylo_Tracer.py Hybrid_Tracer --input_GF_list GF_ID2path.imap --input_Seq_GF_list Seq_GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--target_node N1.txt]
```
### Hybrid_Visualizer
```
Description:
    To visualize hybridization signals, highlighting support from gene tree topologies and D-statistic signals
Required parameter:
    --hyde_out        File containing the result of hyde of Hybrid_Tracer
    --input_sps_tree  A species tree file with Newick format
Optional parameter:
    --node            Node model, stack up all the heatmaps for each monophyletic clade respectively, only the squares in all heatmaps were light, the square after superimposition will be light
Usage:
    Phylo_Tracer.py Hybrid_Visualizer --hyde_out hyde.out  --input_sps_tree sptree.nwk [--node]
```
### HaploFinder
```
Description:
    To distinguish gene conversion by tracing subgenome haplotypes through phylogenomic profiling
Required parameter:
    --input_GF_list      File containing paths to gene tree files, one per line
    --input_imap         File with classification information of species corresponding to genes
    --species_a          Name of species A
    --species_b          Name of species B
    --species_a_gff      GFF file of species A
    --species_b_gff      GFF file of species B
    --species_a_lens     Lens file of species A
    --species_b_lens     Lens file of species B
    --gd_support         GD node support [50-100],and default = 0.0005
    --visual_chr_a       A file containing the chromosome numbers of species A is required to visualize specific chromosome regions
    --visual_chr_b       A file containing the chromosome numbers of species B is required to visualize specific chromosome regions
    --size               The size of each point in the dotplot graph and default = 0.0005
Usage:
    Phylo_Tracer.py HaploFinder --input_GF_list GF.list --input_imap gene2sps.imap --species_a A --species_b B --species_a_gff A.gff --species_b_gff B.gff --species_a_lens A.lens --species_b_lens B.lens  --gd_support 50 [--visual_chr_a chr_a.txt --visual_chr_b chr_b.txt --size 0.0001]
```

## Bug Reports

You can report bugs or request features through our [GitHub Issues page](https://github.com/YiyongZhao/PhyloTracer/issues). If you have any questions, suggestions, or issues, please do not hesitate to contact us.

## Contributing

If you're interested in contributing code or reporting bugs, we welcome your ideas and contributions to improve PhyloTracer! Please check out [Contribution Guidelines](https://docs.github.com/en/issues).

## Version History

Check the [Changelog](https://github.com/YiyongZhao/PhyloTracer/commits/PhyloTracer_v1.0.0) for details on different versions and updates.

## License

PhyloTracer is licensed under the [MIT LICENSE](LICENSE).




