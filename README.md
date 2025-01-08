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

A User-Friendly Toolkit for Comprehensive inference for the Manipulation of Tree Format, Gene Tree Rooting, the Origin and Loss of Gene Duplication Events, Ortholog Retrieval, Phylogenetic Noise Elimination, Gene Tree Topology Summary, Species Hybridization Detection, and Visualization.

---
## Introduction

`PhyloTracer` aims to provide more accurate rooting of gene trees, serving as a foundation for inferring putative orthologous genes. It also includes functions to statistically summarize the topology types for models like ABAB-ABBA, aiding in identifying hybridization signals.

---
## Module features
1. **PhyloTree_CollapseExpand:** To transform a phylogenetic tree with Newick format into a ‘comb’ structure based on a predefined support value threshold. It can also revert the 'comb' structure to the binary tree, allowing it to meet the standard software analysis requirements.
2. **PhyloSupport_Scaler:** To recalibrate support values (bootstrap/posterior probability) for a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational requirements.
3. **BranchLength_NumericConverter:** To convert branch length values from string format to numeric format for a phylogenetic tree, which is critical for quantitative analysis and computational operations.
4. **Phylo_Rooter:** Enables an accuracy method for gene tree rooting and enhancing the downstream evolutionary genomic analysis.
5. **OrthoFilter_LB:** To prune phylogenomic noises from single-copy and multi-copy gene family trees by removing the tips with long branch lengths.
6. **OrthoFilter_Mono:** To prune phylogenomic noise from single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing well-established relationships unless otherwise required for specific objectives.
7. **TreeTopology_Summarizer:** To enumerate the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades.
8. **Tree_Visualizer:** To visualize gene duplication events on phylogenetic trees by marking duplication nodes, adding color-coded labels, integrating gene expression data, and summarizing duplications on the species tree for intuitive evolutionary analysis.
9. **GD_Detector:** To identify gene duplication events by reconciliation of gene family trees to a species tree.
10. **GD_Visualizer:** To visualize gene duplication detection in the context of a species tree.
11. **GD_Loss_Tracker:** To track the gene duplication loss event starting across each node/tip from a specific GD burst event in the species tree.
12. **GD_Loss_Visualizer:** To visualize the summary of gene duplication loss event counts for each node/tip in the context of the species tree.
13. **Ortho_Retriever:** To rapidly infer putative single-copy orthologs by splitting paralogs from large-scale gene family trees across multiple species.
14. **Hybrid_Tracer:** To detect hybridization signals for each potential GD burst event across species tree by using the D-statistic (ABAB-BABA) test.
15. **Hybrid_Visualizer:** To visualize hybridization signals, highlighting support from gene tree topologies and D-statistic signals.
16. **HaploFinder:** Distinguishing ancient genome recombination events including gene conversions and crossovers by tracing subgenome haplotypes through phylogenomic profiling.
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
bash install_package.sh

#Reminder for potential visualization issues: qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found and this application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
#Alternative available platform plugins include: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-#egl, wayland-xcomposite-glx, webgl, xcb. before running PhyloTracer, please execute the following bash command:
export QT_QPA_PLATFORM=linuxfb
```
---
### Install from PyPI with pip:

```bash
pip install PhyloTracer
```
---
## Features
* Incorporating the principles of maximizing the outgroup depth score, minimizing the Robinson-Foulds (RF) distance, reducing the variance in ingroup branch lengths, and maximizing the overlap ratio of gene duplication species enhances the accuracy of root determination.
* Use GD clade to verify hybridization signals between species.
* Introducing the concept of long-branch genes for noise filtration in gene trees.
* Introducing the concept of inserted genes for monophyletic filtering in single-copy gene trees.
  
---
## Installation

### Required dependencies:

* Python 3.0+
* Python modules:
  * ete3
  * HyDe
  * pandas
  * numpy
  * tqdm
  * time
  * pypdf4
  * matplotlib
  * pyqt5
    
Note: PhyloTracer uses basic functions of analysis and visualization of trees from Python framework [ete3](http://etetoolkit.org/) and detects species hybridization signals using ABAB-BABA test by [HyDe](https://github.com/pblischak/HyDe).

---
## Example input files
The following input file should have two columns and be separated by a tab key.
```
------------GF_ID2path.imap--------------------------------------------------------------------------------------------------------
OG_104001  example_data/Phylo_Rooter/OG_104001.treefile   
OG_104002  example_data/Phylo_Rooter/OG_104002.treefile    
OG_104003  example_data/Phylo_Rooter/OG_104003.treefile

------------gene2length.imap---------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  201
ATCG00500.1           1467
Glyma.07G273800.2     3417

------------gene2sps.imap #should sorted by second column with the recommended bash command: sort -k2,2 gene2sps.imap-------------------
AMTR_s00796p00010580  Amborella_trichopoda
ATCG00500.1           Arabidopsis_thaliana
Glyma.07G273800.2     Glycine_max

------------gene2family.imap------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborellaceae
ATCG00500.1           Brassicaceae
Glyma.07G273800.2     Fabaceae

------------gene2order.imap--------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Amborellales
ATCG00500.1           Brassicales
Glyma.07G273800.2     Fabales

------------gene2taxa.imap--------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Angiosperm
ATCG00500.1           Malvids
Glyma.07G273800.2     Fabids

------------gene2clade.imap--------------------------------------------------------------------------------------------------------
AMTR_s00796p00010580  Nitrogen-fixing
ATCG00500.1           Nitrogen-fixing
Glyma.07G273800.2     non-Nitrogen-fixing

------------gene2expression.imap---------------------------------------------------------------------------------------------------
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
    --tree_style          The tree style, 'r' is meaning rectangular, 'c' is meaning circular
Optional parameter:
    --gene_categories     File with taxonomic information for species
    --gene_family         File with family classification information corresponding to genes
    --input_sps_tree      If you provide the --gene_family parameter, you should provide a species tree file with Newick format
    --gene_expression     Gene expression level files
    --visual_gd           Visualize the gd node of gene family trees
Usage:
    Phylo_Tracer.py Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--gene_categories [ gene2order.imap  gene2taxa.imap gene2clade.imap] --keep_branch {1,0} --tree_style {r,c} --gene_family gene2family.imap --input_sps_tree sptree.nwk  --gene_expression gene2expression.imap --visual_gd]
```
### GD_Detector
```
Description:
    To identify gene duplication events by reconciling gene trees and species tree
Required parameter:
    --input_GF_list            File containing paths to gene tree files, one per line
    --input_imap               File with classification information of species corresponding to genes
    --gd_support               GD node support and default = 50
    --subclade_support         The subclade support of GD node and default = 50
    --dup_species_proportion   The proportion of overlapped species from two subclades for a GD event with range [0-1] and default = 0.2
    --dup_species_num          The number of species with species duplications under the GD node
    --input_sps_tree           A species tree file with Newick format
    --deepvar                  Maximum variance of deepth and default = 1
Usage:
    Phylo_Tracer.py GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0.5 --dup_species_num 1 --input_sps_tree sptree.nwk --deepvar 1
```
### GD_Visualizer
```
Description:
    To visualize gene duplication detection results and integrate these findings onto the species tree
Required parameter:
    --input_sps_tree  A numbered species tree file with Newick format
    --gd_result       Result file of GD_Detector
Usage:
    Phylo_Tracer.py GD_Visualizer --input_sps_tree sptree.nwk --gd_result gd_result.txt
```
### GD_Loss_Tracker
```
Description:
    To analyze and summarize gene duplication loss events traverse each node from the species tree for each tips
Required parameter:
    --input_GF_list      File containing paths to gene tree files, one per line
    --input_sps_tree     A species tree file with Newick format
    --input_imap         File with classification information of species corresponding to genes
    --all                If specified, detects gene duplications (GD) and loss events across all nodes
    --start_node         File containing species of node for detecting gene duplications (GD) and losses. This limits the detection to the subtree rooted at the specified node
    --end_species        The species where detection ends. Only events affecting species up to and including the specified species will be detected
Usage:
    Phylo_Tracer.py GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--all --start_node --end_species ]
```
### GD_Loss_Visualizer
```
Description:
    To visualize the summary of gene duplication loss events on the context of species tree
Required parameter:
    --input_sps_tree   A numbered species tree file with Newick format
Usage:
    Phylo_Tracer.py GD_Loss_Visualizer --input_sps_tree numbered_species_tree.nwk
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
    --target_node        File with the species name that maps specific gd node
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




