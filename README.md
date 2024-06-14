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
                                                                                             
   PhyloTracer: A Comprehensive and User-Friendly Toolkit for Evolutionary Genomic Analysis.
                                                                                             
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


# PhyloTracer

A User-Friendly Toolkit for Comprehensive inference for Manipunation of Tree Foramt, Gene Tree Rooting, the Origin and Loss of Gene Duplication Events, Ortholog Retrieval, Phylogenetic Noise Elimination, Gene Tree Topology Summary, Species Hybridization Detection, and Visualization.

## Introduction

`PhyloTracer` aims to provide more accurate rooting of gene trees, serving as a foundation for inferring putative orthologous genes. It also includes functions to statistically summarize the topology types for models like ABAB-ABBA, aiding in the identification of hybridization signals.


## Module features
1. **PhyloTree_CollapseExpand:** Transforms a phylogenetic tree in Newick format into a ‘comb’ structure based on predefined support value threshold. It can also revert this 'comb' structure back to a fully resolved binary tree, allowing dynamic topology adjustments.
2. **PhyloSupport_Scaler:** Recalibrate support values from bootstrap or posterior probability in a phylogenetic tree, scaling them between \[0,1] and \[1,100] ranges for computational compatibility, and vice versa to meet various analytical needs.
3. **BranchLength_NumericConverter:** Converts branch length values of a phylogenetic tree from string to numerical format, critical for quantitative analysis and computational operations.
4. **Phylo_Rooter:** Enhances the accuracy of gene tree rooting, providing a robust framework for phylogenetic inference.
5. **OrthoFilter_LB:** Prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length.
6. **OrthoFilter_Mono:** Prunes phylogenomic noise from both single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing well-established relationships unless otherwise required for specific objectives.
7. **TreeTopology_Summarizer:** Enumerates the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades.
8. **Tree_Visualizer:** Visualizes and integrates gene duplication detection results into the species tree.
9. **GD_Detector:** identification of gene duplication events by reconciliaiton of gene and species trees.
10. **GD_Visualizer:** Visualizes gene duplication detection results and integrates these findings into the species tree.
11. **GD_Loss_Tracker:** Analyzes and summarizes gene duplication loss events across each node from species tree for each tips .
12. **GD_Loss_Visualizer:** Visualizes the summary of gene duplication loss event on the context of speices tree.
13. **Ortho_Retriever:** Infers single-copy putative orthologs by spliting paralogs from large-scale gene family trees across multiple species.
14. **Hybrid_Tracer:** Uses the ABAB-BABA test to detect hybridization signals for each potential GD burst events across species tree detect species hybridization events for .
15. **Hybrid_Visualizer:** Visualizes hybridization signals, highlighting support from gene tree topologies and D-statistic signals.
16. **HaploFinder:** Distinguishes gene conversion by tracing subgenome haplotypes through phylogenomic profiling.

## Input file requirements

The following input files must strictly comply with the requirements, with tabs separating each column; otherwise, an error will occur:

```bash
--GF.txt-------------------------------------------------------------------------------------------------
OG_104001  example_data/Phylo_Rooter/OG_104001.treefile   
OG_104002  example_data/Phylo_Rooter/OG_104002.treefile    
OG_104003  example_data/Phylo_Rooter/OG_104003.treefile   
.   
.   
.
---------------------------------------------------------------------------------------------------------
❗️: this file must be sorted
--imap.txt-----------------------------------------------------------------------------------------------
ACT_0000001  ACT   
ACT_0000002  ACT   
ACT_0000003  ACT   
AMB_0000001  AMB   
AMB_0000002  AMB   
AMB_0000003  AMB   
AQU_0000001  AQU   
AQU_0000002  AQU   
AQU_0000003  AQU   
.   
.   
.
--------------------------------------------------------------------------------------------------------- 
--length.txt---------------------------------------------------------------------------------------------
ACT_0000001  501   
ACT_0000002  267   
ACT_0000003  903   
AMB_0000001  339   
AMB_0000002  756   
AMB_0000003  1275   
AQU_0000001  2733   
AQU_0000002  219   
AQU_0000003  1131   
.   
.   
.

--taxa.txt or or family.txt or order.txt or class.txt-----------------------------------------------------
ACT  Actinidia_chinensis | ACT  Actinidiaceae | ACT  Ericales     | ACT  Magnoliopsida
AMB  Amborella_trichopoda| AMB  Amborellaceae | AMB  Amborellales | AMB  Magnoliopsida
AQU  Aquilegia_coerulea  | AQU  Ranunculaceae | AQU  Ranunculales | AQU  Magnoliopsida
.
.
.
----------------------------------------------------------------------------------------------------------
```

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
    
Note: PhyloTracer use basic funcitons of analysis and visualization of trees from Python framework [ete3](http://etetoolkit.org/) and detect species hybridizaiotn signals using ABAB-BABA test by [HyDe](https://github.com/pblischak/HyDe).

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

### Install from PyPI with pip:

```bash
pip install PhyloTracer
```

## Features
* Incorporating the principles of maximizing the outgroup depth score, minimizing the Robinson-Foulds (RF) distance, reducing the variance in ingroups' branch lengths, and maximizing the overlap ratio of gene duplication species enhances the accuracy of root determination.   
* Introducing the concept of long-branch genes for noise filtration in gene trees.   
* Introducing the concept of inserted genes for monophyletic filtering in single-copy gene trees.   

## Usage
### PhyloTree_CollapseExpand
```bash
Description
    Collapse the branches in the gene tree whose support is equal to a specific value and convert the comb structure into a binary tree structure
Usage
    Phylo_Tracer.py PhyloTree_CollapseExpand  --input_GF_list  GF.txt  --support_value  50
Arguement
    --input_GF_list		File containing paths to gene tree files, one per line.
    --support_value		Nodes with support values less than or equal to this value will be collapsed.
```
### PhyloSupport_Scaler
```bash
Description
    Recalibrate support values from bootstrap or posterior probability in a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational compatibility, and vice versa to meet various analytical needs.
Usage
    Phylo_Tracer.py PhyloSupport_Scaler  --input_GF_list  GF.txt  --scale  1
Arguement
    --input_GF_list    File containing paths to gene tree files, one per line.
    --scale    1 to scale up (0 to 1 to 1 to 100), 0 to scale down (1 to 100 to 0 to 1).
```
### BranchLength_NumericConverter
```bash
Description
    Normalize the branch length of the gene tree, retaining specific decimal places
Usage
    Phylo_Tracer.py  BranchLength_NumericConverter  --input_GF_list  GF.txt  [ --decimal_place  10 ]
Arguement
    --input_GF_list		File containing paths to gene tree files, one per line.
Optional arguement
    --decimal_place		Set how many decimal places to keep(default=10)
```
### Phylo_Rooter
```bash
Description
    Enhances the accuracy of gene tree rooting, providing a robust framework for phylogenetic inference.
Usage
    Phylo_Tracer.py  Phylo_Rooter  --input_GF_list  GF.txt  --input_imap  imap.txt  --input_gene_length length.txt  --input_sps_tree  sptree.nwk 
Arguement
    --input_GF_list		  File containing paths to gene tree files, one per line.
    --input_imap        File with classification information of species corresponding to genes.
    --input_gene_length	File with information corresponding to gene lengths.
    --input_sps_tree		Species tree file.
```
### OrthoFilter_LB
```bash
Description
    Prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length.
Usage
    Phylo_Tracer.py  OrthoFilter_LB  --input_GF_list  GF.txt  --input_taxa  taxa.txt  --long_branch_index 	10  [ --visual ]
Arguement
    --input_GF_list		    File containing paths to gene tree files, one per line.
    --input_taxa			    File with taxonomic information for species.
    --long_branch_index 	Long branch index
Optional arguement
    --visual				      Visualize the results if set
```
### OrthoFilter_Mono
```bash
Description
  Prunes phylogenomic noise from both single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing 	well-established relationships unless otherwise required for specific objectives.
Usage
    Phylo_Tracer.py  OrthoFilter_Mono  --input_GF_list  GF.txt  --input_taxa  taxa.txt  	--long_branch_index 	10  --insert_branch_index  10  [ --visual ]
Arguement
    --input_GF_list		    File containing paths to gene tree files, one per line.
    --input_taxa			    File with taxonomic information for species.
    --long_branch_index 	Long branch index
    --insert_branch_index	Insert_branch_index
Optional arguement
    --visual				      Visualize the results if set

TreeTopology_Summarizer
Description
Enumerates the frequency of both absolute and relative topologies for single-copy gene trees or interested 	predefined clades.
Usage
Phylo_Tracer.py  TreeTopology_Summarizer  --input_GF_list  GF.txt  --input_imap  imap.txt  --outfile 	filename
Arguement
--input_GF_list		File containing paths to gene tree files, one per line.
--input_imap			File with classification information of species corresponding to genes.
--outfile 				Output filename







Tree_Visualizer
Description
Visualizes and integrates gene duplication detection results into the species tree.
Usage
Phylo_Tracer.py Tree_Visualizer  --input_GF_list  GF.txt  --input_imap  imap.txt  [ --gene_categories 	[ taxa.txt  family.txt  order.txt  class.txt ]  --keep_branch {1,0}  --tree_style {r,c}  --gene_family  	gene2family.txt  --input_sps_tree  sptree.nwk  --gene_expression  gene_expression.csv ]
Arguement
--input_GF_list		File containing paths to gene tree files, one per line.
--input_imap			File with classification information of species corresponding to genes.
--keep_branch		1 or 0 indicates whether or not to preserve branch length information
--tree_style			The treestyle, ‘r’ is meaning rectangular, ‘c’ is meaning circular
Optional arguement
--gene_categories	File with taxonomic information for species.
--gene_family		If you want to mark gene families you need to provide this file
--input_sps_tree		If you provide the --gene_family parameter, you must provide the species tree
--gene_expression	Gene expression level files

GD_Detector
Description
identification of gene duplication events by reconciliaiton of gene and species trees.
Usage
Phylo_Tracer.py  GD_Detector  --input_GF_list  GF.txt  --input_imap  imap.txt  --gd_support 	50 	--clade_support  50  --dup_species_radio  0.5  --dup_species_num  2  --input_sps_tree  sptree.nwk
Arguement
--input_GF_list		File containing paths to gene tree files, one per line.
--input_imap			File with classification information of species corresponding to genes.
--gd_support			GD node support [50-100]
--clade_support		The children support of GD node [0-100]
--dup_species_radio	The proportion of species with species duplications under the GD node [0-1]
--dup_species_num	The number of species with species duplications under the GD node
--input_sps_tree		Species tree file

GD_Visualizer
Description
Visualizes gene duplication detection results and integrates these findings into the species tree.
Usage
Phylo_Tracer.py  GD_Visualizer  --input_sps_tree  sptree.nwk  --gd_result  gd_result.txt
Arguement
--input_sps_tree		Species tree file
--gd_result			Result file of GD_Detector




GD_Loss_Tracker
Description
Analyzes and summarizes gene duplication loss events across each node from species tree for each tips .
Usage
Phylo_Tracer.py  GD_Loss_Tracker   --input_GF_list  GF.txt  --input_sps_tree  sptree.nwk  		--output_folder 	filename
Arguement
--input_GF_list		File containing paths to gene tree files, one per line.
--input_sps_tree		Species tree file
--output_folder		Output foldername

GD_Loss_Visualizer
Description
Visualizes the summary of gene duplication loss event on the context of speices tree.
Ortho_Retriever: Infers single-copy putative orthologs by spliting paralogs from large-scale gene family 	trees across multiple species.
Usage
Phylo_Tracer.py  GD_Loss_Visualizer  --input_folder  input_foldername  --output_folder  	output_foldername
Arguement
--input_folder			Input foldername
--output_folder		Output foldername

Ortho_Retriever
Description
Infers single-copy putative orthologs by spliting paralogs from large-scale gene family trees across multiple 	species.
Usage
Phylo_Tracer.py  Ortho_Retriever  --input_GF_list  GF.txt  --input_imap  imap.txt  	--input_gene_length  sptree.nwk  
Arguement
--input_GF_list		File containing paths to gene tree files, one per line.
--input_imap			File with classification information of species corresponding to genes.
--input_gene_length	File with information corresponding to gene lengths.

Hybrid_Tracer
Description
Uses the ABAB-BABA test to detect hybridization signals for each potential GD burst events across species 	tree detect species hybridization events for.
Usage
Phylo_Tracer.py  Hybrid_Tracer  --input_GF_list  GF.txt  --input_Seq_GF_list  Seq_GF.txt 		--input_sps_tree  sptree.nwk  --input_imap  imap.txt
Arguement
--input_GF_list		File containing paths to gene tree files, one per line.
--input_Seq_GF_list	File containing paths to sequence alignment files corresponding to the gene trees.
--input_imap			File with classification information of species corresponding to genes.
--input_sps_tree 		Species tree file



Hybrid_Visualizer
Description
Visualizes hybridization signals, highlighting support from gene tree topologies and D-statistic signals.
HaploFinder
Description
Distinguishes gene conversion by tracing subgenome haplotypes through phylogenomic profiling.

## Bug Reports

You can report bugs or request features through our [GitHub Issues page](https://github.com/YiyongZhao/PhyloTracer/issues). If you have any questions, suggestions, or encounter any issues, please do not hesitate to contact us.

## Contributing

If you're interested in contributing code or reporting bugs, we welcome your ideas and contributions to improve PhyloTracer! Please check out [Contribution Guidelines](https://docs.github.com/en/issues).

## Version History

Check the [Changelog](https://github.com/YiyongZhao/PhyloTracer/commits/PhyloTracer_v1.0.0) for details on different versions and updates.

## License

PhyloTracer is licensed under the [MIT LICENSE](LICENSE).




