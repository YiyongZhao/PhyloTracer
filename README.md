
<div align="center">
  
# <img src="logo/PhyloTracer_logo.png" width="80" height="80" align="center"> PhyloTracer V1.0.0 </div> 

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
    Contacts: Taoli(Taoli@gmail.com); Yiyong Zhao(yiyongzhao1991@gmail.com)                  
                                                                                             
###############################################################################################
```

# PhyloTracer

A User-Friendly Toolkit for Comprehensive inference for Manipunation of Tree Foramt, Gene Tree Rooting, the Origin and Loss of Gene Duplication Events, Ortholog Retrieval, Phylogenetic Noise Elimination, Gene Tree Topology Summary, Species Hybridization Detection, and Visualization.

## Introduction

PhyloTracer aims to provide more accurate rooting of gene trees, serving as a foundation for inferring true orthologous genes. It also includes functionality to statistically summarize the topology types for models like ABAB-ABBA, aiding in the identification of hybridization signals.


## Module features
1. **PhyloTree_CollapseExpand:** This module enables the transformation of a phylogenetic tree in Newick format into a ‘comb’ structure, contingent upon a predefined threshold of support values. Conversely, it can also expand a ‘comb’ structure into a fully resolved binary tree, thus facilitating dynamic topology adjustments based on support value thresholds.
2. **PhyloSupport_Scaler:** This module is designed to recalibrate support values within a phylogenetic tree. It can scale a support value PP from a range of [0,1] to a Bootstrap range of [1,100], or inversely scale down a support value from [0,100] to a range of [0,1], thereby accommodating different analytical requirements.
3. **BranchLength_NumericConverter:** This module converts the branch length values of a phylogenetic tree from a string type to a numerical format. This conversion is essential for subsequent quantitative analyses and computational operations on the tree.
4. **Phylo_Rooter:** This module enhances the accuracy of gene tree rooting, providing a robust framework for phylogenetic inference.
5. **OrthoFilter:** This module is designed to handle the process of single-copy and multi-copy gene families. It utilizes provided labels to filter gene trees, identify potential orthologous genes, and remove long-branch outliers and phylogenomic noise. Additionally, OrthoFilter ensures that known labeled taxa (such as family, order, etc.) form monophyletic groups in the gene tree.
6. **TreeTopology_Summarizer:** Enumerates the frequency  of  both absolute and relative topologies for single-copy gene trees and interested clades.
7. **Tree_Visualizer:** Provides an intuitive visualization for phylogenetic trees, enable labeling tips with multi-layered annotation labeling.
8. **GD_Detector:**  Facilitates the identification of gene duplication events by reconciliaiton of gene and species trees.
9. **GD_Visualizer:** Visualizes gene duplication detection results and integrates these findings into the species tree.
10. **GD_Loss_Tracker:** Analyzes and summarizes gene duplication loss events for each tips across species tree.
11. **GD_Loss_Visualizer:** Presents a visual summary of gene duplication loss event on the context of speices tree.
12. **Ortho_Retriever:** Putative orthologs inferrenec from large-scale gene family trees across numerous species.
15. **Hybrid_Tracer:** Utilizes the ABAB-BABA test to detect species hybridization signals.
16. **Hybrid_Visualizer:** Visualizes hybridization signals, highlighting gene tree topology ratios (ABB+BAA) that support allopolyploidy, D-statistic signals, and GCN-based predictions.
17. **HaploFinder:** Distinguishing Gene Conversion by Tracing Subgenome Haplotypes through Phylogenomic Profiling.
18. 
## Installation

### Required dependencies:

* Python 3.0+
* Python modules:
  * ete3
  * pandas
  * numpy
  * tqdm
  * time
  * pypdf4
  * matplotlib
  * pyqt5

### Clone and Install Environment:

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

## Usage

## Example Scenarios:

1. **Phylo_Rooter:**

```bash
python PhyloTracer.py Phylo_Rooter --input_GF_list GF.txt --input_imap imap.txt --input_sps_tree 30sptree.nwk --input_gene_length length.txt
```

2. **PhyloNoise_Filter:**

```bash
python PhyloTracer.py PhyloNoise_Filter --input_GF_list GF.txt --input_taxa taxa
```
    
3. **TreeTopology_Summarizer:**

```bash
python PhyloTracer.py TreeTopology_Summarizer --input_GF_list GF.txt --input_imap imap.txt
```

4. **Tree_Visualizer:**

```bash
Python PhyloTracer.py Tree_Visualizer --input_GF_list GF.txt --input_imap imap.txt --gene_categories genus order --keep_branch 1 --tree_style r
```

5. **GD_Detector:**

```bash
python PhyloTracer.py GD_Detector --input_GF_list GF.txt --input_imap imap.txt --input_sps_tree 30sptree.nwk --support 50 --dup_species_radio 0.5 --dup_species_num 2
```

6. **Ortho_Retriever:**

```bash
python PhyloTracer.py Ortho_Retriever --input_GF_list GF.txt --input_imap imap.txt --input_sps_tree 30sptree.nwk --input_gene_length length.txt
```

7. **GeneDynamics_Visualizer:**

```bash
python PhyloTracer.py GeneDynamics_Visualizer  --input_sps_tree sptree.nwk --input_summary_tree summary_tree
```
    

## Contributing

If you are interested in contributing code or reporting bugs, please check the [Contribution Guidelines](CONTRIBUTING.md).

## Version History

Check the [Changelog](CHANGELOG.md) for details on different versions and updates.

## License

PhyloTracer is licensed under the [License Name]. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or suggestions, feel free to contact us via [email](mailto:your.email@example.com).
Feel free to customize it further according to your preferences or additional details you want to include.


