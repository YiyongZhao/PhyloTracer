
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
14. **Hybrid_Tracer:** Uses the ABAB-BABA test to detect hybridization signals for each potential GD burst events across species tree
15. detect species hybridization events for .
16. **Hybrid_Visualizer:** Visualizes hybridization signals, highlighting support from gene tree topologies and D-statistic signals.
17. **HaploFinder:** Distinguishes gene conversion by tracing subgenome haplotypes through phylogenomic profiling.
    
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
    
`PhyloTracer`use Python framework ([ete3](http://etetoolkit.org/)) for the analysis and visualization of trees. Hyde ([seqwish](https://github.com/pblischak/HyDe)) was used to detect species hybridizaiotn signals by ABAB-BABA test.

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
pip install `PhyloTracer`
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

Check the Changelog: https://github.com/YiyongZhao/PhyloTracer/commits/PhyloTracer_v1.0.0 for details on different versions and updates.

## License

PhyloTracer is licensed under the [License Name]. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or suggestions, feel free to contact us via [email](mailto:your.email@example.com).
Feel free to customize it further according to your preferences or additional details you want to include.


