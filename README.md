
<div align="center">
  
# <img src="logo/PhyloTracer_logo.png" width="80" height="80" align="center"> PhyloTracer V1.1.1 </div> 

```
###############################################################################################
#                                                                                             #
# ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ████████╗██████╗  █████╗  ██████╗███████╗██████╗  #
# ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗╚══██╔══╝██╔══██╗██╔══██╗██╔════╝██╔════╝██╔══██╗ #
# ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║   ██║   ██████╔╝███████║██║     █████╗  ██████╔╝ #
# ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║   ██║   ██╔══██╗██╔══██║██║     ██╔══╝  ██╔══██╗ #
# ██║     ██║  ██║   ██║   ███████╗╚██████╔╝   ██║   ██║  ██║██║  ██║╚██████╗███████╗██║  ██║ #
# ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝    ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚══════╝╚═╝  ╚═╝ #                            
#                                                                                             #
#    PhyloTracer: A User-Friendly Toolkit for Gene Tree RootingGene Duplication, Expansion    #
#    and Contraction ldentification.Ortholog Retrieval, Phylogenetic Noise Elimination,       #
#    SpeciesHybridization Detection, and Visualization.                                       #
#                                                                                             #
#    Pypi: https://pypi.org/project/PhyloTracer                                               #
#    Github: https://github.com/YiyongZhao/PhyloTracer                                        #
#    Licence: MIT license                                                                     #
#    Release Date: 2023-7                                                                     #
#    Contacts: Taoli(Taoli@gmail.com); Yiyong Zhao(yiyongzhao1991@gmail.com)                  #
#                                                                                             #
###############################################################################################
```


# PhyloTracer

PhyloTracer is a User-Friendly Toolkit for Gene Tree RootingGene Duplication, Expansion and Contraction ldentification.Ortholog Retrieval, Phylogenetic Noise Elimination, SpeciesHybridization Detection, and Visualization.

## Introduction

PhyloTracer aims to provide more accurate rooting of gene trees, serving as a foundation for inferring true orthologous genes. It also includes functionality to statistically summarize the topology types for models like ABAB-ABBA, aiding in the identification of hybridization signals.

## Features

1. **Phylo_Rooter:** Provides more accurate gene tree rooting.
2. **PhyloNoise_Filter:** Filters based on provided labels to identify true orthologous genes.
3. **TreeTopology_Summarizer:** Counts the occurrences of the absolute topology and relative topology of single copy gene tree.
4. **Tree_Visualizer:** Offers an intuitive visualization of phylogenetic trees.
5. **GD_Detector:** Assists in identifying hybridization signals in gene trees.
6. **GD_Visualizer:** Visualize the results of GD_Detector and summarize them on the species tree.
7. **GD_Loss_Tracker:** Explore and summarize the situation of GD_loss.
8. **GD_Loss_Visualizer:** Visualize and summarize the results of GD_Loss Tracker.
9. **Ortho_Retriever:** Splits multi-copy gene trees.
10. **GeneDynamics_Tracker:** Explore and summarize the loss and replication of genes.
11. **GeneDynamics_Visualizer:** According to gene_gain_and_generate visualized PDF files for loss information generation.
12. **Hybrid_Tracer:** Detecting hybrid signals.
13. **Hybrid_Visualizer:** Visualize the results of Hybrid_Tracer.

## Installation

### Requirements:
* Python 3.0+
* Python Modules:
  * ete3
  * pandas
  * numpy
  * tqdm
  * time
  * PyPDF4
  * matplotlib

### Clone and Install:

```bash
git clone https://github.com/YiyongZhao/PhyloTracer.git  
cd PhyloTracer  
python setup.py install
```

### Install from PyPI with pip:

```bash
pip install PhyloTracer
```

## Usage

### Example Scenarios:

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
    python PhyloTracer.py Tree_Visualizer --input_GF_list GF.txt --input_imap imap.txt --gene_categories genus order --keep_branch 1 --tree_style r
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
```

Feel free to customize it further according to your preferences or additional details you want to include.

`

