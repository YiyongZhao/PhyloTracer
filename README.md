
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
#    PhyloTracer: A User-Friendly Toolkit for Gene Tree Rooting, Gene Duplication             #
#    Identification, Ortholog Retrieval, Phylogenetic Noise Elimination, Species              #
#    Hybridization Detection,and Visualization.                                               #
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

PhyloTracer is a tool for gene tree analysis, offering six major functionalities: rooting, monophyly filtering, topology statistics, phylogenetic tree visualization, detection of gene duplication events, and the splitting of multi-copy trees.

## Introduction

PhyloTracer aims to provide more accurate rooting of gene trees, serving as a foundation for inferring true orthologous genes. It also includes the functionality to statistically summarize the topology types for models like ABAB-ABBA, aiding in the identification of hybridization signals.

## Features

1. **Rooting:** Provides more accurate gene tree rooting.
2. **Monophyly Filtering:** Filters based on provided labels to identify true orthologous genes.
3. **Topology Statistics:** Counts the occurrences of different topology types for summarizing ABAB-ABBA models.
4. **Phylogenetic Tree Visualization:** Offers an intuitive visualization of phylogenetic trees.
5. **Duplication Detection (GD):** Assists in identifying hybridization signals in gene trees.
6. **Multi-Copy Tree Splitting:** Splits multi-copy gene trees.

# Installation
***
# Requirements:
* Python 3.0+
* Python Modules:
  * ete3
  * pandas
  * numpy
  * tqdm
  * time
  * PyPDF4
  * matplotlib

# Clone HyDe repository from GitHub
git clone https://github.com/YiyongZhao/PhyloTracer.git
cd PhyloTracer

# Now install PhyloTracer module
python setup.py install

# Test the installation
make test

# Install from PyPI with pip
pip install PhyloTracer

***
To install PhyloTracer, simply download it use :
git clone https://github.com/YiyongZhao/PhyloTracer.git
