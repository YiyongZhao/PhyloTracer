
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


# User Mannual and Guide:
## gene-tree rooting and gene dup mark

PhyloTracer v1.1.1
A bioinformatics software that uses RF distance, branch length, and gene length for rooting

I. Prerequisites and installation
If you are having trouble installing PhyloTracer, please contact the lead developer,taoli, via email to get help.
To install using pip, we strongly recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:

PhyloTracer
IMPORTANT: PhyloTracer is currently only compatible with Python 3

Installation via conda
Create a conda environment under Python 3 and activate the environment

#install pandas, numpy, tqdm, time, and ete3
conda install -c etetoolkit ete3
conda install pandas
conda install numpy
conda install tqdm
conda install time

# deactivate conda environment
deactivate

To install PhyloTracer, simply download it use :git
git clone https://github.com/YiyongZhao/PhyloTracer.git
