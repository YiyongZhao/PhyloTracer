
<div align="center">
  
# <img src="logo/PhyloTools_logo.jpg" width="60" height="60" align="center"> PhyloTools V1.1.1 </div> 

```
############################################################################################################  
#                       ____  _           _     _____           _                                          #
#                      |  _ \| |__  _   _| | __|_   _|__   ___ | |___                                      #
#                      | |_) | '_ \| | | | |/ _ \| |/ _ \ / _ \| / __|                                     #
#                      |  __/| | | | |_| | | (_) | | (_) | (_) | \__ \                                     #
#                      |_|   |_| |_|\__, |_|\___/|_|\___/ \___/|_|___/                                     #
#                                   |___/                                                                  #
#                                                                                                          #
#   PhyloTools: A User-Friendly Toolkit for Gene Tree Rooting, Gene Duplication Identification, Ortholog   #
#   Retrieval, Phylogenetic Noise Elimination, Species Hybridization Detection,and Visualization.          #                            
#                                                                                                          #
#   Pypi: https://pypi.org/project/PhyloTools                                                              #
#   Github: https://github.com/YiyongZhao/PhyloTools                                                       #              
#   Licence: MIT license                                                                                       #
#   Release Date: 2023-7                                                                                   #
#   Please cite: Li et al. 2024, XXXX.                                                                     #
#   Contacts: Taoli(Taoli@gmail.com); Yiyong Zhao(yzhao@bwh.harvard.edu)                                   #
#                                                                                                          #
############################################################################################################
```


# User Mannual and Guide:
## gene-tree rooting and gene dup mark

PhyloTools v1.1.1
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
git clone https://github.com/lt11300/PhyloTracer.git
