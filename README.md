
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

1. **Phylo_Rooting:** Provides more accurate gene tree rooting.
2. **Eliminate_PhyloNoise:** Filters based on provided labels to identify true orthologous genes.
3. **Statistical_Topology:** Counts the occurrences of different topology types for summarizing ABAB-ABBA models.
4. **Tree_Visualization:** Offers an intuitive visualization of phylogenetic trees.
5. **GD_Detector:** Assists in identifying hybridization signals in gene trees.
6. **Ortho_Split:** Splits multi-copy gene trees.

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

## Usage

### Example Scenarios:

1. **Eliminate PhyloNoise:**

    ```bash
    python PhyloTracer.py Eliminate_PhyloNoise --input_GF_list GF_list.txt --input_taxa taxa
    ```

2. **Gene Duplication Detection (GD Detector):**

    ```bash
    python PhyloTracer.py GD_Detector --input_GF_list GF.txt --input_imap imap.txt --input_sps_tree 30sptree.nwk --support 50 --dup_species_radio 0.5 --dup_species_num 2
    ```

3. **Orthologous Gene Splitting:**

    ```bash
    python PhyloTracer.py Ortho_Split --input_GF_list GF.txt --input_imap imap.txt --input_sps_tree 30sptree.nwk --input_gene_length length.txt
    ```

4. **Statistical Topology Analysis:**

    ```bash
    python PhyloTracer.py Statistical_Topology --input_GF_list GF.txt --input_imap imap.txt
    ```

5. **Tree Visualization:**

    ```bash
    python PhyloTracer.py Tree_Visualization --input_GF_list GF.txt --input_imap imap --gene_categories genus order --keep_branch 1 --tree_style r
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

Make sure to replace placeholders like "[License Name]", "[dependency1 dependency2]", "[Changelog link]", "[email]", and others with the actual information for your project.

***
To install PhyloTracer, simply download it use :
git clone https://github.com/YiyongZhao/PhyloTracer.git

`

