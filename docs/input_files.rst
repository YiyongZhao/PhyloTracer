.. include:: links.rst

.. _Input_Files:

PhyloTracer Data Files
======================

The main data files for running PhyloTracer. Examples of these files can be found in the  ``example_data/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/tree/main/example_data/>`__.

Path Map
--------
A tab-delimited file mapping gene family identifiers to their corresponding gene tree file paths. Each line contains a gene family ID followed by the path to its tree file.
**Example:** ``GF_ID2path.imap``

.. code::

  OG_1  example_data/Phylo_Rooter/OG_1.treefile   
  OG_2  example_data/Phylo_Rooter/OG_2.treefile    
  OG_3  example_data/Phylo_Rooter/OG_3.treefile
  .
  .
  .

Gene To Species Map
-------------------

A two-column tab-delimited file mapping each gene identifier to its corresponding species name. This file is required by most PhyloTracer modules.

**Example:** ``gene2sps.imap``

.. code::

  AMTR_s00796p00010580  Amborella_trichopoda
  ATCG00500.1           Arabidopsis_thaliana
  Glyma.07G273800.2     Glycine_max
  .
  .
  .

GeneLength Map
--------------

A two-column tab-delimited file mapping each gene identifier to its sequence length (in base pairs or amino acids). Used by Phylo_Rooter and Ortho_Retriever.

.. code::

  AMTR_s00796p00010580  201
  ATCG00500.1           1467
  Glyma.07G273800.2     3417
  .
  .
  .

Gene To Family Map
------------------

A two-column tab-delimited file mapping each gene identifier to its taxonomic family. Used by Tree_Visualizer for family-level annotations.

**Example:** ``gene2family.imap``

.. code::

  AMTR_s00796p00010580  Amborellaceae
  ATCG00500.1           Brassicaceae
  Glyma.07G273800.2     Fabaceae
  .
  .
  .

Gene To Order Map
-----------------

A two-column tab-delimited file mapping each gene identifier to its taxonomic order. Used by Tree_Visualizer for order-level category annotations.

**Example:** ``gene2order.imap``

.. code::

  AMTR_s00796p00010580  Amborellales
  ATCG00500.1           Brassicales
  Glyma.07G273800.2     Fabales
  .
  .
  .

Gene To Taxa Map
----------------

A two-column tab-delimited file mapping each gene identifier to a higher-level taxonomic group (e.g., Angiosperm, Malvids). Used for taxonomic category annotations in visualization.

**Example:** ``gene2taxa.imap``

.. code::

  AMTR_s00796p00010580  Angiosperm
  ATCG00500.1           Malvids
  Glyma.07G273800.2     Fabids
  .
  .
  .

Gene To Clade Map
-----------------

A two-column tab-delimited file mapping each gene identifier to a predefined clade or lineage label (e.g., Nitrogen-fixing). Used by OrthoFilter_Mono for monophyletic constraint filtering.

**Example:** ``gene2clade.imap``

.. code::

  AMTR_s00796p00010580  Nitrogen-fixing
  ATCG00500.1           Nitrogen-fixing
  Glyma.07G273800.2     non-Nitrogen-fixing
  .
  .
  .

Gene To Expression Map
----------------------

A two-column tab-delimited file mapping each gene identifier to its expression level. Used by Tree_Visualizer for overlaying expression data on gene tree figures.

**Example:** ``gene2expression.imap``

.. code::

  AMTR_s00796p00010580  5.0
  ATCG00500.1           12.0
  Glyma.07G273800.2     0.0
  .
  .
  .

You can add any number of imap files. They will sequentially provide annotations to the right of the gene tips according to the order of input.
