.. include:: links.rst

.. _Input_Files:

PhyloTracer Data Files
======================

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

Path Map
--------
The main data files for running PhyloTracer are (1) the DNA sequences and (2) the
map of sampled individuals to their respective taxa/populations. Examples of these files
can be found in the  ``examples/`` directories within the main PhyloTracer
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

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

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

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

.. code::

  AMTR_s00796p00010580  201
  ATCG00500.1           1467
  Glyma.07G273800.2     3417
  .
  .
  .

Gene To Family Map
------------------

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

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

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

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

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

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

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

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

The main data files for running PhyloTracer. Examples of these files can be found in the  ``examples/`` directories within the main PhyloTracer folder. Below we describe each file in more detail. If you installed PhyloTracer from PyPI and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/YiyongZhao/PhyloTracer/examples/>`__.

**Example:** ``gene2expression.imap``

.. code::

  AMTR_s00796p00010580  5.0
  ATCG00500.1           12.0
  Glyma.07G273800.2     0.0
  .
  .
  .

You can add any number of imap files. They will sequentially provide annotations to the right of the gene tips according to the order of input.
