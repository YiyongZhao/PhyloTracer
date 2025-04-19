.. PhyloTracer documentation master file, created by
   sphinx-quickstart on Wed Jan  8 15:57:49 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PhyloTracer: A Versatile Toolkit for Comparative Genomics and Phylogenomics Analysis
====================================================================================

**PhyloTracer** is a user-friendly toolkit designed for comprehensive phylogenetic analysis, including tree format manipulation, gene tree rooting, detection of gene duplication origins and losses, ortholog retrieval, phylogenetic noise elimination, gene tree topology summarization, species hybridization detection, and visualization. By enabling more accurate rooting of gene trees, PhyloTracer lays a solid foundation for inferring putative orthologous genes. Additionally, it provides tools to statistically summarize topology types, such as ABAB-ABBA models, facilitating the identification of hybridization signals.

Features
--------

- Incorporating the principles of maximizing the outgroup depth score, minimizing the Robinson-Foulds (RF) distance, reducing the variance in ingroup branch lengths, and maximizing the overlap ratio of gene duplication species enhances the accuracy of root determination.
- Use GD clade to verify hybridization signals between species.
- Introducing the concept of long-branch genes for noise filtration in gene trees.
- Introducing the concept of inserted genes for monophyletic filtering in single-copy gene trees.

Documentation
=============

.. toctree::
   :maxdepth: 1

   installation.rst
   input_files.rst
   analyze.rst


