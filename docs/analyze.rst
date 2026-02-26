.. include:: links.rst

.. _usage:

Usage
=====

This section provides detailed usage instructions for the various modules and features of **PhyloTracer**. Each command includes a description, required parameters, optional parameters (if any), and an example of how to use it.

Recommended workflow order for tree-processing modules is root-first:
``Phylo_Rooter -> PhyloTree_CollapseExpand -> PhyloSupport_Scaler -> BranchLength_NumericConverter -> downstream modules``.


1. PhyloTree_CollapseExpand
---------------------------

Description:
    To transform a phylogenetic tree in Newick format into a 'comb' structure based on a predefined support value threshold. It can also revert this `comb` structure to a fully resolved binary tree, allowing dynamic topology adjustments.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--support_value``   Nodes whose support is less than or equal to this value will be collapsed (default is 50).

Optional Parameters:
    - ``--revert``          Revert the 'comb' structure to a fully resolved binary tree.

Example:
    .. code-block:: bash

        PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50 [--revert]


2. PhyloSupport_Scaler
-----------------------

Description:
    To recalibrate support value from bootstrap or posterior probability in a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational compatibility, and vice versa to meet various analytical needs.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--scale_to``        Specify '1' to scale support values from 1-100 to 0-1, or '100' to scale from 0-1 to 1-100.

Example:
    .. code-block:: bash

        PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 1


3. BranchLength_NumericConverter
---------------------------------

Description:
    To convert branch length values of a phylogenetic tree from string to numerical format.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.

Optional Parameters:
    - ``--decimal_place``   Specify the number of decimal places for branch lengths (default is 10).

Example:
    .. code-block:: bash

        PhyloTracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap [--decimal_place 10]


4. Phylo_Rooter
---------------

Description:
    Enables an accurate method for gene tree rooting and enhancing the downstream evolutionary genomic analysis.

Required Parameters:
    - ``--input_GF_list``       File containing paths to gene tree files, one per line.
    - ``--input_imap``          File with species classification information corresponding to genes.
    - ``--input_gene_length``   File with gene length information.
    - ``--input_sps_tree``      A species tree file in Newick format.

Example:
    .. code-block:: bash

        PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap --input_sps_tree sptree.nwk


5. OrthoFilter_LB
-----------------

Description:
    To prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length.

Required Parameters:
    - ``--input_GF_list``           File containing paths to gene tree files, one per line.
    - ``--input_imap``              File with species classification information corresponding to genes.
    - ``--rrbr_cutoff``  RRBR cutoff based on root-to-tip distance (default is 5).
    - ``--srbr_cutoff``  SRBR cutoff based on sister-relative branch ratio (default is 2.5).

Optional Parameters:
    - ``--visual``                  Visualize the results of gene family trees before and after removing long branches.

Example:
    .. code-block:: bash

        PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --rrbr_cutoff 5 --srbr_cutoff 2.5 [--visual]


6. OrthoFilter_Mono
-------------------

Description:
    To prune phylogenomic noise from both single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing well-established relationships unless otherwise required for specific objectives.

Required Parameters:
    - ``--input_GF_list``           File containing paths to gene tree files, one per line.
    - ``--input_taxa``              File with taxonomic information for species.
    - ``--input_imap``              File with species classification information corresponding to genes.
    - ``--input_sps_tree``          Species tree file in Newick format.

Optional Parameters:
    - ``--purity_cutoff``           Target purity for dominant lineage (default is 0.95).
    - ``--max_remove_fraction``     Maximum fraction of tips allowed to be removed (default is 0.5).
    - ``--visual``                  Visualize results before and after filtering.

Example:
    .. code-block:: bash

        PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_taxa gene2clade.imap --input_sps_tree sptree.nwk [--purity_cutoff 0.95 --max_remove_fraction 0.5 --visual]


7. TreeTopology_Summarizer
---------------------------

Description:
    To enumerate and visualize the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--input_imap``      File with species classification information corresponding to genes.

Optional Parameters:
    - ``--visual_top``      Number of top-ranked topologies to visualize (default is 10).

Example:
    .. code-block:: bash

        PhyloTracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--visual_top 10]


8. Tree_Visualizer
-------------------

Description:
    To mark tips of gene trees with provided tags, identify GD nodes, and integrate gene duplication results onto the species tree.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--input_imap``      File with species classification information corresponding to genes.
    - ``--keep_branch``     Whether to preserve branch length information (1 or 0).
    - ``--tree_style``      Tree style ('r' for rectangular, 'c' for circular).

Optional Parameters:
    - ``--gene_categories`` Taxonomic information for species.
    - ``--gene_family``     File with family classification information corresponding to genes.
    - ``--input_sps_tree``  Species tree file in Newick format (required with --gene_family).
    - ``--gene_expression`` Gene expression level files.
    - ``--visual_gd``       Visualize GD nodes of gene family trees.

Example:
    .. code-block:: bash

        PhyloTracer Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--gene_categories gene2order.imap --keep_branch 1 --tree_style r --gene_family gene2family.imap --input_sps_tree sptree.nwk --gene_expression gene2expression.imap --visual_gd]


9. GD_Detector
---------------

Description:
    To identify gene duplication events by reconciling gene trees with a species tree.

Required Parameters:
    - ``--input_GF_list``           File containing paths to gene tree files, one per line.
    - ``--input_imap``              File with species classification information corresponding to genes.
    - ``--gd_support``              GD node support (default is 50).
    - ``--subclade_support``        Subclade support of GD node (default is 50).
    - ``--dup_species_proportion``  Proportion of overlapping species for a GD event (default is 0.2).
    - ``--dup_species_num``         Number of species duplications under the GD node.
    - ``--input_sps_tree``          Species tree file in Newick format.
    - ``--deepvar``                 Maximum variance of depth (default is 1).

Optional Parameters:
    - ``--gdtype_mode``             GD type assignment mode: ``relaxed`` (default) or ``strict``.

Example:
    .. code-block:: bash

        PhyloTracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0.5 --dup_species_num 1 --input_sps_tree sptree.nwk --deepvar 1 [--gdtype_mode relaxed]


10. GD_Visualizer
------------------

Description:
    To visualize gene duplication detection results and integrate findings onto the species tree.

Required Parameters:
    - ``--input_sps_tree``  A numbered species tree file in Newick format.
    - ``--gd_result``       Result file from GD_Detector.
    - ``--input_imap``      File with species classification information corresponding to genes.

Example:
    .. code-block:: bash

        PhyloTracer GD_Visualizer --input_sps_tree sptree.nwk --gd_result gd_result.txt --input_imap gene2sps.imap


11. GD_Loss_Tracker
--------------------

Description:
    Analyze and summarize gene duplication loss events across nodes and tips in the species tree.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--input_sps_tree``  Species tree file in Newick format.
    - ``--input_imap``      File with species classification information corresponding to genes.

Optional Parameters:
    - ``--target_species``            Only count loss paths ending in this species (e.g., Arabidopsis_thaliana). Can be used multiple times.
    - ``--mrca_node``                 Only count loss paths passing through the MRCA of two species. Format: ``SpeciesA,SpeciesB`` (comma-separated, no space). Can be used multiple times.
    - ``--include_unobserved_species`` If set, species unobserved in a gene family are still classified by left/right presence instead of labeled as missing_data.

Example:
    .. code-block:: bash

        PhyloTracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--target_species Arabidopsis_thaliana --mrca_node SpeciesA,SpeciesB --include_unobserved_species]


12. GD_Loss_Visualizer
-----------------------

Description:
    To visualize the summary of gene duplication loss events on the context of species tree.

Required Parameters:
    - ``--gd_loss_result``  Detailed table generated by GD_Loss_Tracker (gd_loss_summary.txt).
    - ``--input_sps_tree``  A numbered species tree file in Newick format.

Example:
    .. code-block:: bash

        PhyloTracer GD_Loss_Visualizer --gd_loss_result gd_loss_summary.txt --input_sps_tree numbered_species_tree.nwk


13. Ortho_Retriever
--------------------

Description:
    To infer single-copy putative orthologs by splitting paralogs from large-scale gene family trees for multiple species.

Required Parameters:
    - ``--input_GF_list``       File containing paths to gene tree files, one per line.
    - ``--input_imap``          File with species classification information corresponding to genes.
    - ``--input_gene_length``   File with gene length information.

Example:
    .. code-block:: bash

        PhyloTracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap


14. Hybrid_Tracer
------------------

Description:
    To use the ABAB-BABA test to detect hybridization signals for each potential GD burst event across species tree detect species hybridization events for.

Required Parameters:
    - ``--input_GF_list``       File containing paths to gene tree files, one per line.
    - ``--input_Seq_GF_list``   File containing paths to sequence alignment files corresponding to the gene trees.
    - ``--input_imap``          File with species classification information corresponding to genes.
    - ``--input_sps_tree``      A species tree file in Newick format.

Optional Parameters:
    - ``--mrca_node``           Restrict Hybrid_Tracer to the MRCA of two species. Format: ``SpeciesA,SpeciesB`` (comma-separated, no space).
    - ``--split_groups``        Number of partitions for HYDE batch processing (default is 1).

Example:
    .. code-block:: bash

        PhyloTracer Hybrid_Tracer --input_GF_list GF_ID2path.imap --input_Seq_GF_list Seq_GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--mrca_node SpeciesA,SpeciesB --split_groups 1]


15. Hybrid_Visualizer
----------------------

Description:
    To visualize hybridization signals, highlighting support from gene tree topologies and D-statistic signals.

Required Parameters:
    - ``--hyde_out``        File containing the result of HyDe from Hybrid_Tracer.
    - ``--input_sps_tree``  A species tree file in Newick format.

Optional Parameters:
    - ``--node``            Node model, stack up all the heatmaps for each monophyletic clade respectively, only the squares in all heatmaps were light, the square after superimposition will be light.

Example:
    .. code-block:: bash

        PhyloTracer Hybrid_Visualizer --hyde_out hyde.out --input_sps_tree sptree.nwk [--node]


16. HaploFinder
----------------

Description:
   To distinguish gene conversion by tracing subgenome haplotypes through phylogenomic profiling.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--input_imap``      File with species classification information corresponding to genes.
    - ``--species_a``       Name of species A.
    - ``--species_b``       Name of species B.
    - ``--species_a_gff``   GFF file of species A.
    - ``--species_b_gff``   GFF file of species B.
    - ``--species_a_lens``  Lens file of species A.
    - ``--species_b_lens``  Lens file of species B.
    - ``--gd_support``      GD node support (default is 50).

Optional Parameters:
    - ``--mode``            Run mode: ``haplofinder`` (default) for GD analysis, ``split`` for FASTA partitioning by color labels.
    - ``--pair_support``    Minimum support of ortholog/speciation pair nodes (default is 50).
    - ``--visual_chr_a``    File containing chromosome numbers of species A for visualization.
    - ``--visual_chr_b``    File containing chromosome numbers of species B for visualization.
    - ``--size``            Size of points in the dotplot graph (default is 0.0005).

Example:
    .. code-block:: bash

        PhyloTracer HaploFinder --input_GF_list GF.list --input_imap gene2sps.imap --species_a A --species_b B --species_a_gff A.gff --species_b_gff B.gff --species_a_lens A.lens --species_b_lens B.lens --gd_support 50 [--mode haplofinder --pair_support 50 --visual_chr_a chr_a.txt --visual_chr_b chr_b.txt --size 0.0001]


For additional modules and detailed usage examples, refer to the relevant sections in this documentation.
