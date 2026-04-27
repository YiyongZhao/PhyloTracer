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

Optional Parameters:
    - ``--support_value``   Nodes whose support is strictly less than this value will be collapsed (default is 50).
    - ``--revert``          Revert the 'comb' structure to a fully resolved binary tree.
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - ``collapse_expand_tree/``  Directory containing the collapsed (or reverted) gene trees in Newick format.

Example:
    .. code-block:: bash

        PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap [--support_value 50] [--revert] [--output_dir DIR]


2. PhyloSupport_Scaler
-----------------------

Description:
    To recalibrate support value from bootstrap or posterior probability in a phylogenetic tree, scaling them between [0,1] and [1,100] ranges for computational compatibility, and vice versa to meet various analytical needs.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--scale_to``        Specify '1' to scale support values from 1-100 to 0-1, or '100' to scale from 0-1 to 1-100.

Optional Parameters:
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - ``support_scaler_tree/``  Directory containing gene trees with rescaled support values.

Example:
    .. code-block:: bash

        PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 1 [--output_dir DIR]


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
    Roots gene trees using a six-metric composite scoring framework guided by the species tree.
    All six metrics are normalized to [0, 1] (higher = better) via direction-aware min-max
    normalization before weighting, ensuring each metric contributes according to its assigned
    weight regardless of raw scale.

    Scoring metrics (all mapped to higher-is-better after normalization):

    - ``OD`` (Outgroup Depth, *lower raw value = better*): topological distance from species-tree
      root to the mapped outgroup clade. A smaller value means a more basal outgroup placement.
    - ``BLV`` (Branch Length Variance, *lower = better*): absolute difference in tip-to-root
      branch-length variance between the two root child-clades. Smaller values indicate more
      balanced, clock-like branch lengths after rooting.
    - ``GD`` (Gene Duplication count, *lower = better*): number of GD events inferred under the
      candidate root. Parsimony favours fewer duplications.
    - ``SO`` (Species Overlap at largest GD node, *higher = better*): Jaccard overlap of species
      sets in the two child-clades of the largest GD node. High overlap indicates a true
      duplication rather than a rooting artefact.
    - ``GDC`` (GD Consistency, *higher = better*): mean of (size_symmetry x Jaccard) across all
      GD nodes. Measures how consistently the GD nodes are biologically coherent.
    - ``MulRF`` (Multi-copy RF distance, *lower = better*): normalised Robinson-Foulds distance
      between the rooted gene tree and the species tree. Measures topological concordance.

Required Parameters:
    - ``--input_GF_list``       Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree per line.
    - ``--input_imap``          Two-column mapping file (gene_id<TAB>species_name).
    - ``--input_sps_tree``      Species tree file in Newick format.

Optional Parameters:
    - ``--weights OD BLV GD SO GDC RF``
                                Six space-separated float weights applied when
                                ``--weight-strategy empirical``; must sum to 1.0.
                                Default: ``0.30 0.10 0.30 0.10 0.10 0.10``.
    - ``--weight-strategy {empirical,entropy}``
                                Scoring weight strategy (default: ``empirical``).

                                - ``empirical``: uses the biologically-informed prior weights from
                                  ``--weights``; recommended when the number of candidate roots is
                                  small (< 10) where entropy estimates are unreliable.
                                - ``entropy``: Entropy Weight Method (EWM) -- derives weights
                                  automatically from the candidate distribution; assigns higher
                                  weight to metrics with greater variance among candidates.
                                  Recommended when many diverse candidates exist.
                                Note: when ``empirical`` is used and ``--weights`` is not provided,
                                the default weights are ``0.30 0.10 0.30 0.10 0.10 0.10``.
    - ``--output_dir``          Output directory (default: current working directory).

Output:
    - ``rooted_trees/``         Directory of rooted gene trees in Newick format.
    - ``stat_matrix.csv``       Per-tree scoring table with columns:
                                Tree, score, deep (OD), var (BLV), GD, species_overlap (SO),
                                gd_consistency (GDC), RF (MulRF).

Example:
    .. code-block:: bash

        # Default (empirical weights)
        PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk

        PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --weights 0.30 0.10 0.30 0.10 0.10 0.10

        # Entropy (data-driven) weighting
        PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --weight-strategy entropy


5. MulRF_Distance
------------------

Description:
    To compute species-level MulRF topological conflict distances between gene trees.
    Supports two modes:

    - **Mode 1** – Gene Tree vs Gene Tree: pairwise comparisons within one gene family list (self-comparisons excluded).
    - **Mode 2** – Gene Tree vs Species Tree: each gene tree is compared against a reference species tree.

    This is a topology-distance metric (not a sequence/genetic distance). Practical uses include quantifying gene-tree/species-tree discordance, comparing evolutionary pattern similarity among gene families, and supporting rooting decisions in Phylo_Rooter.

Required Parameters:
    - ``--mode``            Comparison mode: ``1`` = Gene Tree vs Gene Tree, ``2`` = Gene Tree vs Species Tree.
    - ``--input_GF_list``   Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line.
    - ``--input_imap``      Two-column mapping file (gene_id<TAB>species_name).

Optional Parameters:
    - ``--input_sps_tree``  Species tree file in Newick format (required when ``--mode 2``).
    - ``--output``          Output TSV filename (default: ``mulrf_distance.tsv``).

Output:
    A tab-separated file (``mulrf_distance.tsv``) with columns:

    - ``tre_id_1``, ``tre_id_2`` – tree identifiers for the compared pair.
    - ``gene_tree_1_leaf_count``, ``gene_tree_2_leaf_count`` – leaf counts.
    - ``gene_tree_1_species_count``, ``gene_tree_2_species_count`` – species counts.
    - ``shared_species_count`` – number of shared species between the two trees.
    - ``mulrf_distance`` – raw MulRF topological conflict distance.
    - ``maximum_possible_mulrf_distance`` – theoretical maximum for this pair.
    - ``normalized_mulrf_distance`` – ``mulrf_distance / maximum_possible_mulrf_distance``.
    - ``shared_species_bipartition_count`` – shared species-level bipartitions.
    - ``gene_tree_1_only_bipartition_count``, ``gene_tree_2_only_bipartition_count`` – bipartitions unique to each tree.

Example:
    .. code-block:: bash

        PhyloTracer MulRF_Distance --mode 1 --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--output mulrf_mode1.tsv]
        PhyloTracer MulRF_Distance --mode 2 --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk [--output mulrf_mode2.tsv]


6. OrthoFilter_LB
-----------------

Description:
    To prune phylogenomic noises from both single-copy and multi-copy gene family trees by removing the tips with long branch length.

Required Parameters:
    - ``--input_GF_list``           File containing paths to gene tree files, one per line.
    - ``--input_imap``              File with species classification information corresponding to genes.
    - ``--rrbr_cutoff``  RRBR cutoff based on root-to-tip distance (default is 5).
    - ``--srbr_cutoff``  SRBR cutoff based on sister-relative branch ratio (default is 2.5).

Optional Parameters:
    - ``--lb_mode``                 Long-branch decision mode: ``or`` (remove when RRBR OR SRBR exceeds cutoff) or ``and`` (remove only when both exceed cutoff). Default is ``or``.
    - ``--visual``                  Visualize the results of gene family trees before and after removing long branches.
    - ``--output_dir``              Output directory (default: current working directory).

Output:
    - ``orthofilter_lb/pruned_tree/``   Directory containing pruned gene trees.
    - ``orthofilter_lb/delete_gene/``    Directory containing lists of removed genes per family.

Example:
    .. code-block:: bash

        PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --rrbr_cutoff 5 --srbr_cutoff 2.5 [--lb_mode or] [--visual] [--output_dir DIR]


7. OrthoFilter_Mono
-------------------

Description:
    To prune phylogenomic noise from both single-copy and multi-copy gene family trees. It removes outliers and paralogs based on predefined taxonomic constraints (e.g., ensuring members from taxa such as families or orders form monophyletic groups). Caution: Groupings should be selected with care, prioritizing well-established relationships unless otherwise required for specific objectives.

Required Parameters:
    - ``--input_GF_list``           File containing paths to gene tree files, one per line.
    - ``--input_taxa``              Two-column mapping file (``species_id<TAB>clade_or_lineage_label``).
    - ``--input_imap``              File with species classification information corresponding to genes.
    - ``--input_sps_tree``          Species tree file in Newick format.

Optional Parameters:
    - ``--purity_cutoff``           Target purity for dominant lineage (default is 0.95).
    - ``--max_remove_fraction``     Maximum fraction of tips allowed to be removed (default is 0.5).
    - ``--visual``                  Visualize results before and after filtering.
    - ``--output_dir``              Output directory (default: current working directory).

Output:
    - ``orthofilter_mono/pruned_tree/``   Directory containing pruned gene trees.
    - ``orthofilter_mono/delete_gene/``    Directory containing lists of removed genes per family.

Example:
    .. code-block:: bash

        PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_taxa Clade.imap --input_sps_tree sptree.nwk [--purity_cutoff 0.95 --max_remove_fraction 0.5 --visual] [--output_dir DIR]


8. TreeTopology_Summarizer
---------------------------

Description:
    To enumerate and visualize the frequency of both absolute and relative topologies for single-copy gene trees or interested predefined clades.

Required Parameters:
    - ``--input_GF_list``   Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree per line.
    - ``--input_imap``      Two-column mapping file (gene_id<TAB>species_name).

Optional Parameters:
    - ``--visual_top``      Number of top-ranked topologies to visualize (default is 10).
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - ``absolute_topology_summary.txt``   Tab-delimited table listing each unique absolute topology and its frequency count.
    - ``relative_topology_summary.txt``   Tab-delimited table listing each unique relative topology and its frequency count.
    - Merged vector PDF files visualizing the top-ranked topologies.

Example:
    .. code-block:: bash

        PhyloTracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--visual_top 10] [--output_dir DIR]


9. Tree_Visualizer
-------------------

Description:
    To mark tips of gene trees with provided tags, identify GD nodes, and integrate gene duplication results onto the species tree.

Required Parameters:
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--input_imap``      File with species classification information corresponding to genes.

Optional Parameters:
    - ``--keep_branch``     Whether to preserve branch length information (``1`` = yes, ``0`` = no).
    - ``--tree_style``      Tree style (``r`` for rectangular, ``c`` for circular; default is ``r``).
    - ``--gene_categories`` One or more two-column files in format ``gene_id<TAB>category_label``; each file is one annotation layer (e.g., family/order/clade).
    - ``--input_sps_tree``  Species tree file in Newick format.
    - ``--heatmap_matrix``  Gene-associated numeric matrix file (recommended: .txt/.tsv tab-delimited; also supports .csv/.xls/.xlsx), genes as row index. Values should be in the range [0, 100]; normalize your data before input if needed.
    - ``--visual_gd``       Visualize GD nodes of gene family trees.
    - ``--gd_support``      Minimum support of a GD candidate node used by ``--visual_gd`` (range: 0-100, default is 50).
    - ``--subclade_support``  Minimum support required in GD child subclades used by ``--visual_gd`` (range: 0-100, default is 0).
    - ``--dup_species_proportion``  Minimum overlap ratio of duplicated species between GD child clades used by ``--visual_gd`` (range: 0-1, default is 0.2).
    - ``--dup_species_num``  Minimum number of overlapping duplicated species under a GD node used by ``--visual_gd`` (default is 2).
    - ``--deepvar``         Maximum tolerated depth-variance score used by ``--visual_gd`` (default is 1).
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - Per-gene-family annotated tree PDFs with tip labels, category color strips, and optional heatmap columns.
    - If ``--visual_gd`` is set, GD nodes are highlighted on the gene tree figures.
    - If ``--input_sps_tree`` is provided, a species-tree PDF annotated with family-level GD counts.

Example:
    .. code-block:: bash

        PhyloTracer Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap [--gene_categories Family.imap Order.imap Clade.imap --keep_branch 1 --tree_style r --input_sps_tree sptree.nwk --heatmap_matrix heatmap_matrix.txt --visual_gd --gd_support 50 --subclade_support 0 --dup_species_proportion 0.2 --dup_species_num 2 --deepvar 1] [--output_dir DIR]


10. GD_Detector
---------------

Description:
    To identify gene duplication events by reconciling gene trees with a species tree.

Required Parameters:
    - ``--input_GF_list``           File containing paths to gene tree files, one per line.
    - ``--input_imap``              File with species classification information corresponding to genes.
    - ``--gd_support``              GD node support (default is 50).
    - ``--subclade_support``        Subclade support of GD node (default is 0).
    - ``--dup_species_proportion``  Proportion of overlapping species for a GD event (default is 0.2).
    - ``--dup_species_num``         Minimum number of overlapping duplicated species under a GD node (default is 2).
    - ``--input_sps_tree``          Species tree file in Newick format.
    - ``--deepvar``                 Maximum variance of depth (default is 1).

Optional Parameters:
    - ``--gdtype_mode``             GD type assignment mode: ``relaxed`` (default) or ``strict``.
    - ``--output_dir``              Output directory (default: current working directory).

Output:
    - ``gd_result_relaxed.txt`` (or ``gd_result_strict.txt``)  Tab-delimited GD event table.
    - ``gd_type_relaxed.tsv`` (or ``gd_type_strict.tsv``)      GD type assignment table.
    - ``numed_sptree.nwk``                                       Numbered species tree for use with GD_Visualizer.

Example:
    .. code-block:: bash

        PhyloTracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 0 --dup_species_proportion 0.2 --dup_species_num 2 --input_sps_tree sptree.nwk --deepvar 1 [--gdtype_mode relaxed] [--output_dir DIR]


11. GD_Visualizer
------------------

Description:
    To visualize gene duplication detection results and integrate findings onto the species tree.

Required Parameters:
    - ``--input_sps_tree``  A numbered species tree file in Newick format (use ``numed_sptree.nwk`` from GD_Detector).
    - ``--gd_result``       Result file from GD_Detector.
    - ``--input_imap``      File with species classification information corresponding to genes.

Optional Parameters:
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - Species tree PDF(s) annotated with GD event counts on each node/branch.

Example:
    .. code-block:: bash

        PhyloTracer GD_Visualizer --input_sps_tree numed_sptree.nwk --gd_result gd_result_relaxed.txt --input_imap gene2sps.imap [--output_dir DIR]


12. GD_Loss_Tracker
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
    - ``--node_count_mode``           Node counting mode for ``path_count_*`` transition statistics: ``parsimony`` (default) collapses shared descendant losses to one ancestral count when supported; ``accumulate`` keeps all repeated descendant transitions.
    - ``--output_dir``                Output directory (default: current working directory).

Output:
    - ``gd_loss.csv``                 Standard GD-loss table. Columns begin with ``Tree ID``, ``GD ID``, ``GD Burst Node``, and ``Loss Node``; the loss node is defined as the first node where copy number decreases along ``Loss Path``. The table also contains per-species A/B copy assignments, presence flags, loss type, and path-derived transition summaries.
    - ``gd_loss_count_summary.txt``   Aggregated loss counts per species-tree node/tip.
    - ``gd_loss.xlsx``                Excel workbook combining the summary tables.

Example:
    .. code-block:: bash

        PhyloTracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--target_species Arabidopsis_thaliana --mrca_node SpeciesA,SpeciesB --include_unobserved_species --node_count_mode parsimony] [--output_dir DIR]


13. GD_Loss_Visualizer
-----------------------

Description:
    To visualize the summary of gene duplication loss events on the context of species tree.

Required Parameters:
    - ``--gd_loss_result``  CSV table generated by GD_Loss_Tracker (gd_loss.csv).
    - ``--input_sps_tree``  A numbered species tree file in Newick format (required for species tree layout).

Optional Parameters:
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - Species tree PDF annotated with per-node and per-tip loss counts (colored by loss category: 2-2, 2-1, 2-0).
    - Example outputs may be stored in paired ``parsimony/`` and
      ``accumulate/`` folders so each PDF stays next to the exact
      ``gd_loss.csv`` and ``gd_loss_node_summary.tsv`` used to render it.

Example:
    .. code-block:: bash

        PhyloTracer GD_Loss_Visualizer --gd_loss_result gd_loss.csv --input_sps_tree numed_sptree.nwk [--output_dir DIR]


14. Ortho_Retriever
--------------------

Description:
    To infer single-copy putative orthologs by splitting paralogs from large-scale gene family trees for multiple species. Optional synteny support can be applied, and strict sister-clade outgroup genes can also be attached to the written trees.

Required Parameters:
    - ``--input_GF_list``       File containing paths to gene tree files, one per line.
    - ``--input_imap``          File with species classification information corresponding to genes.
    - ``--input_gene_length``   File with gene length information.

Optional Parameters:
    - ``--input_synteny_blocks``  Optional raw synteny block file used to refine candidate ortholog sets.
    - ``--add_outgroup``          Attach strict sister-clade outgroup genes to output trees. Outgroup species must not overlap with ingroup species and must be single-copy within the candidate sister clade.
    - ``--output_dir``          Output directory (default: current working directory).

Output:
    - ``ortho_retriever_summary.txt``  Summary of ortholog extraction results per gene family.
    - ``ortholog_trees.tsv``           Tab-delimited list of extracted single-copy ortholog subtrees.
    - ``ortholog_synteny_report.tsv``  Written only when ``--input_synteny_blocks`` is used.
    - ``ortholog_outgroup_report.tsv`` Written only when ``--add_outgroup`` is used; records selected outgroup genes and whether each tree received a valid strict outgroup.

Example:
    .. code-block:: bash

        PhyloTracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap [--input_synteny_blocks collinearity] [--add_outgroup] [--output_dir DIR]


15. Hybrid_Tracer
------------------

Description:
    To detect hybridization signals from GD-derived gene sets using coalescent-based phylogenetic invariants (HyDe) on species-tree mapped GD nodes.

Required Parameters:
    - ``--input_GF_list``       File containing paths to gene tree files, one per line.
    - ``--input_Seq_GF_list``   File containing paths to sequence alignment files corresponding to the gene trees.
    - ``--input_imap``          File with species classification information corresponding to genes.
    - ``--input_sps_tree``      A species tree file in Newick format.

Optional Parameters:
    - ``--mrca_node``           Restrict Hybrid_Tracer to the MRCA of two species. Format: ``SpeciesA,SpeciesB`` (comma-separated, no space).
    - ``--split_groups``        Number of partitions for HYDE batch processing (default is 1).
    - ``--output_dir``          Output directory (default: current working directory).

Output:
    - ``hyde_out.txt``            Full HyDe output table with columns: P1, Hybrid, P2, Zscore, Pvalue, Gamma, and site-pattern counts (AAAA, AAAB, AABA, AABB, AABC, ABAA, ABAB, ABAC, ABBA, BAAA, ABBC, CABC, BACA, BCAA, ABCD).
    - ``hyde_filtered_out.txt``   Filtered HyDe results (significant tests only).

Example:
    .. code-block:: bash

        PhyloTracer Hybrid_Tracer --input_GF_list GF_ID2path.imap --input_Seq_GF_list Seq_GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap [--mrca_node SpeciesA,SpeciesB --split_groups 1] [--output_dir DIR]

    Outputs: ``hyde_out.txt``, ``hyde_filtered_out.txt``, and ``hyde_summary.txt``.


16. Hybrid_Visualizer
----------------------

Description:
    To visualize hybridization signals, highlighting hybridization proportions (gamma) from HyDe analysis.

Required Parameters:
    - ``--hyde_out``        File containing the result of HyDe from Hybrid_Tracer.
    - ``--input_sps_tree``  A species tree file in Newick format.

Optional Parameters:
    - ``--node``            Node model, stack up all the heatmaps for each monophyletic clade respectively, only the squares in all heatmaps were light, the square after superimposition will be light.
    - ``--output_dir``      Output directory (default: current working directory).

Output:
    - Heatmap PDF(s) showing hybridization gamma values across species pairs, annotated on the species tree.

Example:
    .. code-block:: bash

        PhyloTracer Hybrid_Visualizer --hyde_out hyde.out --input_sps_tree sptree.nwk [--node] [--output_dir DIR]

    Combined figure legend: red = focal hybrid/clade, blue = target internal node label in ``--node`` mode, yellow = γ values, white = tested hybridization combinations.


17. HaploFinder
----------------

Description:
   To distinguish gene conversion by tracing subgenome haplotypes through phylogenomic profiling.

Required Parameters (haplofinder mode):
    - ``--input_GF_list``   File containing paths to gene tree files, one per line.
    - ``--input_imap``      File with species classification information corresponding to genes.
    - ``--input_sps_tree``  Species tree file in Newick format.
    - ``--species_a``       Name of species A.
    - ``--species_b``       Name of species B.
    - ``--species_a_gff``   GFF file of species A.
    - ``--species_b_gff``   GFF file of species B.
    - ``--species_a_lens``  Lens file of species A.
    - ``--species_b_lens``  Lens file of species B.

Optional Parameters (haplofinder mode):
    - ``--gd_support``      Minimum support of GD nodes used for pair extraction (range: 0-100, default is 50).
    - ``--pair_support``    Minimum support of ortholog/speciation pair nodes (range: 0-100, default is 50).
    - ``--visual_chr_a``    File containing chromosome numbers of species A for visualization.
    - ``--visual_chr_b``    File containing chromosome numbers of species B for visualization.
    - ``--size``            Size of points in the dotplot graph (default is 0.0005).
    - ``--output_dir``      Output directory (default: current working directory).

Split mode:
    When ``--mode split`` is used, HaploFinder partitions FASTA sequences by subgenome color labels.

    Required parameters for split mode:
        - ``--input_GF_list``   Tab-delimited mapping file (GF_ID<TAB>gene_tree_path).
        - ``--input_imap``      Two-column mapping file (gene_id<TAB>species_name).
        - ``--input_fasta``     Input FASTA file (.fa/.fasta).
        - ``--cluster_file``    Split-mode cluster metadata file.
        - ``--hyb_sps``         Hybrid species name used for subgenome assignment.
        - ``--parental_sps``    Parental species names (quoted, space-separated string).
        - ``--species_b_gff``   Genome annotation file for species B.

    Optional parameters (both modes):
        - ``--mode``            Run mode: ``haplofinder`` (default) or ``split``.
        - ``--output_dir``      Output directory (default: current working directory).

Output (haplofinder mode):
    - ``gd_pairs_dotplot.pdf``    Dotplot PDF showing GD-pair chromosomal positions.
    - ``gd_pairs_dotplot.png``    Dotplot PNG version.
    - ``color_label.txt``         Color-label assignments for each gene.

Output (split mode):
    - ``haplofinder_split/split_assignment.tsv``        Subgenome assignment table.
    - ``haplofinder_split/split_subgenome_A.fasta``     FASTA for subgenome A.
    - ``haplofinder_split/split_subgenome_B.fasta``     FASTA for subgenome B.
    - ``haplofinder_split/split_subgenome_unknown.fasta``  FASTA for unassigned genes.
    - ``haplofinder_split/split_summary.txt``            Summary statistics.

Example:
    .. code-block:: bash

        PhyloTracer HaploFinder --input_GF_list GF.list --input_imap gene2sps.imap --input_sps_tree sptree.nwk --species_a A --species_b B --species_a_gff A.gff --species_b_gff B.gff --species_a_lens A.lens --species_b_lens B.lens [--gd_support 50 --pair_support 50 --visual_chr_a chr_a.txt --visual_chr_b chr_b.txt --size 0.0001] [--output_dir DIR]

        PhyloTracer HaploFinder --mode split --input_GF_list gf.txt --input_imap gene2sps.imap --input_fasta proteins.fa --cluster_file cluster.tsv --hyb_sps Hybrid --parental_sps "P1 P2" --species_b_gff B.gff [--output_dir DIR]


For additional modules and detailed usage examples, refer to the relevant sections in this documentation.
