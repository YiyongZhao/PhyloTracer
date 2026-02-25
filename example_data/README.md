# 2) 依次跑 1-16 模块（示例数据）
# 1_Phylo_Rooter
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/1_Phylo_Rooter
PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap --input_sps_tree sptree.nwk

# 2_PhyloTree_CollapseExpand
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/2_PhyloTree_CollapseExpand
PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50

# 3_PhyloSupport_Scaler
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/3_PhyloSupport_Scaler
PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 1

# 4_BranchLength_NumericConverter
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/4_BranchLength_NumericConverter
PhyloTracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap --decimal_place 10

# 5_OrthoFilter_LB
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/5_OrthoFilter_LB
PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --absolute_branch_length 5 --relative_branch_length 2.5 --visual

# 6_OrthoFilter_Mono
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/6_OrthoFilter_Mono
PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_taxa gene2clade.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --purity_cutoff 0.95 --max_remove_fraction 0.5 --visual

# 7_TreeTopology_Summarizer
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/7_TreeTopology_Summarizer
PhyloTracer TreeTopology_Summarizer --input_GF_list gf.txt --input_imap imap --visual_top 10

# 8_Tree_Visualizer（我实际跑的是20棵代表树）
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/8_Tree_Visualizer
PhyloTracer Tree_Visualizer --input_GF_list GF_ID2path.visual20.imap --input_imap gene2sps.imap --gene_categories gene2family.imap gene2order.imap gene2clade.imap --keep_branch 1 --tree_style r --gene_family gene2family.imap --input_sps_tree sptree.nwk --gene_expression expression.csv --visual_gd

# 9_GD_Detector
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/9_GD_Detector
PhyloTracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --gd_support 50 --subclade_support 50 --dup_species_proportion 0 --dup_species_num 2 --input_sps_tree sptree.nwk --deepvar 1 --gdtype_mode relaxed

# 10_GD_Visualizer
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/10_GD_Visualizer
PhyloTracer GD_Visualizer --input_sps_tree sptree.nwk --gd_result gd_result.txt --input_imap /Users/apple/Documents/GitHub/PhyloTracer/example_data/9_GD_Detector/gene2sps.imap

# 11_GD_Loss_Tracker
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/11_GD_Loss_Tracker
PhyloTracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree /Users/apple/Documents/GitHub/PhyloTracer/example_data/9_GD_Detector/numed_sptree.nwk --input_imap gene2sps.imap --node_count_mode nonaccumulate

# 12_GD_Loss_Visualizer
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/12_GD_Loss_Visualizer
PhyloTracer GD_Loss_Visualizer --input_sps_tree numed_sptree.nwk --gd_loss_result gd_loss_summary.txt

# 13_Ortho_Retriever
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/13_Ortho_Retriever
PhyloTracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap

# 14_Hybrid_Tracer
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/14_Hybrid_Tracer
PhyloTracer Hybrid_Tracer --input_GF_list gf.txt --input_Seq_GF_list gf_aln.txt --input_sps_tree sptree.nwk --input_imap imap.list --split_groups 2

# 15_Hybrid_Visualizer（新接口，不用 --hyde_filtered_out）
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/15_Hybrid_Visualizer
PhyloTracer Hybrid_Visualizer --hyde_out hyde_out.txt --input_sps_tree sptree.nwk --node

# 16_HaploFinder
cd /Users/apple/Documents/GitHub/PhyloTracer/example_data/16_Haplofinder
PhyloTracer HaploFinder --mode haplofinder --input_GF_list gf.txt --input_imap imap.txt --input_sps_tree sptree.nwk --species_a ard --species_b ari --species_a_gff ard.gff --species_b_gff arh.gff --species_a_lens ard.lens --species_b_lens arh.lens --gd_support 50 --pair_support 50
