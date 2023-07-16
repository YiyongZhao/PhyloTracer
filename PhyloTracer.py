import sys, textwrap
import argparse
import time
from Ortho_Split import *
from Phylo_Rooting import *
from Tree_Visualization import *
##################################################################
print(textwrap.dedent("""\
 __                  __   ___  __             __   __   __
|__) |__| \ / |    /  \    |   |__)    /_\  /      |__  |__)
|      |  |  |   |__ \__/   |   |   \  /    \ \__   |__  |  \
                                             
"""))
print('############################################################\n\
PhyloTracer v1.1.1\n\
A bioinformatics software that uses RF distance, branch length, and gene length for rooting\n')

if sys.version_info.major==2:
    print('You are using Python 2. Please upgrade to Python 3. PhyloRoot quit now...')
    quit() 

# create a parameter parser
parser = argparse.ArgumentParser(description='PhyloTracer is a bioinformatics utility package for rooting and labeling duplicate genes in gene trees.')

# add parameters
# Tree_visualization
#usage:python PhyloTracer.py --Tree_visualization --input_GF_list GF --input_imap imap --gene_categories order family --keep_branch no --tree_style r
Tree_visualization_parser = parser.add_argument_group('Tree_visualization function arguments',description='>usage:python PhyloTracer.py --Tree_visualization --input_GF_list GF --input_imap imap --gene_categories order family --keep_branch no --tree_style r')
Tree_visualization_parser.add_argument('--Tree_visualization', action='store_true', help='Tree_visualization')
Tree_visualization_parser.add_argument('--tree_style', metavar='style', default='r', help='Tree style are "c" (ircular) or "r"(ectangular) (default: rectangular)')
Tree_visualization_parser.add_argument('--gene_categories', metavar='category', nargs='+', help='Gene category information')
Tree_visualization_parser.add_argument('--keep_branch', metavar='style', default='yes',help='Whether to preserve branch length information')

# Phylo_Rooting
#usage:python PhyloTracer.py --Phylo_Rooting --input_GF_list GF --input_imap imap --input_gene_length length --input_sps_tree sptree
Phylo_Rooting_parser = parser.add_argument_group('Phylo_Rooting function arguments',description='>usage:python PhyloTracer.py --Phylo_Rooting --input_GF_list GF --input_imap imap --input_gene_length length --input_sps_tree sptree')
Phylo_Rooting_parser.add_argument('--Phylo_Rooting', action='store_true', help='Phylo_Rooting')
Phylo_Rooting_parser.add_argument('--input_sps_tree', metavar='file', help='Input species tree file')


#Ortho_Split
#>python PhyloTracer.py --Ortho_Split --input_GF_list GF --input_imap imap --input_gene_length length --out_file_name test.txt
Ortho_Split_parser = parser.add_argument_group('Ortho_Split function arguments',description='>usage:python PhyloTracer.py --Ortho_Split --input_GF_list GF --input_imap imap --input_gene_length length --out_file_name test.txt')
Ortho_Split_parser.add_argument('--Ortho_Split', action='store_true', help='Ortho_Split')
Ortho_Split_parser.add_argument('--out_file_name', help='output file name')

# public parameters
parser.add_argument('--input_GF_list', metavar='file', help='Input gene tree list')
parser.add_argument('--input_imap', metavar='file', help='Input imap file')
parser.add_argument('--input_gene_length', metavar='file', help='Input gene length list')


# analyse command line parameters
args = parser.parse_args()

def main():
    # Perform the corresponding functions according to the parameters
    if args.Tree_visualization:
        # excute the Tree_visualization function
        if args.input_GF_list and args.input_imap and args.gene_categories and args.tree_style and args.keep_branch :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            tree_style = args.tree_style
            gene_categories = args.gene_categories
            keep_branch=args.keep_branch
            gene_category_list = [read_and_return_dict(i) for i in gene_categories]
            gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer(input_imap)
            tre_dic=read_and_return_dict(input_GF_list)
            view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic, gene_category_list,tree_style,keep_branch)
            end_time = time.time()
            execution_time = end_time - start_time
            print("程序执行时间：", execution_time, "秒")

    if args.Ortho_Split:
        # excute the Ortho_Split function
        if args.input_GF_list and args.input_imap and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_gene_length =args.input_gene_length
            out_file_name=args.out_file_name
            gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer(input_imap)
            tre_dic=read_and_return_dict(input_GF_list)
            len_dic=read_and_return_dict(input_gene_length)
            renamed_len_dic=rename_len_dic(len_dic,gene2new_named_gene_dic)
            # The code here calls the root function, using the parameters input_GF_list, input_imap, input_sps_tree, and input_gene_length.
            split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic,out_file_name)
            end_time = time.time()
            execution_time = end_time - start_time
            print("程序执行时间：", execution_time, "秒")
    
    if args.Phylo_Rooting:
        # excute the Phylo_Rooting function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            input_gene_length =args.input_gene_length
            gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer(input_imap)
            len_dic=read_and_return_dict(input_gene_length)
            renamed_len_dic=rename_len_dic(len_dic,gene2new_named_gene_dic)
            sptree=PhyloTree(input_sps_tree)
            tre_dic=read_and_return_dict(input_GF_list)   
            root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic,sptree)
            end_time = time.time()
            execution_time = end_time - start_time
            print("程序执行时间：", execution_time, "秒")

if __name__ == "__main__":
    main()
