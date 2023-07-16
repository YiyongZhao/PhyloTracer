import sys, textwrap
import argparse
import time
from Ortho_Split import *
from Phylo_Rooting import *
from Tree_Visualization import *
##################################################################
print(textwrap.dedent("""
###############################################################################################
#                                                                                             #
#                    ____  __          __    ______                                           #
#                   / __ \/ /_  __  __/ /___/_  __/________ _________  _____                  #
#                  / /_/ / __ \/ / / / / __ \/ / / ___/ __ `/ ___/ _ \/ ___/                  #
#                 / ____/ / / / /_/ / / /_/ / / / /  / /_/ / /__/  __/ /                      #
#                /_/   /_/ /_/\__, /_/\____/_/ /_/   \__,_/\___/\___/_/                       #
#                            /____/                                                           #
#                                                                                             #
#                                                                                             #        
#                                   PhyloTracer v1.1.1                                        #
#     A User-Friendly Toolkit for Gene Tree Rooting, Gene Duplication Identification,         #
#     Ortholog Retrieval, Species Hybridization Detection,and Visualization.                  #  
#                                                                                             #
#     Contacts: Taoli(l948777439@gmail.com); Yiyong Zhao(yiyongzhao1991@gmail.com)            #
#     Licence: GPL-3.0                                                                        #
#     Release Date: 2023-7                                                                    #
#                                                                                             #
###############################################################################################
"""))

if sys.version_info.major==2:
    print('You are using Python 2. Please upgrade to Python 3. PhyloRoot quit now...')
    quit() 

# create a parameter parser
import argparse

class CustomUsageFormatter(argparse.HelpFormatter):
    def _format_usage(self, usage, actions, groups, prefix):
        usage = super()._format_usage(usage, actions, groups, prefix)
        return usage.replace("[", "[-").replace("]", "]")

def setup_parser():
    parser = argparse.ArgumentParser(description='PhyloTracer is a bioinformatics utility package for rooting and labeling duplicate genes in gene trees.', formatter_class=CustomUsageFormatter)
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Tree_visualization command
    tree_visualization_parser = subparsers.add_parser('Tree_visualization', help='Perform tree visualization')
    tree_visualization_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
    tree_visualization_parser.add_argument('--input_imap', metavar='file', required=True, help='Input imap file')
    tree_visualization_parser.add_argument('--gene_categories', metavar='category', nargs='+', required=True, help='Gene category information')
    tree_visualization_parser.add_argument('--keep_branch', action='store_true', help='Whether to preserve branch length information')
    tree_visualization_parser.add_argument('--tree_style', metavar='style', default='rectangular', choices=['rectangular', 'circular'], help='Tree style: "rectangular" or "circular" (default: rectangular)')

    # Phylo_Rooting command
    phylo_rooting_parser = subparsers.add_parser('Phylo_Rooting', help='Perform phylo rooting')
    phylo_rooting_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
    phylo_rooting_parser.add_argument('--input_imap', metavar='file', required=True, help='Input imap file')
    phylo_rooting_parser.add_argument('--input_gene_length', metavar='file', required=True, help='Input gene length list')
    phylo_rooting_parser.add_argument('--input_sps_tree', metavar='file', required=True, help='Input species tree file')

    # Ortho_Split command
    ortho_split_parser = subparsers.add_parser('Ortho_Split', help='Perform ortho split')
    ortho_split_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
    ortho_split_parser.add_argument('--input_imap', metavar='file', required=True, help='Input imap file')
    ortho_split_parser.add_argument('--input_gene_length', metavar='file', required=True, help='Input gene length list')


    return parser

# Analyze command line parameters
parser = setup_parser()
args = parser.parse_args()

def main():
    # Perform the corresponding functions according to the parameters
    if args.command == 'Tree_visualization':
        # Execute the Tree_visualization function
        if args.input_GF_list and args.input_imap and args.gene_categories and args.tree_style and args.keep_branch:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            tree_style = args.tree_style
            gene_categories = args.gene_categories
            keep_branch = args.keep_branch
            gene_category_list = [read_and_return_dict(i) for i in gene_categories]
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Tree_visualization command are missing.")

    elif args.command == 'Phylo_Rooting':
        # Execute the Phylo_Rooting function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            input_gene_length = args.input_gene_length
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            sptree = PhyloTree(input_sps_tree)
            tre_dic = read_and_return_dict(input_GF_list)
            root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic, sptree)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Phylo_Rooting command are missing.")
        
    elif args.command == 'Ortho_Split':
        # Execute the Ortho_Split function
        if args.input_GF_list and args.input_imap and args.input_gene_length :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_gene_length = args.input_gene_length
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic, renamed_len_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Ortho_Split command are missing.")

    else:
        print("Invalid command.")


if __name__ == "__main__":
    main()
