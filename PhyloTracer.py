import sys, textwrap
import argparse
import time
from Ortho_Split import *
from Phylo_Rooting import *
from Tree_Visualization import *
from Eliminate_PhyloNoise import *
from GD_Detector import *
from Statistical_Topology import *

from __init__ import *
##################################################################
print(textwrap.dedent("""
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
#    Please cite: Li et al. 2024, XXXX.                                                       #
#    Contacts: Taoli(Taoli@gmail.com); Yiyong Zhao(yzhao@bwh.harvard.edu)                     #
#                                                                                             #
###############################################################################################
"""))

if sys.version_info.major==2:
    print('You are using Python 2. Please upgrade to Python 3. PhyloRoot quit now...')
    quit() 

# create a parameter parser
class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action(self, action):
        if action.dest == 'command':
            # Override the format of the subparsers
            choices = self._metavar_formatter(action, action.choices)
            return f"available programs:\n{''.join(choices)}\n"
        return super()._format_action(action)

parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter, add_help=False)
subparsers = parser.add_subparsers(dest='command', help='available programs:')

# Tree_visualization command
tree_visualization_parser = subparsers.add_parser('Tree_Visualization', help='Tree Visualization help')
tree_visualization_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
tree_visualization_parser.add_argument('--input_imap', metavar='file', required=True, help='Input imap file')
tree_visualization_parser.add_argument('--gene_categories', metavar='file', nargs='+',  help='Gene category information')
tree_visualization_parser.add_argument('--keep_branch', type=int,  choices=[1, 0],help='[1/0] you can only input 1 or 0 Whether to preserve branch length information')
tree_visualization_parser.add_argument('--tree_style',  choices=['r', 'c'],default='r', help='Tree style: [r/c] (rectangular) or (circular) (default: rectangular)')

# Phylo_Rooting command
phylo_rooting_parser = subparsers.add_parser('Phylo_Rooting', help='Phylo Rooting help')
phylo_rooting_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
phylo_rooting_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
phylo_rooting_parser.add_argument('--input_gene_length', metavar='file',  help='Input gene length list')
phylo_rooting_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Ortho_Split command
ortho_split_parser = subparsers.add_parser('Ortho_Split', help='Ortho Split help')
ortho_split_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
ortho_split_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
ortho_split_parser.add_argument('--input_gene_length', metavar='file',  required=True, help='Input gene length list')
ortho_split_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Statistical_Topology command
statistical_Topology_parser = subparsers.add_parser('Statistical_Topology', help='Statistical_Topology help')
statistical_Topology_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
statistical_Topology_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')

# GD_Detector command
gd_detector_parser = subparsers.add_parser('GD_Detector', help='GD_Detector help')
gd_detector_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
gd_detector_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
gd_detector_parser.add_argument('--support', type=int,required=True, help='GD node support [50-100]')
gd_detector_parser.add_argument('--dup_species_radio', type=float ,required=True,help='The proportion of species with species duplications under the GD node [0-1]')
gd_detector_parser.add_argument('--dup_species_num', type=int ,required=True,help='The number of species with species duplications under the GD node')
gd_detector_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Eliminate_PhyloNoise command
eliminate_phyloNoise_parser = subparsers.add_parser('Eliminate_PhyloNoise', help='Eliminate_PhyloNoise help')
eliminate_phyloNoise_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
eliminate_phyloNoise_parser.add_argument('--input_taxa', metavar='file',  required=True, help='Input taxa file')

# Gene_Gain_And_Loss_Visualization command
Gene_Gain_And_Loss_Visualization_parser = subparsers.add_parser('Gene_Gain_And_Loss_Visualization', help='Gene_Gain_And_Loss_Visualization help')
Gene_Gain_And_Loss_Visualization_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')
Gene_Gain_And_Loss_Visualization_parser.add_argument('--input_summary_tree', metavar='file',  required=True, help='Input summary_species tree file')

parser.add_argument('-h', '--help', action='store_true', help=argparse.SUPPRESS)
# Analyze command line parameters


args = parser.parse_args()


def main():
    # Perform the corresponding functions according to the parameters
    if args.command == 'Tree_Visualization':
        # Execute the Tree_visualization function
        if args.input_GF_list and args.input_imap :
            start_time = time.time()
            os.makedirs(os.path.join(os.getcwd(), "pdf_result"))
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            tree_style = args.tree_style
            gene_categories = args.gene_categories
            keep_branch = args.keep_branch
            gene_category_list = [read_and_return_dict(i) for i in gene_categories]
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch,new_named_gene2gene_dic)
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
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            input_gene_length = args.input_gene_length
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            sptree = PhyloTree(input_sps_tree)
            sptree=rename_species_tree(sptree,voucher2taxa_dic)
            split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic,sptree,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Ortho_Split command are missing.")
            
    elif args.command == 'Statistical_Topology':
        # Execute the Statistical_Topology function
        if args.input_GF_list and args.input_imap :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            statistical_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Statistical_Topology command are missing.")


    elif args.command == 'GD_Detector':
        # Execute the GD_Detector function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.support and args.dup_species_radio and args.dup_species_num :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            support=args.support
            dup_species_percent = args.dup_species_radio
            dup_species_num = args.dup_species_num
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer(input_imap)
            sptree=PhyloTree(input_sps_tree)
            sptree=rename_species_tree(sptree, voucher2taxa_dic)
            sptree=num_tre_node(sptree)
            empty_count_dic= {node.name: 0 for node in sptree.traverse()}
            sps_tree=sptree.copy()
            tre_dic = read_and_return_dict(input_GF_list)
            empty_count_dic=batch_gfs_traverse(tre_dic, support, empty_count_dic,sptree,gene2new_named_gene_dic) 
            mark_sptree(sptree,empty_count_dic,voucher2taxa_dic)
            filename = 'result.txt'
            write_gene_duplication_events(filename, tre_dic, support,dup_species_percent, dup_species_num,sps_tree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for GD_Detector command are missing.")
            
    elif args.command == 'Eliminate_PhyloNoise':
        # Execute the Eliminate_PhyloNoise function
        if args.input_GF_list and args.input_taxa  :
            start_time = time.time()
            os.makedirs(os.path.join(os.getcwd(), "pruned_tree"))
            input_GF_list = args.input_GF_list
            input_taxa=args.input_taxa
            tre_dic = read_and_return_dict(input_GF_list)
            taxa_dic=read_and_return_dict(input_taxa)
            prune_main()
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Eliminate_PhyloNoise command are missing.")

    elif args.command == 'Gene_Gain_And_Loss_Visualization':
        # Execute the Gene_Gain_And_Loss_Visualization function
        if args.input_sps_tree and args.input_summary_tree :
            start_time = time.time()
            input_sps_tree = args.input_sps_tree
            input_summary_tree=args.input_summary_tree
            dic=get_gain_and_loss_dic(input_summary_tree)
            t=Tree(input_sps_tree)
            t.ladderize()
            t.sort_descendants("support")
            t.sort_descendants()
            ts=create_tree_style()
            mark_sptree(t,dic)
            t.render(file_name='gene_gain_loss.pdf',tree_style=ts)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Eliminate_PhyloNoise command are missing.")
            
    else:
        print("Usage: python PhyloTracer.py  [-h]  {Tree_Visualization, Phylo_Rooting, Ortho_Split, Statistical_Topology, GD_Detector, Eliminate_PhyloNoise}")
        print()
        print("optional arguments:")
        print('  -h, --help            show this help message and exit')
        print()
        print('available programs::')
        print('  {Tree_Visualization, Phylo_Rooting, Ortho_Split, Statistical_Topology, GD_Detector, Eliminate_PhyloNoise}')


if __name__ == "__main__":
    main()

