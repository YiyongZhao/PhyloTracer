import sys, textwrap
import argparse
import time
from Phylo_Rooter import *
from PhyloNoise_Filter import *
from Tree_Visualizer import *
from GD_Detector import *
from GD_Visualizer import *
from Ortho_Retriever import *
from GeneDynamics_Tracker import *
from GeneDynamics_Visualizer import *
#from Hybrid_Tracer import *
#from Hybrid_Visualizer import *


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

# Tree_Visualizer command
Tree_Visualizer_parser = subparsers.add_parser('Tree_Visualizer', help='Tree_Visualizer help')
Tree_Visualizer_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
Tree_Visualizer_parser.add_argument('--input_imap', metavar='file', required=True, help='Input imap file')
Tree_Visualizer_parser.add_argument('--gene_categories', metavar='file', nargs='+',  help='Gene category information')
Tree_Visualizer_parser.add_argument('--keep_branch', type=int,  choices=[1, 0],help='[1/0] you can only input 1 or 0 Whether to preserve branch length information')
Tree_Visualizer_parser.add_argument('--tree_style',  choices=['r', 'c'],default='r', help='Tree style: [r/c] (rectangular) or (circular) (default: rectangular)')

# Phylo_Rooter command
Phylo_Rooter_parser = subparsers.add_parser('Phylo_Rooter', help='Phylo_Rooter help')
Phylo_Rooter_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
Phylo_Rooter_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
Phylo_Rooter_parser.add_argument('--input_gene_length', metavar='file',  help='Input gene length list')
Phylo_Rooter_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Ortho_Retriever command
Ortho_Retriever_parser = subparsers.add_parser('Ortho_Retriever', help='Ortho_Retriever help')
Ortho_Retriever_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
Ortho_Retriever_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
Ortho_Retriever_parser.add_argument('--input_gene_length', metavar='file',  required=True, help='Input gene length list')
Ortho_Retriever_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Statistical_Topology command
statistical_Topology_parser = subparsers.add_parser('Statistical_Topology', help='Statistical_Topology help')
statistical_Topology_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
statistical_Topology_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')

# GD_Detector command
GD_Detector_parser = subparsers.add_parser('GD_Detector', help='GD_Detector help')
GD_Detector_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
GD_Detector_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
GD_Detector_parser.add_argument('--support', type=int,required=True, help='GD node support [50-100]')
GD_Detector_parser.add_argument('--dup_species_radio', type=float ,required=True,help='The proportion of species with species duplications under the GD node [0-1]')
GD_Detector_parser.add_argument('--dup_species_num', type=int ,required=True,help='The number of species with species duplications under the GD node')
GD_Detector_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# PhyloNoise_Filter command
PhyloNoise_Filter_parser = subparsers.add_parser('PhyloNoise_Filter', help='PhyloNoise_Filter help')
PhyloNoise_Filter_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
PhyloNoise_Filter_parser.add_argument('--input_taxa', metavar='file',  required=True, help='Input taxa file')

# Gene_Gain_And_Loss_Visualization command
Gene_Gain_And_Loss_Visualization_parser = subparsers.add_parser('Gene_Gain_And_Loss_Visualization', help='Gene_Gain_And_Loss_Visualization help')
Gene_Gain_And_Loss_Visualization_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')
Gene_Gain_And_Loss_Visualization_parser.add_argument('--input_summary_tree', metavar='file',  required=True, help='Input summary_species tree file')

parser.add_argument('-h', '--help', action='store_true', help=argparse.SUPPRESS)
# Analyze command line parameters


args = parser.parse_args()


def main():
    # Perform the corresponding functions according to the parameters
    if args.command == 'Tree_Visualizer':
        # Execute the Tree_Visualizer function
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
            if args.gene_family and args.input_sps_tree:
                input_gene2fam = args.gene_family
                gene2fam = read_and_return_dict(input_gene2fam)
                input_sps_tree = args.input_sps_tree
                sptree = Tree(input_sps_tree)
                mark_gene_to_sptree_main(tre_dic,gene_category_list,sptree,gene2fam)
                view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch,new_named_gene2gene_dic,gene2fam)
            else:
                view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch,new_named_gene2gene_dic,gene2fam)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Tree_Visualizer command are missing.")

    elif args.command == 'Phylo_Rooter':
        # Execute the Phylo_Rooter function
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
            root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic, sptree,voucher2taxa_dic)
            end_time = time.time()
            end_time = time.time()
            execution_time = end_time - start_time
            print("Program execution time:", execution_time, "s")
        else:
            print("Required arguments for Phylo_Rooter command are missing.")
        
    elif args.command == 'Ortho_Retriever':
        # Execute the Ortho_Retriever function
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
            print("Required arguments for Ortho_Retriever command are missing.")
            
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
            
    elif args.command == 'PhyloNoise_Filter':
        # Execute the PhyloNoise_Filter function
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
            print("Required arguments for PhyloNoise_Filter command are missing.")

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
            print("Required arguments for Gene_Gain_And_Loss_Visualization command are missing.")
            
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

