import sys, textwrap
import argparse
import time
from PhyloTree_CollapseExpand import *
from PhyloSupport_Scaler import *
from BranchLength_NumericConverter import *
from Phylo_Rooter import *
from OrthoFilter import *
from TreeTopology_Summarizer import *
from Tree_Visualizer import *
from GD_Detector import *
from GD_Visualizer import *
from GD_Loss_Tracker import *
from GD_Loss_Visualizer import *
from Ortho_Retriever import *
from Gene_GainLoss_Finder import *
from Gene_GainLoss_Visualizer import *
from Hybrid_Tracer_ABAB_BABA import *
from Hybrid_Tracer_GCN import *
from Hybrid_Visualizer import *
from Phylo_Collapse import *
from Phylo_Collapse_Visualizer import*



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


# PhyloTree_CollapseExpand command
PhyloTree_CollapseExpand_parser = subparsers.add_parser('PhyloTree_CollapseExpand', help='PhyloTree_CollapseExpand help')
PhyloTree_CollapseExpand_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
PhyloTree_CollapseExpand_parser.add_argument('--support_value', type=int,required=True, help='If the support of the node is less than this value, it will be folded')

# PhyloSupport_Scaler command
PhyloSupport_Scaler_parser = subparsers.add_parser('PhyloSupport_Scaler', help='PhyloSupport_Scaler help')
PhyloSupport_Scaler_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
PhyloSupport_Scaler_parser.add_argument('--scale', type=str,  choices=['1', '0'],help='[1/0] you can only input 1 or 0 Whether to scale branch support information')

# BranchLength_NumericConverter command
BranchLength_NumericConverter_parser = subparsers.add_parser('BranchLength_NumericConverter', help='BranchLength_NumericConverter help')
BranchLength_NumericConverter_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
BranchLength_NumericConverter_parser.add_argument('--decimal_place', type=int, help='Set how many decimal places to keep')

# Phylo_Rooter command
Phylo_Rooter_parser = subparsers.add_parser('Phylo_Rooter', help='Phylo_Rooter help')
Phylo_Rooter_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
Phylo_Rooter_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
Phylo_Rooter_parser.add_argument('--input_gene_length', metavar='file',  help='Input gene length list')
Phylo_Rooter_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')
Phylo_Rooter_parser.add_argument('--threads', type=int, default=8, help='Thread Number')

# OrthoFilter command
OrthoFilter_parser = subparsers.add_parser('OrthoFilter', help='OrthoFilter help')
OrthoFilter_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
OrthoFilter_parser.add_argument('--input_taxa', metavar='file',  required=True, help='Input taxa file')
OrthoFilter_parser.add_argument('--long_branch_index', type=int, default=5, required=True, help='Long branch index')
OrthoFilter_parser.add_argument('--insert_branch_index', type=int, default=5, required=True, help='Insert_branch_index')

# TreeTopology_Summarizer command
TreeTopology_Summarizer_parser = subparsers.add_parser('TreeTopology_Summarizer', help='TreeTopology_Summarizer help')
TreeTopology_Summarizer_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
TreeTopology_Summarizer_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
TreeTopology_Summarizer_parser.add_argument('--outfile', metavar='file',  required=True, help='Out filename')

# Tree_Visualizer command
Tree_Visualizer_parser = subparsers.add_parser('Tree_Visualizer', help='Tree_Visualizer help')
Tree_Visualizer_parser.add_argument('--input_GF_list', metavar='file', required=True, help='Input gene tree list')
Tree_Visualizer_parser.add_argument('--input_imap', metavar='file', required=True, help='Input imap file')
Tree_Visualizer_parser.add_argument('--gene_categories', metavar='file', nargs='+',  help='Gene category information')
Tree_Visualizer_parser.add_argument('--keep_branch', type=str,  choices=['1', '0'],help='[1/0] you can only input 1 or 0 Whether to preserve branch length information')
Tree_Visualizer_parser.add_argument('--tree_style',  choices=['r', 'c'],default='r', help='Tree style: [r/c] (rectangular) or (circular) (default: rectangular)')
Tree_Visualizer_parser.add_argument('--gene_family', metavar='file',  required=False, help='Input species tree file')
Tree_Visualizer_parser.add_argument('--input_sps_tree', metavar='file',  required=False, help='Input species tree file')
Tree_Visualizer_parser.add_argument('--gene_expression', metavar='file',  required=False, help='gene_expression')

# GD_Detector command
GD_Detector_parser = subparsers.add_parser('GD_Detector', help='GD_Detector help')
GD_Detector_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
GD_Detector_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
GD_Detector_parser.add_argument('--gd_support', type=int,required=True, help='GD node support [50-100]')
GD_Detector_parser.add_argument('--clade_support', type=int,required=True, help='The children support of GD node [50-100]')
GD_Detector_parser.add_argument('--dup_species_radio', type=float ,required=True,help='The proportion of species with species duplications under the GD node [0-1]')
GD_Detector_parser.add_argument('--dup_species_num', type=int ,required=True,help='The number of species with species duplications under the GD node')
GD_Detector_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Ortho_Retriever command
Ortho_Retriever_parser = subparsers.add_parser('Ortho_Retriever', help='Ortho_Retriever help')
Ortho_Retriever_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
Ortho_Retriever_parser.add_argument('--input_imap', metavar='file',  required=True, help='Input imap file')
Ortho_Retriever_parser.add_argument('--input_gene_length', metavar='file',  required=True, help='Input gene length list')
Ortho_Retriever_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Phylo_Collapse command
Phylo_Collapse_parser = subparsers.add_parser('Phylo_Collapse', help='Phylo_Collapse help')
Phylo_Collapse_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='Input gene tree list')
Phylo_Collapse_parser.add_argument('--input_taxa', metavar='file',  required=True, help='Input taxa file')
Phylo_Collapse_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='Input species tree file')

# Phylo_Collapse_Visualizer command
Phylo_Collapse_Visualizer_parser = subparsers.add_parser('Phylo_Collapse_Visualizer', help='Phylo_Collapse_Visualizer help')

parser.add_argument('-h', '--help', action='store_true', help=argparse.SUPPRESS)
# Analyze command line parameters


args = parser.parse_args()

def format_time(seconds):
    days = seconds // (24 * 3600)
    hours = (seconds % (24 * 3600)) // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return f"{int(days)} d, {int(hours)} h, {int(minutes)} m, {seconds:.2f} s"

def main():
    if args.command == 'PhyloTree_CollapseExpand':
        # Execute the PhyloTree_CollapseExpand function
        if args.input_GF_list and args.support_value:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            support_value = args.support_value
            tre_dic = read_and_return_dict(input_GF_list)
            collapse_expand_main(tre_dic, support_value)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for PhyloTree_CollapseExpand command are missing.")


    elif args.command == 'PhyloSupport_Scaler':
        # Execute the PhyloSupport_Scaler function
        if args.input_GF_list and args.scale :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            scale = args.scale
            tre_dic = read_and_return_dict(input_GF_list)
            support_scaler_main(tre_dic,scale)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for PhyloSupport_Scaler command are missing.")


    elif args.command == 'BranchLength_NumericConverter':
        # Execute the BranchLength_NumericConverter function
        if args.input_GF_list:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            decimal_place=args.decimal_place
            tre_dic = read_and_return_dict(input_GF_list)
            branch_length_numeric_converter_main(tre_dic,decimal_place)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for BranchLength_NumericConverter command are missing.")


    elif args.command == 'Phylo_Rooter':
        # Execute the Phylo_Rooter function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            input_gene_length = args.input_gene_length
            num_threads=args.threads
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic,taxa2voucher_dic = gene_id_transfer(input_imap)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            sptree = PhyloTree(input_sps_tree)
            renamed_sptree=rename_input_tre(sptree,taxa2voucher_dic)
            tre_dic = read_and_return_dict(input_GF_list)
            root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic, renamed_sptree,voucher2taxa_dic,num_threads)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Phylo_Rooter command are missing.")
        
    elif args.command == 'OrthoFilter':
        # Execute the PhyloNoise_Filter function
        if args.input_GF_list and args.input_taxa and args.long_branch_index and args.insert_branch_index:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_taxa=args.input_taxa
            long_brancch_index=args.long_branch_index
            insert_branch_index=args.insert_branch_index
            tre_dic = read_and_return_dict(input_GF_list)
            taxa_dic=read_and_return_dict(input_taxa)
            prune_main(tre_dic,taxa_dic,long_brancch_index,insert_branch_index)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for OrthoFilter command are missing.")


    elif args.command == 'TreeTopology_Summarizer':
        # Execute the TreeTopology_Summarizer function
        if args.input_GF_list and args.input_imap :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            outfile=args.outfile
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            statistical_main(tre_dic,outfile,gene2new_named_gene_dic,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for TreeTopology_Summarizer command are missing.")


    # Perform the corresponding functions according to the parameters
    elif args.command == 'Tree_Visualizer':
        # Execute the Tree_Visualizer function
        if args.input_GF_list and args.input_imap :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            tree_style = args.tree_style
            gene_categories = args.gene_categories
            keep_branch = args.keep_branch
            gene_category_list = [read_and_return_dict(i) for i in gene_categories]
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            if args.gene_family and args.input_sps_tree:
                input_gene2fam = args.gene_family
                gene2fam = read_and_return_dict(input_gene2fam)
                input_sps_tree = args.input_sps_tree
                sptree = Tree(input_sps_tree)
                mark_gene_to_sptree_main(tre_dic,gene_category_list,sptree,gene2fam)
                view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch,new_named_gene2gene_dic,gene2fam)
            if args.gene_expression:
                gene2fam=None
                df=pd.read_excel(args.gene_expression, index_col=0) 
                view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch,new_named_gene2gene_dic,gene2fam,df)
            else:
                gene2fam=None
                view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, gene_category_list, tree_style, keep_branch,new_named_gene2gene_dic,gene2fam)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Tree_Visualizer command are missing.")

    
    elif args.command == 'GD_Detector':
        # Execute the GD_Detector function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.gd_support and args.clade_support and args.dup_species_radio and args.dup_species_num :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            gd_support=args.gd_support
            clade_support=args.clade_support
            dup_species_percent = args.dup_species_radio
            dup_species_num = args.dup_species_num
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            sptree=PhyloTree(args.input_sps_tree)
            num_tre_node(sptree)
            renamed_sptree=rename_species_tree(sptree, voucher2taxa_dic)
            tre_dic = read_and_return_dict(input_GF_list)
            filename = 'result.txt'
            write_gd_result(filename, tre_dic, gd_support,clade_support,dup_species_percent, dup_species_num,renamed_sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for GD_Detector command are missing.")


    elif args.command == 'Ortho_Retriever':
        # Execute the Ortho_Retriever function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            input_gene_length = args.input_gene_length
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            sptree = PhyloTree(input_sps_tree)
            sptree=rename_species_tree(sptree,voucher2taxa_dic)
            split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic,sptree,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Ortho_Retriever command are missing.")
            

    elif args.command == 'Phylo_Collapse':
        # Execute the PhyloNoise_Filter function
        if args.input_GF_list and args.input_taxa and args.input_sps_tree :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_taxa=args.input_taxa
            sptree=Tree(args.input_sps_tree)
            tre_dic = read_and_return_dict(input_GF_list)
            taxa_dic=read_and_return_dict(input_taxa)
            collapse_main(tre_dic,taxa_dic,sptree)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Phylo_Collapse command are missing.")


    elif args.command == 'Phylo_Collapse_Visualizer':
        start_time = time.time()
        if not os.path.exists('output/collapse_tree/'):
            print("Error: 'output/collapse_tree' folder does not exist in the current directory. Please perform Phylo_Collapse processing first")
        else:
            tre_dic = {i.split('.')[0]:'output/collapse_tree/'+i for i in os.listdir('output/collapse_tree/')}
            taxa_dic=read_and_return_dict('node2taxa.txt')
            c_color_dic=get_taxa_to_color_dict(taxa_dic)
            collapse_visual_main(tre_dic,c_color_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        #else:
            #print("Required arguments for Phylo_Collapse_Visualizer command are missing.")
     
    else:
        print("Usage: python PhyloTracer.py  [-h]  {GeneDynamics_Tracker, Hybrid_Visualizer, Ortho_Retriever, GD_Visualizer, GeneDynamics_Visualizer, Hybrid_Tracer, OrthoFilter_Multi, TreeTopology_Summarizer, Tree_Visualizer, GD_Loss_Tracker, Phylo_Collapse, GD_Detector, OrthoFilter_Single, Phylo_Rooter, Phylo_Tracer, GD_Loss_Visualizer, Phylo_Collapse_Visualizer}")
        print()
        print("optional arguments:")
        print('  -h, --help            show this help message and exit')
        print()
        print('available programs::')
        print('  {GeneDynamics_Tracker, Hybrid_Visualizer, Ortho_Retriever, GD_Visualizer, GeneDynamics_Visualizer, Hybrid_Tracer, OrthoFilter_Multi, TreeTopology_Summarizer, Tree_Visualizer, GD_Loss_Tracker, Phylo_Collapse, GD_Detector, OrthoFilter_Single, Phylo_Rooter, Phylo_Tracer, GD_Loss_Visualizer, Phylo_Collapse_Visualizer}')


if __name__ == "__main__":
    main()

