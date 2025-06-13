import sys, textwrap
import argparse
import time
from PhyloTree_CollapseExpand import *
from PhyloSupport_Scaler import *
from BranchLength_NumericConverter import *
from Phylo_Rooter import *
from OrthoFilter_LB import *
from OrthoFilter_Mono import *
from TreeTopology_Summarizer import *
from Tree_Visualizer import *
from GD_Detector import *
from GD_Visualizer import *
from GD_Loss_Tracker import *
from GD_Loss_Visualizer import *
from Ortho_Retriever import *
from Hybrid_Tracer import *
from Hybrid_Visualizer import *
from HaploFinder import *



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

# Check Python version
if sys.version_info.major == 2:
    print('You are using Python 2. Please upgrade to Python 3. PhyloTracer quit now...')
    sys.exit()

# Custom help formatter for better formatting of subcommands and section titles
class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action(self, action):
        if action.dest == 'command' and action.choices:
            # Format the list of available subcommands
            choices = ', '.join(action.choices)
            return f"Available programs:\n  {choices}\n"
        return super()._format_action(action)

    def add_usage(self, usage, actions, groups, prefix=None):
        # Customize the "usage" section prefix
        if prefix is None:
            prefix = "Usage: "  # Change "usage" to "Usage"
        super().add_usage(usage, actions, groups, prefix)

    def start_section(self, heading):
        # Customize section titles
        if heading == "options":
            heading = "Options"  # Change "options" to "Options"
        super().start_section(heading)

# Create the main parser
parser = argparse.ArgumentParser(prog='PhyloTracer',description='A toolkit for phylogenetic tree analysis and visualization.',formatter_class=argparse.RawTextHelpFormatter,add_help=False)

# Add subparsers for commands
subparsers = parser.add_subparsers(dest='command',help='Available programs for phylogenetic tree processing')



# PhyloTree_CollapseExpand command
PhyloTree_CollapseExpand_parser = subparsers.add_parser('PhyloTree_CollapseExpand', help='Transform a phylogenetic tree into a "comb" structure based on a support value or revert it to a binary tree',formatter_class=CustomHelpFormatter)
PhyloTree_CollapseExpand_parser.add_argument('--input_GF_list', metavar='FILE', required=True, help='Path to a file containing gene tree file paths, one per line')
PhyloTree_CollapseExpand_parser.add_argument('--support_value', type=int, required=True, help='Support value threshold: nodes with support <= this value will be transformed')
PhyloTree_CollapseExpand_parser.add_argument('--revert', action='store_true', help='Revert the "comb" structure back to a fully resolved binary tree')

# PhyloSupport_Scaler command
PhyloSupport_Scaler_parser = subparsers.add_parser('PhyloSupport_Scaler', help='Scale support values of gene trees', formatter_class=CustomHelpFormatter)
PhyloSupport_Scaler_parser.add_argument('--input_GF_list', metavar='FILE', required=True, help='Path to a file containing gene tree file paths, one per line')
PhyloSupport_Scaler_parser.add_argument('--scale_to', type=str, choices=['1', '100'], help='Scale support values: "1" to scale from 1-100 to 0-1, or "100" to scale from 0-1 to 1-100')

# BranchLength_NumericConverter command
BranchLength_NumericConverter_parser = subparsers.add_parser('BranchLength_NumericConverter', help='Convert branch length values to a specified number of decimal places',formatter_class=CustomHelpFormatter)
BranchLength_NumericConverter_parser.add_argument('--input_GF_list', metavar='file', required=True, help='File containing paths to gene tree files, one per line')
BranchLength_NumericConverter_parser.add_argument('--decimal_place', type=int, help='Return the branch length values to 10 decimal places and default = 10')

# Phylo_Rooter command
Phylo_Rooter_parser = subparsers.add_parser('Phylo_Rooter', help='Phylo_Rooter help',formatter_class=CustomHelpFormatter)
Phylo_Rooter_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
Phylo_Rooter_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
Phylo_Rooter_parser.add_argument('--input_gene_length', metavar='file',  help='File with information corresponding to gene length')
Phylo_Rooter_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='A species tree file with Newick format')

# OrthoFilter_LB command
OrthoFilter_LB_parser = subparsers.add_parser('OrthoFilter_LB', help='OrthoFilter_LB help',formatter_class=CustomHelpFormatter)
OrthoFilter_LB_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
OrthoFilter_LB_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
OrthoFilter_LB_parser.add_argument('--absolute_branch_length', type=int, default=5, required=True, help='Absolute branch length multiplier and default = 5')
OrthoFilter_LB_parser.add_argument('--relative_branch_length', type=float, default=2.5, required=True, help='Relative branch length multiplier and default = 2.5')
OrthoFilter_LB_parser.add_argument('--visual', action='store_true', help='Visualize the results of gene family trees before and after removing long branches')

# OrthoFilter_Mono command
OrthoFilter_Mono_parser = subparsers.add_parser('OrthoFilter_Mono', help='OrthoFilter_Mono help',formatter_class=CustomHelpFormatter)
OrthoFilter_Mono_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
OrthoFilter_Mono_parser.add_argument('--input_taxa', metavar='file',  required=True, help='Input taxa file')
OrthoFilter_Mono_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
OrthoFilter_Mono_parser.add_argument('--branch_length_multiples', type=int, default=5, required=True, help='Branch_length_multiples and default = 10')
OrthoFilter_Mono_parser.add_argument('--insert_branch_index', type=int, default=5, required=True, help='Insert_branch_index and default = 10')
OrthoFilter_Mono_parser.add_argument('--visual', action='store_true', help='Visualize the results of gene family trees before and after removing long branches')

# TreeTopology_Summarizer command
TreeTopology_Summarizer_parser = subparsers.add_parser('TreeTopology_Summarizer', help='TreeTopology_Summarizer help',formatter_class=CustomHelpFormatter)
TreeTopology_Summarizer_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
TreeTopology_Summarizer_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
TreeTopology_Summarizer_parser.add_argument('--visual_top', type=int,  required=False, default=10,help='Visualize the result of provide number,default=10')

# Tree_Visualizer command
Tree_Visualizer_parser = subparsers.add_parser('Tree_Visualizer', help='Tree_Visualizer help',formatter_class=CustomHelpFormatter)
Tree_Visualizer_parser.add_argument('--input_GF_list', metavar='file', required=True, help='File containing paths to gene tree files, one per line')
Tree_Visualizer_parser.add_argument('--input_imap', metavar='file', required=True, help='File with classification information of species corresponding to genes')
Tree_Visualizer_parser.add_argument('--gene_categories', metavar='file', nargs='+',  help='File with taxonomic information for species')
Tree_Visualizer_parser.add_argument('--keep_branch', type=str,  choices=['1', '0'],help='1 or 0 indicates whether or not to preserve branch length information')
Tree_Visualizer_parser.add_argument('--tree_style',  choices=['r', 'c'],default='r', help='The tree style, r is meaning rectangular, c is meaning circular')
Tree_Visualizer_parser.add_argument('--gene_family', metavar='file',  required=False, help='If you want to mark gene families you need to provide this file')
Tree_Visualizer_parser.add_argument('--input_sps_tree', metavar='file',  required=False, help='A species tree file with Newick format')
Tree_Visualizer_parser.add_argument('--gene_expression', metavar='file',  required=False, help='Gene expression level files')
Tree_Visualizer_parser.add_argument('--visual_gd', action='store_true', help='Visualize the gd node of gene family trees')

# GD_Detector command
GD_Detector_parser = subparsers.add_parser('GD_Detector', help='GD_Detector help',formatter_class=CustomHelpFormatter)
GD_Detector_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
GD_Detector_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
GD_Detector_parser.add_argument('--gd_support', type=int,required=True, help='GD node support [50-100]')
GD_Detector_parser.add_argument('--subclade_support', type=int,required=True, help='The subclade support of GD node [0-100]')
GD_Detector_parser.add_argument('--dup_species_proportion', type=float ,required=True,help='The proportion of overlappped species from two subclade for a GD event with range [0-1] and default = 0.2')
GD_Detector_parser.add_argument('--dup_species_num', type=int ,required=True,help='The number of species with species duplications under the GD node')
GD_Detector_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='A species tree file with Newick format')
GD_Detector_parser.add_argument('--deepvar',type=int,required=True, help='Maximum variance of deepth and default = 1')

# GD_Visualizer command
GD_Visualizer_parser = subparsers.add_parser('GD_Visualizer', help='GD_Visualizer help',formatter_class=CustomHelpFormatter)
GD_Visualizer_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='A numbered species tree file with Newick format')
GD_Visualizer_parser.add_argument('--gd_result', metavar='file',  required=True, help='Result file of GD_Detecto')
GD_Visualizer_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')

# GD_Loss_Tracker command
GD_Loss_Tracker_parser = subparsers.add_parser('GD_Loss_Tracker', help='GD_Loss_Tracker help',formatter_class=CustomHelpFormatter)
GD_Loss_Tracker_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
GD_Loss_Tracker_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='A species tree file with Newick format')
GD_Loss_Tracker_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')


# GD_Loss_Visualizer command
GD_Loss_Visualizer_parser = subparsers.add_parser('GD_Loss_Visualizer', help='GD_Loss_Visualizer help',formatter_class=CustomHelpFormatter)
GD_Loss_Visualizer_parser.add_argument('--gd_loss_result', metavar='file',  required=True, help='Result file of gd loss count summary of GD_Loss_Tracker')
GD_Loss_Visualizer_parser.add_argument('--input_sps_tree', metavar='file',  required=False, help='A numbered species tree file with Newick format')

# Ortho_Retriever command
Ortho_Retriever_parser = subparsers.add_parser('Ortho_Retriever', help='Ortho_Retriever help',formatter_class=CustomHelpFormatter)
Ortho_Retriever_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
Ortho_Retriever_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
Ortho_Retriever_parser.add_argument('--input_gene_length', metavar='file',  required=True, help='File with information corresponding to gene length')

# Hybrid_Tracer
Hybrid_Tracer_parser = subparsers.add_parser('Hybrid_Tracer', help='Hybrid_Tracer help',formatter_class=CustomHelpFormatter)
Hybrid_Tracer_parser.add_argument('--input_GF_list', metavar='file',  required=True, help='File containing paths to gene tree files, one per line')
Hybrid_Tracer_parser.add_argument('--input_Seq_GF_list', metavar='file',  required=True, help='File containing paths to sequence alignment files corresponding to the gene trees')
Hybrid_Tracer_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='A species tree file with Newick format')
Hybrid_Tracer_parser.add_argument('--input_imap', metavar='file',  required=True, help='File with classification information of species corresponding to genes')
Hybrid_Tracer_parser.add_argument('--target_node',  metavar='file', required=False, help='Specific node to process. Use "all" to process all gd_names.')
Hybrid_Tracer_parser.add_argument('--split_gd', action='store_true', help='Split the gd to run hyde')

# Hybrid_Visualizer
Hybrid_Visualizer_parser = subparsers.add_parser('Hybrid_Visualizer', help='Hybrid_Visualizer help',formatter_class=CustomHelpFormatter)
Hybrid_Visualizer_parser.add_argument('--hyde_out', metavar='file',  required=True, help='File containing result of hyde')
Hybrid_Visualizer_parser.add_argument('--input_sps_tree', metavar='file',  required=True, help='A species tree file with Newick format')
Hybrid_Visualizer_parser.add_argument('--node', action="store_true", default=False, help="Node model, stack up all the heatmaps for each monophyletic clade respectively, only the squares in all heatmaps were light, the square after superimposition will be light")

#HaploFinder
HaploFinder = subparsers.add_parser('HaploFinder', help='HaploFinder help',formatter_class=CustomHelpFormatter)
HaploFinder.add_argument('--input_GF_list', metavar='FILE', required=True, help='File containing paths to gene tree files, one per line file')
HaploFinder.add_argument('--input_imap', metavar='FILE', required=True, help='File with classification information of species corresponding to genes')
HaploFinder.add_argument('--species_a', type=str, required=True, help='Name of species A')
HaploFinder.add_argument('--species_b', type=str, required=True, help='Name of species B')
HaploFinder.add_argument('--species_a_gff', metavar='FILE', required=True, help='GFF file of species A')
HaploFinder.add_argument('--species_b_gff', metavar='FILE', required=True, help='GFF file of species B')
HaploFinder.add_argument('--species_a_lens', metavar='FILE', required=True, help='Lens file of species A')
HaploFinder.add_argument('--species_b_lens', metavar='FILE', required=True, help='Lens file of species B')
HaploFinder.add_argument('--visual_chr_a', metavar='FILE', required=False, help='A file containing the chromosome numbers of species a')
HaploFinder.add_argument('--visual_chr_b', metavar='FILE', required=False, help='A file containing the chromosome numbers of species b')
HaploFinder.add_argument('--gd_support', type=int,required=True, default=50,help='GD node support [50-100]')
HaploFinder.add_argument('--pair_support', type=int,required=True, default=50,help='gene pair support [50-100]')
HaploFinder.add_argument('--size', type=float, required=False, help='The size of each point in the dopolot graph and default = 0.0005')
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
            collapse_expand_main(tre_dic, support_value,revert=args.revert)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for PhyloTree_CollapseExpand command are missing.")


    elif args.command == 'PhyloSupport_Scaler':
        # Execute the PhyloSupport_Scaler function
        if args.input_GF_list and args.scale_to :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            scale = args.scale_to
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
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic,taxa2voucher_dic = gene_id_transfer(input_imap)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            sptree = PhyloTree(input_sps_tree)
            renamed_sptree=rename_input_tre(sptree,taxa2voucher_dic)
            tre_dic = read_and_return_dict(input_GF_list)
            root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic, renamed_sptree,voucher2taxa_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Phylo_Rooter command are missing.")
        
    elif args.command == 'OrthoFilter_LB':
        # Execute the OrthoFilter_LB function
        if args.input_GF_list and args.input_imap and args.absolute_branch_length and args.relative_branch_length :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            absolute_branch_length=args.absolute_branch_length-1
            relative_branch_length=args.relative_branch_length-1
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic,taxa2voucher_dic = gene_id_transfer(input_imap)
            
            tre_dic = read_and_return_dict(input_GF_list)
            prune_main_LB(tre_dic,voucher2taxa_dic,gene2new_named_gene_dic, new_named_gene2gene_dic,absolute_branch_length,relative_branch_length,visual=args.visual)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for OrthoFilter_LB command are missing.")


    elif args.command == 'OrthoFilter_Mono':
        # Execute the OrthoFilter_Mono function
        if args.input_GF_list and args.input_taxa and args.branch_length_multiples and args.insert_branch_index:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_taxa=args.input_taxa
            input_imap = args.input_imap
            long_branch_index=args.branch_length_multiples
            insert_branch_index=args.insert_branch_index
            gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic,taxa2voucher_dic = gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            taxa_dic=read_and_return_dict(input_taxa)
            prune_main_Mono(tre_dic,taxa_dic,long_branch_index,insert_branch_index,new_named_gene2gene_dic,gene2new_named_gene_dic,visual=args.visual)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for OrthoFilter_Mono command are missing.")


    elif args.command == 'TreeTopology_Summarizer':
        # Execute the TreeTopology_Summarizer function
        if args.input_GF_list and args.input_imap :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            top_n=args.visual_top if args.visual_top else None
            outfile='topology'
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            statistical_main(tre_dic,outfile,gene2new_named_gene_dic,voucher2taxa_dic,top_n)
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
            species_category_list = [read_and_return_dict(i) for i in gene_categories]
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            
            tre_dic = read_and_return_dict(input_GF_list)
            gene2fam = None
            df = None
            
            if args.gene_family and args.input_sps_tree:
                gene2fam = read_and_return_dict(args.gene_family)
                gene2sps = read_and_return_dict(input_imap)
                sptree = Tree(args.input_sps_tree)
                sptree1=rename_input_tre(sptree,taxa2voucher_dic)
                mark_gene_to_sptree_main(tre_dic, species_category_list, sptree1, gene2fam, gene2sps,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic)
            
            if args.gene_expression:
                file_extension = os.path.splitext(args.gene_expression)[1]  # 获取文件扩展名
                
                if file_extension == '.xlsx' or file_extension == '.xls':
                    df = pd.read_excel(args.gene_expression, index_col=0)
                elif file_extension == '.csv':
                    df = pd.read_csv(args.gene_expression, index_col=0)
                else:
                    raise ValueError("Unsupported file format. Please provide an Excel or CSV file.")

            view_main(tre_dic, gene2new_named_gene_dic, voucher2taxa_dic, species_category_list, tree_style, keep_branch, new_named_gene2gene_dic, gene2fam, df,visual=args.visual_gd)

            # 计算并格式化程序执行时间
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Tree_Visualizer command are missing.")

    
    elif args.command == 'GD_Detector':
        # Execute the GD_Detector function
        if args.input_GF_list and args.input_imap and args.input_sps_tree and args.gd_support is not None and args.subclade_support is not None and args.dup_species_proportion is not None and args.dup_species_num is not None and args.deepvar is not None:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            gd_support=args.gd_support
            subclade_support=args.subclade_support
            dup_species_percent = args.dup_species_proportion
            dup_species_num = args.dup_species_num
            deep_var=args.deepvar
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            sptree=PhyloTree(args.input_sps_tree)
            num_tre_node(sptree)
            sptree.write(outfile='numed_sptree.nwk',format=1)
            renamed_sptree=rename_input_tre(sptree, taxa2voucher_dic)
            tre_dic = read_and_return_dict(input_GF_list)
            filename = 'gd_result.txt'
            write_gd_result(filename, tre_dic, gd_support,subclade_support,dup_species_percent, dup_species_num,renamed_sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,deep_var)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for GD_Detector command are missing.")


    elif args.command == 'GD_Visualizer':
        # Execute the GD_Visualizer function
        if args.input_sps_tree and args.gd_result and args.input_imap :
            start_time = time.time()
            sptree=Tree(args.input_sps_tree,format=1)
            gd_result = args.gd_result
            taxa=read_and_return_dict(args.input_imap)
            gd_visualizer_main(sptree,gd_result,taxa)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for GD_Visualizer command are missing.")


    elif args.command == 'GD_Loss_Tracker':
        # Execute the GD_Loss_Tracker function
        if args.input_GF_list and args.input_sps_tree and args.input_imap :
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_sps_tree = args.input_sps_tree
            # dir_path = os.path.join(os.getcwd(), "gd_loss/")
            # if os.path.exists(dir_path):
            #     shutil.rmtree(dir_path)
            # os.makedirs(dir_path)

            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            sptree=PhyloTree(input_sps_tree)
            num_sptree(sptree)

            tre_dic=read_and_return_dict(input_GF_list)

            get_path_str_num_dic(tre_dic,sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic)
            
            

            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for GD_Loss_Tracker command are missing.")

    elif args.command == 'GD_Loss_Visualizer':
        # Execute the GD_Loss_Visualizer function
        if args.input_sps_tree and args.gd_loss_result :
        # if args.input_folder and  args.output_folder :
            start_time = time.time()
            

            input_sps_tree=args.input_sps_tree
            sptree=Tree(input_sps_tree,format=1)
            # result=process_gd_loss_summary(args.gd_loss_gf_result)
            
            #visualizer_sptree(result,sptree)
            

            a=parse_file_to_dic(args.gd_loss_result)
            visualizer_sptree(a,sptree)

            #generate_plt(input_dir,out_dir)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for GD_Loss_Visualizer command are missing.")

    elif args.command == 'Ortho_Retriever':
        # Execute the Ortho_Retriever function
        if args.input_GF_list and args.input_imap  and args.input_gene_length:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_imap = args.input_imap
            input_gene_length = args.input_gene_length
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            tre_dic = read_and_return_dict(input_GF_list)
            len_dic = read_and_return_dict(input_gene_length)
            renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
            split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Ortho_Retriever command are missing.")


    elif args.command == 'Hybrid_Tracer':
        # Execute the Hybrid_Tracer function
        if args.input_GF_list and args.input_Seq_GF_list  and args.input_sps_tree and args.input_imap:
            start_time = time.time()
            input_GF_list = args.input_GF_list
            input_Seq_GF_list = args.input_Seq_GF_list
            input_sps_tree = args.input_sps_tree
            input_imap= args.input_imap
            sptree=read_phylo_tree(input_sps_tree)
            target_node_name = process_start_node(args.target_node, sptree) if args.target_node else None
            print(target_node_name)
            tre_dic = read_and_return_dict(input_GF_list)
            seq_path_dic = read_and_return_dict(input_Seq_GF_list)
            gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(input_imap)
            
            num_tre_node(sptree)
            rename_sptree=rename_input_tre(sptree,taxa2voucher_dic)
            sptree.write(outfile='numed_sptree.nwk',format=1)

            
            if target_node_name:
                print(f"Target node determined: {target_node_name}")
            else:
                print("No target_species_file provided, processing all gd.")
            hyde_main(tre_dic,seq_path_dic,rename_sptree,gene2new_named_gene_dic,voucher2taxa_dic,taxa2voucher_dic,new_named_gene2gene_dic,target_node=target_node_name,split_gd=args.split_gd)
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Hybrid_Tracer command are missing.")

    elif args.command == 'Hybrid_Visualizer':
        # Execute the Hybrid_Visualizer function
        if args.hyde_out  and args.input_sps_tree :
            start_time = time.time()
            hyde_out=args.hyde_out
            input_sps_tree = args.input_sps_tree
            sptree=read_tree(input_sps_tree)
            if args.node:
                hyde_visual_node_main(hyde_out,sptree) 
            else:
                hyde_visual_leaf_main(hyde_out,sptree)
            
            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for Hybrid_Visualizer command are missing.")

    elif args.command == 'HaploFinder':
        required_args = [args.input_GF_list, args.input_imap, args.species_a, args.species_b,
                     args.species_a_gff, args.species_b_gff, args.species_a_lens, args.species_b_lens,args.gd_support,args.pair_support]
    
        if all(required_args):
            start_time = time.time()
            
            process_gd_pairs = process_gd_result(args.input_GF_list, args.input_imap, args.species_a, args.species_b,args.gd_support,args.pair_support)
            
            # GFF and lens variables
            gff1, gff2 = args.species_a_gff, args.species_b_gff
            lens1, lens2 = args.species_a_lens, args.species_b_lens
            spe1, spe2 = args.species_a, args.species_b

            target_chr1 = args.visual_chr_a
            target_chr2 = args.visual_chr_b
            size = args.size if args.size else 0.001

            generate_dotplot(gff1, gff2, lens1, lens2, process_gd_pairs, spe1, spe2, 'gd_pairs', target_chr1, target_chr2, size)

            end_time = time.time()
            execution_time = end_time - start_time
            formatted_time = format_time(execution_time)
            print("Program execution time:", formatted_time)
        else:
            print("Required arguments for HaploFinder command are missing.")
     
    else:
        print("Usage: python PhyloTracer.py  [-h]  {BranchLength_NumericConverter, GD_Detector, GD_Loss_Tracker, GD_Loss_Visualizer, GD_Visualizer, HaploFinder, Hybrid_Tracer, Hybrid_Visualizer, OrthoFilter_LB, OrthoFilter_Mono, Ortho_Retriever, PhyloSupport_Scaler, PhyloTree_CollapseExpand, Phylo_Rooter, TreeTopology_Summarizer, Tree_Visualizer}")
        print()
        print("Optional arguments:")
        print('  -h, --help            show this help message and exit')
        print()
        print('Available programs:')
        print('  {BranchLength_NumericConverter, GD_Detector, GD_Loss_Tracker, GD_Loss_Visualizer, GD_Visualizer, HaploFinder, Hybrid_Tracer, Hybrid_Visualizer, OrthoFilter_LB, OrthoFilter_Mono, Ortho_Retriever, PhyloSupport_Scaler, PhyloTree_CollapseExpand, Phylo_Rooter, TreeTopology_Summarizer, Tree_Visualizer}')


if __name__ == "__main__":
    main()

