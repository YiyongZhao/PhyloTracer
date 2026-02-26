"""
Command-line entry point for the PhyloTracer toolkit.

This module defines subcommands, CLI arguments, and command dispatch for
phylogenetic tree processing, orthology filtering, duplication detection,
hybridization analyses, and visualization workflows.
"""

import argparse
import logging
import os
import sys
import tempfile
import textwrap
import time

logger = logging.getLogger(__name__)

import pandas as pd
from ete3 import PhyloTree, Tree
from tqdm import tqdm

from phylotracer import (
    gene_id_transfer,
    num_sptree,
    read_and_return_dict,
    read_phylo_tree,
    read_tree,
    rename_input_tre,
)
from phylotracer.BranchLength_NumericConverter import *
from phylotracer.GD_Detector import *
from phylotracer.GD_Loss_Tracker import *
from phylotracer.GD_Loss_Visualizer import *
from phylotracer.GD_Visualizer import *
from phylotracer.HaploFinder import *
from phylotracer.Hybrid_Tracer import *
from phylotracer.Hybrid_Visualizer import *
from phylotracer.Ortho_Retriever import *
from phylotracer.OrthoFilter_LB import *
from phylotracer.OrthoFilter_Mono import *
from phylotracer.Phylo_Rooter import *
from phylotracer.PhyloSupport_Scaler import *
from phylotracer.PhyloTree_CollapseExpand import *
from phylotracer.Tree_Visualizer import *
from phylotracer.TreeTopology_Summarizer import *

# Keep CLI logger stable after wildcard imports that may export `logger`.
logger = logging.getLogger(__name__)

##################################################################

BANNER = textwrap.dedent("""
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
#    Contacts: Taoli(Taoli@gmail.com); Yiyong Zhao(yiyong.zhao@yale.edu)                     #
#                                                                                             #
###############################################################################################
""")

# Custom help formatter for better formatting of subcommands and section titles


class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    """Unified CLI help formatter with stable alignment across subcommands."""

    def __init__(self, prog):
        super().__init__(prog, max_help_position=34, width=120)

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


def bounded_int(min_value: int = None, max_value: int = None):
    """Create an argparse int type with inclusive bounds."""
    def _validator(value: str) -> int:
        parsed = int(value)
        if min_value is not None and parsed < min_value:
            raise argparse.ArgumentTypeError(f"must be >= {min_value}")
        if max_value is not None and parsed > max_value:
            raise argparse.ArgumentTypeError(f"must be <= {max_value}")
        return parsed
    return _validator


def bounded_float(min_value: float = None, max_value: float = None):
    """Create an argparse float type with inclusive bounds."""
    def _validator(value: str) -> float:
        parsed = float(value)
        if min_value is not None and parsed < min_value:
            raise argparse.ArgumentTypeError(f"must be >= {min_value}")
        if max_value is not None and parsed > max_value:
            raise argparse.ArgumentTypeError(f"must be <= {max_value}")
        return parsed
    return _validator


# Create the main parser
parser = argparse.ArgumentParser(
    prog='PhyloTracer',
    description='A comprehensive toolkit for phylogenetic tree analysis, gene family evolution inference, and visualization, supporting tree transformation, rooting, topology summarization, duplication-loss analysis, and hybridization detection.',
    formatter_class=CustomHelpFormatter,
    add_help=False,
)

# Add subparsers for commands
subparsers = parser.add_subparsers(dest='command', help='Available subcommands')


# PhyloTree_CollapseExpand command
PhyloTree_CollapseExpand_parser = subparsers.add_parser('PhyloTree_CollapseExpand', help='Collapse low-support branches or restore resolved binary trees', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer PhyloTree_CollapseExpand --input_GF_list GF_ID2path.imap --support_value 50')
PhyloTree_CollapseExpand_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
PhyloTree_CollapseExpand_parser.add_argument('--support_value', metavar='INT', type=bounded_int(0, 100), required=True, help='Node support cutoff used for collapsing internal branches, default = 50')
PhyloTree_CollapseExpand_parser.add_argument('--revert', action='store_true', help='If set, expand previously collapsed comb structures back to binary form, default = False')
PhyloTree_CollapseExpand_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# PhyloSupport_Scaler command
PhyloSupport_Scaler_parser = subparsers.add_parser('PhyloSupport_Scaler', help='Rescale branch support values for all gene trees', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer PhyloSupport_Scaler --input_GF_list GF_ID2path.imap --scale_to 100')
PhyloSupport_Scaler_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
PhyloSupport_Scaler_parser.add_argument('--scale_to', metavar='1|100', choices=['1', '100'], required=True, help='Target support scale: "1" for [0,1], "100" for [0,100]')
PhyloSupport_Scaler_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# BranchLength_NumericConverter command
BranchLength_NumericConverter_parser = subparsers.add_parser('BranchLength_NumericConverter', help='Standardize branch-length precision for all gene trees', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer BranchLength_NumericConverter --input_GF_list GF_ID2path.imap')
BranchLength_NumericConverter_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
BranchLength_NumericConverter_parser.add_argument('--decimal_place', metavar='INT', type=bounded_int(0), default=10, help='Number of decimal places to keep for branch lengths, default = 10')
BranchLength_NumericConverter_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# Phylo_Rooter command
Phylo_Rooter_parser = subparsers.add_parser('Phylo_Rooter', help='Root gene trees using species-tree guidance and gene length', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer Phylo_Rooter --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap --input_sps_tree sptree.nwk')
Phylo_Rooter_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
Phylo_Rooter_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
Phylo_Rooter_parser.add_argument('--input_gene_length', metavar='GENE_LENGTH_LIST', required=True, help='Two-column mapping file (gene_id<TAB>gene_length)')
Phylo_Rooter_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Species tree file in Newick format')
Phylo_Rooter_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# OrthoFilter_LB command
OrthoFilter_LB_parser = subparsers.add_parser('OrthoFilter_LB', help='Remove long-branch outliers from gene trees', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer OrthoFilter_LB --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --rrbr_cutoff 5 --srbr_cutoff 2.5')
OrthoFilter_LB_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
OrthoFilter_LB_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
OrthoFilter_LB_parser.add_argument('--rrbr_cutoff', metavar='FLOAT', type=bounded_float(0.0), default=5, required=True, help='RRBR cutoff based on root-to-tip distance, default = 5')
OrthoFilter_LB_parser.add_argument('--srbr_cutoff', metavar='FLOAT', type=bounded_float(0.0), default=2.5, required=True, help='SRBR cutoff based on sister-relative branch ratio, default = 2.5')
OrthoFilter_LB_parser.add_argument('--visual', action='store_true', help='If set, export before/after tree visualization PDFs, default = False')
OrthoFilter_LB_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# OrthoFilter_Mono command
OrthoFilter_Mono_parser = subparsers.add_parser('OrthoFilter_Mono', help='Prune alien lineages inside dominant clades', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer OrthoFilter_Mono --input_GF_list GF_ID2path.imap --input_taxa gene2clade.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk')
OrthoFilter_Mono_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
OrthoFilter_Mono_parser.add_argument('--input_taxa', metavar='TAXA_LIST', required=True, help='Two-column mapping file (gene_id<TAB>clade_or_lineage_label)')
OrthoFilter_Mono_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
OrthoFilter_Mono_parser.add_argument('--purity_cutoff', metavar='FLOAT', type=bounded_float(0.0, 1.0), default=0.95, help='Target purity for dominant lineage, default = 0.95')
OrthoFilter_Mono_parser.add_argument('--max_remove_fraction', metavar='FLOAT', type=bounded_float(0.0, 1.0), default=0.5, help='Maximum fraction of tips allowed to be removed, default = 0.5')
OrthoFilter_Mono_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Species tree file in Newick format')
OrthoFilter_Mono_parser.add_argument('--visual', action='store_true', help='If set, export before/after pruning visualization PDFs, default = False')
OrthoFilter_Mono_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# TreeTopology_Summarizer command
TreeTopology_Summarizer_parser = subparsers.add_parser('TreeTopology_Summarizer', help='Summarize and rank gene-tree topologies', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer TreeTopology_Summarizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap')
TreeTopology_Summarizer_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
TreeTopology_Summarizer_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
TreeTopology_Summarizer_parser.add_argument('--visual_top', metavar='INT', type=bounded_int(1), required=False, default=10, help='Number of top-ranked topologies to visualize, default = 10')
TreeTopology_Summarizer_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# Tree_Visualizer command
Tree_Visualizer_parser = subparsers.add_parser('Tree_Visualizer', help='Render gene-tree figures with optional metadata overlays', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer Tree_Visualizer --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap')
Tree_Visualizer_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
Tree_Visualizer_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
Tree_Visualizer_parser.add_argument('--gene_categories', metavar='TAXA_LIST', nargs='+', help='One or more two-column files (gene_id<TAB>category_label) for color annotations')
Tree_Visualizer_parser.add_argument('--keep_branch', metavar='0|1', choices=['0', '1'], help='Whether to preserve branch lengths in plotting: 1=yes, 0=no')
Tree_Visualizer_parser.add_argument('--tree_style', metavar='r|c', choices=['r', 'c'], default='r', help='Tree layout style: r=rectangular, c=circular')
Tree_Visualizer_parser.add_argument('--gene_family', metavar='GENE_FAMILY', help='Two-column mapping file (gene_id<TAB>family_label)')
Tree_Visualizer_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', help='Species tree file in Newick format (required with --gene_family)')
Tree_Visualizer_parser.add_argument('--gene_expression', metavar='EXPRESSION_FILE', help='Gene expression matrix file (.csv/.xls/.xlsx), genes as row index')
Tree_Visualizer_parser.add_argument('--visual_gd', action='store_true', help='If set, overlay predicted GD nodes on gene-tree figures, default = False')
Tree_Visualizer_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# GD_Detector command
GD_Detector_parser = subparsers.add_parser('GD_Detector', help='Detect gene-duplication events and classify GD types', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer GD_Detector --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --gd_support 50 --subclade_support 50 --dup_species_proportion 0 --dup_species_num 2 --deepvar 1')
GD_Detector_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
GD_Detector_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
GD_Detector_parser.add_argument('--gd_support', metavar='INT', type=bounded_int(0, 100), required=True, default=50, help='Minimum support of a GD candidate node (accepted range: 0-100; typical: 50-100), default = 50')
GD_Detector_parser.add_argument('--subclade_support', metavar='INT', type=bounded_int(0, 100), required=True, default=0, help='Minimum support required in GD child subclades (accepted range: 0-100), default = 0')
GD_Detector_parser.add_argument('--dup_species_proportion', metavar='FLOAT', type=bounded_float(0.0, 1.0), required=True, default=0.2, help='Minimum overlap ratio of duplicated species between the two GD child clades (range: 0-1), default = 0.2')
GD_Detector_parser.add_argument('--dup_species_num', metavar='INT', type=bounded_int(1), required=True, default=2, help='Minimum number of overlapping duplicated species under a GD node, default = 2')
GD_Detector_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Species tree file in Newick format')
GD_Detector_parser.add_argument('--deepvar', metavar='INT', type=bounded_int(0), required=True, default=1, help='Maximum tolerated depth-variance score for GD screening, default = 1')
GD_Detector_parser.add_argument('--gdtype_mode', choices=['relaxed', 'strict'], default='relaxed', help='GD type mode: relaxed (species overlap only) or strict (overlap + depth constraint), default = relaxed')
GD_Detector_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')
# GD_Visualizer command
GD_Visualizer_parser = subparsers.add_parser('GD_Visualizer', help='Visualize GD counts or GD-type summaries on a species tree', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer GD_Visualizer --input_sps_tree numed_sptree.nwk --gd_result gd_result.txt --input_imap gene2sps.imap')
GD_Visualizer_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Numbered species tree file in Newick format')
GD_Visualizer_parser.add_argument('--gd_result', metavar='GD_RESULT', required=True, help='GD result table produced by GD_Detector')
GD_Visualizer_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
GD_Visualizer_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# GD_Loss_Tracker command
GD_Loss_Tracker_parser = subparsers.add_parser('GD_Loss_Tracker', help='Track inferred post-GD loss paths across species tree branches', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer GD_Loss_Tracker --input_GF_list GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap')
GD_Loss_Tracker_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
GD_Loss_Tracker_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Species tree file in Newick format')
GD_Loss_Tracker_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
GD_Loss_Tracker_parser.add_argument('--target_species', metavar='SP', action='append', default=None, help='Only count loss paths ending in this species (e.g., Arabidopsis_thaliana). Can be used multiple times.')
GD_Loss_Tracker_parser.add_argument('--mrca_node', metavar='SP1,SP2', action='append', default=None, help='Only count loss paths passing through the MRCA of SP1 and SP2. Format: SpeciesA,SpeciesB (comma-separated, no space). Can be used multiple times.')
GD_Loss_Tracker_parser.add_argument('--include_unobserved_species', action='store_true', help='Classification policy for species absent from the current gene family, default = False. If set, treat unobserved species as classifiable (2-2/2-1/2-0) based on left/right presence; if not set, label them as missing_data. This flag does not change loss_path copy-state values (0/1/2) or GD event detection.')
GD_Loss_Tracker_parser.add_argument('--node_count_mode', choices=['nonaccumulate', 'accumulate'], default='nonaccumulate', help='Node counting mode for path_count_* transition statistics, default = nonaccumulate. nonaccumulate: within one GD event, repeated transitions on the same internal node are counted once; accumulate: keep all repeated transitions. This flag does not change the (0/1/2) copy-state numbers shown in loss_path.')
GD_Loss_Tracker_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# GD_Loss_Visualizer command
GD_Loss_Visualizer_parser = subparsers.add_parser('GD_Loss_Visualizer', help='Visualize GD-loss summary results on species tree topology', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer GD_Loss_Visualizer --gd_loss_result gd_loss_summary.txt --input_sps_tree numed_sptree.nwk')
GD_Loss_Visualizer_parser.add_argument('--gd_loss_result', metavar='GD_LOSS_RESULT', required=True, help='Detailed table generated by GD_Loss_Tracker (gd_loss_summary.txt)')
GD_Loss_Visualizer_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=False, help='Numbered species tree file in Newick format')
GD_Loss_Visualizer_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# Ortho_Retriever command
Ortho_Retriever_parser = subparsers.add_parser('Ortho_Retriever', help='Retrieve ortholog sets from rooted gene trees', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer Ortho_Retriever --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_gene_length gene2length.imap')
Ortho_Retriever_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
Ortho_Retriever_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
Ortho_Retriever_parser.add_argument('--input_gene_length', metavar='GENE_LENGTH_LIST', required=True, help='Two-column mapping file (gene_id<TAB>gene_length)')
Ortho_Retriever_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# Hybrid_Tracer
Hybrid_Tracer_parser = subparsers.add_parser('Hybrid_Tracer', help='Prepare HYDE input from GD-supported candidate triplets', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer Hybrid_Tracer --input_GF_list GF_ID2path.imap --input_Seq_GF_list Seq_GF_ID2path.imap --input_sps_tree sptree.nwk --input_imap gene2sps.imap')
Hybrid_Tracer_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree path per line')
Hybrid_Tracer_parser.add_argument('--input_Seq_GF_list', metavar='ALIGNMENT_LIST', required=True, help='Tab-delimited mapping file (GF_ID<TAB>alignment_path), aligned with --input_GF_list IDs')
Hybrid_Tracer_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Species tree file in Newick format')
Hybrid_Tracer_parser.add_argument('--input_imap', metavar='IMAP', required=True, help='Two-column mapping file (gene_id<TAB>species_name)')
Hybrid_Tracer_parser.add_argument('--mrca_node', metavar='SP1,SP2', action='append', default=None, help='Restrict Hybrid_Tracer to the MRCA of SP1 and SP2. Format: SpeciesA,SpeciesB (comma-separated, no space). If multiple are provided, only the first valid pair is used.')
Hybrid_Tracer_parser.add_argument('--split_groups', type=bounded_int(1), required=False, default=1, help='Number of partitions for HYDE batch processing, default = 1')
Hybrid_Tracer_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# Hybrid_Visualizer
Hybrid_Visualizer_parser = subparsers.add_parser('Hybrid_Visualizer', help='Visualize HYDE hybridization signals', formatter_class=CustomHelpFormatter, epilog='Example:\n  PhyloTracer Hybrid_Visualizer --hyde_out hyde_out.txt --input_sps_tree sptree.nwk')
Hybrid_Visualizer_parser.add_argument('--hyde_out', metavar='HYDE_OUT', required=True, help='HYDE output table file')
Hybrid_Visualizer_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=True, help='Species tree file in Newick format')
Hybrid_Visualizer_parser.add_argument('--node', action="store_true", default=False, help="Use node-mode heatmaps (monophyletic clade stacking) instead of leaf-mode output")
Hybrid_Visualizer_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')

# HaploFinder
haplofinder_parser = subparsers.add_parser('HaploFinder', help='Detect and visualize haplotype-level GD signals; also supports FASTA split mode', formatter_class=CustomHelpFormatter, epilog='Examples:\n  PhyloTracer HaploFinder --mode haplofinder --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_sps_tree sptree.nwk --species_a A --species_b B --species_a_gff A.gff --species_b_gff B.gff --species_a_lens A.lens --species_b_lens B.lens\n  PhyloTracer HaploFinder --mode split --input_GF_list GF_ID2path.imap --input_imap gene2sps.imap --input_fasta proteins.fa --cluster_file cluster.tsv --hyb_sps Hybrid --parental_sps \"P1 P2\" --species_b_gff B.gff')
haplofinder_parser.add_argument('--input_GF_list', metavar='GENE_TREE_LIST', required=False, default=None, help='Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); required in haplofinder mode')
haplofinder_parser.add_argument('--input_imap', metavar='IMAP', required=False, default=None, help='Two-column mapping file (gene_id<TAB>species_name); required in both haplofinder and split modes')
haplofinder_parser.add_argument('--input_sps_tree', metavar='NEWICK_TREE', required=False, help='Species tree file in Newick format; required in haplofinder mode')
haplofinder_parser.add_argument('--species_a', metavar='SPECIES_A', type=str, required=False, help='Name of species A (required for haplofinder mode)')
haplofinder_parser.add_argument('--species_b', metavar='SPECIES_B', type=str, required=False, help='Name of species B (required for haplofinder mode)')
haplofinder_parser.add_argument('--species_a_gff', metavar='GFF', required=False, help='Genome annotation file for species A in GFF/GTF-compatible format')
haplofinder_parser.add_argument('--species_b_gff', metavar='GFF', required=False, help='Genome annotation file for species B in GFF/GTF-compatible format')
haplofinder_parser.add_argument('--species_a_lens', metavar='LENS', required=False, help='Chromosome-length file for species A (chr<TAB>length)')
haplofinder_parser.add_argument('--species_b_lens', metavar='LENS', required=False, help='Chromosome-length file for species B (chr<TAB>length)')
haplofinder_parser.add_argument('--visual_chr_a', metavar='CHR_LIST', required=False, help='Optional chromosome list file for species A visualization subset')
haplofinder_parser.add_argument('--visual_chr_b', metavar='CHR_LIST', required=False, help='Optional chromosome list file for species B visualization subset')
haplofinder_parser.add_argument('--gd_support', metavar='INT', type=bounded_int(0, 100), required=False, default=50, help='Minimum support of GD nodes used for pair extraction (accepted range: 0-100, default = 50)')
haplofinder_parser.add_argument('--pair_support', metavar='INT', type=bounded_int(0, 100), required=False, default=50, help='Minimum support of ortholog/speciation pair nodes (accepted range: 0-100, default = 50)')
haplofinder_parser.add_argument('--size', metavar='FLOAT', type=bounded_float(0.0), required=False, default=0.0005, help='Point size in dotplot rendering (positive float, default = 0.0005)')
haplofinder_parser.add_argument('--mode', metavar='MODE', type=str, choices=['haplofinder', 'split'], default='haplofinder', help='Run mode: "haplofinder" for GD analysis, "split" for FASTA partitioning by color labels, default = haplofinder')
# Split mode specific arguments
haplofinder_parser.add_argument('--hyb_sps', metavar='HYBRID_SPECIES', type=str, required=False, help='Hybrid species name used for subgenome assignment (required in split mode)')
haplofinder_parser.add_argument('--parental_sps', metavar='PARENTAL_SPECIES', type=str, required=False, help='Parental species names used for split-mode assignment; provide as a single quoted, space-separated string')
haplofinder_parser.add_argument('--input_fasta', metavar='FASTA_FILE', required=False, help='Input FASTA file (.fa/.fasta), required in split mode')
haplofinder_parser.add_argument('--cluster_file', metavar='CLUSTER_FILE', required=False, help='Split-mode cluster metadata file (legacy compatibility field; currently required by CLI checks)')
haplofinder_parser.add_argument('--output_dir', metavar='DIR', default=None, help='Output directory (default: current working directory)')
parser.add_argument('-h', '--help', action='store_true', help=argparse.SUPPRESS)
# Analyze command line parameters

_parser = parser
CURRENT_COMMAND = None


def format_time(seconds):
    days = seconds // (24 * 3600)
    hours = (seconds % (24 * 3600)) // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return f"{int(days)} d, {int(hours)} h, {int(minutes)} m, {seconds:.2f} s"


def report_execution_time(start_time: float) -> None:
    """Log formatted runtime based on a given start timestamp."""
    runtime_logger_name = (
        f"phylotracer.{CURRENT_COMMAND}" if CURRENT_COMMAND else __name__
    )
    runtime_logger = logging.getLogger(runtime_logger_name)
    runtime_logger.info(
        "Program execution time: %s",
        format_time(time.time() - start_time),
    )


def handle_phylo_tree_collapse_expand(cli_args):
    if cli_args.input_GF_list and cli_args.support_value is not None:
        start_time = time.time()
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        collapse_expand_main(tre_dic, cli_args.support_value, revert=cli_args.revert)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for PhyloTree_CollapseExpand command are missing.")


def handle_phylo_support_scaler(cli_args):
    if cli_args.input_GF_list and cli_args.scale_to:
        start_time = time.time()
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        support_scaler_main(tre_dic, cli_args.scale_to)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for PhyloSupport_Scaler command are missing.")


def handle_branch_length_numeric_converter(cli_args):
    if cli_args.input_GF_list:
        start_time = time.time()
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        branch_length_numeric_converter_main(tre_dic, cli_args.decimal_place)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for BranchLength_NumericConverter command are missing.")


def handle_phylo_rooter(cli_args):
    if cli_args.input_GF_list and cli_args.input_imap and cli_args.input_sps_tree and cli_args.input_gene_length:
        start_time = time.time()
        gene2new_named_gene_dic, new_named_gene2gene_dic, _, taxa2voucher_dic = gene_id_transfer(cli_args.input_imap)
        len_dic = read_and_return_dict(cli_args.input_gene_length)
        renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
        sptree = PhyloTree(cli_args.input_sps_tree)
        renamed_sptree = rename_input_tre(sptree, taxa2voucher_dic)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic, renamed_sptree)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for Phylo_Rooter command are missing.")


def handle_ortho_filter_lb(cli_args):
    if (
        cli_args.input_GF_list
        and cli_args.input_imap
        and cli_args.rrbr_cutoff is not None
        and cli_args.srbr_cutoff is not None
    ):
        start_time = time.time()
        rrbr_cutoff = cli_args.rrbr_cutoff
        srbr_cutoff = cli_args.srbr_cutoff
        gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, _ = gene_id_transfer(cli_args.input_imap)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        prune_main_LB(
            tre_dic,
            voucher2taxa_dic,
            gene2new_named_gene_dic,
            new_named_gene2gene_dic,
            rrbr_cutoff,
            srbr_cutoff,
            visual=cli_args.visual,
        )
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for OrthoFilter_LB command are missing.")


def handle_ortho_filter_mono(cli_args):
    if (
        cli_args.input_GF_list
        and cli_args.input_imap
        and cli_args.input_taxa
        and cli_args.input_sps_tree
        and cli_args.purity_cutoff is not None
        and cli_args.max_remove_fraction is not None
    ):
        start_time = time.time()
        gene2new_named_gene_dic, new_named_gene2gene_dic, _, taxa2voucher_dic = gene_id_transfer(cli_args.input_imap)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        taxa_dic = read_and_return_dict(cli_args.input_taxa)
        sptree = PhyloTree(cli_args.input_sps_tree)
        renamed_sptree = rename_input_tre(sptree, taxa2voucher_dic)
        prune_main_Mono(
            tre_dic,
            taxa_dic,
            renamed_sptree,
            cli_args.purity_cutoff,
            cli_args.max_remove_fraction,
            new_named_gene2gene_dic,
            gene2new_named_gene_dic,
            visual=cli_args.visual,
        )
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for OrthoFilter_Mono command are missing.")


def handle_tree_topology_summarizer(cli_args):
    if cli_args.input_GF_list and cli_args.input_imap:
        start_time = time.time()
        gene2new_named_gene_dic, _, voucher2taxa_dic, _ = gene_id_transfer(cli_args.input_imap)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        statistical_main(tre_dic, 'topology', gene2new_named_gene_dic, voucher2taxa_dic, cli_args.visual_top)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for TreeTopology_Summarizer command are missing.")


def handle_tree_visualizer(cli_args):
    if cli_args.input_GF_list and cli_args.input_imap:
        start_time = time.time()
        if cli_args.visual_gd and not cli_args.input_sps_tree:
            logger.error("Tree_Visualizer: --input_sps_tree is required when --visual_gd is enabled.")
            return
        species_category_list = []
        if cli_args.gene_categories is not None:
            species_category_list = [read_and_return_dict(i) for i in cli_args.gene_categories]

        gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(cli_args.input_imap)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        gene2fam = None
        df = None
        sptree1 = None

        if cli_args.gene_family and not cli_args.input_sps_tree:
            logger.error("Tree_Visualizer: --input_sps_tree is required when --gene_family is provided.")
            return

        if cli_args.input_sps_tree:
            sptree = Tree(cli_args.input_sps_tree)
            sptree1 = rename_input_tre(sptree, taxa2voucher_dic)

        if cli_args.gene_family:
            gene2fam = read_and_return_dict(cli_args.gene_family)
            gene2sps = read_and_return_dict(cli_args.input_imap)
            mark_gene_to_sptree_main(
                tre_dic,
                species_category_list,
                sptree1,
                gene2fam,
                gene2sps,
                gene2new_named_gene_dic,
                new_named_gene2gene_dic,
                voucher2taxa_dic,
            )

        if cli_args.gene_expression:
            file_extension = os.path.splitext(cli_args.gene_expression)[1]
            if file_extension in ('.xlsx', '.xls'):
                df = pd.read_excel(cli_args.gene_expression, index_col=0)
            elif file_extension == '.csv':
                df = pd.read_csv(cli_args.gene_expression, index_col=0)
            else:
                raise ValueError("Unsupported file format. Please provide an Excel or CSV file.")

        view_main(
            tre_dic,
            sptree1,
            gene2new_named_gene_dic,
            voucher2taxa_dic,
            species_category_list,
            cli_args.tree_style,
            cli_args.keep_branch,
            new_named_gene2gene_dic,
            gene2fam,
            df,
            visual=cli_args.visual_gd,
        )
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for Tree_Visualizer command are missing.")


def handle_gd_detector(cli_args):
    if (
        cli_args.input_GF_list
        and cli_args.input_imap
        and cli_args.input_sps_tree
        and cli_args.gd_support is not None
        and cli_args.subclade_support is not None
        and cli_args.dup_species_proportion is not None
        and cli_args.dup_species_num is not None
        and cli_args.deepvar is not None
    ):
        start_time = time.time()
        gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(cli_args.input_imap)
        sptree = PhyloTree(cli_args.input_sps_tree)
        num_sptree(sptree)
        sptree.write(outfile='numed_sptree.nwk', format=1, format_root_node=True)
        renamed_sptree = rename_input_tre(sptree, taxa2voucher_dic)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        output_file = f"gd_result_{cli_args.gdtype_mode}.txt"
        write_gene_duplication_results(
            output_file,
            tre_dic,
            cli_args.gd_support,
            cli_args.subclade_support,
            cli_args.dup_species_proportion,
            cli_args.dup_species_num,
            renamed_sptree,
            gene2new_named_gene_dic,
            new_named_gene2gene_dic,
            voucher2taxa_dic,
            cli_args.deepvar,
            gdtype_mode=cli_args.gdtype_mode,
        )
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for GD_Detector command are missing.")


def handle_gd_visualizer(cli_args):
    if cli_args.input_sps_tree and cli_args.gd_result and cli_args.input_imap:
        start_time = time.time()
        sptree = read_tree(cli_args.input_sps_tree)
        taxa = read_and_return_dict(cli_args.input_imap)
        gd_visualizer_main(sptree, cli_args.gd_result, taxa)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for GD_Visualizer command are missing.")


def handle_gd_loss_tracker(cli_args):
    if cli_args.input_GF_list and cli_args.input_sps_tree and cli_args.input_imap:
        start_time = time.time()
        gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(cli_args.input_imap)
        sptree = PhyloTree(cli_args.input_sps_tree)
        num_sptree(sptree)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)

        allowed_gd_species_sets = set()
        if cli_args.mrca_node:
            for pair_str in cli_args.mrca_node:
                if ',' not in pair_str:
                    logger.warning("Invalid mrca_node format '%s'. Skipping.", pair_str)
                    continue
                sp1, sp2 = [x.strip() for x in pair_str.split(',', 1)]
                try:
                    mrca_clade = sptree.get_common_ancestor(sp1, sp2)
                    species_under_mrca = frozenset(mrca_clade.get_leaf_names())
                    allowed_gd_species_sets.add(species_under_mrca)
                    logger.info("GD event restricted to EXACT MRCA of %s and %s", sp1, sp2)
                    logger.info("   -> Covers species: %s", sorted(species_under_mrca))
                except Exception as exc:
                    logger.error("Cannot find MRCA for %s,%s: %s", sp1, sp2, exc)
                    logger.error("Available species: %s", sorted([leaf.name for leaf in sptree.get_leaves()]))

        get_path_str_num_dic(
            tre_dic,
            sptree,
            gene2new_named_gene_dic,
            new_named_gene2gene_dic,
            voucher2taxa_dic,
            taxa2voucher_dic,
            target_species_list=cli_args.target_species,
            allowed_gd_species_sets=allowed_gd_species_sets,
            include_unobserved_species=cli_args.include_unobserved_species,
            node_count_mode=cli_args.node_count_mode,
        )
        parse_text_to_excel('gd_loss_count_summary.txt')
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for GD_Loss_Tracker command are missing.")


def handle_gd_loss_visualizer(cli_args):
    if cli_args.input_sps_tree and cli_args.gd_loss_result:
        start_time = time.time()
        sptree = Tree(cli_args.input_sps_tree, format=1)
        visualizer_sptree(cli_args.gd_loss_result, sptree)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for GD_Loss_Visualizer command are missing.")


def handle_ortho_retriever(cli_args):
    if cli_args.input_GF_list and cli_args.input_imap and cli_args.input_gene_length:
        start_time = time.time()
        gene2new_named_gene_dic, new_named_gene2gene_dic, _, _ = gene_id_transfer(cli_args.input_imap)
        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        len_dic = read_and_return_dict(cli_args.input_gene_length)
        renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
        split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic, renamed_len_dic)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for Ortho_Retriever command are missing.")


def handle_hybrid_tracer(cli_args):
    if cli_args.input_GF_list and cli_args.input_Seq_GF_list and cli_args.input_sps_tree and cli_args.input_imap:
        start_time = time.time()
        sptree = read_phylo_tree(cli_args.input_sps_tree)
        num_sptree(sptree)
        target_node_name = None
        if cli_args.mrca_node:
            for pair_str in cli_args.mrca_node:
                if ',' not in pair_str:
                    logger.warning("Invalid mrca_node format '%s'. Skipping.", pair_str)
                    continue
                sp1, sp2 = [x.strip() for x in pair_str.split(',', 1)]
                try:
                    mrca_clade = sptree.get_common_ancestor(sp1, sp2)
                    target_node_name = mrca_clade.name
                    logger.info("Hybrid_Tracer target MRCA selected: %s,%s -> %s", sp1, sp2, target_node_name)
                    break
                except Exception as exc:
                    logger.error("Cannot find MRCA for %s,%s: %s", sp1, sp2, exc)
                    logger.error("Available species: %s", sorted([leaf.name for leaf in sptree.get_leaves()]))

        if target_node_name:
            logger.info("Target node determined: %s", target_node_name)
        else:
            logger.info("No valid --mrca_node provided, processing all gd.")

        tre_dic = read_and_return_dict(cli_args.input_GF_list)
        seq_path_dic = read_and_return_dict(cli_args.input_Seq_GF_list)
        gene2new_named_gene_dic, _, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(cli_args.input_imap)
        rename_sptree = rename_input_tre(sptree, taxa2voucher_dic)

        sptree.write(outfile='numed_sptree.nwk', format=1)

        hyde_main(
            tre_dic,
            seq_path_dic,
            rename_sptree,
            gene2new_named_gene_dic,
            voucher2taxa_dic,
            target_node=target_node_name,
            gd_group=cli_args.split_groups,
        )
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for Hybrid_Tracer command are missing.")


def handle_hybrid_visualizer(cli_args):
    if cli_args.hyde_out and cli_args.input_sps_tree:
        start_time = time.time()
        sptree = read_tree(cli_args.input_sps_tree)
        if cli_args.node:
            hyde_visual_node_main(cli_args.hyde_out, sptree)
        else:
            hyde_visual_leaf_main(cli_args.hyde_out, sptree)
        report_execution_time(start_time)
    else:
        logger.error("Required arguments for Hybrid_Visualizer command are missing.")


def handle_haplofinder(cli_args):
    if cli_args.mode == 'split':
        if (
            cli_args.input_GF_list and cli_args.input_fasta and cli_args.input_imap
            and cli_args.cluster_file and cli_args.hyb_sps and cli_args.parental_sps
            and cli_args.species_b_gff
        ):
            start_time = time.time()
            parental_sps = cli_args.parental_sps.split() if cli_args.parental_sps else []
            split_sequences(
                cli_args.input_GF_list,
                cli_args.input_imap,
                cli_args.hyb_sps,
                parental_sps,
                cli_args.species_b_gff,
                cli_args.input_fasta,
                cli_args.cluster_file,
            )
            report_execution_time(start_time)
        else:
            logger.error("Required arguments for split mode are missing: --input_GF_list, --input_imap, --input_fasta, --cluster_file, --hyb_sps, --parental_sps, --species_b_gff")
    else:
        required_args = [
            cli_args.input_GF_list,
            cli_args.input_imap,
            cli_args.input_sps_tree,
            cli_args.species_a,
            cli_args.species_b,
            cli_args.species_a_gff,
            cli_args.species_b_gff,
            cli_args.species_a_lens,
            cli_args.species_b_lens,
            cli_args.gd_support,
            cli_args.pair_support,
        ]
        if all(arg is not None for arg in required_args):
            start_time = time.time()
            process_gd_pairs = process_gd_result(
                cli_args.input_GF_list,
                cli_args.input_imap,
                cli_args.input_sps_tree,
                cli_args.species_a,
                cli_args.species_b,
                cli_args.gd_support,
                cli_args.pair_support,
            )
            size = cli_args.size if cli_args.size is not None else 0.001
            generate_dotplot(
                cli_args.species_a_gff,
                cli_args.species_b_gff,
                cli_args.species_a_lens,
                cli_args.species_b_lens,
                process_gd_pairs,
                cli_args.species_a,
                cli_args.species_b,
                'gd_pairs',
                cli_args.visual_chr_a,
                cli_args.visual_chr_b,
                size,
            )
            report_execution_time(start_time)
        else:
            logger.error("Required arguments for HaploFinder command are missing.")


def print_cli_usage():
    print("Usage: python PhyloTracer.py  [-h]  {BranchLength_NumericConverter, GD_Detector, GD_Loss_Tracker, GD_Loss_Visualizer, GD_Visualizer, HaploFinder, Hybrid_Tracer, Hybrid_Visualizer, OrthoFilter_LB, OrthoFilter_Mono, Ortho_Retriever, PhyloSupport_Scaler, PhyloTree_CollapseExpand, Phylo_Rooter, TreeTopology_Summarizer, Tree_Visualizer}")
    print()
    print("Optional arguments:")
    print('  -h, --help            show this help message and exit')
    print()
    print('Available programs:')
    print('  {BranchLength_NumericConverter, GD_Detector, GD_Loss_Tracker, GD_Loss_Visualizer, GD_Visualizer, HaploFinder, Hybrid_Tracer, Hybrid_Visualizer, OrthoFilter_LB, OrthoFilter_Mono, Ortho_Retriever, PhyloSupport_Scaler, PhyloTree_CollapseExpand, Phylo_Rooter, TreeTopology_Summarizer, Tree_Visualizer}')


COMMAND_HANDLERS = {
    'PhyloTree_CollapseExpand': handle_phylo_tree_collapse_expand,
    'PhyloSupport_Scaler': handle_phylo_support_scaler,
    'BranchLength_NumericConverter': handle_branch_length_numeric_converter,
    'Phylo_Rooter': handle_phylo_rooter,
    'OrthoFilter_LB': handle_ortho_filter_lb,
    'OrthoFilter_Mono': handle_ortho_filter_mono,
    'TreeTopology_Summarizer': handle_tree_topology_summarizer,
    'Tree_Visualizer': handle_tree_visualizer,
    'GD_Detector': handle_gd_detector,
    'GD_Visualizer': handle_gd_visualizer,
    'GD_Loss_Tracker': handle_gd_loss_tracker,
    'GD_Loss_Visualizer': handle_gd_loss_visualizer,
    'Ortho_Retriever': handle_ortho_retriever,
    'Hybrid_Tracer': handle_hybrid_tracer,
    'Hybrid_Visualizer': handle_hybrid_visualizer,
    'HaploFinder': handle_haplofinder,
}

DEFAULT_OUTPUT_DIRS = {
    "PhyloTree_CollapseExpand": "collapse_expand_tree",
    "PhyloSupport_Scaler": "support_scaler_tree",
    "BranchLength_NumericConverter": "converter_tree",
    "Phylo_Rooter": "rooted_trees",
    "OrthoFilter_LB": "orthofilter_lb",
    "OrthoFilter_Mono": "orthofilter_mono",
    "TreeTopology_Summarizer": "tree_topology_summarizer",
    "Tree_Visualizer": "tree_visualizer",
    "GD_Detector": "gd_detector",
    "GD_Visualizer": "gd_visualizer",
    "GD_Loss_Tracker": "gd_loss_tracker",
    "GD_Loss_Visualizer": "gd_loss_visualizer",
    "Ortho_Retriever": "ortho_retriever",
    "Hybrid_Tracer": "hybrid_tracer",
    "Hybrid_Visualizer": "hybrid_visualizer",
    "HaploFinder": "haplofinder",
}


def _abspath_if_exists_str(path_value):
    """Return absolute path for non-empty string values."""
    if isinstance(path_value, str) and path_value:
        return os.path.abspath(path_value)
    return path_value


def _rewrite_mapping_file_with_abs_paths(mapping_file: str) -> str:
    """
    Rewrite a two-column mapping file so the 2nd column path is absolute.

    This is used for GF/Seq-GF mapping files whose second column stores
    file-system paths. Relative paths are resolved against the mapping file
    directory.
    """
    mapping_file = os.path.abspath(mapping_file)
    map_dir = os.path.dirname(mapping_file)
    fd, tmp_path = tempfile.mkstemp(prefix="phylotracer_absmap_", suffix=".imap")
    os.close(fd)
    with open(mapping_file, "r", encoding="utf-8") as fin, open(
        tmp_path, "w", encoding="utf-8"
    ) as fout:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line.strip():
                fout.write(raw)
                continue
            parts = line.split()
            if len(parts) < 2:
                fout.write(raw)
                continue
            key = parts[0]
            value = parts[1]
            if not os.path.isabs(value):
                value = os.path.abspath(os.path.join(map_dir, value))
            fout.write(f"{key}\t{value}\n")
    return tmp_path


def _normalize_input_paths(args):
    """Normalize CLI input paths before switching to output_dir."""
    # Paths passed directly as CLI values
    path_args = [
        "input_imap",
        "input_gene_length",
        "input_taxa",
        "input_sps_tree",
        "gd_result",
        "gd_loss_result",
        "hyde_out",
        "gene_family",
        "gene_expression",
        "input_fasta",
        "cluster_file",
        "species_a_gff",
        "species_b_gff",
        "species_a_lens",
        "species_b_lens",
        "visual_chr_a",
        "visual_chr_b",
    ]
    for attr in path_args:
        if hasattr(args, attr):
            setattr(args, attr, _abspath_if_exists_str(getattr(args, attr)))

    # List-like path args
    if hasattr(args, "gene_categories") and args.gene_categories:
        args.gene_categories = [os.path.abspath(p) for p in args.gene_categories]

    # Mapping files with path-bearing 2nd column
    if hasattr(args, "input_GF_list") and args.input_GF_list:
        args.input_GF_list = _rewrite_mapping_file_with_abs_paths(args.input_GF_list)
    if hasattr(args, "input_Seq_GF_list") and args.input_Seq_GF_list:
        args.input_Seq_GF_list = _rewrite_mapping_file_with_abs_paths(args.input_Seq_GF_list)

    # Output directory: if not provided, use command-specific default.
    if hasattr(args, "output_dir"):
        if not args.output_dir:
            args.output_dir = DEFAULT_OUTPUT_DIRS.get(args.command, None)
        if args.output_dir:
            args.output_dir = os.path.abspath(args.output_dir)


def main():
    global CURRENT_COMMAND
    if sys.version_info.major == 2:
        print('You are using Python 2. Please upgrade to Python 3. PhyloTracer quit now...')
        sys.exit()
    print(BANNER)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    args = _parser.parse_args()
    _normalize_input_paths(args)
    handler = COMMAND_HANDLERS.get(args.command)
    if handler is None:
        print_cli_usage()
        return
    CURRENT_COMMAND = args.command
    if hasattr(args, 'output_dir') and args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        os.chdir(args.output_dir)
    handler(args)


if __name__ == "__main__":
    main()
