from __init__ import *
import re

def process_gd_result(gd_file: str) -> list:
    """
    Process the gene duplication result file and extract relevant data.

    Args:
        gd_file (str): Path to the gene duplication result file.

    Returns:
        list: A list of tuples containing gene duplication ID and level.
    """
    gds = []
    with open(gd_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            gds.append((parts[1], parts[5]))
    return gds
	
def get_count_dic(gene_duplications: list) -> dict:
    """
    Count unique gene duplication events by level

    Args:
        gene_duplications (list): List of tuples containing (gene_dup_id, level)

    Returns:
        dict: Dictionary with levels as keys and counts as values
    """
    count_dict = {}
    seen_ids = set()

    for dup_id, level in gene_duplications:
        if dup_id not in seen_ids:
            count_dict[level] = count_dict.get(level, 0) + 1
            seen_ids.add(dup_id)

    return count_dict


def mark_sptree(sptree: object, count_dic: dict, taxa: dict) -> object:
    """
    Marks a species tree with gene duplication counts and taxa labels.

    Args:
        sptree (object): The species tree object to be marked.
        count_dic (dict): A dictionary mapping node names to gene duplication counts.
        taxa (dict): A dictionary mapping leaf names to taxa labels.

    Returns:
        object: The marked species tree object.
    """
    sptree.ladderize()
    sptree.sort_descendants("support")

    ts = TreeStyle()
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = 'black'
    ts.branch_vertical_margin = -1
    ts.legend_position = 1  
    ts.legend.add_face(TextFace("Legend:", fsize=8, bold=True), column=0)
    ts.legend.add_face(TextFace("Red numbers: Gene duplication events", fsize=8, fgcolor="red"), column=0)
    ts.legend.add_face(TextFace("Blue numbers: Node identifiers", fsize=8, fgcolor="blue"), column=0)

    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        node.set_style(nstyle)

        num = str(count_dic.get(node.name, 0))
        node.add_face(TextFace(num, fsize=5, fgcolor="red"), column=0, position="branch-top")

        if re.match(r'^N\d+$', node.name):
            node.add_face(TextFace(node.name+' ', fsize=5, fgcolor="blue"), column=0, position="branch-bottom")

    for leaf in sptree:
        if leaf.name in taxa:
            leaf.name = taxa[leaf.name]

    realign_branch_length(sptree)
    rejust_root_dist(sptree)

    return sptree.render('phylotracer_gd_visualizer.pdf', w=210, units="mm", tree_style=ts)

def gd_visualizer_main(sptree, gd_result, taxa):
    """
    Main function to visualize gene duplication events on a species tree.

    Args:
        sptree (object): The species tree object to be visualized.
        gd_result (str): Path to the gene duplication result file.
        taxa (dict): A dictionary mapping leaf names to taxa labels.
    """
    gds = process_gd_result(gd_result)
    count_dic = get_count_dic(gds)
    mark_sptree(sptree, count_dic, taxa)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python GD_Visualizer.py <species_tree_file> <gd_result_file> <taxa_file>")
        sys.exit(1)
    
    species_tree_file = sys.argv[1]
    gd_result_file = sys.argv[2]
    taxa_file = sys.argv[3]
    
    sptree = load_species_tree(species_tree_file)
    taxa_data = load_taxa(taxa_file)
    
    gd_visualizer_main(sptree, gd_result_file, taxa_data)


