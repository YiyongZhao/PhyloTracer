from __init__ import *
from PIL import Image, ImageDraw,ImageFont
import os
import math

def get_only_sps_tree(Phylo_t: object) -> object:
    """
    Creates a new tree where leaf names are simplified to only include the species code.

    Args:
        Phylo_t (object): The input phylogenetic tree object.

    Returns:
        object: A new tree object with simplified leaf names.
    """
    Phylo_t_c = Phylo_t.copy()
    for node in Phylo_t_c:
        node.name = node.name.split('_')[0]
    return Phylo_t_c
    
def get_max_tree(trees: list[object]) -> object:
    """
    Finds the tree with the maximum number of leaves from a list of trees.

    Args:
        trees (list[object]): A list of tree objects.

    Returns:
        object | None: The tree object with the maximum number of leaves, or None if the list is empty.
    """
    max_leaf_count = 0
    max_leaf_count_tree = None

    for tree in trees:
        leaf_count = len(tree.get_leaves())
        if leaf_count > max_leaf_count:
            max_leaf_count = leaf_count
            max_leaf_count_tree = tree

    return max_leaf_count_tree

def group_trees_by_topology_with_ids(trees_with_ids: list[tuple], result_dict: dict) -> list[tuple]:
    """
    Groups a list of trees with IDs based on their topology using Robinson-Foulds distance.

    Args:
        trees_with_ids (list[tuple]): A list of (tree_id, tree_object) tuples to be grouped.
        result_dict (dict): A dictionary to store the grouped trees,
                           where keys are tree Newick strings and values are lists of (tree_id, tree) tuples.

    Returns:
        list[tuple]: A list of (tree_id, tree) tuples that did not group with the first tree.
    """
    if len(trees_with_ids) == 1:
        tree_id, tree = trees_with_ids[0]
        result_dict[tree.write(format=9)] = [trees_with_ids[0]]
        return []
    else:
        same_rf_trees = []
        different_trees = []
        first_tree_id, first_tree = trees_with_ids[0]
        
        for tree_id, tree in trees_with_ids:
            rf = first_tree.robinson_foulds(tree)[0]
            if rf == 0:
                same_rf_trees.append((tree_id, tree))
            else:
                different_trees.append((tree_id, tree))

        # 选择叶子数最多的树作为代表
        max_tree_tuple = max(same_rf_trees, key=lambda x: len(x[1].get_leaves()))
        result_dict[max_tree_tuple[1].write(format=9)] = same_rf_trees
        return different_trees

def process_tree_with_ids(trees_with_ids: list[tuple], result_dict: dict):
    """
    Recursively processes a list of trees with IDs to group them by topology.

    Args:
        trees_with_ids (list[tuple]): A list of (tree_id, tree_object) tuples to be processed.
        result_dict (dict): A dictionary to store the grouped trees.
    """
    remaining_trees = trees_with_ids
    while len(remaining_trees) >= 1:
        remaining_trees = group_trees_by_topology_with_ids(remaining_trees, result_dict)

def write_relative_summary_with_ids(outfile: str, dic: dict, voucher2taxa_dic: dict, top_n: int = None):
    """
    Writes a summary of relative tree topologies with tree IDs to a file and optionally visualizes top N.
    Each unique topology is represented by only one line in the output.

    Args:
        outfile (str): The base name for the output file.
        dic (dict): A dictionary where keys are tree Newick strings and values are lists of (tree_id, tree) tuples.
        voucher2taxa_dic (dict): A dictionary mapping voucher names to taxa labels for renaming.
        top_n (int | None): The number of top topologies to visualize, or None to skip visualization.
    """
    sorted_dict = dict(sorted(dic.items(), key=lambda item: len(item[1]), reverse=True))
    
    if top_n:
        output_path = f'merge_relative_top{top_n}.png'
        # 需要修改visualize_top_trees函数来处理新的数据结构
        # visualize_top_trees(sorted_dict, output_path, voucher2taxa_dic, top_n)
    
    with open(f'relative_{outfile}.txt', 'w') as file:
        file.write('Topology_ID\tTopology_num\tTopology\tTree_IDs\n')
        for index, clade in enumerate(sorted_dict.keys()):
            trees_with_ids = sorted_dict[clade]
            
            # 选择该拓扑结构中叶子数最多的树作为代表
            representative_tree_tuple = max(trees_with_ids, key=lambda x: len(x[1].get_leaves()))
            _, representative_tree = representative_tree_tuple
            
            # 重命名树的叶子节点
            tree_renamed = rename_input_tre(representative_tree, voucher2taxa_dic)
            tree_str = tree_renamed.write(format=9)
            
            # 收集所有属于该拓扑结构的树的ID
            tree_ids = [tree_id for tree_id, _ in trees_with_ids]
            tree_ids_str = ','.join(tree_ids)
            
            # 输出：拓扑ID、拓扑数量、拓扑结构、树ID列表
            file.write(f'{index}\t{len(trees_with_ids)}\t{tree_str}\t{tree_ids_str}\n')

def write_absolute_summary(outfile: str, trees_with_ids: list[tuple], voucher2taxa_dic: dict, top_n: int = None):
    """
    Groups trees by absolute topology (after folding same taxa), writes summary to file, and optionally visualizes top N topologies.

    Args:
        outfile (str): Output file base name.
        trees_with_ids (list): List of tuples (tree_id, tree_object) to be grouped.
        voucher2taxa_dic (dict): Mapping for renaming tree leaves.
        top_n (int | None): Number of top topologies to visualize, or None to skip visualization.
    """
    topolo_dic = {}
    for tree_id, tree in trees_with_ids:
        
        tree_str = tree.write(format=9)
        if tree_str in topolo_dic:
            topolo_dic[tree_str].append((tree_id, tree))
        else:
            topolo_dic[tree_str] = [(tree_id, tree)]

    # Sort by number of trees per topology, then by topology string length
    sorted_dict = dict(sorted(topolo_dic.items(), key=lambda item: (len(item[1]), len(item[0])), reverse=True))

    if top_n:
        output_path = f'merge_absolutely_top{top_n}.png'
        # 需要相应修改visualize_top_trees函数来处理新的数据结构
        visualize_top_trees(sorted_dict, output_path, voucher2taxa_dic, top_n)

    with open(f'absolute_{outfile}.txt', 'w') as file:
        file.write('Topology_ID\tTopology_num\tTopology\tTree_IDs\n')
        for num, k in enumerate(sorted_dict):
            t1 = Tree(k)
            t2 = rename_input_tre(t1, voucher2taxa_dic)
            t1_str = t2.write(format=9)
            # 提取所有tree id
            tree_ids = [tree_id for tree_id, _ in sorted_dict[k]]
            tree_ids_str = ','.join(tree_ids)
            file.write(f'{num}\t{len(sorted_dict[k])}\t{t1_str}\t{tree_ids_str}\n')

def visualize_top_trees(tree_count_dict: dict, output_path: str, voucher2taxa_dic: dict, top_n: int = 10):
    """
    Visualizes the top N tree topologies by rendering each as an image and combining them into a grid.

    Args:
        tree_count_dict (dict): Dictionary where keys are tree Newick strings and values are lists of trees.
        output_path (str): Path to save the combined image.
        voucher2taxa_dic (dict): Mapping for renaming tree leaves.
        top_n (int, optional): Number of top topologies to visualize. Defaults to 10.
    """
    cols = math.ceil(math.sqrt(top_n))
    rows = math.ceil(top_n / cols)
    sorted_trees = list(tree_count_dict.items())[:top_n]
    temp_dir = "temp_trees"
    os.makedirs(temp_dir, exist_ok=True)
    tree_images = []
    for i, (tree_str, count) in enumerate(sorted_trees):
        tree0 = read_tree(tree_str)
        tree = rename_input_tre(tree0, voucher2taxa_dic)
        tree.ladderize()
        tree.resolve_polytomy(recursive=True)
        tree.sort_descendants("support")

        realign_branch_length(tree)
        rejust_root_dist(tree)

        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0
        nstyle["hz_line_type"] = 0
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        for node in tree.traverse():
            node.set_style(nstyle)
        ts = TreeStyle()
        ts.scale = 10
        ts.show_leaf_name = True
        ts.title.add_face(TextFace(f" Count: {len(count)}"), column=1)
        ts.show_scale = False
        tree_path = os.path.join(temp_dir, f"tree_{i + 1}.png")
        tree.render(tree_path, h=3200, tree_style=ts)
        tree_images.append(tree_path)
    if tree_images:
        with Image.open(tree_images[0]) as img:
            single_image_width, single_image_height = img.size
    else:
        single_image_width, single_image_height = 1600, 1600
    final_width = cols * single_image_width
    final_height = rows * single_image_height
    combined_image = Image.new("RGB", (final_width, final_height), "white")
    for idx, img_path in enumerate(tree_images):
        with Image.open(img_path) as img:
            row, col = divmod(idx, cols)
            x_offset = col * single_image_width
            y_offset = row * single_image_height
            combined_image.paste(img, (x_offset, y_offset))
    combined_image.save(output_path)
    for img_path in tree_images:
        os.remove(img_path)
    os.rmdir(temp_dir)

def statistical_main(
    tre_dic: dict,
    outfile: str,
    gene2new_named_gene_dic: dict,
    voucher2taxa_dic: dict,
    top_n: int
) -> None:
    """
    Main statistical analysis function for summarizing tree topologies.

    Args:
        tre_dic (dict): Dictionary mapping tree identifiers to tree file paths or objects.
        outfile (str): Base name for output files.
        gene2new_named_gene_dic (dict): Mapping from gene IDs to renamed gene IDs.
        voucher2taxa_dic (dict): Mapping from voucher names to taxon labels.
        top_n (int): Number of top topologies to visualize and summarize.

    Returns:
        None
    """
    trees_with_ids = []  # 改为保存(tree_id, tree_object)的列表
    for k, v in tre_dic.items():
        t = read_tree(v)
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        t1.sort_descendants()
       
        t2 = get_only_sps_tree(t1)
        if len(get_species_set(t2))==1:
            continue
        # Check if it's a single-copy gene tree, exit if not
        if len(get_species_set(t2)) != len(t2):
            print(f"Error: Gene tree {k} is not a single-copy gene tree!")
            print(f"Number of species: {len(get_species_set(t2))}, Number of genes: {len(t2)}")
            print("This program requires single-copy gene trees as input. Program terminated.")
            exit(1)
        trees_with_ids.append((k, t2))  # 保存tree id和tree对象的元组
    
    write_absolute_summary(outfile, trees_with_ids, voucher2taxa_dic, top_n)
    
    # 修改相对拓扑结构统计部分
    dic_with_ids = {}
    process_tree_with_ids(trees_with_ids, dic_with_ids)
    write_relative_summary_with_ids(outfile, dic_with_ids, voucher2taxa_dic, top_n)


if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt')   
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    outfile='result'
    statistical_main(tre_dic,outfile,gene2new_named_gene_dic,new_named_gene2gene_dic)
