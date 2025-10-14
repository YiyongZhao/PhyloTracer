import matplotlib.pyplot as plt
import matplotlib.colors as colors
import re
import string
from __init__ import *


def realign_branch_length(Phylo_t1: object) -> object:
    """
    Realign the branch lengths of a phylogenetic tree for visualization purposes.
    The function ladderizes the tree, resolves polytomies, sorts descendants by support,
    and adjusts branch lengths so that the tree is visually balanced.

    Args:
        Phylo_t1 (object): The input phylogenetic tree object.

    Returns:
        object: The phylogenetic tree object with realigned branch lengths.
    """
    Phylo_t1.ladderize()
    Phylo_t1.resolve_polytomy(recursive=True)
    Phylo_t1.sort_descendants("support")
    max_depth = get_max_deepth(Phylo_t1)
    for node in Phylo_t1.traverse():
        if not node.is_root():
            node.dist = 1
            degree = node.get_distance(node.get_tree_root()) + 1
            depth = get_max_deepth(node)
            node.dist = max_depth - depth - degree
    clade_up = Phylo_t1.get_children()[0]
    clade_down = Phylo_t1.get_children()[1]
    difference = abs(get_max_deepth(clade_up) - get_max_deepth(clade_down)) + 1
    clade_up.dist += difference
    clade_down.dist += difference
    return Phylo_t1

def rejust_root_dist(sptree: object) -> None:
    """
    Adjust the branch lengths of the root's immediate children to balance the tree for visualization.
    The longer clade is assigned a minimal distance, while the shorter clade's distance is set based on its depth.

    Args:
        sptree (object): The input phylogenetic tree object (rooted).

    Returns:
        None
    """
    clade_up = sptree.get_children()[0]
    clade_down = sptree.get_children()[1]
    if len(clade_up) > len(clade_down):
        clade_up.dist = 1
        if clade_down.is_leaf():
            clade_down.dist = get_max_deepth(sptree) - 1
        else:
            clade_down.dist = get_max_deepth(sptree) - get_max_deepth(clade_down)
    else:
        clade_down.dist = 1
        if clade_up.is_leaf():
            clade_up.dist = get_max_deepth(sptree) - 1
        else:
            clade_up.dist = get_max_deepth(sptree) - get_max_deepth(clade_up)

    return sptree
#############################################################################################################
def dup_nodeids_from_numbered_gfs(phylo_t: object) -> tuple:
    """
    Root the tree if necessary, number all nodes, and find duplication nodes.
    Returns a new tree object with numbered nodes and a list of duplication node names.

    Args:
        phylo_t (object): The input phylogenetic tree object.

    Returns:
        tuple: (numbered_tree (object), dup_node_name_list (list))
    """
    if not is_rooted(phylo_t):
        phylo_t = root_tre_with_midpoint_outgroup(phylo_t)
    phylo_t = num_tre_node(phylo_t)
    dup_node_name_list = find_dup_node(phylo_t)
    phylo_t1 = Tree(phylo_t.write())
    num_tre_node(phylo_t1)
    return phylo_t1, dup_node_name_list

def create_tree_style(tree_style: str, tree_id: str, visual: bool = False) -> object:
    """
    Create and configure a TreeStyle object for tree visualization.

    Args:
        tree_style (str): The display mode for the tree (e.g., 'circular', 'rectangular').
        tree_id (str): The identifier for the tree, used in the title.
        visual (bool, optional): Whether to add duplication event legends. Defaults to False.

    Returns:
        object: Configured TreeStyle object.
    """
    ts = TreeStyle()
    ts.title.add_face(TextFace("★", fgcolor="white", ftype='Arial'), column=0)
    ts.title.add_face(TextFace(tree_id, ftype='Arial'), column=1)
    if visual:
        ts.title.add_face(TextFace("★", fgcolor="red", ftype='Arial'), column=0)
        ts.title.add_face(TextFace("Interspecific gene duplication event", ftype='Arial'), column=1)
        ts.title.add_face(TextFace("★", fgcolor="blue", ftype='Arial'), column=0)
        ts.title.add_face(TextFace("Intraspecific gene duplication event", ftype='Arial'), column=1)
    ts.mode = tree_style
    ts.scale = 20
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = 'black'
    ts.branch_vertical_margin = -1

    
    return ts

def set_node_style(node: object, dup_node_name_list: list, visual: bool = False) -> object:
    """
    Set the visual style for a tree node, highlighting duplication nodes if present.

    Args:
        node (object): The tree node to style.
        dup_node_name_list (list): List of duplication node names.
        visual (bool, optional): Whether to apply duplication event coloring. Defaults to False.

    Returns:
        object: Configured NodeStyle object.
    """
    nstyle = NodeStyle()
    splist = set(get_species_list(node))
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    nstyle["vt_line_type"] = 0
    nstyle["hz_line_type"] = 0
    nstyle["size"] = 0
    nstyle["shape"] = "circle"
    nstyle["fgcolor"] = "black"
    

    if visual:
        if node.name in dup_node_name_list and len(splist) == 1:
            node.add_face(TextFace("★", fsize=7, fgcolor="blue",ftype='Arial'), column=1, position="branch-top")
        elif node.name in dup_node_name_list and len(splist) != 1:
            node.add_face(TextFace("★", fsize=7, fgcolor="red",ftype='Arial'), column=1, position="branch-top")
    
    node.set_style(nstyle)

    
def get_treestyle(Phylo_t:object,tree_style:str,tre_ID:str,visual:bool=False)->object:
    Phylo_t1, dup_node_list = dup_nodeids_from_numbered_gfs(Phylo_t)
    ts = create_tree_style(tree_style,tre_ID,visual)
    dup_node_name_list=[node.name for node in dup_node_list]
    for node in Phylo_t1.traverse():
        set_node_style(node, dup_node_name_list,visual)
    
    return Phylo_t1, ts


####################################################################################################################
def get_color_dict(dictionary: dict) -> dict:
    """
    Generate a color mapping dictionary for unique values in the input dictionary using a rainbow colormap.
    Each key in the input dictionary is mapped to a string combining its value and the corresponding color hex code.

    Args:
        dictionary (dict): Input dictionary with keys to be colored and values as categories.

    Returns:
        dict: A dictionary mapping each key to a string formatted as 'category@#hexcolor'.
    """
    colormap = plt.get_cmap("rainbow")
    unique_values = set(dictionary.values())
    colors_list = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict = dict(zip(unique_values, colors_list))
    sps_color_list = {k: f"{v}@{color_dict.get(v)}" for k, v in dictionary.items() if v in color_dict}
    return sps_color_list


def generate_color_dict(gene_categories: list[dict]) -> list[dict]:
    """
    Generate a list of color mapping dictionaries for a list of gene category dictionaries.

    Args:
        gene_categories (list[dict]): List of dictionaries, each mapping gene identifiers to categories.

    Returns:
        list[dict]: List of color mapping dictionaries for each input dictionary.
    """
    return [get_color_dict(category_dict) for category_dict in gene_categories]


def generate_string(index: int) -> str:
    """
    Generate a string label based on the given index, using uppercase and lowercase English letters.
    For index >= 52, returns a string representation of the index.

    Args:
        index (int): The index to convert to a string label.

    Returns:
        str: The generated string label.
    """
    letters = list(string.ascii_uppercase) + list(string.ascii_lowercase)
    if index < len(letters):
        return "@" + letters[index]
    else:
        first_letter = letters[(index - len(letters)) // 52]
        second_letter = letters[(index - len(letters)) % 52]
        return "@" + first_letter + second_letter

def get_new_sorted_dict(gene2fam: dict) -> dict:
    """
    Generate a sorted dictionary mapping each unique family to its color and a unique suffix.
    The color is extracted from the color dictionary, and the suffix is generated by index.

    Args:
        gene2fam (dict): Dictionary mapping gene identifiers to family names.

    Returns:
        dict: Sorted dictionary mapping family names to color hex codes with unique suffixes.
    """
    uniq_fam = set(get_color_dict(gene2fam).values())
    fam2color = {i.split('@')[0]: i.split('@')[-1] for i in uniq_fam}
    sorted_dict = dict(sorted(fam2color.items(), key=lambda x: x[0], reverse=False))
    for index, key in enumerate(sorted_dict.keys()):
        suffix = generate_string(index)
        sorted_dict[key] += suffix
    return sorted_dict


def fuzzy_match(search_string: str, key: str) -> re.Match:
    """
    Perform a fuzzy (regex-based) match between the search string and the key.

    Args:
        search_string (str): The regex pattern to search for.
        key (str): The string to be searched.

    Returns:
        re.Match: The match object if a match is found, otherwise None.
    """
    return re.search(search_string, key)

def tips_mark(
    Phylo_t1: object,
    voucher2taxa_dic: dict,
    color_dicts: list,
    sps_color_dict: dict,
    tre_ID,
    ts,
    new_named_gene2gene_dic: dict,
    dir_path,
    gene_color_dict=None,
    df=None
) -> object:
    """
    Annotate the phylogenetic tree with colored faces for species and gene information.
    This function adds visual marks to tree nodes based on species and gene categories, using color dictionaries.

    Args:
        Phylo_t1 (object): The phylogenetic tree object to annotate.
        voucher2taxa_dic (dict): Mapping from voucher names to taxa labels.
        color_dicts (list): List of color mapping dictionaries for gene categories.
        sps_color_dict (dict): Color mapping dictionary for species.
        tre_ID: Tree identifier.
        ts: TreeStyle object for visualization.
        new_named_gene2gene_dic (dict): Mapping from renamed gene IDs to original gene IDs.
        dir_path: Directory path for output or resources.
        gene_color_dict (dict, optional): Color mapping dictionary for genes. Defaults to None.
        df (optional): DataFrame or additional data for annotation. Defaults to None.

    Returns:
        object: The annotated phylogenetic tree object.
    """
    def add_face_to_node(node: object, face: object, column: int, position: str = "aligned") -> None:
        """
        Add a face (visual annotation) to a tree node if it has not already been added at the specified column and position.

        Args:
            node (object): The tree node to which the face will be added.
            face (object): The face object to add (e.g., TextFace).
            column (int): The column index for the face.
            position (str): The position for the face (default is "aligned").
        """
        if (node, column, position) not in faces_added:
            node.add_face(face, column=column, position=position)
            faces_added.add((node, column, position))

    def generate_face_mark(node: object, species: str, column: int, color_dict: dict) -> None:
        """
        Generate and add a colored face mark for a species on a tree node.
        If the species is in the color dictionary, use its color; otherwise, use white.

        Args:
            node (object): The tree node to annotate.
            species (str): The species label for the node.
            column (int): The column index for the face.
            color_dict (dict): Mapping from species to color annotation strings.
        """
        if (node, column, "aligned") not in faces_added:
            if species in color_dict:
                color = color_dict[species].split('@')[-1]
                face = TextFace("  ▐" + '  ' + color_dict[species].split('@')[0], fgcolor=color, ftype='Arial')
                add_face_to_node(node, face, column, position="aligned")
            else:
                color = 'white'
                face = TextFace("  ▐" + '  ' + ''*len(species), fgcolor=color, ftype='Arial')
                add_face_to_node(node, face, column, position="aligned")


    def add_species_face(node: object, gene: str, species: str, sps_color_dict: dict) -> None:
        """
        Add colored faces for gene and species information to a tree node.
        The gene face is colored and italicized, and the species face uses the species color.

        Args:
            node (object): The tree node to annotate.
            gene (str): The gene label for the node.
            species (str): The species label for the node.
            sps_color_dict (dict): Mapping from species to color annotation strings.
        """
        if species in sps_color_dict:
            color = sps_color_dict[species].split('@')[-1]
            gene_face = TextFace(' ' + gene, fgcolor=color, ftype='Arial', fstyle='italic')
            node.add_face(gene_face, column=-1)
        if species in sps_color_dict:
            color = sps_color_dict[species].split('@')[-1]
            species_face = TextFace("  ▐" + '  ' + sps_color_dict[species].split('@')[0], fgcolor=color, ftype='Arial', fstyle='italic')
            add_face_to_node(node, species_face, 0, position="aligned")

    def add_gene_face(node: object, gene: str, column: int) -> None:
        """
        Add a gene label face to a tree node at the specified column.

        Args:
            node (object): The tree node to annotate.
            gene (str): The gene label for the node.
            column (int): The column index for the face.
        """
        matched_key = None
        matched_value = None
        for key in gene_color_dict:
            if fuzzy_match(key, gene):
                matched_key = key
                matched_value = gene_color_dict.get(key)
                break
        if matched_value:
            color = matched_value.split('@')[-1]
            face5 = TextFace("  ▐" + '  ' +  matched_value.split('@')[0], fgcolor=color,ftype='Arial')
            add_face_to_node(node, face5, column, position="aligned")
        

    def get_color(value: float) -> str:
        """
        Map a numeric value to a specific color hex code for heatmap visualization.
    
        Args:
            value (float): The numeric value to be mapped to a color. Can be NaN.
    
        Returns:
            str: The hex color code corresponding to the value. Returns 'white' if value is NaN.
        """
        if np.isnan(value): 
            return 'white'
        else:
            if 0 <= value <= 5:
                return '#006599'
            elif 5 < value <= 10:
                return '#408ca6'
            elif 10 < value <= 15:
                return '#7fb2b2'
            elif 15 < value <= 20:
                return '#bfd9bf'
            elif 20 < value <= 25:
                return '#ffffcc'
            elif 25 < value <= 30:
                return '#f7deab'
            elif 30 < value <= 35:
                return '#eebc88'
            elif 35 < value <= 40:
                return '#e69966'
            elif 40 < value <= 45:
                return '#dc7845'
            elif 45 < value <= 50:
                return '#d55623'
            else:
                return '#cc3300'
        
    def add_heat_map_to_node(tree:object, df:pd.DataFrame, new_named_gene2gene_dic: dict, start_col: int) -> None:
        """
        Add a heatmap to each node of the tree based on values from a DataFrame.
    
        Args:
            tree: The phylogenetic tree object whose nodes will be annotated.
            df: pandas.DataFrame containing the heatmap values, indexed by gene names.
            new_named_gene2gene_dic (dict): Mapping from node names to gene names.
            start_col (int): The starting column index for the heatmap faces.
    
        Returns:
            None
        """
        columns = df.columns.tolist()
        for node in tree:
            gene = new_named_gene2gene_dic[node.name]
            matched_key = None
            for key in df.index:
                if fuzzy_match(key, gene):
                    matched_key = key
                    break
            for ind, col_name in enumerate(columns):
                col_idx = start_col + ind
                if node.is_leaf() and matched_key in df.index:
                    color = get_color(df[col_name][matched_key])
                else:
                    color = '#F5F5F5'
                face = RectFace(width=10, height=10, fgcolor=color, bgcolor=color)
                node.add_face(face, column=col_idx, position='aligned')
                    
    def add_header_to_tree(ts: any, df: pd.DataFrame, new_start: int) -> None:
        """
        Add a rotated header to the tree visualization, labeling each column with the DataFrame's column names.
    
        Args:
            ts (any): The TreeStyle or visualization object to which the header will be added.
            df (pd.DataFrame): DataFrame whose columns will be used as header labels.
            new_start (int): The starting column index for header placement.
    
        Returns:
            None
        """
        labels = df.columns.to_list()
        for ind, i in enumerate(labels):
            face = TextFace(' ' + i, fgcolor='black', ftype='Arial', fsize=9)
            face.rotation = -90
            face.vt_align = 2
            ts.aligned_header.add_face(face, new_start + ind)
    
    
    def add_color_bar(ts: any) -> None:
        """
        Add a color bar legend to the tree visualization, showing the mapping from value ranges to colors.
    
        Args:
            ts (any): The TreeStyle or visualization object to which the color bar will be added.
    
        Returns:
            None
        """
        ts.legend.add_face(TextFace(' '), column=0)
        bar_face = TextFace('Color Bar ', ftype='Arial')
        ts.legend.add_face(bar_face, column=0)
        cols = ['#006599', '#408ca6', '#7fb2b2', '#bfd9bf', '#ffffcc', '#f7deab', '#eebc88', '#e69966', '#dc7845', '#d55623', '#cc3300']
        bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        col_dic = dict(zip(bounds, cols))
        n = 1
        for k, v in col_dic.items():
            colorbar_face = RectFace(width=20, height=20, fgcolor=v, bgcolor=v)
            ts.legend.add_face(TextFace(' ' + str(k)), column=n)
            ts.legend.add_face(colorbar_face, column=n)
            n += 1
        ts.legend_position = 2

    
    
    faces_added = set()  # used to track the added faces
    

    for node in Phylo_t1.traverse():
        if node.is_leaf():
            gene = new_named_gene2gene_dic[node.name]
            species = voucher2taxa_dic[node.name.split("_")[0]]
            rename_species=node.name.split("_")[0]
            add_species_face(node, gene,rename_species,sps_color_dict)
            column = 1
            for color_dict in color_dicts:
                generate_face_mark(node, gene, column, color_dict)
                column += 1

            if gene_color_dict is not None and gene in gene_color_dict:
                add_gene_face(node, gene, column)
            else:
                face_placeholder = TextFace("  ▐", fgcolor='white', ftype='Arial')
                add_face_to_node(node, face_placeholder, column, position="aligned")
            column += 1
    if df is not None:
        add_heat_map_to_node(Phylo_t1,df,new_named_gene2gene_dic,column)
        add_header_to_tree(ts,df,column)
        add_color_bar(ts)

    
    return Phylo_t1.render(dir_path+tre_ID+'.pdf',w=210, units="mm",tree_style=ts)


def get_matched_value(gene: str, gene2fam: dict) -> tuple:
    """
    Fuzzy match a gene name to the keys in gene2fam and return the matched key and its family value.

    Args:
        gene (str): The gene name to match.
        gene2fam (dict): Dictionary mapping gene names to family names.

    Returns:
        tuple: (matched_gene, matched_family) if found, otherwise (None, None).
    """
    for key, value in gene2fam.items():
        if fuzzy_match(gene, key):
            return key, value
    return None, None


def get_fam_dic(t: any, gene2fam: dict, new_named_gene2gene_dic: dict) -> dict:
    """
    Build a dictionary mapping each gene family to a list of node names belonging to that family.

    Args:
        t (any): The phylogenetic tree object to iterate over its nodes.
        gene2fam (dict): Dictionary mapping gene names to family names.
        new_named_gene2gene_dic (dict): Mapping from node names to gene names.

    Returns:
        dict: Dictionary mapping family names to lists of node names.
    """
    fam_dic = {}
    for i in t:
        gene = new_named_gene2gene_dic[i.name]
        match_gene, match_family = get_matched_value(gene, gene2fam)
        if match_family and match_family in gene2fam.values():
            fam_dic.setdefault(match_family, []).append(i.name)
    return fam_dic

def find_combinations(my_list: list) -> list:
    """
    Generate all unique pairwise combinations from a list.

    Args:
        my_list (list): The input list from which to generate combinations.

    Returns:
        list: List of tuple pairs, each representing a unique combination.
    """
    combinations = []
    for i in range(len(my_list)):
        for j in range(i+1, len(my_list)):
            combinations.append((my_list[i], my_list[j]))
    return combinations


def get_dup_family_dic(t: any, gene2fam: dict, new_named_gene2gene_dic: dict) -> dict:
    """
    Build a dictionary mapping each gene family to the set of node names representing duplication events.

    Args:
        t (any): The phylogenetic tree object.
        gene2fam (dict): Dictionary mapping gene names to family names.
        new_named_gene2gene_dic (dict): Mapping from node names to gene names.

    Returns:
        dict: Dictionary mapping family names to sets of node names (duplication nodes).
    """
    fam_node_dic = {}
    fam_dic = get_fam_dic(t, gene2fam, new_named_gene2gene_dic)
    for k, v in fam_dic.items():
        nodes = set()
        com = find_combinations(v)
        for i in com:
            clade = t.get_common_ancestor(i)
            nodes.add(clade.name)
        fam_node_dic[k] = nodes
    return fam_node_dic


def mapping_sptree(t: any, sptree: any, sp_node_dic: dict, gene2fam: dict, new_named_gene2gene_dic: dict) -> None:
    """
    Map duplication nodes from the gene tree to the species tree based on shared species.

    Args:
        t (any): The gene phylogenetic tree object.
        sptree (any): The species phylogenetic tree object.
        sp_node_dic (dict): Dictionary to store mapping from species tree nodes to gene families.
        gene2fam (dict): Dictionary mapping gene names to family names.
        new_named_gene2gene_dic (dict): Mapping from node names to gene names.

    Returns:
        None
    """
    fam_node_dic = get_dup_family_dic(t, gene2fam, new_named_gene2gene_dic)
    for k, v in fam_node_dic.items():
        for i in v:
            clade = t & i
            sps = get_species_list(clade)
            uniq_sps = set(sps)
            clade2sptree = sptree.get_common_ancestor(uniq_sps)
            sp_node_dic.setdefault(clade2sptree.name, set()).add(k)


def get_sptree_style(sorted_dict: dict) -> any:
    """
    Create and configure a TreeStyle object for species tree visualization with legend.

    Args:
        sorted_dict (dict): Dictionary mapping family names to color and label strings.

    Returns:
        any: Configured TreeStyle object for visualization.
    """
    ts = TreeStyle()
    ts.legend_position = 1
    ts.mode = 'r'
    ts.scale = 30
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = True
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = 'black'
    ts.branch_vertical_margin = -1
    for k, v in sorted_dict.items():
        ts.legend.add_face(TextFace(v.split('@')[1] + ' ' + k, fsize=20, fgcolor=v.split('@')[0]), column=0)
    return ts

def create_species_mappings(dict_list: list[dict[str, str]]) -> list[dict[str, str]]:
    """
    Create mappings from species to family, order, and clade based on gene-level dictionaries.

    Args:
        dict_list (list[dict[str, str]]):
            A list containing four dictionaries: gene2sps, gene2family, gene2order, gene2clade.
            Each dictionary maps gene names to their respective species, family, order, or clade.

    Returns:
        list[dict[str, str]]: A list containing three dictionaries:
            - sps2family: species to family mapping
            - sps2order: species to order mapping
            - sps2clade: species to clade mapping
    """
    gene2sps, gene2family, gene2order, gene2clade = dict_list
    sps2family = {}
    sps2order = {}
    sps2clade = {}
    for gene, sps in gene2sps.items():
        if gene in gene2family:
            sps2family[sps] = gene2family[gene]
        if gene in gene2order:
            sps2order[sps] = gene2order[gene]
        if gene in gene2clade:
            sps2clade[sps] = gene2clade[gene]
    return [sps2family, sps2order, sps2clade]

def mark_gene_to_sptree(
    sptree,
    sp_node_dic: dict,
    gene_categories: list,
    sorted_dict: dict,
    gene2sps: dict,
    voucher2taxa_dic: dict
) -> None:
    """
    Mark gene categories on the species tree by adding colored faces to nodes.

    Args:
        sptree: The species tree object (e.g., ete3.Tree).
        sp_node_dic (dict): Mapping from species names to tree nodes.
        gene_categories (list): List of gene category dictionaries.
        sorted_dict (dict): Dictionary for sorting or grouping genes.
        gene2sps (dict): Mapping from gene names to species names.
        voucher2taxa_dic (dict): Mapping from voucher codes to taxa names.

    Returns:
        None
    """
    gene_categories_1 = gene_categories.copy()
    gene_categories_1.insert(0, gene2sps)
    gene_categories_2 = create_species_mappings(gene_categories_1)
    color_dicts = generate_color_dict(gene_categories_2)
    faces_added = set()  # used to track the added faces
    def add_face_to_node(node, face, column, position="aligned"):
        if (node, column, position) not in faces_added:
            node.add_face(face, column=column, position=position)
            faces_added.add((color_dict[node.name].split('@')[0], column, position))
    def generate_face_mark(node, species, column, color_dict):
        if (color_dict[node.name].split('@')[0], column, "aligned") not in faces_added:
            if species in color_dict:
                color = color_dict[species].split('@')[-1]
                face = TextFace("   ▐" + '  ' + color_dict[species].split('@')[0], fgcolor=color,ftype='Arial', fstyle='italic')
                add_face_to_node(node, face, column, position="aligned")
        else:
            color = color_dict[species].split('@')[-1]
            face = TextFace("   ▐" + '  ' , fgcolor=color,ftype='Arial', fstyle='italic')
            add_face_to_node(node, face, column, position="aligned")
                
    for i in sptree.traverse():
        if not i.is_leaf():
            for k, v in sp_node_dic.items():
                if i.name == k:
                    n = len(v)
                    for index, value in enumerate(sorted(v)):
                        position = 'branch-top' if n == 1 or index < n / 2 else 'branch-bottom'
                        column = index if n == 1 or index < n / 2 else index - n / 2
                        
                        i.add_face(TextFace(sorted_dict[value].split('@')[1], fgcolor=sorted_dict[value].split('@')[0]), column=column, position=position)
        else:
            column = 1
            species=i.name
            for color_dict in color_dicts:
                generate_face_mark(i, species, column, color_dict)
                column += 1


def rename_sptree(sptree) -> None:
    """
    Set the style for each node in the species tree for visualization.

    Args:
        sptree: The species tree object (e.g., ete3.Tree).

    Returns:
        None
    """
    for i in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0
        nstyle["hz_line_type"] = 0
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        i.set_style(nstyle)
        #if i.name in sp:
        #    i.name=sp[i.name]

def view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic,gene_categories,tree_style,keep_branch,new_named_gene2gene_dic,gene2fam=None,df=None,visual:bool=False):
    dir_path = os.path.join(os.getcwd(), "tree_visualizer/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    color_dicts = generate_color_dict(gene_categories)
    sps_color_dict=get_color_dict(voucher2taxa_dic)
    gene_color_dict = None  # 保证变量总是有定义
    if gene2fam is not None:
        gene_color_dict=get_color_dict(gene2fam)
    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID,tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        Phylo_t0=read_phylo_tree(tre_path)
        Phylo_t0=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        Phylo_t1,ts=get_treestyle(Phylo_t0,tree_style,tre_ID,visual)
        Phylo_t1.ladderize()
        Phylo_t1.resolve_polytomy(recursive=True)
        Phylo_t1.sort_descendants("support")
        if keep_branch !='1' :
            realign_branch_length(Phylo_t1)
            rejust_root_dist(Phylo_t1)

        tips_mark(Phylo_t1,voucher2taxa_dic,color_dicts,sps_color_dict,tre_ID,ts,new_named_gene2gene_dic,dir_path,gene_color_dict,df)
        pbar.update(1)
    pbar.close()

def mark_gene_to_sptree_main(tre_dic,gene_categories,sptree,gene2fam,gene2sps,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic):
    sorted_dict=get_new_sorted_dict(gene2fam)
    num_tre_node(sptree)
    sp_node_dic={}
    for k,v in tre_dic.items():
        t=read_phylo_tree(v)
        num_tre_node(t)
        t1=rename_input_tre(t,gene2new_named_gene_dic)
        mapping_sptree(t1,sptree,sp_node_dic,gene2fam,new_named_gene2gene_dic)

    sp2=rename_input_tre(sptree,voucher2taxa_dic)

    mark_gene_to_sptree(sp2,sp_node_dic,gene_categories,sorted_dict,gene2sps,voucher2taxa_dic)

    rename_sptree(sp2)
    ts=get_sptree_style(sorted_dict)
    realign_branch_length(sp2)
    clade_up=sp2.get_children()[0]
    clade_down=sp2.get_children()[1]
    clade_up.dist=1
    clade_down.dist=get_max_deepth(clade_up)
    sp2.render(file_name='genefamily_map2_sptree.pdf',tree_style=ts)

####################################################################################
if __name__ == "__main__":
    os.makedirs(os.path.join(os.getcwd(), "pdf_result"))
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    #gene2fam=read_and_return_dict('gene2fam')
    sp2order=read_and_return_dict('order')
    sp2family=read_and_return_dict('genus')
    tre_dic=read_and_return_dict('GF.txt')
    gene_categories=[]
    gene_categories.append(sp2family)
    gene_categories.append(sp2order)
    tree_style='r'
    keep_branch =1 
    view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic,gene_categories,tree_style,keep_branch,new_named_gene2gene_dic)
    
    
