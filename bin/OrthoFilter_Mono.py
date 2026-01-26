from __init__ import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from PyPDF4 import PdfFileReader, PdfFileWriter
from BranchLength_NumericConverter import trans_branch_length,write_tree_to_newick


def rename_input_single_tre(Phylo_t: object, taxa_dic: dict, new_named_gene2gene_dic: dict) -> object:
    """
    Rename the leaf nodes of a phylogenetic tree by mapping gene names to species names and appending the original node name.
    Args:
        Phylo_t (object): The phylogenetic tree object to be renamed.
        taxa_dic (dict): A dictionary mapping gene names to species names.
        new_named_gene2gene_dic (dict): A dictionary mapping node names to gene names.
    Returns:
        object: The phylogenetic tree object with renamed leaf nodes.
    """
    for node in Phylo_t:
        gene = new_named_gene2gene_dic.get(node.name, node.name)
        species = taxa_dic.get(gene, None)
        if species:
            node.name = species + '_' + node.name
    return Phylo_t

def rename_output_tre(Phylo_t: object, new_named_gene2gene_dic: dict) -> object:
    """
    Restore the original gene names for the leaf nodes of a phylogenetic tree after processing.
    Args:
        Phylo_t (object): The phylogenetic tree object whose leaf nodes will be renamed.
        new_named_gene2gene_dic (dict): A dictionary mapping processed node names back to original gene names.
    Returns:
        object: The phylogenetic tree object with restored leaf node names.
    """
    for node in Phylo_t:
        parts = node.name.split("_")
        new_str = "_".join(parts[1:]) if len(parts) > 1 else node.name
        node.name = new_named_gene2gene_dic.get(new_str, new_str)
    return Phylo_t

def is_multi_copy_tree(Phylo_t: object) -> bool:
    """
    Determine whether the phylogenetic tree is a multi-copy gene tree.
    Returns True if the number of leaf nodes does not equal the number of unique species.
    Args:
        Phylo_t (object): The phylogenetic tree object to check.
    Returns:
        bool: True if the tree is a multi-copy gene tree, False otherwise.
    """
    leafs = Phylo_t.get_leaf_names()
    uniq_species = get_species_set(Phylo_t)
    if len(leafs) != len(uniq_species):
        return True
    return False

def set_style(Phylo_t:object, color_dict:dict, new_named_gene2gene_dic:dict) -> object:
    """
    Set the visual style for each node in the phylogenetic tree for visualization purposes.
    For each node, applies a default style (black circle, size 0). For leaf nodes, adds a colored and italicized label based on species and gene name.
    
    Args:
        Phylo_t (object): The phylogenetic tree object to be styled.
        color_dict (dict): A dictionary mapping gene names to color strings (e.g., 'geneA': '#FF0000').
        new_named_gene2gene_dic (dict): A dictionary mapping processed node names to original gene names.
    
    Returns:
        object: The styled phylogenetic tree object, ready for visualization.
    """
    for node in Phylo_t.traverse():
        nstyle=NodeStyle()
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        node.set_style(nstyle)
        if  node.is_leaf():
            parts = node.name.split("_")
            species_name=parts[0]
            new_str = "_".join(parts[1:]) if len(parts) > 1 else node.name
            gene=new_named_gene2gene_dic.get(new_str, new_str)
            color_info = color_dict.get(gene, None)
            if color_info:
                color = color_info.split('*')[-1]
                face = TextFace(species_name+'_'+gene, fgcolor=color,fstyle='italic')
                node.add_face(face, column=0)
    return Phylo_t

def get_color_dict(dictory:dict) -> dict:
    """
    Generate a color dictionary for unique values in the input dictionary, mapping each value to a unique color string.
    
    Args:
        dictory (dict): Input dictionary whose values will be assigned colors.
    
    Returns:
        dict: A dictionary mapping each key to a color string (e.g., 'geneA': '#FF0000').
    """
    colormap = plt.get_cmap("gist_rainbow")
    unique_values=set(dictory.values())
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_dict = {k: v + '*' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}
    return sps_color_dict

def get_single_clades(Phylo_t:object, empty_set:set) -> None:
    """
    Recursively find and collect all clades (subtrees) in the phylogenetic tree that contain only a single species.
    
    Args:
        Phylo_t (object): The phylogenetic tree object to search.
        empty_set (set): A set to store clades containing only one species.
    
    Returns:
        None
    """
    if calculate_species_num(Phylo_t)==1:
        empty_set.add(Phylo_t)
        return
    for child in Phylo_t.get_children():
        get_single_clades(child,empty_set)                        

def get_node_single_taxa_dict(Phylo_t:object) -> dict:
    """
    Generate a dictionary mapping each single-taxa (species) to the list of clades (subtrees) containing only that species.
    
    Args:
        Phylo_t (object): The phylogenetic tree object to search.
    
    Returns:
        dict: A dictionary mapping species name to a list of clades containing only that species.
    """
    def get_single_taxa_caldes(empty_set):
        single_taxa_dict = {}
        for clade in empty_set:
            single_taxa = get_species_list(clade)[0]
            single_taxa_dict.setdefault(single_taxa, []).append(clade)
        return single_taxa_dict   
    empty_set=set()
    get_single_clades(Phylo_t,empty_set)
    single_taxa_dict =get_single_taxa_caldes(empty_set)      
    sorted_dict = dict(sorted(single_taxa_dict.items(), key=lambda item: len(item[1]), reverse=True))
    return sorted_dict                         

def is_single_tree(Phylo_t:object) -> bool:
    """
    Determine whether the given phylogenetic tree contains only a single species.
    
    Args:
        Phylo_t (object): The phylogenetic tree object to check.
    
    Returns:
        bool: True if the tree contains only one species, False otherwise.
    """
    single_taxa_dict = get_node_single_taxa_dict(Phylo_t)
    # if 'basal angiosperms'  in single_taxa_dict.keys():
    #     single_taxa_dict.pop('basal angiosperms')
    # if 'Outgroup' in single_taxa_dict.keys():
    #     single_taxa_dict.pop('Outgroup')
    def check_single(dictionary):
        """
        Check if all values in the dictionary are the same, indicating a single unique value.
        
        Args:
            dictionary (dict): The dictionary to check.
        
        Returns:
            bool: True if all values are the same, False otherwise.
        """
        for value in dictionary.values():
            if len(value) != 1:
                return False
        return True

    return check_single(single_taxa_dict)  

def merge_pdfs_side_by_side(file1: str, file2: str, output_file: str) -> None:
    """
    Merge two PDF files side by side into a single output PDF file.
    
    Args:
        file1 (str): Path to the first PDF file.
        file2 (str): Path to the second PDF file.
        output_file (str): Path to the output merged PDF file.
    
    Returns:
        None
    """
    pdf_writer = PdfFileWriter()
    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        pdf1 = PdfFileReader(f1)
        pdf2 = PdfFileReader(f2)
        page1 = pdf1.getPage(0)
        page2 = pdf2.getPage(0)
        new_page_width = page1.mediaBox.getWidth() + page2.mediaBox.getWidth()
        new_page_height = max(page1.mediaBox.getHeight(), page2.mediaBox.getHeight())
        new_page = pdf_writer.addBlankPage(new_page_width, new_page_height)
        new_page.mergeTranslatedPage(page1, 0, 0)
        new_page.mergeTranslatedPage(page2, page1.mediaBox.getWidth(), 0)
        with open(output_file, 'wb') as output:
            pdf_writer.write(output)

def calculate_avg_length(node: object) -> int:
    """
    Calculate the average branch length for all descendant branches of the given node.
    
    Args:
        node (object): The node whose descendant branch lengths are to be averaged.
    
    Returns:
        int: The average branch length of all descendant branches.
    """
    total_length = 0.0
    leaf_count = 0
    for leaf in node.iter_leaves():
        total_length += node.get_distance(leaf)
        leaf_count += 1
    if leaf_count > 0:
        avg_length = total_length / leaf_count
    else:
        avg_length = 0.0
    return avg_length+(node.dist/leaf_count)
    
def calculate_branch_length_relative_score(node: object) -> int:
    """
    Calculate a relative score based on the branch lengths of the given node and its descendants.
    
    Args:
        node (object): The node for which to calculate the relative branch length score.
    
    Returns:
        int: The calculated relative score based on branch lengths.
    """
    if not node.is_leaf():
        avg_length=calculate_avg_length(node)
        sister = node.get_sisters()[0] if not node.is_root() and node.get_sisters() else None
        if sister and not sister.is_leaf():
            sister_length = calculate_avg_length(sister)
        else:
            sister_length = sister.dist if sister else 0.0
        return (avg_length / sister_length) if sister_length != 0 else 0.0
    else:
        branch_length=node.dist
        sister = node.get_sisters()[0] if not node.is_root() and node.get_sisters() else None
        if sister and not sister.is_leaf():
            sister_length = calculate_avg_length(sister)
        else:
            sister_length = sister.dist if sister else 0.0
        return (branch_length / sister_length) if sister_length != 0 else 0.0

def calculate_insertion_index(node):
    """
    Calculate the insertion index for a given node, typically used for determining the position to insert a node in a tree structure.
    
    Args:
        node (object): The node for which to calculate the insertion index.
    
    Returns:
        int: The calculated insertion index.
    """
    insertion_index = 1.0 
    current_node = node
    while current_node.up:
        current_tips = len(current_node.get_leaves())
        parent_node = current_node.up
        parent_tips = len(parent_node.get_leaves())
        
        if parent_tips == 0:
            insertion_index = 0.0  
            break
        current_insertion_index = len(node.get_leaves())/ current_tips
        insertion_index = min(insertion_index, current_insertion_index)
        
        current_node = parent_node  
    
    return insertion_index

  
def prune_single(Phylo_t):
    """
    Prune (remove) all single-species clades from the given phylogenetic tree.
    
    Args:
        Phylo_t (object): The phylogenetic tree object to be pruned.
    
    Returns:
        object: The pruned phylogenetic tree object.
    """
    rm_list=[]
    single_taxa_dict=get_node_single_taxa_dict(Phylo_t)
    for k,v in single_taxa_dict.items():
        if len(v) ==1 :
            pass
        elif len(v)==2:
            if len(v[0])>len(v[1]):
                leafs=v[1].get_leaf_names()
                total_leafs=Phylo_t.get_leaf_names()
                diff = [a for a in total_leafs if a not in set(leafs)]
                Phylo_t.prune(diff,preserve_branch_length=True)

            else:
                leafs=v[0].get_leaf_names()
                total_leafs=Phylo_t.get_leaf_names()
                diff = [a for a in total_leafs if a not in set(leafs)]
                Phylo_t.prune(diff,preserve_branch_length=True)
        else:
            result=[]
            for i in v :
                insertion_index=calculate_insertion_index(i)
                result.append(insertion_index)
            min_score=min(result)
            d=[]
            for j in v :
                insertion_index1=calculate_insertion_index(j)
                if insertion_index1 ==min_score:
                    branch_length_relative_score=calculate_branch_length_relative_score(j)
                    d.append(branch_length_relative_score)
            min_score1=min(d)
            for h in v :
                branch_length_relative_score1=calculate_branch_length_relative_score(h)
                if branch_length_relative_score1 ==min_score1:
                    rm_list.append(h.name)
                    
    for i in rm_list :
        clade=Phylo_t&i
        if clade.is_leaf():
            clade.delete()
        else:
            leafs=clade.get_leaf_names()
            total_leafs=Phylo_t.get_leaf_names()
            diff = [a for a in total_leafs if a not in set(leafs)]
            Phylo_t.prune(diff,preserve_branch_length=True)

def get_root_relative_branch_ratio(leaf:object, avg_length:int) -> int:
    """
    Calculate the relative branch length ratio of a leaf node compared to the average branch length.

    Args:
        leaf (object): The leaf node whose branch length is to be compared.
        avg_length (int): The average branch length among all descendants.

    Returns:
        int: The relative ratio of the leaf's branch length to the average branch length.
    """
    branch_length=leaf.dist
    return (branch_length-avg_length)/avg_length

def get_sister_relative_branch_ratio(node: object, sister: object) -> int:
    """
    Calculate the relative branch length ratio between a node and its sister node.

    Args:
        node (object): The node whose branch length is to be compared.
        sister (object): The sister node for comparison.

    Returns:
        int: The relative ratio of the node's branch length to its sister's branch length.
    """
    if node.is_leaf():
        branch_length = node.dist
    else:
        branch_length = calculate_avg_length(node) if node else 0.0

    if sister and not sister.is_leaf():
        sister_length = calculate_avg_length(sister)
    else:
        sister_length = sister.dist if sister else 0.0

    if sister_length != 0:
        return (branch_length - sister_length) / sister_length
    else:
        return 0.0

    
def calculate_insertion_depth(clade: object, node: object) -> int:
    """
    Calculate the insertion depth of a node within a given clade (subtree).

    Args:
        clade (object): The clade (subtree) in which the node is inserted.
        node (object): The node whose insertion depth is to be calculated.

    Returns:
        int: The topological distance from the clade root to the node.
    """
    if clade is None:
        return 0
    return clade.get_distance(node, topology_only=True)

def calculate_insertion_coverage(clade: object, node: object) -> float:
    """
    Calculate the coverage ratio of a node within a given clade (subtree).

    Args:
        clade (object): The clade (subtree) containing the node.
        node (object): The node whose coverage is to be calculated.

    Returns:
        float: The ratio of the number of leaves in the node to the number of leaves in the clade.
    """
    if clade is None or node is None:
        return 0.0  
    return len(node) / len(clade)

def calculate_insertion(clade:object, node:object)->int:
    """
    Calculate the insertion score of a node within a clade, based on depth and coverage.

    Args:
        clade (object): The clade (subtree) in which the node is inserted.
        node (object): The node to be evaluated.

    Returns:
        int: The insertion score, calculated as depth divided by coverage.
    """
    depth = calculate_insertion_depth(clade, node)
    coverage = calculate_insertion_coverage(clade, node)
    return depth / coverage if coverage != 0 else 0

def get_target_clade(clade:object)->object:
    """
    Find the target clade containing exactly two species, traversing from the given clade up to the root.

    Args:
        clade (object): The starting clade (subtree) for the search.

    Returns:
        object: The ancestor clade containing exactly two species, or None if not found.
    """
    node2root = clade.get_ancestors()
    target_clade = None
    for ancestor in node2root:
        if len(get_species_set(ancestor)) > 2:
            break
        elif len(get_species_set(ancestor)) == 2:
            target_clade = ancestor      
    return target_clade

def is_ancestor_sister_same(node:object)->bool:
    """
    Determine whether the sister of the node and the sister of its ancestor have the same species name.

    Args:
        node (object): The node to be checked.

    Returns:
        bool: True if the sister of the node and the sister of its ancestor have the same species name, otherwise False.
    """
    sps_list = get_species_list(node)
    sps = sps_list[0] if sps_list else None
    sisters = node.get_sisters()
    sis = sisters[0] if sisters else None
    sis_name_list = get_species_list(sis) if sis else []
    sis_name = sis_name_list[0] if sis_name_list else None
    if sps is None or sis_name is None or sps == sis_name:
        return False
    an = node.up
    if not an or an.is_root():
        return False
    an_sisters = an.get_sisters()
    an_si = an_sisters[0] if an_sisters else None
    an_name_list = get_species_list(an_si) if an_si else []
    an_name = an_name_list[0] if an_name_list else None
    return sis_name == an_name

def get_tips_avg_length(Phylo_t:object)->int:
    """
    Calculate the average distance from the root to all leaf nodes (tips) in the phylogenetic tree.

    Args:
        Phylo_t (object): The phylogenetic tree object.

    Returns:
        int: The average distance from the root to all tips.
    """
    total = sum([Phylo_t.get_distance(leaf) for leaf in Phylo_t])
    count = len(Phylo_t)
    return (total / count) if count else 0

def get_node_avg_length(Phylo_t:object)->int:
    """
    Calculate the average branch length of all nodes in the phylogenetic tree.

    Args:
        Phylo_t (object): The phylogenetic tree object.

    Returns:
        int: The average branch length of all nodes in the tree.
    """
    node_length=[]
    node_num=0
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            node_num+=1
            node2root=Phylo_t.get_distance(node)
            node_length.append(node2root)
    return (sum(node_length)/node_num) if node_num else 0

def _resolve_original_gene_name(leaf_name: str, new_named_gene2gene_dic: dict) -> str:
    parts = leaf_name.split('_')
    candidate = '_'.join(parts[1:]) if len(parts) > 1 else leaf_name
    return new_named_gene2gene_dic.get(candidate, new_named_gene2gene_dic.get(leaf_name, leaf_name))


def remove_long_gene(Phylo_t:object, long_branch_index:int, outfile:str, tre_ID:str, new_named_gene2gene_dic:dict) -> object:
    """
    Remove genes with long branches from the phylogenetic tree based on distance thresholds.

    Args:
        Phylo_t (object): The phylogenetic tree object.
        long_branch_index (int): The threshold for identifying long branches.
        outfile (str): Output file handle for writing removal information.
        tre_ID (str): Tree identifier for output records.
        new_named_gene2gene_dic (dict): Mapping from processed node names to original gene names.

    Returns:
        object: The pruned phylogenetic tree object with long branches removed.
    """
    Phylo_t1 = Phylo_t.copy()
    remove_gene_set = set()
    tips_avg_length = get_tips_avg_length(Phylo_t1)
    node_avg_length = get_node_avg_length(Phylo_t1)
    while True:
        is_modified = False 
        
        for leaf in Phylo_t1:
            sister = leaf.get_sisters()[0] if not leaf.is_root() and leaf.get_sisters() else None
            distance = Phylo_t1.get_distance(leaf)
            distance2root_radio = abs(distance / tips_avg_length) if tips_avg_length else 0
            leaf2sister_radio = abs(get_sister_relative_branch_ratio(leaf, sister))
            
            if distance2root_radio > long_branch_index or leaf2sister_radio > long_branch_index:
                gene_name = _resolve_original_gene_name(leaf.name, new_named_gene2gene_dic)
                outfile.write(tre_ID + '\t' + '*' + '\t' + gene_name + '\t' + str(distance2root_radio) + '\t' + str(leaf2sister_radio) + '\t' + '\n')
                remove_gene_set.add(leaf.name)
                is_modified = True  
            else:
                gene_name = _resolve_original_gene_name(leaf.name, new_named_gene2gene_dic)
                outfile.write(tre_ID + '\t' + '\t' + '\t' + gene_name + '\t' + str(distance2root_radio) + '\t' + str(leaf2sister_radio) + '\t' + '\n')


        if not is_modified:
            break
            
        total_leafs_set = set(Phylo_t1.get_leaf_names())
        diff = total_leafs_set - remove_gene_set
        Phylo_t1.prune(diff, preserve_branch_length=True)
    
    return Phylo_t1


def remove_insert_gene(Phylo_t:object,long_branch_index:int,outfile:str,tre_ID:str,new_named_gene2gene_dic:dict)->object:
    """
    Remove inserted genes from the phylogenetic tree based on insertion criteria.

    Args:
        Phylo_t (object): The phylogenetic tree object.
        long_branch_index (int): The threshold for identifying insertions.
        outfile (str): Output file handle for writing removal information.
        tre_ID (str): Tree identifier for output records.
        new_named_gene2gene_dic (dict): Mapping from processed node names to original gene names.

    Returns:
        object: The pruned phylogenetic tree object with insertions removed.
    """
    Phylo_t1 = Phylo_t.copy()    
    taxa_clade = set()
    get_single_clades(Phylo_t1, taxa_clade)
    remove_gene_set = set()

    for clade in taxa_clade:
        temp_name= '_'.join(clade.name.split('_')[1:])
        if is_ancestor_sister_same(clade):
            target_clade = get_target_clade(clade)
            if target_clade==None:
                continue
            else:
                index_num = calculate_insertion_depth(target_clade, clade)
                index_over = calculate_insertion_coverage(target_clade, clade)
                insert=calculate_insertion(target_clade, clade)
                if clade.is_leaf():
                    gene_name = _resolve_original_gene_name(clade.name, new_named_gene2gene_dic)
                    outfile.write(tre_ID+'\t'+'@'+'\t'+gene_name+'\t'+str(index_num)+'\t'+str(index_over)+'\t'+str(insert)+'\n')
                else:
                    outfile.write(tre_ID+'\t'+'@'+'\t'+clade.name+'\t'+str(index_num)+'\t'+str(index_over)+'\t'+str(insert)+'\n')
                remove_gene_set.add(clade.name)              
    total_leafs_set=set(Phylo_t1.get_leaf_names())
    diff=total_leafs_set - remove_gene_set
    Phylo_t1.prune(diff,preserve_branch_length=True)
    return Phylo_t1

def generate_pdf_before(tre_ID,Phylo_t)->None:
    """
    Generate a PDF visualization of the phylogenetic tree before pruning or modification.

    Args:
        tre_ID (str): Tree identifier used in the PDF file name.
        Phylo_t (object): The phylogenetic tree object to visualize.

    Returns:
        None
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(tre_ID + '_before', fsize=10), column=0)
    Phylo_t.render(file_name=tre_ID + '_before.pdf', tree_style=ts)

def generate_pdf_after(tre_ID,Phylo_t)->None:
    """
    Generate a PDF visualization of the phylogenetic tree after pruning or modification.

    Args:
        tre_ID (str): Tree identifier used in the PDF file name.
        Phylo_t (object): The phylogenetic tree object to visualize.

    Returns:
        None
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(tre_ID + '_after', fsize=10), column=0)
    Phylo_t.render(file_name=tre_ID + '_after.pdf', tree_style=ts)

def prune_main_Mono(tre_dic:dict, taxa_dic:dict, long_branch_index:int, insert_branch_index:int, new_named_gene2gene_dic:dict, gene2new_named_gene_dic:dict, visual:bool=False) -> None:
    """
    Main function for pruning phylogenetic trees in mono-copy ortholog analysis, including long branch and insertion branch removal, and optional visualization.

    Args:
        tre_dic (dict): Dictionary mapping tree IDs to tree file paths.
        taxa_dic (dict): Mapping from gene names to species names.
        long_branch_index (int): Threshold for identifying long branches.
        insert_branch_index (int): Threshold for identifying insertion branches.
        new_named_gene2gene_dic (dict): Mapping from processed node names to original gene names.
        gene2new_named_gene_dic (dict): Mapping from original gene names to processed node names.
        visual (bool, optional): Whether to generate and merge PDF visualizations. Defaults to False.

    Returns:
        None
    """
    color_dic = get_color_dict(taxa_dic)
    dir_path1 = os.path.join(os.getcwd(), "orthofilter_mono/pruned_tree/")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    if visual:
        dir_path2 = os.path.join(os.getcwd(), "orthofilter_mono/pruned_tree_pdf/")
        if os.path.exists(dir_path2):
            shutil.rmtree(dir_path2)
        os.makedirs(dir_path2)
    dir_path3 = os.path.join(os.getcwd(), "orthofilter_mono/long_branch_gene/")
    if os.path.exists(dir_path3):
        shutil.rmtree(dir_path3)
    os.makedirs(dir_path3)

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t0 = Tree(tre_path)
        t0.ladderize()
        t0.resolve_polytomy(recursive=True)
        t0.sort_descendants("support")
        t=rename_input_tre(t0,gene2new_named_gene_dic)
        num_tre_node(t)
        
        if is_multi_copy_tree(t):
            o = open(os.path.join(dir_path3, tre_ID + '_delete_gene.txt'), 'w')
            o.write('tre_ID' + '\t' + 'long_branch_label' + '\t' + 'gene' + '\t' +'root_relative_branch_ratio' + '\t' + 'sister_relative_branch_ratio' + '\n')
            rename_input_single_tre(t, taxa_dic,new_named_gene2gene_dic)
            
            if  visual:
                set_style(t, color_dic,new_named_gene2gene_dic)
                generate_pdf_before(tre_ID,t)
            t1=t
            # t1 = remove_long_gene(t, long_branch_index, o, tre_ID,new_named_gene2gene_dic)
            o.write('\n')
            if len(get_species_set(t1)) !=1:
                o.write('tre_ID' + '\t' + 'insert_branch_label' + '\t' + 'gene' + '\t' + 'insertion_depth' + '\t' + 'insertion_coverage' + '\t' + 'calculate_insertion' + '\n')
                t2 = remove_insert_gene(t1, insert_branch_index, o, tre_ID,new_named_gene2gene_dic)
                o.close()
                if  visual:
                    generate_pdf_after(tre_ID,t2)
                    merge_pdfs_side_by_side(tre_ID + '_before.pdf', tre_ID + '_after.pdf', os.path.join(dir_path2, tre_ID + '.pdf'))
                    os.remove(tre_ID + '_before.pdf')
                    os.remove(tre_ID + '_after.pdf')
                t3 = rename_output_tre(t2,new_named_gene2gene_dic)
                tree_str = trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)
            else:
                if  visual:
                    generate_pdf_after(tre_ID,t1)
                    merge_pdfs_side_by_side(tre_ID + '_before.pdf', tre_ID + '_after.pdf', os.path.join(dir_path2, tre_ID + '.pdf'))
                    os.remove(tre_ID + '_before.pdf')
                    os.remove(tre_ID + '_after.pdf')
                t3 = rename_output_tre(t1,new_named_gene2gene_dic)
                tree_str = trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)

            
        else:
            o = open(os.path.join(dir_path3, tre_ID + '_delete_gene.txt'), 'w')
            o.write('tre_ID' + '\t' + 'long_branch_label' + '\t' + 'gene' + '\t' +'root_relative_branch_ratio' + '\t' + 'sister_relative_branch_ratio' + '\n')
            rename_input_single_tre(t, taxa_dic,new_named_gene2gene_dic)
            if  visual:
                set_style(t, color_dic,new_named_gene2gene_dic)
                generate_pdf_before(tre_ID,t)
            t1=t
            # t1 = remove_long_gene(t, long_branch_index, o, tre_ID,new_named_gene2gene_dic)
            o.write('\n')
            if len(get_species_set(t1)) !=1:
                o.write('tre_ID' + '\t' + 'insert_branch_label' + '\t' + 'gene' + '\t' + 'insertion_depth' + '\t' + 'insertion_coverage' + '\t' + 'calculate_insertion' + '\n')
                t2=remove_insert_gene(t1,insert_branch_index,o,tre_ID,new_named_gene2gene_dic)
            
                while not is_single_tree(t2):
                    prune_single(t2)
                    pass
                if  visual:
                    generate_pdf_after(tre_ID,t2)
                    merge_pdfs_side_by_side(tre_ID+'_before.pdf', tre_ID+'_after.pdf', os.path.join(dir_path2, tre_ID + '.pdf'))
                    os.remove(tre_ID+'_before.pdf')
                    os.remove(tre_ID+'_after.pdf')

                t3=rename_output_tre(t2,new_named_gene2gene_dic)
                tree_str= trans_branch_length(t3)
                write_tree_to_newick(tree_str,tre_ID,dir_path1)
            else:
                if  visual:
                    generate_pdf_after(tre_ID,t1)
                    merge_pdfs_side_by_side(tre_ID + '_before.pdf', tre_ID + '_after.pdf', os.path.join(dir_path2, tre_ID + '.pdf'))
                    os.remove(tre_ID + '_before.pdf')
                    os.remove(tre_ID + '_after.pdf')
                t3 = rename_output_tre(t1,new_named_gene2gene_dic)
                tree_str =trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)
            
            o.close()
        pbar.update(1)
    pbar.close()

def prune_mono_copy_trees(tre_dic:dict, taxa_dic:dict, long_branch_index:int, insert_branch_index:int, new_named_gene2gene_dic:dict, gene2new_named_gene_dic:dict, visual:bool=False) -> None:
    return prune_main_Mono(tre_dic, taxa_dic, long_branch_index, insert_branch_index, new_named_gene2gene_dic, gene2new_named_gene_dic, visual)
    
if __name__ == "__main__":
    os.makedirs(os.path.join(os.getcwd(), "pruned_tree"))
    taxa_dic=read_and_return_dict('taxa.txt')
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    long_brancch_index=5
    insert_branch_index=5
    prune_main_Mono(tre_dic,taxa_dic,long_brancch_index,insert_branch_index)
