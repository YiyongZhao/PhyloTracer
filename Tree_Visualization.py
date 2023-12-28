from ete3 import NodeStyle, Tree, TreeStyle,TextFace
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from __init__ import *

def realign_branch_length(Phylo_t1:object)->object:
    Phylo_t1.ladderize()
    Phylo_t1.resolve_polytomy(recursive=True)
    Phylo_t1.sort_descendants("support")
    max_deep=get_max_deepth(Phylo_t1)
    for node in Phylo_t1.traverse():
        if not node.is_root():
            node.dist=1
            degree=node.get_distance(node.get_tree_root()) + 1
            deep=get_max_deepth(node)
            node.dist=max_deep-deep-degree
    clade_up=Phylo_t1.get_children()[0]
    clade_down=Phylo_t1.get_children()[1]
    difference=abs(get_max_deepth(clade_up)-get_max_deepth(clade_down))+1
    clade_up.dist=clade_up.dist+difference  
    clade_down.dist=clade_down.dist+difference   
    
    return Phylo_t1

#############################################################################################################
def Dup_NodeIDs_from_Numbered_GFs(Phylo_t:object)->object:
    if not isRoot_Judger(Phylo_t):
        Phylo_t = root_tre_with_midpoint_outgroup(Phylo_t)
    Phylo_t = num_tre_node(Phylo_t)
    dup_node_name_list = find_dup_node(Phylo_t)
    Phylo_t1 = Tree(Phylo_t.write())
    num_tre_node(Phylo_t1)
    
    return Phylo_t1, dup_node_name_list

def create_tree_style(tree_style,tre_ID):
    ts = TreeStyle()
    ts.legend.add_face(TextFace("★", fsize=20, fgcolor="red"), column=0)
    ts.legend.add_face(TextFace("Interspecific gene duplication event", fsize=20), column=1)
    ts.legend.add_face(TextFace("★", fsize=20, fgcolor="blue"), column=0)
    ts.legend.add_face(TextFace("Intraspecific gene duplication event", fsize=20), column=1)
    ts.title.add_face(TextFace(tre_ID, fsize=30), column=2)
    
    ts.legend_position = 1
    
    ts.mode = tree_style
    ts.scale = 50
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.extra_branch_line_type =0
    ts.extra_branch_line_color='black'
    ts.branch_vertical_margin = -1
	
    
    return ts

def set_node_style(node:object, dup_node_name_list:list):
    nstyle = NodeStyle()
    splist = set(get_species_list(node))
    if node.name in dup_node_name_list and len(splist) == 1:
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0 
        nstyle["hz_line_type"] = 0 
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        node.add_face(TextFace("★", fsize=7, fgcolor="blue"), column=1, position="branch-top")
    elif node.name in dup_node_name_list and len(splist) != 1:
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0 
        nstyle["hz_line_type"] = 0 
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        node.add_face(TextFace("★", fsize=7, fgcolor="red"), column=1, position="branch-top")
    else:
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0 
        nstyle["hz_line_type"] = 0 
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
	node.set_style(nstyle)
    
def get_treestyle(Phylo_t:object,tree_style:str,tre_ID:str)->object:
    Phylo_t1, dup_node_name_list = Dup_NodeIDs_from_Numbered_GFs(Phylo_t)
    ts = create_tree_style(tree_style,tre_ID)

    for node in Phylo_t1.traverse():
        set_node_style(node, dup_node_name_list)
    
    return Phylo_t1, ts


####################################################################################################################
def get_color_dict(dictory:dict)->dict:
    colormap = plt.get_cmap("gist_rainbow")
    # color dictionary
    unique_values=set(dictory.values())   
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_list = {k: v + '*' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}
    
    return sps_color_list

def generate_color_dict(gene_categories:list)->list:
    return [get_color_dict(i) for i in gene_categories]

def fuzzy_match(search_string, key):
    return re.search(search_string, key)

def tips_mark(Phylo_t1:object,voucher2taxa_dic:dict,gene_categories:list,tre_ID,ts,new_named_gene2gene_dic:dict)->object:
    sps_color_dict=get_color_dict(voucher2taxa_dic)
    if 'gene2fam'  in globals():
        gene_color_dict=get_color_dict(gene2fam)
    color_dicts = generate_color_dict(gene_categories)
    faces_added = set()  # used to track the added faces

    def add_face_to_node(node, face, column, position="aligned"):
        if (node, column, position) not in faces_added:
            node.add_face(face, column=column, position=position)
            faces_added.add((node, column, position))

    def generate_face_mark(node, species, column, color_dict):
        if (node, column, "aligned") not in faces_added:
            if species in color_dict:
                color = color_dict[species].split('*')[-1]
                face = TextFace("   ▐" + '  ' + color_dict[species].split('*')[0], fgcolor=color,ftype='Arial', fstyle='italic')
                add_face_to_node(node, face, column, position="aligned")

    def add_species_face(node, species):
        if species in sps_color_dict:
            color = sps_color_dict[species].split('*')[-1]
            face = TextFace(' ' + gene, fgcolor=color,ftype='Arial',fstyle='italic')
            node.add_face(face, column=-1)

        if species in sps_color_dict:
            color = sps_color_dict[species].split('*')[-1]
            face4 = TextFace("   ▐" + '  ' + sps_color_dict[species].split('*')[0], fgcolor=color, ftype='Arial',fstyle='italic')
            add_face_to_node(node, face4, 0, position="aligned")

    def add_gene_face(node, gene):
        matched_key = None
        matched_value = None
        for key in gene_color_dict:
            if fuzzy_match(key, gene):
                matched_key = key
                matched_value = gene_color_dict.get(key)
                break
        if matched_value:
            color = matched_value.split('*')[-1]
            face5 = TextFace("   ▐" + '  ' +  matched_value.split('*')[0], fgcolor=color, ftype='Arial',fstyle='italic')
            add_face_to_node(node, face5, column, position="aligned")

    
    for node in Phylo_t1.traverse():
        if node.is_leaf():
            gene = new_named_gene2gene_dic[node.name]
            species = voucher2taxa_dic[node.name.split("_")[0]]
            rename_species=node.name.split("_")[0]
            add_species_face(node, rename_species)
            column = 1
            for color_dict in color_dicts:
                generate_face_mark(node, species, column, color_dict)
                column += 1
            if 'gene2fam' in globals() and  gene in gene_color_dict:
                add_gene_face(node, gene)
            else:
                pass

    return Phylo_t1.render(file_name='pdf_result/'+str(tre_ID)+'.PDF',tree_style=ts)

def view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic,gene_categories,tree_style,keep_branch,new_named_gene2gene_dic):
    for tre_ID,tre_path in tre_dic.items():
        Phylo_t0=read_phylo_tree(tre_path)
        Phylo_t0=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        Phylo_t1,ts=get_treestyle(Phylo_t0,tree_style,tre_ID)
        Phylo_t1.ladderize()
        Phylo_t1.resolve_polytomy(recursive=True)
        Phylo_t1.sort_descendants("support")
        if keep_branch !=1 :
            realign_branch_length(Phylo_t1)
        tips_mark(Phylo_t1,voucher2taxa_dic,gene_categories,tre_ID,ts,new_named_gene2gene_dic)
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
	
