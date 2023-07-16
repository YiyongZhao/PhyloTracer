from ete3 import NodeStyle, Tree, TreeStyle,TextFace
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from Base import *

def find_dup_node(Phylo_t:object)->list:#After searching for duplication events, the list of node names where duplication events occurred is as follows:
    events = Phylo_t.get_descendant_evol_events()
    dup_node_name_list = []
    for ev in events:
        if ev.etype == "D":
            i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
            events_node_name_list = i.split(',')
            common_ancestor_node_name = Phylo_t.get_common_ancestor(events_node_name_list)
            dup_node_name_list.append(common_ancestor_node_name.name)
    return dup_node_name_list

def realign_branch_length(Phylo_t1:object)->object:
    max_deep=get_max_deepth(Phylo_t1)
    
    for node in Phylo_t1.traverse():
        if not node.is_leaf():
            node.dist=1
        else:
            ancestors_list=node.get_ancestors()#list→All ancestor nodes of the current node.
            current_deepth=len(ancestors_list)
            node.dist=max_deep-current_deepth
    
    return Phylo_t1

#############################################################################################################
def Dup_NodeIDs_from_Numbered_GFs(Phylo_t:object)->object:
    if not isRoot_Judger(Phylo_t):
        Phylo_t = root_tre_with_midpoint_outgroup(Phylo_t)
    Phylo_t = num_tre_node(Phylo_t)
    dup_node_name_list = find_dup_node(Phylo_t)
    Phylo_t1 = Tree(Phylo_t.write())
    num_tre_node(Phylo_t1)
    Phylo_t1.ladderize()
    Phylo_t1.sort_descendants("support")
    
    return Phylo_t1, dup_node_name_list

def create_tree_style(tree_style,tre_ID):
    ts = TreeStyle()
    ts.legend.add_face(TextFace("★", fsize=20, fgcolor="blue"), column=0)
    ts.legend.add_face(TextFace("Intraspecific duplication", fsize=20), column=1)
    ts.legend.add_face(TextFace("★", fsize=20, fgcolor="red"), column=0)
    ts.legend.add_face(TextFace("Interspecific duplication", fsize=20), column=1)
    ts.title.add_face(TextFace(tre_ID, fsize=40), column=2)
    ts.legend_position = 1
    
    ts.mode = tree_style
    ts.orientation = 0
    ts.branch_vertical_margin = 0
    ts.scale = 50
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = -1
    
    return ts

def set_node_style(node:object, dup_node_name_list:list):
    nstyle = NodeStyle()
    splist = set(get_species_list(node))
    if node.name in dup_node_name_list and len(splist) == 1:
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "blue"
        node.add_face(TextFace("★", fsize=10, fgcolor="blue"), column=1, position="branch-top")
    elif node.name in dup_node_name_list and len(splist) != 1:
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "red"
        node.add_face(TextFace("★", fsize=10, fgcolor="red"), column=1, position="branch-top")
    else:
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
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
    sps_color_list = {k: v + '-' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}
    
    return sps_color_list

def generate_color_dict(gene_categories:list)->list:
    return [get_color_dict(i) for i in gene_categories]

def tips_mark(Phylo_t1:object,voucher2taxa_dic:dict,gene_categories:list,tre_ID,ts)->object:
    sps_color_dict=get_color_dict(voucher2taxa_dic)
    color_dicts = generate_color_dict(gene_categories)
    faces_added = set()  # used to track the added faces

    def add_face_to_node(node, color_dict:dict,species: str, column: int, position="aligned"):
        if (node, column, position) not in faces_added and species in color_dict:
            color = color_dict[species].split('-')[1]
            face = TextFace("   ▐" + '  ' + color_dict[species].split('-')[0], fgcolor=color)
            node.add_face(face, column=column, position=position)
            faces_added.add((node, column, position))

    def generate_face_mark(node, column:int, species:str, color_dict:dict):
        if (node, column, "aligned") not in faces_added:
            if species in color_dict:
                add_face_to_node(node, color_dict,species, column, position="aligned")
                

    for node in Phylo_t1.traverse():
        if node.is_leaf():
            # Get species names
            re_named_species = node.name.split("_")[0]
            gene = node.name[4:]
            species = voucher2taxa_dic[re_named_species]

            # Set the color for this species name
            if re_named_species in sps_color_dict:
                color = sps_color_dict[re_named_species].split('-')[1]
                face = TextFace(' ' + gene, fgcolor=color)
                node.add_face(face, column=-1)

            if re_named_species in sps_color_dict:
                color = sps_color_dict[re_named_species].split('-')[1]
                face4 = TextFace("   ▐" + '  ' + species, fgcolor=color)
                node.add_face(face4, column=0, position="aligned")

            column = 1
            for color_dict in color_dicts:
                generate_face_mark(node, column, species,color_dict)
                column += 1
    return Phylo_t1.render(str(tre_ID)+'.PDF',tree_style=ts)

def view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic,gene_categories,tree_style,keep_branch):
    for tre_ID,tre_path in tre_dic.items():
        Phylo_t0=read_tree(tre_path)
        rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        Phylo_t1,ts=get_treestyle(Phylo_t0,tree_style,tre_ID)
        if keep_branch !='yes' :
            realign_branch_length(Phylo_t1)
        tips_mark(Phylo_t1,voucher2taxa_dic,gene_categories,tre_ID,ts)
####################################################################################
if __name__ == "__main__":
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    gene2fam=read_and_return_dict('gene2fam')
    sp2order=read_and_return_dict('sp2oeder')
    sp2family=read_and_return_dict('sp2family')
    tre_dic=read_and_return_dict('GF_list')
    gene_categories=[]
    gene_categories.append(sp2family)
    gene_categories.append(sp2order)
    tree_style='r'
    view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic,gene_categories,tree_style,keep_branch)
	


