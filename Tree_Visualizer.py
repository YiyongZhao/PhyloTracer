import matplotlib.pyplot as plt
import matplotlib.colors as colors
import re
import string
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

def rejust_root_dist(sptree):
    clade_up=sptree.get_children()[0]
    clade_down=sptree.get_children()[1]
    if len(clade_up)>len(clade_down):
        clade_up.dist=1 
        if clade_down.is_leaf():
            clade_down.dist=get_max_deepth(sptree)-1
        else:
            clade_down.dist=get_max_deepth(sptree)-get_max_deepth(clade_down)
    else:
        clade_down.dist=1 
        if clade_up.is_leaf():
            clade_up.dist=get_max_deepth(sptree)-1
        else:
            clade_up.dist=get_max_deepth(sptree)-get_max_deepth(clade_up)

    return sptree
#############################################################################################################
def Dup_NodeIDs_from_Numbered_GFs(Phylo_t:object)->object:
    if not is_rooted(Phylo_t):
        Phylo_t = root_tre_with_midpoint_outgroup(Phylo_t)
    Phylo_t = num_tre_node(Phylo_t)
    dup_node_name_list = find_dup_node(Phylo_t)
    Phylo_t1 = Tree(Phylo_t.write())
    num_tre_node(Phylo_t1)
    
    return Phylo_t1, dup_node_name_list

def create_tree_style(tree_style,tre_ID):
    ts = TreeStyle()
    ts.title.add_face(TextFace(' ', fsize=20), column=1)
    ts.title.add_face(TextFace(tre_ID, fsize=20), column=0)
    ts.title.add_face(TextFace("★", fsize=20, fgcolor="red"), column=0)
    ts.title.add_face(TextFace("Interspecific gene duplication event", fsize=20), column=1)
    ts.title.add_face(TextFace("★", fsize=20, fgcolor="blue"), column=0)
    ts.title.add_face(TextFace("Intraspecific gene duplication event", fsize=20), column=1)
    
    ts.mode = tree_style
    ts.scale = 20
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
    colormap = plt.get_cmap("rainbow")
    # color dictionary
    unique_values=set(dictory.values())   
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_list = {k: v + '@' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}
    
    return sps_color_list

def generate_color_dict(gene_categories:list)->list:
    return [get_color_dict(i) for i in gene_categories]

def generate_string(index):
    letters = list(string.ascii_uppercase) + list(string.ascii_lowercase)
    if index < len(letters):
        return "@" + letters[index]
    else:
        first_letter = letters[(index - len(letters)) // 52]
        second_letter = letters[(index - len(letters)) % 52]
        return "@" + first_letter + second_letter

def get_new_sorted_dict(gene2fam):
    uniq_fam=set(get_color_dict(gene2fam).values())
    fam2color={i.split('@')[0]:i.split('@')[-1] for i in uniq_fam}
    sorted_dict = dict(sorted(fam2color.items(), key=lambda x: x[0], reverse=False))
    for index, key in enumerate(sorted_dict.keys()):
        suffix = generate_string(index)
        sorted_dict[key] += suffix
    return sorted_dict

def fuzzy_match(search_string, key):
    return re.search(search_string, key)

def tips_mark(Phylo_t1:object,voucher2taxa_dic:dict,gene_categories:list,tre_ID,ts,new_named_gene2gene_dic:dict,dir_path,gene2fam=None,df=None)->object:
    sps_color_dict=get_color_dict(voucher2taxa_dic)
    if gene2fam is not None:
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
                color = color_dict[species].split('@')[-1]
                face = TextFace("   ▐" + '  ' + color_dict[species].split('@')[0], fgcolor=color,ftype='Arial')
                add_face_to_node(node, face, column, position="aligned")

    def add_species_face(node, species):
        if species in sps_color_dict:
            color = sps_color_dict[species].split('@')[-1]
            face = TextFace(' ' + only_gene, fgcolor=color,ftype='Arial', fstyle='italic')
            node.add_face(face, column=-1)

        if species in sps_color_dict:
            color = sps_color_dict[species].split('@')[-1]
            face4 = TextFace("   ▐" + '  ' + sps_color_dict[species].split('@')[0], fgcolor=color,ftype='Arial', fstyle='italic')
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
            color = matched_value.split('@')[-1]
            face5 = TextFace("   ▐" + '  ' +  matched_value.split('@')[0], fgcolor=color,ftype='Arial', fstyle='italic')
            add_face_to_node(node, face5, column, position="aligned")

    
    for node in Phylo_t1.traverse():
        if node.is_leaf():
            gene = new_named_gene2gene_dic[node.name]
            species = voucher2taxa_dic[node.name.split("_")[0]]
            rename_species=node.name.split("_")[0]
            only_gene='_'.join(gene.split('_')[1:])
            add_species_face(node, rename_species)
            column = 1
            for color_dict in color_dicts:
                generate_face_mark(node, species, column, color_dict)
                column += 1
            if gene2fam is not None and  gene[4:] in gene_color_dict:
                add_gene_face(node, gene[4:])
            else:
                pass

    def get_color(value):
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
        
    def add_heat_map_to_node(tree,df,new_named_gene2gene_dic,gene_categories):
        for node in tree:
            gene = new_named_gene2gene_dic[node.name]
            only_gene='_'.join(gene.split('_')[1:])
            if node.is_leaf() and only_gene in df.index:
                n=len(gene_categories)+1
                for i in df.columns.tolist():
                    color = get_color(df[i][only_gene])
                    face=RectFace(width=10, height=10, fgcolor=color,bgcolor=color)
                    node.add_face(face, column=n, position='aligned')
                    n+=1
                    
    def add_header_to_tree(ts,df,gene_categories):
        labels = df.columns.to_list()
        n=len(gene_categories)+1
        for i in labels:
            face=TextFace(' '+i,fgcolor='black',fsize=9)
            face.rotation=-90
            face.vt_align =2
            ts.aligned_header.add_face(face,n)
            n+=1
            
    def add_color_bar(ts):
        ts.legend.add_face(TextFace(' '), column=0)
        bar_face=TextFace('Color Bar Title ')
        ts.legend.add_face(bar_face, column=0)
        cols=['#006599', '#408ca6', '#7fb2b2', '#bfd9bf', '#ffffcc','#f7deab', '#eebc88', '#e69966', '#dc7845', '#d55623', '#cc3300']
        bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        col_dic=dict(zip(bounds,cols))
        n=1
        for k,v in col_dic.items():
            colorbar_face = RectFace(width=20, height=20, fgcolor=v,bgcolor=v)
            ts.legend.add_face(TextFace(' '+str(k)), column=n)
            ts.legend.add_face(colorbar_face, column=n)
            n+=1
            
        ts.legend_position=2

    if df is not None:
        add_heat_map_to_node(Phylo_t1,df,new_named_gene2gene_dic,gene_categories)
        add_header_to_tree(ts,df,gene_categories)
        add_color_bar(ts)
    return Phylo_t1.render(dir_path+tre_ID+'.pdf',tree_style=ts)


def get_matched_value(gene, gene2fam):
    for key, value in gene2fam.items():
        if fuzzy_match(key, gene):
            return key, value
    return None, None

    
def get_fam_dic(t, gene2fam):
    fam_dic = {}
    for i in t:
        gene = i.name[4:]
        match_gene, match_family = get_matched_value(gene, gene2fam)
        if match_family and match_family in gene2fam.values():
            fam_dic.setdefault(match_family, []).append(i.name)
    return fam_dic

def find_combinations(my_list):
    combinations = []
    for i in range(len(my_list)):
        for j in range(i+1, len(my_list)):
            combinations.append((my_list[i], my_list[j]))
    return combinations

def get_dup_family_dic(t,gene2fam):
    fam_node_dic={}
    fam_dic=get_fam_dic(t,gene2fam)
    for k,v in fam_dic.items():
        nodes=set()
        com=find_combinations(v)
        for i in com:
            clade=t.get_common_ancestor(i)
            nodes.add(clade.name)
        fam_node_dic[k]=nodes
    return fam_node_dic

def mapping_sptree(t,sptree,sp_node_dic,gene2fam):
    
    fam_node_dic=get_dup_family_dic(t,gene2fam)
    for k,v in fam_node_dic.items():
        for i in v :
            clade=t&i
            sps=get_species_list(clade)
            uniq_sps=set(sps)
            clade2sptree=sptree.get_common_ancestor(uniq_sps)
            sp_node_dic.setdefault(clade2sptree.name, set()).add(k)
            
def get_sptree_style(sorted_dict):
    ts = TreeStyle()
    ts.legend_position = 1 
    ts.mode ='r'
    ts.scale = 30
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = True
    #ts.show_branch_support = True
    ts.extra_branch_line_type =0
    ts.extra_branch_line_color='black'
    ts.branch_vertical_margin=-1
    for k,v in sorted_dict.items():
        ts.legend.add_face(TextFace(v.split('@')[1]+' '+k,fsize=20,fgcolor=v.split('@')[0]), column=0)
    return ts

def mark_gene_to_sptree(sptree, sp_node_dic,gene_categories,sorted_dict):
    color_dicts = generate_color_dict(gene_categories)
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


def rename_sptree(sptree):
    for i in sptree.traverse() :
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

def view_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic,gene_categories,tree_style,keep_branch,new_named_gene2gene_dic,gene2fam=None,df=None):
    dir_path = os.path.join(os.getcwd(), "output/pdf_result/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID,tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        Phylo_t0=read_phylo_tree(tre_path)
        Phylo_t0=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        Phylo_t1,ts=get_treestyle(Phylo_t0,tree_style,tre_ID)
        Phylo_t1.ladderize()
        Phylo_t1.resolve_polytomy(recursive=True)
        Phylo_t1.sort_descendants("support")
        if keep_branch !='1' :
            realign_branch_length(Phylo_t1)
            rejust_root_dist(Phylo_t1)
        tips_mark(Phylo_t1,voucher2taxa_dic,gene_categories,tre_ID,ts,new_named_gene2gene_dic,dir_path,gene2fam,df)
        pbar.update(1)
    pbar.close()

def mark_gene_to_sptree_main(tre_dic,gene_categories,sptree,gene2fam):
    sorted_dict=get_new_sorted_dict(gene2fam)
    num_tre_node(sptree)
    sp_node_dic={}
    for k,v in tre_dic.items():
        t=read_phylo_tree(v)
        num_tre_node(t)
        mapping_sptree(t,sptree,sp_node_dic,gene2fam)
    mark_gene_to_sptree(sptree,sp_node_dic,gene_categories,sorted_dict)
    rename_sptree(sptree)
    ts=get_sptree_style(sorted_dict)
    realign_branch_length(sptree)
    clade_up=sptree.get_children()[0]
    clade_down=sptree.get_children()[1]
    clade_up.dist=1 
    clade_down.dist=get_max_deepth(clade_up)-1
    sptree.render(file_name='424species_branch7.pdf',tree_style=ts)

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
    
	
