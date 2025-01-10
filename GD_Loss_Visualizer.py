from __init__ import *
import matplotlib.pyplot as plt
import re 
from PIL import Image
import math

def extract_numbers_in_parentheses(string):
    pattern = r"\((\d+)\)"
    matches = re.findall(pattern, string)
    return matches

def get_gd_node_to_species_loss_dic(sp_dic):
    loss_dic = {'No-duplicate loss': [], 'One-duplicate loss': [], 'Two-duplicate loss': []}

    for k, v in sp_dic.items():
        s = extract_numbers_in_parentheses(k)
        s_set = set(s)
        
        if s_set == {'2'}:
            loss_dic['No-duplicate loss'].append(v)
        elif s_set == {'2', '1'} or s_set == {'2', '1', '0'} or s_set == {'1', '0'}:
            loss_dic['One-duplicate loss'].append(v)
        else:
            loss_dic['Two-duplicate loss'].append(v)
    
    for key in loss_dic:
        filtered_values = []
        for val in loss_dic[key]:
            if isinstance(val, str) and not val.isdigit():
                continue
            filtered_values.append(int(val) if isinstance(val, str) and val.isdigit() else float(val))
        loss_dic[key] = filtered_values

    sum_loss_dic = {key: sum(value) for key, value in loss_dic.items()}
    
    return sum_loss_dic

def visualizer_sum_loss_dic(sum_loss_dic, sps, gd_id):
    keys = list(sum_loss_dic.keys())
    values = list(sum_loss_dic.values())
    color = 'lightblue'

    plt.figure(figsize=(10, 6))  

    bars = plt.bar(keys, values, color=color)

    plt.ylabel('GD num', fontsize=14) 
    plt.title('Count of ' + sps + ' Species under Node ' + gd_id, fontsize=16)  
    plt.yticks(fontsize=12)  

    for key, value in zip(keys, values):
        plt.text(key, value, str(value), ha='center', va='bottom', fontsize=12)

    plt.tight_layout()  
    plt.savefig(f'{gd_id}_{sps}.png')
    plt.cla()
    plt.close("all")  
    
def generate_plt(full_path):
    new_dic = read_and_return_dict(full_path)
 
    sum_loss_dic = get_gd_node_to_species_loss_dic(new_dic)

    first_key = list(new_dic.keys())[1]
    
    gd_id = first_key.split('->')[0].split('(')[0]
    sps = first_key.split('->')[-1].split('(')[0]
    
    visualizer_sum_loss_dic(sum_loss_dic, sps, gd_id)

def process_gd_loss_summary(file):
    element_counts = {}
    
    with open(file, 'r') as f:
        lines = f.readlines()

        for line in lines[2:]:  
            line = line.strip()
            if not line:  
                continue
            path = line.split('\t')[0]
            input_string = line.split('\t')[1]
            treeid= [og.strip().strip("'") for og in input_string.split(',')]
            elements = path.split('->')
            for i in elements:
                if i in element_counts:
                    element_counts[i].update(treeid)  
                else:
                    element_counts[i] = set(treeid)  
    new_dic={}
    for key in element_counts:
        new_dic[key] = len(element_counts[key])

    return new_dic

def transform_dict(original_dict):
    new_dict = {}

    for key, value in original_dict.items():
        gene_name, number = key.split('(') 
        number = number.strip(')') 

        if gene_name not in new_dict:
            new_dict[gene_name] = {}
        
        if number not in new_dict[gene_name]:
            new_dict[gene_name][number] = 0
        
        new_dict[gene_name][number] += value  

    return new_dict

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

def visualizer_sptree(result,sptree):
    new_dict=transform_dict(result)
    sptree.ladderize()
    sptree.sort_descendants("support")

    ts = TreeStyle()
    ts.title.add_face(TextFace("Green color : No-duplicate loss",fsize=5), column=1)
    ts.title.add_face(TextFace("Blue color : One-duplicate loss", fsize=5),column=1)
    ts.title.add_face(TextFace("Red color : Two-duplicate loss", fsize=5),column=1)
    ts.title.add_face(TextFace("The number represents the number of statistical GDs", fsize=5),column=1)

    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = True
    ts.show_branch_support = False
    ts.extra_branch_line_type =0
    ts.extra_branch_line_color='black'
    ts.branch_vertical_margin = -1

    for node in sptree.traverse():

        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        node.set_style(nstyle)
        if node.name in new_dict:
            loss_dic=new_dict[node.name]
            for k,v in loss_dic.items():
                if k=='2':
                    color='green'
                    node.add_face(TextFace(v, fsize=5, fgcolor=color), column=0, position="branch-top")
                elif k=='1':
                    color='blue'
                    node.add_face(TextFace(v, fsize=5, fgcolor=color), column=1, position="branch-top")
                else:
                    color='red'
                    node.add_face(TextFace(v, fsize=5, fgcolor=color), column=2, position="branch-top")
    # for leaf in sptree:
    #     if leaf.name in taxa:
    #         leaf.name=taxa[leaf.name]
            # node.add_face(TextFace(taxa[leaf.name], fsize=5, fgcolor='black',ftype='Arial'), column=0, position="aligned")
    # realign_branch_length(sptree)
    # rejust_root_dist(sptree)
    sptree.render('gd_loss_visualizer.PDF',w=210, units="mm",tree_style=ts)
if __name__ == "__main__":
    out='outfile'
    input_folder='input'
    os.makedirs(out, exist_ok=True)
    generate_plt(input_folder,out)

