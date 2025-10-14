from hmac import new
from __init__ import *
import matplotlib.pyplot as plt
import re 
from PIL import Image
import math
from ete3 import PieChartFace,CircleFace


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

def parse_file_to_dic(filepath):

    result_dic = {}
    gd_id_set=set()
    with open(filepath) as f:
        for line in f:
            cols = line.strip().split('\t')
    
            gdid = cols[1]
            
               
            if len(cols) < 8:
                continue

            if gdid not in gd_id_set:
                gdid_species_set = set()
                path_strs = cols[7].split('@')
                for path_str in path_strs:
                    for node in path_str.split('->'):
                        if '(' in node and node.endswith(')'):
                            sp, num = node[:-1].split('(')
                            num = int(num)
                            if sp not in result_dic:
                                result_dic[sp] = {2:0, 1:0, 0:0}
                            gdid_species_set.add((sp, int(num)))
                for sp, num in gdid_species_set:
                    result_dic[sp][num] += 1

                gd_id_set.add(gdid)
            else:
                continue

    
    return result_dic


def build_legend(ts):
  
    ts.title.add_face(TextFace(' ',fsize=10), column=0)
    ts.title.add_face(CircleFace(5, color='green'), column=1)
    ts.title.add_face(TextFace(' Two copies retained ',fsize=10), column=2)
    ts.title.add_face(CircleFace(5, color='blue'), column=3)
    ts.title.add_face(TextFace(' One copy retained ',fsize=10), column=4)
    ts.title.add_face(CircleFace(5, color='red'), column=5)
    ts.title.add_face(TextFace(' No copy retained ',fsize=10), column=6)
    ts.legend.add_face(TextFace(" Legend: (Green, Blue, Red)", fsize=9, fgcolor="gray"), column=0)
    ts.legend_position = 3

def visualizer_sptree(result_dic,sptree):
    new_dic={}
    for sp, count_dic in result_dic.items():
        total = sum(count_dic.values())
        if total == 0:
            continue
        a=[]
        for v in count_dic.values():
            b= int(v/ total * 100)
            a.append(b)
        diff = 100 - sum(a)
        if diff != 0:
            
            idx = a.index(max(a))
            a[idx] += diff
        new_dic[sp]=a

    sptree.ladderize()
    sptree.sort_descendants("support")

    for node in sptree.traverse():
        node_id = node.name
        if node_id in new_dic:
            pie_values = new_dic[node_id]
            total = sum(pie_values)

           
            pie_face = PieChartFace(pie_values, width=10, height=10, colors=['green', 'blue', 'red'])
            node.add_face(pie_face, column=0, position="float")

            labels = [str(va) for va in result_dic[node_id].values()]
            
            label_text = "(" + ",".join(labels) + ")"
            label_face = TextFace(label_text, fsize=5)
            node.add_face(label_face, column=-1, position="float")
            node.add_face(TextFace(' ', fsize=5), column=-1, position="float")
            

        nstyle = NodeStyle()
        nstyle["size"] = 0
        node.set_style(nstyle)

    
    ts = TreeStyle()
    ts.scale = 60
    

    build_legend(ts)
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

    realign_branch_length(sptree)
    rejust_root_dist(sptree)
    sptree.render('gd_loss_visualizer.PDF',w=210, units="mm",tree_style=ts)


if __name__ == "__main__":
    out='outfile'
    input_folder='input'
    os.makedirs(out, exist_ok=True)
    generate_plt(input_folder,out)

