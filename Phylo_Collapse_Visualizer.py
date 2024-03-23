import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib 
matplotlib.use('Agg') 
from __init__ import *


def create_tree_style(tre_ID):
    ts = TreeStyle()
    ts.title.add_face(TextFace(tre_ID, fsize=30), column=2)
    
    ts.legend_position = 1
    
    ts.mode = 'r'
    ts.scale = 20
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.extra_branch_line_type =0
    ts.extra_branch_line_color='black'
    ts.branch_vertical_margin = -1
    
    
    return ts


def set_style_node(t,c_color_dict):
    nstyle = NodeStyle()
    for i in t.traverse():
        nstyle["vt_line_width"] =0
        nstyle["hz_line_width"] = 0
        nstyle["vt_line_type"] = 0 
        nstyle["hz_line_type"] = 0 
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        i.set_style(nstyle)
        if  i.is_leaf():
            v=i.name.split('_')[0]
            if v in c_color_dict:
                color1 = c_color_dict[v]
                face = TextFace(''+i.name, fgcolor=color1,fstyle='italic')
                i.add_face(face, column=0)
    return t

def get_taxa_to_color_dict(dictory:dict)->dict:
    colormap = plt.get_cmap("gist_rainbow")
    # color dictionary\n",
    unique_values=set(dictory.values())
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_list = {v: color_dict.get(v) for k, v in dictory.items() if v in color_dict}

    return sps_color_list 

def collapse_visual_main(tre_dic,c_color_dic):
    dir_path1 = os.path.join(os.getcwd(), "collapse_tree_pdf")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    for k,v in tre_dic.items():
        t=Tree(v)
        t.sort_descendants("support")
        ts=create_tree_style(k)
        t1=set_style_node(t,c_color_dic)
        t1.render(file_name=os.path.join(dir_path1, k + '.pdf'),tree_style=ts)

