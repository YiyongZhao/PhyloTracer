from __init__ import *
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from ete3 import RectFace

def realign_branch_length(Phylo_t1:object)->object:
    Phylo_t1.resolve_polytomy(recursive=True)
    Phylo_t1.sort_descendants("support")
    max_deep=get_max_deepth(Phylo_t1)
    nstyle = NodeStyle()
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    nstyle["vt_line_type"] = 0 
    nstyle["hz_line_type"] = 0 
    nstyle["size"] = 0
    nstyle["shape"] = "circle"
    nstyle["fgcolor"] = "black"
    for node in Phylo_t1.traverse():
        node.set_style(nstyle)
        
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
    Phylo_t2=rejust_root_dist(Phylo_t1)
    
    return Phylo_t2

def get_color_dict(dictory:dict)->dict:
    colormap = plt.get_cmap("rainbow")
    # color dictionary
    unique_values=set(dictory.values())   
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_list = {k: v + '@' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}
    
    return sps_color_list

def create_tree_style():
    ts = TreeStyle()
    ts.mode = 'r'
    ts.show_scale = False
    ts.scale=20
    ts.show_border = False
    ts.margin_bottom = 10
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.extra_branch_line_type =0
    ts.extra_branch_line_color='black'
    ts.branch_vertical_margin = -1
    
    
    return ts

def add_label_to_tree(tree,color_dic):
    
    for i in tree:
        sps=i.name.split('_')[0]
        if sps in color_dic:
            color = color_dic[sps].split('@')[-1]
            face=TextFace(' ' + '_'.join(i.name.split('_')[1:]),fgcolor='black')
            face1=TextFace(' ' + color_dic[sps].split('@')[0],fgcolor=color)
            i.add_face(face, column=0, position='aligned')
            i.add_face(face1, column=1, position='aligned')
            
    realign_branch_length(tree)
    
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
        
def add_heat_map_to_node(tree,df):
    for node in tree:
        if node.is_leaf() and node.name in df.index:
            n=2
            for i in df.columns.tolist():
                color = get_color(df[i][node.name])
                
                face=RectFace(width=10, height=10, fgcolor=color,bgcolor=color)
                node.add_face(face, column=n, position='aligned')
                n+=1
                
def add_header_to_tree(ts,df):
    labels = df.columns.to_list()
    n=2
    for i in labels:
        face=TextFace(i,fgcolor='black',fsize=9)
        face.rotation=-90
        face.vt_align =2
        ts.aligned_header.add_face(face,n)
        n+=1
        
def add_color_bar(ts):
    bar_face=TextFace('  Color Bar Title  ')
    ts.title.add_face(bar_face, column=0)
    cols=['#006599', '#408ca6', '#7fb2b2', '#bfd9bf', '#ffffcc','#f7deab', '#eebc88', '#e69966', '#dc7845', '#d55623', '#cc3300']
    bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    col_dic=dict(zip(bounds,cols))
    n=1
    for k,v in col_dic.items():
        colorbar_face = RectFace(width=20, height=20, fgcolor=v,bgcolor=v)
        ts.title.add_face(TextFace(' '+str(k)), column=n)
        ts.title.add_face(colorbar_face, column=n)
        n+=1

def gene_expression_main(tre_dic,expression_dic,color_dic):
	dir_path = os.path.join(os.getcwd(), "output/gene_expression_pdf_result/")
	if os.path.exists(dir_path):
		shutil.rmtree(dir_path)
	os.makedirs(dir_path)
	pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
	for tre_ID,tre_path in tre_dic.items():
		pbar.set_description(f"Processing {tre_ID}")
		t=read_tree(tre_path)
		t.sort_descendants('support')
		df = pd.read_excel(expression_dic[tre_ID], index_col=0)
		ts=create_tree_style()
		add_label_to_tree(t,color_dic)
		add_heat_map_to_node(t,df)
		add_header_to_tree(ts,df)
		add_color_bar(ts)
		t.render(dir_path+tre_ID+'_gene_expression.pdf',tree_style=ts,w=210, units="mm")

		pbar.update(1)
	pbar.close()

if __name__ == "__main__":
	os.makedirs(os.path.join(os.getcwd(), "geneexpression_pdf_result"))
	label=read_and_return_dict('order.txt')
	color_dic=get_color_dict(label)
	tre_dic=read_and_return_dict('gf')
	expression_dic=read_and_return_dict('expression_gf')
	gene_expression_main(tre_dic,expression_dic,color_dic)