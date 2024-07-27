from __init__ import *
import matplotlib
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from PIL import Image
from collections import defaultdict
import seaborn as sns


def merge_files(directory, suffix, output_file):
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for filename in os.listdir(directory):
            if filename.endswith(suffix):
                file_path = os.path.join(directory, filename)
                with open(file_path, 'r', encoding='utf-8') as infile:
                    for line in infile:
                        outfile.write(line)

def parse_hyde_out(filename):
    lst=[]
    with open(filename,'r') as f:
        for i in f :
            i1=i.strip().split('\t')
            if i1[0]=='P1':
                continue
            elif i1[5]=='nan':
                continue
            elif 0<float(i1[5])<1: 
                tup=[i1[0],i1[1],i1[2],i1[5]]
                lst.append(tup)
    return lst

def calculate_three_tup(hyde_out_lst):
    count_dic=defaultdict(list)
    for i in hyde_out_lst :
        key='-'.join(i[:3])
        count_dic[key].append(i)
    return count_dic

def get_hybrid_dic(summary_dic):
    hybrid_dic={}
    for k,v in summary_dic.items():
        hybrid_tup=k.split('-')
        if hybrid_tup[1] in hybrid_dic :
            hybrid_dic[hybrid_tup[1]].append(k)
        else:
            hybrid_dic[hybrid_tup[1]]=[k]
    return hybrid_dic

def calculate_gamma(lst):#一个三元组统计的list
    gamma=0
    for i in lst:
        if i[3]=='nan':
            continue
        elif i[3] =='-inf':
            continue
        else:
            num=float(i[3])
            gamma+=num 
    return gamma/len(lst)
    
def generate_tree_leaf(t, hybrid_sps,filename):
    nstyle = NodeStyle()
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    nstyle["vt_line_type"] = 0 
    nstyle["hz_line_type"] = 0 
    nstyle["size"] = 0
    nstyle["shape"] = "circle"
    nstyle["fgcolor"] = "black"
    
    leaf_names = t.get_leaf_names()
    max_l = max(len(name) for name in leaf_names)
    
    for i in t.traverse():
        i.set_style(nstyle)
        if i.is_leaf():
            leaf_name = i.name.ljust(max_l, '·')  # 使用填充确保字符长度相同
            fgcolor = 'red' if i.name == hybrid_sps else 'black'
            i.add_face(TextFace(leaf_name, fgcolor=fgcolor), column=1, position='aligned')
            
    ts = TreeStyle()
    ts.scale = 10
    ts.show_leaf_name = False
    ts.show_scale = False
    realign_branch_length(t)
    rejust_root_dist(t)
    
    t.render(file_name=filename+"_img_faces.png", h=3200, tree_style=ts)


    
def generate_tree_node(t,node,filename):
    nodes=[i.name for i in node.traverse()]
    max_l=max([len(j )for j in t.get_leaf_names()])
    def set_node_style(n, color):  
        nstyle = NodeStyle()  
        nstyle["vt_line_width"] = 1  
        nstyle["hz_line_width"] = 1  
        nstyle["vt_line_type"] = 0  
        nstyle["hz_line_type"] = 0  
        nstyle["size"] = 0  
        nstyle["shape"] = "circle"  
        nstyle["vt_line_color"] = color
        nstyle["hz_line_color"] = color
        n.set_style(nstyle) 
        
    for i in t.traverse():
        if i.name ==node.name:
            i.add_face(TextFace(i.name, fsize=4, fgcolor="blue"), column=1, position="branch-top")
        if i.name in nodes:  
            set_node_style(i, "red")  
        else:  
            set_node_style(i, "black")
            
    for leaf in t:
        if leaf.name in node.get_leaf_names():
            leaf.add_face(TextFace(leaf.name.ljust(max_l, '·') , fgcolor='red', fsize=10), column=1, position='aligned')
        else:
            leaf.add_face(TextFace(leaf.name.ljust(max_l, '·') , fgcolor='black', fsize=10), column=1, position='aligned')
            
    ts = TreeStyle()
    ts.scale = 10
    ts.show_leaf_name=False
    
    ts.show_scale = False
    
    realign_branch_length(t)
    rejust_root_dist(t)
    
    t.render(file_name=filename+"_img_faces.png",h=3200,tree_style=ts)

def from_summary_get_hyb_to_date(summary_dic, filter_dic,a_hyb_to_three_tup_list, leafs):
    tup_df = pd.DataFrame(index=leafs, columns=leafs, data=0)
    gamma_df = pd.DataFrame(index=leafs, columns=leafs, data=0)
    filtered_gamma_df=pd.DataFrame(index=leafs, columns=leafs, data=0)
    for i in a_hyb_to_three_tup_list:
        tup = i.split('-')
        p1,hyb,p2=tup
        same_tups = summary_dic[i]
        num = len(same_tups)
        gamma = calculate_gamma(same_tups)
        
        gamma_3=round(gamma, 3)
        tup_df.loc[p1, p2] = num
        gamma_df.loc[p1, p2] = gamma_3
        filtered_gamma_df[p1][p2]=len(filter_dic[i])
        
    return tup_df,gamma_df,filtered_gamma_df
        

def hyde_visual_cmap():
    cmap = []
    color = 1
    lucency = 0
    for each in range(5000):
        cmap.insert(each, np.array([0,0,color,lucency]))
        color = color - (1/5000)
        lucency = lucency + (1/5000)
    color = 0
    lucency = 1
    for each in range(5001, 9999):
        cmap.insert(each, np.array([color,0,0,lucency]))
        color = color + (1/5000)
        lucency = lucency - (1/5000)
    newcmp = ListedColormap(cmap)
    return newcmp

def get_necmap():
    def create_lightened_cmap(rgb_color, reverse=False):
        colors = [rgb_color]
        for i in range(50):  
            lightened_color = tuple((np.array(rgb_color) + (1 - np.array(rgb_color)) * (i / 50)).tolist())
            colors.append(lightened_color)
        cmap_values = np.linspace(0, 1, len(colors))
        if reverse:
            colors = colors[::-1]
        lightened_cmap = LinearSegmentedColormap.from_list('lightened_cmap', list(zip(cmap_values, colors)))
        return lightened_cmap
    
    rgb_color1 = (33/255, 205/255, 67/255)
    rgb_color2 = (33/255, 171/255, 205/255)


    lightened_cmap1 = create_lightened_cmap(rgb_color1,reverse=True)
    lightened_cmap2 = create_lightened_cmap(rgb_color2)

    combined_cmap = LinearSegmentedColormap.from_list('combined_cmap', np.vstack((lightened_cmap1(np.linspace(0, 1, 256)), lightened_cmap2(np.linspace(0, 1, 256)))))

    return combined_cmap

def create_hot_map_node(summary_dic,filter_dic,hyb_dic,node, t,filename):
    def average_dataframes(summary_tup_df): 
        tup_result_df = pd.DataFrame(0, index=summary_tup_df[0].index, columns=summary_tup_df[0].columns)  
        for df in summary_tup_df:  
            tup_result_df += df  
        tup_result_df = tup_result_df.div(len(summary_tup_df)) 
        return tup_result_df
    
    node_s=node.get_leaf_names()
    leafs=t.get_leaf_names()
    df_lst=[]
    for i in node_s:
        a_hyb_to_three_tup_list=hyb_dic[i]
        tup_df,gamma_df,filtered_gamma_df=from_summary_get_hyb_to_date(summary_dic,filter_dic,a_hyb_to_three_tup_list,leafs)
        df_lst.append((tup_df,gamma_df,filtered_gamma_df))

    summary_tup_df=[d[0] for d in df_lst] 
    summary_gamma_df=[d[1] for d in df_lst] 
    summary_filter_df=[d[2] for d in df_lst] 
    
    tup_result_df = average_dataframes(summary_tup_df).astype(int)
    gamma_result_df=average_dataframes(summary_gamma_df)
    filter_result_df=average_dataframes(summary_filter_df)
    
    for leaf in node_s:
        tup_result_df.loc[:, leaf] = 0  
        tup_result_df.loc[leaf] = 0 
        gamma_result_df.loc[:, leaf] = 0  
        gamma_result_df.loc[leaf] = 0 
    
    mask_upper = np.triu(np.ones_like(gamma_result_df, dtype=bool), k=1)
    tup_result_df[mask_upper] = 0
    gamma_result_df[mask_upper] = 0

    tup_annot = tup_result_df.astype(str).where(tup_result_df != 0, other='')
    gamma_annot = gamma_result_df.applymap(lambda x: f"{x:.3f}" if x != 0 else '')

    fig = plt.figure(figsize=(30, 30))

    border_width = 0.001
    ax_size = [0 + border_width, 0 + border_width, 1 - 2 * border_width, 1 - 2 * border_width]
    ax = fig.add_axes(ax_size)
    
    
    newcmp = hyde_visual_cmap()
    
    sns.heatmap(tup_result_df, annot=tup_annot, fmt="", cmap='Greys', ax=ax,
                annot_kws={'color':'#cd4321','size': 60,'va':'bottom'},
                xticklabels=False, yticklabels=False, cbar=False, linewidths=1.5, linecolor='black',square=True)

    # 第二个 heatmap，显示非零数值
    sns.heatmap(gamma_result_df, annot=gamma_annot, fmt="", cmap='Greys', ax=ax,
                annot_kws={'color':'#cd21ab','size': 60,'va':'top'},
                xticklabels=False, yticklabels=False, cbar=False, linewidths=1.5, linecolor='black',square=True)

    # 第三个 heatmap，不显示标记
    newcmp = hyde_visual_cmap()
    
    sns.heatmap(gamma_result_df, annot=False, cmap=newcmp, ax=ax,  
                xticklabels=False, yticklabels=False, cbar=False, cbar_kws={'shrink':0.6},
                vmin=0, vmax=1, linewidths=1.5, linecolor='black',square=True)
    
    m = ax.imshow(gamma_df, norm=colors.Normalize(vmin=0, vmax=1), cmap=newcmp)  
    position=fig.add_axes([0.9, 0.2, 0.05, 0.7])
    cbar = plt.colorbar(m, cax=position)
    cbar.ax.tick_params(labelsize=40)  
    cbar.set_label('Heatmap of Y Values', fontsize=40)
    plt.savefig(filename+"_hotmap.png", dpi=200)
    plt.cla()
    plt.close("all")
        
    
def create_hot_map_leaf(summary_dic, a_hyb_to_three_tup_list, t,filename,b1):
    sp = t.get_leaf_names()
    tup_df = pd.DataFrame(index=sp, columns=sp, data=0)
    gamma_df = pd.DataFrame(index=sp, columns=sp, data=0)
    filtered_gamma_df=pd.DataFrame(index=sp, columns=sp, data=0)
    for i in a_hyb_to_three_tup_list:
        tup = i.split('-')
        p1,hyb,p2=tup
        if p1 not in sp or p2 not in sp:
            print(f"Warning: {p1} or {p2} is not in sp.")
            continue  # 跳过不在 sp 中的组合
        same_tups = summary_dic[i]
        num = len(same_tups)
        gamma = calculate_gamma(same_tups)
        
        tup_df.loc[p1, p2] = num
        gamma_df.loc[p1, p2] = round(gamma, 3)
        filtered_gamma_df.loc[p1, p2]=len(b1[i])
        
    fig = plt.figure(figsize=(30, 30))
    
    border_width = 0.001
    ax_size = [0 + border_width, 0 + border_width, 1 - 2 * border_width, 1 - 2 * border_width]
    ax = fig.add_axes(ax_size)
    
    mask_upper = np.triu(np.ones_like(gamma_df, dtype=bool), k=1)
    gamma_df[mask_upper] = 0
    tup_df[mask_upper] = 0

    tup_annot = tup_df.astype(str).where(tup_df != 0, other='')
    gamma_annot = gamma_df.applymap(lambda x: f"{x:.3f}" if x != 0 else '')
    # 第一个 heatmap，显示非零数值
    sns.heatmap(tup_df, annot=tup_annot, fmt="", cmap='Greys', ax=ax,
                annot_kws={'color':'#cd4321','size': 60,'va':'bottom'},
                xticklabels=False, yticklabels=False, cbar=False, linewidths=1.5, linecolor='black',square=True)

    # 第二个 heatmap，显示非零数值
    sns.heatmap(gamma_df, annot=gamma_annot, fmt="", cmap='Greys', ax=ax,
                annot_kws={'color':'#cd21ab','size': 60,'va':'top'},
                xticklabels=False, yticklabels=False, cbar=False, linewidths=1.5, linecolor='black',square=True)

    # 第三个 heatmap，不显示标记
    newcmp = hyde_visual_cmap()
    
    sns.heatmap(gamma_df, annot=False, cmap=newcmp, ax=ax,  
                xticklabels=False, yticklabels=False, cbar=False, cbar_kws={'shrink':0.6},
                vmin=0, vmax=1, linewidths=1.5, linecolor='black',square=True)
    
    m = ax.imshow(gamma_df, norm=colors.Normalize(vmin=0, vmax=1), cmap=newcmp)  
    position=fig.add_axes([0.88, 0.2, 0.05, 0.7])
    cbar = plt.colorbar(m, cax=position)
    cbar.ax.tick_params(labelsize=40)  
    cbar.set_label('Heatmap of y Values', fontsize=50,labelpad=10)
    

    plt.savefig(filename+"_hotmap.png", dpi=300)
    plt.cla()
    plt.close("all")

#     plt.show()
    

def combine_fig(hybrid_species):
    treepic = Image.open(f"{hybrid_species}_img_faces.png")
    treepic_size = treepic.size
    combine_fig_size = treepic_size[0] + treepic_size[1]

    combine = Image.new("RGB", (combine_fig_size, combine_fig_size), "#FFFFFF")

    combine.paste(treepic, (40,40))
    treepic_rotate = treepic.rotate(90, expand = 1) #treepic旋转90度后粘贴
    combine.paste(treepic_rotate, (treepic_size[0] - 40, int(treepic_size[1]) - 40))

    hotpic = Image.open(f"{hybrid_species}_hotmap.png")
    hotpic.thumbnail((treepic_size[1], treepic_size[1]))
    combine.paste(hotpic, (treepic_size[0] - 40, 40))

    combine.save(hybrid_species + ".png")


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

def hyde_visual_leaf_main(hybrid_tracer,sptree):
    temp_directory = hybrid_tracer
    merge_files(temp_directory, '-out.txt', 'merged_out.txt')
    merge_files(temp_directory, '-filtered.txt', 'merged_filtered.txt')

    out1=parse_hyde_out('merged_out.txt')
    out2=parse_hyde_out('merged_filtered.txt')
    result1=calculate_three_tup(out1)
    result2=calculate_three_tup(out2)
    hybrid_dic1=get_hybrid_dic(result1)
    hybrid_dic2=get_hybrid_dic(result2)
    for k,v in hybrid_dic1.items():
        print(f'{k} is processing')
        t1=sptree.copy()
        #t2=rename_input_tre(t1,taxa_dic)
        generate_tree_leaf(t1,k,k)
        create_hot_map_leaf(result1,v,t1,k,result2)
        combine_fig(k)
        os.remove(f"{k}_hotmap.png")
        os.remove(f'{k}_img_faces.png')

def hyde_visual_node_main(hybrid_tracer,sptree):
    temp_directory = hybrid_tracer
    merge_files(temp_directory, '-out.txt', 'merged_out.txt')
    merge_files(temp_directory, '-filtered.txt', 'merged_filtered.txt')

    out1=parse_hyde_out('merged_out.txt')
    out2=parse_hyde_out('merged_filtered.txt')
    result1=calculate_three_tup(out1)
    result2=calculate_three_tup(out2)
    hybrid_dic1=get_hybrid_dic(result1)
    hybrid_dic2=get_hybrid_dic(result2)

    num_tre_node(sptree)
    nodes=[i for i in sptree.traverse() if not i.is_leaf()]

    for node in nodes:
        if node.is_root():
            continue
        else:
            print(f'{node.name} is processing')
            t1=sptree.copy()
            generate_tree_node(t1,node,node.name)
            create_hot_map_node(result1,result2,hybrid_dic1,node,t1,node.name)
            combine_fig(node.name)
            os.remove(f"{node.name}_hotmap.png")
            os.remove(f'{node.name}_img_faces.png')
if __name__ == "__main__":
    sptree=Tree('sptree.nwk')
    hyde_visual_leaf_main(sptree)
    hyde_visual_node_main(sptree) 

