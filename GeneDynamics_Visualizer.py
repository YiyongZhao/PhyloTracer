from __init__ import *

def create_tree_style():
    ts = TreeStyle()
    ts.legend.add_face(TextFace('gene_gain_loss', fsize=20, fgcolor="black"), column=0)
    ts.legend_position = 1
    ts.mode ='r'
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = True
    ts.show_branch_support = True
    ts.extra_branch_line_type =0
    ts.extra_branch_line_color='black'
  
    return ts
  
def get_gain_and_loss_dic(summary_tree):
    with open(summary_tree, 'r') as file:
        tree_data = file.read()
    dic = {}
    for line in tree_data.split('\n'):
        line=line.replace('|','')
        line=line.replace(' ','')
        if line.strip() != '':
            node, value = line.split('[')
            dic[node.strip()] = '[' + value.strip()
    return dic

def mark_sptree(sptree,dic):
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
        if i.name in dic:
            v = dic[i.name].replace('[', '').replace(']', '')
            com1=TextFace('[',fgcolor='black')
            i.add_face(com1,column=1,position='branch-top')

            com2=TextFace(v.split(',')[0],fgcolor='blue')
            i.add_face(com2,column=2,position='branch-top')
            com3=TextFace(',',fgcolor='black')
            i.add_face(com3,column=3,position='branch-top')

            com4=TextFace(v.split(',')[1],fgcolor='red')
            i.add_face(com4,column=4,position='branch-top')
            #com5=TextFace(',',fgcolor='black')
            #i.add_face(com5,column=5,position='branch-top')

            #com6=TextFace(v.split(',')[2],fgcolor='black')
            #i.add_face(com6,column=6,position='branch-top')
            com7=TextFace(']',fgcolor='black')
            i.add_face(com7,column=7,position='branch-top')

def new_num_tre_node(Phylo_t:object)->object:#Numbering the nodes in the tree.
    i = 0
    for node in Phylo_t.traverse('postorder'):
        if not node.is_leaf():
            node.name = "node_" + str(i)
            i += 1
    return Phylo_t
    
if __name__ == "__main__":
    dic=get_gain_and_loss_dic('example_data/Gene_Gain_And_Loss_Visualization/input/summary_tree')
    t=Tree('example_data/Gene_Gain_And_Loss_Visualization/input/sptree.nwk')
    new_num_tre_node(t)
    t.ladderize()
    t.sort_descendants("support")
    t.sort_descendants()
    ts=create_tree_style()
    mark_sptree(t,dic)
    t.render(file_name='gene_gain_loss.pdf',tree_style=ts)
