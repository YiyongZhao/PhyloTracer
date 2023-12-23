from __init__ import *

def create_tree_style():
    ts = TreeStyle()
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

def mark_sptree(sptree):
    for i in t.traverse():
        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0 
        nstyle["hz_line_type"] = 0 
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        i.set_style(nstyle)
        if i.name in d:
            v = d[i.name].replace('[', '').replace(']', '')
            com1=TextFace('[',fgcolor='black')
            i.add_face(com1,column=1,position='branch-top')

            com2=TextFace(v.split(',')[0],fgcolor='blue')
            i.add_face(com2,column=2,position='branch-top')
            com3=TextFace(',',fgcolor='black')
            i.add_face(com3,column=3,position='branch-top')

            com4=TextFace(v.split(',')[1],fgcolor='red')
            i.add_face(com4,column=4,position='branch-top')
            com5=TextFace(',',fgcolor='black')
            i.add_face(com5,column=5,position='branch-top')

            com6=TextFace(v.split(',')[2],fgcolor='black')
            i.add_face(com6,column=6,position='branch-top')
            com7=TextFace(']',fgcolor='black')
            i.add_face(com7,column=7,position='branch-top')

if __name__ == "__main__":
    dic=get_gain_and_loss_dic('/Volumes/Elements/file/q_out/summary.tree')
    t=Tree('/Users/apple/Desktop/jupyter_notebook/qiao/28sps_root.nwk')
    t.ladderize()
    t.sort_descendants("support")
    t.sort_descendants()
    ts=create_tree_style()
    mark_sptree(t,dic)
    t.render(file_name='gene_gain_loss.pdf',tree_style=ts)
