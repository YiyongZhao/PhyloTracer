from __init__ import *

def get_empty_count_dict(sptree: object) -> dict:
    empty_count_dic = {node.name: 0 for node in sptree.traverse()}
    return empty_count_dic

def judge_support(support,support_value):
    if  support <=1 and 0.5 <=support_value <=1:
        if support>support_value:
            return True
        else:
            return False
        
    elif support <=1 and 50 <= support_value <=100:
        support_value=support_value/100
        if support>support_value:
            return True
        else:
            return False
    elif support > 1 and 0.5 <=support_value <=1:
        support_value=support_value*100
        if support>support_value:
            return True
        else:
            return False
    elif support > 1 and 50 <=support_value <=100:
        if support>support_value:
            return True
        else:
            return False
    
def batch_gfs_traverse(tre_dic: dict, support_value: int, empty_count_dic: dict,sptree:object) -> dict:
    for tre_path in tre_dic.values():
        Phylo_t0 = read_tree(tre_path)
        Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
        
        tre_ParaL,GF_leaves_S =find_tre_dup(Phylo_t1)
        for dup in tre_ParaL:
            leafs=dup.replace('<=>',',')
            leaf_list=leafs.split(',')
            ancestor=Phylo_t1.get_common_ancestor(leaf_list)
            if judge_support(ancestor.support,support_value):
                sps_list=set([leaf[0:3] for leaf in leaf_list])
                if len(sps_list)>1:
                    parent=sptree.get_common_ancestor(sps_list)
                    empty_count_dic[parent.name]+=1
                else:
                    single_sps=list(sps_list)[0]
                    parent=sptree&single_sps
                    empty_count_dic[parent.name]+=1
    return empty_count_dic

def mark_sptree(sptree:object,empty_count_dic:dict)->object:
    sptree.ladderize()
    sptree.sort_descendants("support")
    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        node.set_style(nstyle)
        if not node.is_leaf() and node.name in empty_count_dic:
            leafs=node.get_leaf_names()
            num_list=[empty_count_dic[leaf] for leaf in leafs]
            num_list.append(empty_count_dic[node.name])
            num=str(sum(num_list))
            node.add_face(TextFace(num+' GD', fsize=5, fgcolor="red"), column=0, position="branch-top")
            
    return sptree.render('species_tree_GD_Detector.PDF')
    
if __name__ == "__main__":
    sptree=PhyloTree('F:/a/96tree/7sp.nwk')
    num_tre_node(sptree)
    tre_dic=read_and_return_dict('F:/a/96tree/96.txt')
    empty_count_dic=get_empty_count_dict(sptree)
    support_value=50
    empty_count_dic=batch_gfs_traverse(tre_dic, support_value, empty_count_dic,sptree) 
    mark_sptree(sptree,empty_count_dic)

