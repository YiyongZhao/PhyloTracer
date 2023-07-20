from __init__ import *

def get_nodes(Phylo_t: object) -> list:
    nodes = [node for node in Phylo_t.traverse() if not node.is_leaf()]
    return nodes

def Empty_matrix_creation(Leaves: list) -> pd.DataFrame:
    df_empty = pd.DataFrame(0, columns=Leaves, index=Leaves)
    return df_empty

def get_empty_count_dict(Leaves: list) -> dict:
    empty_count_dic = {}
    for L1 in Leaves:
        for L2 in Leaves:
            empty_count_dic[L1 + ":" + L2] = 0
            empty_count_dic[L2 + ":" + L1] = 0
    return empty_count_dic

def Batch_GFs_traverse(tre_dic: dict, support_value: int, empty_count_dic: dict,imap:dict) -> dict:
    for tre_path in tre_dic.values():
        Phylo_t0 = read_tree(tre_path)
        Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
        Nodes = get_nodes(Phylo_t1)
        for index, Node in enumerate(Nodes):
            Node_confidence = Node.support

            sub_clade1, sub_clade2 = Node.children
            
            Leaves_sub_clade1 = sub_clade1.get_leaf_names()
            Leaves_sub_clade2 = sub_clade2.get_leaf_names()
            
            Taxon_sub_clade1 = {imap[leaf] for leaf in Leaves_sub_clade1}
            Taxon_sub_clade2 = {imap[leaf] for leaf in Leaves_sub_clade2}

            if Taxon_sub_clade1.isdisjoint(Taxon_sub_clade2):
                if Node_confidence is not None and Node_confidence >= support_value:
                    for L1 in Leaves_sub_clade1:
                        for L2 in Leaves_sub_clade2:
                            empty_count_dic[imap[L1] + ":" + imap[L2]] += 1
                            empty_count_dic[imap[L2] + ":" + imap[L1]] += 1
    return empty_count_dic

def sps_pair_count_dict(sptree: object, tre_dic: dict, imap: dict) -> dict:
    num_tre_node(sptree)
    Leaves = sptree.get_leaf_names()
    empty_count_dic = get_empty_count_dict(Leaves)
    df_empty = Empty_matrix_creation(Leaves)
    dic_pair_count = Batch_GFs_traverse(tre_dic, 50, empty_count_dic, imap)

    sps_pair_count_dic = {}
    for node in sptree.traverse():
        if not node.is_leaf():
            leafs = node.get_leaf_names()
            node_pair_list = get_empty_count_dict(leafs)

            num_list = []
            for j in dic_pair_count.keys():
                if j in node_pair_list:
                    num_list.append(dic_pair_count[j])
            total = sum(num_list) / 2

            sps_pair_count_dic[node.name] = total

    return sps_pair_count_dic

def Summary_count(dic_pair_count: dict, empty_count_dic: pd.DataFrame) -> pd.DataFrame:
    for key, value in dic_pair_count.items():
        L1, L2 = key.split(":")
        df_empty.loc[L1, L2] = value
    return df_empty

def sps_pair_count_dict(sptree: object, tre_dic: dict, imap: dict) -> dict:
    num_tre_node(sptree)
    Leaves = sptree.get_leaf_names()
    empty_count_dic = get_empty_count_dict(Leaves)
    df_empty = Empty_matrix_creation(Leaves)
    dic_pair_count = Batch_GFs_traverse(tre_dic, 50, empty_count_dic, imap)

    sps_pair_count_dic = {}
    for node in sptree.traverse():
        if not node.is_leaf():
            leafs = node.get_leaf_names()
            node_pair_list = get_empty_count_dict(leafs)

            num_list = []
            for j in dic_pair_count.keys():
                if j in node_pair_list:
                    num_list.append(dic_pair_count[j])
            total = sum(num_list) / 2

            sps_pair_count_dic[node.name] = total

    return sps_pair_count_dic

def mark_sptree(sptree:object,sps_pair_count_dic:dict)->object:
    sptree.ladderize()
    sptree.sort_descendants("support")
    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        node.set_style(nstyle)
        if not node.is_leaf() and node.name in a:
            num=str(sps_pair_count_dic[node.name])
            node.add_face(TextFace(num+' GD', fsize=5, fgcolor="red"), column=0, position="branch-top")
            
    return sptree.render('species_tree_GD_Detector.PDF')
if __name__ == "__main__":
    sptree=PhyloTree('F:/a/96tree/7sp.nwk')
    tre_dic=read_and_return_dict('F:/a/96tree/96.txt')
    imap=read_and_return_dict('F:/a/96tree/gene2taxa.list')
    sps_pair_count_dic=sps_pair_count_dict(sptree,tre_dic,imap)
    mark_sptree(sptree,sps_pair_count_dic)

