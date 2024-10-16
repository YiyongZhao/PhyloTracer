from __init__ import *
def fold_same_taxa(t):
    t1=t.copy()
    for node in t1.traverse():
        if not node.is_leaf():
            taxa=get_species_list(node)
            if len(set(taxa))==1:
                node.name=taxa[0]
                for child in node.get_children():
                    child.detach()  
        else:
            node.name=node.name.split('_')[0]
    return t1

def get_only_sps_tree(Phylo_t):
    Phylo_t_c=Phylo_t.copy()
    for node in Phylo_t_c:
        node.name=node.name.split('_')[0]
    return Phylo_t_c
    
def get_max_tree(trees):
    max_leaf_count = 0
    max_leaf_count_tree = None

    for i in trees:
        leaf_count = len(i.get_leaves())
        if leaf_count > max_leaf_count:
            max_leaf_count = leaf_count
            max_leaf_count_tree = i

    return max_leaf_count_tree

def process(trees, result_dict):
    if len(trees) == 1:
        result_dict[trees[0].write(format=9)] = [trees[0]]
        return []
    else:
        same_rf_trees = []
        different_trees = []
        first_tree = trees[0]
        for i in trees:
                
            rf = first_tree.robinson_foulds(i)[0]
            if rf == 0:
                same_rf_trees.append(i)
            else:
                different_trees.append(i)
           
            
        max_tree = get_max_tree(same_rf_trees)
        #result_dict[max_tree.write(format=9)] = len(same_rf_trees)
        result_dict[max_tree.write(format=9)] = same_rf_trees
            
        #o.write(max_tree.write(format=9)+'\t'+str(len(same_rf_trees))+'\n')
        return different_trees
      
def process_tree(trees,result_dict):
    remaining_trees = trees
    while len(remaining_trees) >= 1:
        remaining_trees = process(remaining_trees, result_dict)

def get_summary_result(outfile,dic,new_named_gene2gene_dic):
    with open(f'relative_{outfile}.txt', 'w') as file:
        file.write(f'topology_id\ttopology_num\ttips_num\ttopology\n')
        for index,clade in enumerate(dic.keys()):
            tres=dic[clade]
            tres1=sorted(tres,key=lambda x: len(x), reverse=True)
            for tree in tres1:
                tree1=rename_input_tre(tree,new_named_gene2gene_dic)
                tree_str=tree1.write(format=9)
                file.write(str(index)+'\t'+str(len(dic[clade]))+'\t'+str(len(tree))+'\t'+tree_str+'\n')

def get_absolutely_result(outfile,trees,voucher2taxa_dic):
    topolo_dic={}
    for i in trees:
        t=fold_same_taxa(i)
        t_str=t.write(format=9)
        if t_str in topolo_dic:
            topolo_dic[t_str]+=1
        else:
            topolo_dic[t_str]=1
            
    sorted_dict = dict(sorted(topolo_dic.items(), key=lambda x: len(x[0]), reverse=True))
    with open(f'absolute_{outfile}.txt', 'w') as file:
        file.write(f'topology_id\ttopology_num\ttopology\n')
        for num,k in enumerate(sorted_dict):
            t1=Tree(k)
            t2=rename_input_tre(t1,voucher2taxa_dic)
            t1_str=t2.write(format=9)
            file.write(f'{num}\t{sorted_dict[k]}\t{t1_str}\n')

def statistical_main(tre_dic,outfile,gene2new_named_gene_dic,voucher2taxa_dic):
    only_sptrees=[]
    sctrees=[]
    for k,v in tre_dic.items():
        t=read_tree(v)
        t1=rename_input_tre(t,gene2new_named_gene_dic)
        t1.sort_descendants()
        t2=get_only_sps_tree(t1)
        if len(get_species_set(t2))==len(t2):
            sctrees.append(t2)
        only_sptrees.append(t2)
        
    get_absolutely_result(outfile,only_sptrees,voucher2taxa_dic)
    
    dic={}
    process_tree(sctrees,dic)
    get_summary_result(outfile,dic,voucher2taxa_dic)
    

if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt')   
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    outfile='result'
    statistical_main(tre_dic,outfile,gene2new_named_gene_dic,new_named_gene2gene_dic)
