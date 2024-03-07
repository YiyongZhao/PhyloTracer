from __init__ import *


def rename_input_tree(Phylo_t:object, gene2new_named_gene_dic:dict) -> object:
    Phylo_t1=Phylo_t.copy()
    for node in Phylo_t1 :
        sps=node.name.split('_')[0]
        if sps in gene2new_named_gene_dic:
            node.name=gene2new_named_gene_dic[sps]+'_'+sps
    return Phylo_t1

def fold_same_taxa_same_sp(t):
    t1=t.copy()
    for node in t1.traverse():
        if not node.is_leaf():
            taxa=get_species_list(node)
            if len(set(taxa))==1:#同一个科或者目
                sps=set([j.split('_')[1] for j in node.get_leaf_names()])
                if len(sps)==1:#同一个物种
                    node.name=taxa[0]+'_'+list(sps)[0]#+'_'+str(len(taxa))
                    #node.add_feature('node2taxa',str(len(taxa)))
                    
                    for child in node.get_children():
                        child.detach()  
    return t1


def fold_same_taxa(t,sptree,node2taxa_dic):
    t1=t.copy()
    
    for node in t1.traverse():
        if not node.is_leaf():
            taxa=get_species_set(node)
            if len(taxa)==1:#指同一个taxa的枝
                sps=set([j.split('_')[1] for j in node.get_leaf_names()])#得到物种
                com=get_mapp(sptree,sps)#得到这个枝到物种树的映射枝的node编号
                map_list=get_mapp_set(node,sptree)#得到这个枝下所有枝到物种树的映射枝的list
                if com in node2taxa_dic:#指这个node枝可以代表这个taxa
                    if com in map_list:#这个taxa枝发生gd,把这个taxa折叠为taxa_1,taxa_2
                        child1=node.get_children()[0]
                        child2=node.get_children()[1]
                        child1.name=list(taxa)[0]+'_1'
                        child2.name=list(taxa)[0]+'_2'
    
                        if len(child1) !=1:
                            for child in child1.get_children():
                                child.detach() 
                        if len(child2) !=1:
                            for child in child2.get_children():
                                child.detach()
                            
                    else:#这个taxa枝没有发生gd 如过这个taxa枝的上下分枝到物种树到映射枝在同一个路径 折叠为taxa_partial,不在同一个路径折叠为taxa_partial_1,taxa_partial_2
                        child1=node.get_children()[0]
                        child2=node.get_children()[1]#得到上下分枝
                        
                        sps1=set([j.split('_')[1] for j in child1.get_leaf_names()])
                        com1=sptree.get_common_ancestor(sps1)#上分枝到物种树到映射枝
                        sps2=set([j.split('_')[1] for j in child2.get_leaf_names()])
                        com2=sptree.get_common_ancestor(sps2)#下分枝到物种树到映射枝
                        
                        path1=com1.get_ancestors()#路径
                        if com2.name in path1:
                            node.name=list(taxa)[0]+'_partial'
                            for child in node.get_children():
                                child.detach() 
                        else:
                            child1.name=list(taxa)[0]+'_partial_1'
                            child2.name=list(taxa)[0]+'_partial_2'

                            if len(child1) !=1:
                                for child in child1.get_children():
                                    child.detach() 
                            if len(child2) !=1:
                                for child in child2.get_children():
                                    child.detach()
                else:
                    node.name=list(taxa)[0]+'_partial'
                    for child in node.get_children():
                        child.detach()
            else:
                continue
    
                    
        else:
            if node.name.split('_')[1] in node2taxa_dic:
                node.name=node.name.split('_')[0]
            
                      
    return t1

        
def get_mapp(sptree, sps_set):
    if len(sps_set)!=1:
        com = sptree.get_common_ancestor(sps_set)
        return com.name
    else:
        com=sptree&list(sps_set)[0]
        return com.name

def get_mapp_set(t, sptree):#t是改过名字的树tips id 为 taxa+'_'+species
    map_sptree_set = set()
    for i in t.traverse("postorder"):
        if not i.is_leaf() and len(i) != len(t):
            sps_set=set([j.split('_')[1] for j in i.get_leaf_names()])
            map_node_name = get_mapp(sptree, sps_set)
            map_sptree_set.add(map_node_name)
    return map_sptree_set


def get_single_clades(Phylo_t:object,empty_set:set):
    if len(get_species_set(Phylo_t))==1:
        empty_set.add(Phylo_t)
        return
    for i in Phylo_t.get_children():
        get_single_clades(i,empty_set)
        
def get_spree_node_to_taxa(sptree:object)->dict:
    taxa=set()
    get_single_clades(sptree,taxa)
    node2taxa_dic={}
    taxa2sps_dic={}
    o=open('node2taxa.txt','w')
    for i in taxa:
        sps=[leaf.split('_')[1] for leaf in i.get_leaf_names()]
        if i.is_leaf():
            node2taxa_dic[i.name.split('_')[1]]=i.name.split('_')[0]
            o.write(i.name.split('_')[1]+'\t'+i.name.split('_')[0]+'\n')
            taxa2sps_dic[list(get_species_set(i))[0]]=sps
        else:
            node2taxa_dic[i.name]=list(get_species_set(i))[0]
            o.write(i.name+'\t'+list(get_species_set(i))[0]+'\n')
            
            taxa2sps_dic[list(get_species_set(i))[0]]=sps
    o.close()
    return node2taxa_dic,taxa2sps_dic
   

def collapse_main(tre_dic,taxa_dic,sptree):
    dir_path1 = os.path.join(os.getcwd(), "collapse_tree")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    named_sptree=rename_input_tree(sptree,taxa_dic)
    num_tre_node(named_sptree)
    node2taxa_dic,taxa2sps_dic=get_spree_node_to_taxa(named_sptree)

    #color_dic=get_color_dict(taxa_dic)
    for k,v in tre_dic.items():
        t=Tree(v)

        num_tre_node(t)
        t1=rename_input_tree(t,taxa_dic)
        t2=fold_same_taxa_same_sp(t1)
        t3=fold_same_taxa(t2,sptree,node2taxa_dic)
        t3.sort_descendants()
        t3.write(outfile=os.path.join(dir_path1, k + '.nwk'),format=0,dist_formatter='%.10f')


if __name__ == "__main__":
    os.makedirs(os.path.join(os.getcwd(), "pruned_tree"))
    taxa_dic=read_and_return_dict('taxa.txt')
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    sptree=Tree('sptree.nwk')
    collapse_main(tre_dic,taxa_dic,sptree)
   

