from __init__ import *
import os

def find_dup_node(Phylo_t:object)->list:#After searching for duplication events, the list of node names where duplication events occurred is as follows:
    events = Phylo_t.get_descendant_evol_events()
    dup_node_name_list = []
    for ev in events:
        if ev.etype == "D":
            i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
            events_node_name_list = i.split(',')
            common_ancestor_node_name = Phylo_t.get_common_ancestor(events_node_name_list)
            dup_node_name_list.append(common_ancestor_node_name.name)
    return dup_node_name_list

def get_empty_count_dict(sptree: object) -> dict:
    empty_count_dic={}
    for node in sptree.traverse():
        if not node.is_leaf():
            clade=node.write(format=9)
            empty_count_dic[clade]=0
    return empty_count_dic

def switch_tree(phylo_tree):
    
    down_child = phylo_tree.get_children()[1]
    up_child = phylo_tree.get_children()[0]
    
    new_tree = Tree()
    new_tree.add_child(down_child)
    new_tree.add_child(up_child)
    
    return new_tree

def get_only_sps_tree(Phylo_t):
    for node in Phylo_t:
        node.name=node.name[0:3]
    return Phylo_t

def folding_tree(Phylo_t):
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            sps=get_species_list(node)
            if len(set(sps))==1:
                node.name=sps[0]
                for child in node.get_children():
                    child.detach() 
    return Phylo_t

def statistical_clade(clade_dic,Phylo_t):
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            clade=node.write(format=9)
            if clade in clade_dic :
                clade_dic[clade]+=1
            else:
                clade_dic[clade]=1
def rejust_clade_dic(clade_dic):
    new_dic={}
    for k,v in clade_dic.items():
        t=Tree(k)
        t1=switch_tree(t)
        t2=t1.write(format=9)
        if k and t2 not in new_dic :
            new_dic[k]=v
        if t2 in new_dic :
            new_dic[t2]=new_dic[t2]+v
            
    return new_dic

def get_duplicated_species_num(node):
    species_lst=get_species_list(node)
    duplicates = []
    unique_elements = set()

    for element in species_lst:
        if element in unique_elements:
            duplicates.append(element)
        else:
            unique_elements.add(element)

    duplicate_count = len(set(duplicates))

    return duplicate_count

def get_intraspecies_and_interspecies_node_list(dup_node_name_list,Phylo_t0):
    intraspecies_node_list=[]
    interspecies_node_list=[]
    for node in dup_node_name_list:
        clade=Phylo_t0&node
        sps=get_species_list(clade)
        if len(set(sps))==1 :
            intraspecies_node_list.append(node)
        else:
            interspecies_node_list.append(node)
    return intraspecies_node_list,interspecies_node_list
    
def get_target_node_list(dup_node_name_list,Phylo_t0,dup_percent,dup_species_num):
    target_node_list=[]
    for node in dup_node_name_list:
        clade=Phylo_t0&node
        sps=get_species_list(clade)
        if get_duplicated_species_num(clade)/len(set(sps)) >=dup_percent and get_duplicated_species_num(clade)>=dup_species_num:
            target_node_list.append(node)
    return target_node_list


def get_summary_GD_node(tre_dic,dup_percent,dup_species_num):
  outfile=open('Summary_GD.txt','w')
  outfile.write('tre_ID'+'\t'+'Basic_gene_duplication_event'+'\t'+'Intraspecific_gene_duplication_event'+'\t'+'Interspecific_gene_duplication_event'+'\t'+'Targetspecific_gene_duplication_event'+'\n')
  
  for tre_ID,tre_path in tre_dic.items():
      outfile.write(str(tre_ID)+'\t')
      Phylo_t0 = read_phylot_tree(tre_path)
      num_tre_node(Phylo_t0)
      dup_node_name_list=find_dup_node(Phylo_t0)
      intraspecies_node_list,interspecies_node_list=get_intraspecies_and_interspecies_node_list(dup_node_name_list,Phylo_t0)
      target_node_list=get_target_node_list(dup_node_name_list,Phylo_t0,dup_percent,dup_species_num)
      Phylo_t0.write(outfile='num_tree/'+str(tre_ID)+'.nwk',format=1)
      if dup_node_name_list != None:
          outfile.write('-'.join(dup_node_name_list)+'\t')
      else:
          outfile.write('\n')
      if intraspecies_node_list != None:
          outfile.write('-'.join(intraspecies_node_list)+'\t')
      else:
          outfile.write('\n')
      if interspecies_node_list != None:
          outfile.write('-'.join(interspecies_node_list)+'\t')
      else:
          outfile.write('\n')
      if target_node_list != None:
          outfile.write('-'.join(target_node_list)+'\t')
      else:
          outfile.write('\n')
      
  outfile.close()
  
        

if __name__ == "__main__":
  os.makedirs(os.path.join(os.getcwd(), "num_tree"))
  tre_dic=read_and_return_dict('GF_list.txt')   
  clade_dic={}
  for tre_path in tre_dic.values():
      Phylo_t0 = read_tree(tre_path)
      Phylo_t1=get_only_sps_tree(Phylo_t0)
      Phylo_t2=folding_tree(Phylo_t1)
      statistical_clade(clade_dic,Phylo_t2)
  new_dic=rejust_clade_dic(clade_dic)
  with open ('Statistical_calde.txt','w') as f :
    for k,v in new_dic.items():
        f.write(k+'\t'+str(v)+'\n')

  dup_percent=0.5
  dup_species_num=2
  get_summary_GD_node(tre_dic,dup_percent,dup_species_num)
      
  
    




