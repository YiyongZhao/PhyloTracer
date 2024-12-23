from __init__ import *
from collections import Counter
def calculate_depth(node1, node2):
    if node1 in node2.iter_ancestors() or node2 in node1.iter_ancestors():
        distance = node1.get_distance(node2, topology_only=True)+2
        return distance
    common_ancestor = node1.get_common_ancestor(node2)
    return abs(node1.get_distance(common_ancestor, topology_only=True) - \
           node2.get_distance(common_ancestor, topology_only=True))

def mapp_gene_tree_to_species(sp_set,sptree):
    if len(sp_set) !=1:
        clade=sptree.get_common_ancestor(sp_set)
    else:
        clade=sptree&list(sp_set)[0]

    return clade

def are_sister_supports_greater_than_num(sister1, sister2, clade_support):
    support1 = sister1.support if hasattr(sister1, 'support') else 0
    support2 = sister2.support if hasattr(sister2, 'support') else 0
    return support1 > clade_support and support2 > clade_support

def find_dup_node(Phylo_t: object, sptree:object,gd_support: int = 50,clade_support:int=0,dup_species_num:int=2,dup_species_percent:int=0,deepvar:int=1) -> list:
    dup_node_name_list = []
    events = Phylo_t.get_descendant_evol_events()
    for ev in events:
        if ev.etype == "D":
            i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
            events_node_name_list = i.split(',')
            common_ancestor_node = Phylo_t.get_common_ancestor(events_node_name_list)
            child1,child2=common_ancestor_node.get_children()
            sp_set=get_species_set(common_ancestor_node)
            mapp_sp_node=mapp_gene_tree_to_species(sp_set,sptree)
            common_ancestor_node.add_feature('map',mapp_sp_node.name)

            if judge_support(common_ancestor_node.support,gd_support):

                child1, child2 = common_ancestor_node.get_children()
                c1=mapp_gene_tree_to_species(get_species_set(child1),sptree)
                c2=mapp_gene_tree_to_species(get_species_set(child2),sptree)
                #dup_sps=sps_dup_num(get_species_list(common_ancestor_node),get_species_set(common_ancestor_node))
                #dup_percent=dup_sps/len(get_species_set(common_ancestor_node))
                #print(sptree.get_distance(c1,c2, topology_only=True))
                if sptree.get_distance(c1,c2, topology_only=True) <=deepvar:
                #if dup_sps>=dup_species_num and dup_percent>=dup_species_percent:
                #if are_sister_supports_greater_than_num(child1,child2,clade_support):
                    dup_node_name_list.append(common_ancestor_node)
    return dup_node_name_list

def write_gd_result(filename, tre_dic, gd_support,clade_support,dup_species_percent, dup_species_num,sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,deepvar):
    with open(filename,'w') as file:
        file.write('#tree_ID'+'\t'+'gd_id'+'\t'+'gd_support'+'\t'+'gene1'+'\t'+'gene2'+'\t'+'level'+'\t'+'species'+'\t'+'GD_dup_sps'+'\t'+'dup_ratio'+'\t'+'gd_type'+'\t'+'comment'+'\n') 
        gd_num=1
        gd_type_dic={}
        for tre_ID, tre_path in tre_dic.items():
            Phylo_t0 = read_phylo_tree(tre_path)
            Phylo_t0=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
            #Phylo_t1=root_tre_with_midpoint_outgroup(Phylo_t0)
            num_tre_node(Phylo_t0)
            dup_node_name_list = find_dup_node(Phylo_t0,sptree,gd_support,clade_support,dup_species_num,dup_species_percent,deepvar)
            #Phylo_t0.write(outfile='num_tree/' + str(tre_ID) + '.nwk', format=1)
            # if tre_ID=='OG_104386_20':
            # for i in dup_node_name_list:
                
            #     clade=i
            #     parent=clade.map
            #     if parent=='N3':
            #         s=rename_input_tre(clade,new_named_gene2gene_dic)
            #         print(s)
            #         # child1, child2 = clade.get_children()
            #         # c1=mapp_gene_tree_to_species(get_species_set(child1),sptree)
            #         # c2=mapp_gene_tree_to_species(get_species_set(child2),sptree)
            #         # if c1.get_distance(c2, topology_only=True) ==2:
            #         s=rename_input_tre(clade,new_named_gene2gene_dic)
            #         print(tre_ID)
            #         print(s)

            for i in dup_node_name_list:
                sp_set=get_species_set(i)
                mapp_sp_node=mapp_gene_tree_to_species(sp_set,sptree)
                clade=i
                parent=clade.map
                child1,child2=clade.get_children()
                dup_sps=sps_dup_num(get_species_list(clade),get_species_set(clade))
                dup_percent=dup_sps/len(get_species_set(clade))
                model=get_model(clade,sptree)
                gene_pair1=gene_pair(clade)
                null='null'
                if not mapp_sp_node.is_leaf():
                    if mapp_sp_node.name in gd_type_dic:
                        gd_type_dic[mapp_sp_node.name].append(model)
                    else:
                        gd_type_dic[mapp_sp_node.name]=[model]

                
                for j in gene_pair1:
                    file.write(str(tre_ID)+'\t'+str(gd_num)+'\t')
                    a,b=j.split('-')
                    if a.split('_')[0] ==b.split('_')[0]:
                        c=a.split('_')[0]
                    else:
                        if a=='null':
                            c=b.split('_')[0]
                        else:
                            c=a.split('_')[0]
        
                    file.write(str(clade.support)+'\t'+new_named_gene2gene_dic.get(a, null)+'\t'+new_named_gene2gene_dic.get(b, null)+'\t'+voucher2taxa_dic.get(parent,parent)+'\t'+voucher2taxa_dic[c]+'\t'+str(dup_sps)+'\t'+str(round(dup_percent,2))+'\t'+model+'\t'+'-'+'\t''\n')
                gd_num+=1
            
        merged_data=count_elements_in_lists(gd_type_dic)
        df = pd.DataFrame.from_dict(merged_data, orient='index').fillna(0)

        df.to_csv('gd_type.csv')

def count_elements_in_lists(data):
    def merge_and_filter_types(data, merge_map):
        merged_data = {}
        for key, counter in data.items():
            new_counter = Counter()
            for t, count in counter.items():
                merged = False
                for new_type, types_to_merge in merge_map.items():
                    if t in types_to_merge:
                        new_counter[new_type] += count
                        merged = True
                        break
                if not merged and t == 'AB<=>AB':
                    new_counter['ABAB'] += count
            merged_data[key] = new_counter
        return merged_data

    counted_data = {key: Counter(value) for key, value in data.items()}
    merge_map = {
    'ABB': ['AB<=>B', 'B<=>AB', 'XB<=>AB', 'AB<=>XB'],
    'AAB': ['AB<=>A', 'A<=>AB', 'AX<=>AB', 'AB<=>AX']
    }

    result = merge_and_filter_types(counted_data, merge_map)
    return result

def get_model(clade, sptree):
    sps = get_species_list(clade)
    sps_clade = sptree.get_common_ancestor(set(sps))
    
    sps_clade_a = sps_clade.get_children()[0] if sps_clade.get_children() else None
    sps_clade_b = sps_clade.get_children()[1] if sps_clade.get_children() and len(sps_clade.get_children()) > 1 else None
    
    if sps_clade_a.is_leaf():
        sps_clade_a.add_feature('label', 'Aa')
    else:
        sps_clade_a_1 = sps_clade_a.get_children()[0] if sps_clade_a.get_children() else None
        sps_clade_a_2 = sps_clade_a.get_children()[1] if sps_clade_a.get_children() and len(sps_clade_a.get_children()) > 1 else None
        
        for leaf in sps_clade_a_1:
            leaf.add_feature('label', 'A')
        for leaf in sps_clade_a_2:
            leaf.add_feature('label', 'a')

    if sps_clade_b.is_leaf():
        sps_clade_b.add_feature('label', 'Bb')
    else:
        sps_clade_b_1 = sps_clade_b.get_children()[0] if sps_clade_b.get_children() else None
        sps_clade_b_2 = sps_clade_b.get_children()[1] if sps_clade_b.get_children() and len(sps_clade_b.get_children()) > 1 else None
        
        for leaf in sps_clade_b_1:
            leaf.add_feature('label', 'B')
        for leaf in sps_clade_b_2:
            leaf.add_feature('label', 'b')

    for j in clade:
        species = j.name.split('_')[0]
        clade1 = sps_clade & species
        if clade1:
            j.add_feature('label', clade1.label)

    up_clade = ''
    for j in clade.get_children()[0]:
        up_clade += j.label
    up_clade = up_clade + '<=>'
    for j in clade.get_children()[1]:
        up_clade += j.label
    clade_up = set(up_clade.split('<=>')[0])
    clade_down = set(up_clade.split('<=>')[1])
    clade_up_1 = ''.join(clade_up)
    clade_up_1_1 = ''.join(sorted(clade_up_1, key=lambda x: (x.lower(), x.isupper())))

    clade_down_1 = ''.join(clade_down)
    clade_down_1_1 = ''.join(sorted(clade_down_1, key=lambda x: (x.lower(), x.isupper())))
    clade_model = clade_up_1_1 + '<=>' + clade_down_1_1
    
    def process_string(s):
        result = []
        i = 0
        while i < len(s):
            if i < len(s) - 1 and ((s[i] == 'A' and s[i+1] == 'a') or (s[i] == 'a' and s[i+1] == 'A')):
                result.append('A')
                i += 2
            elif i < len(s) - 1 and ((s[i] == 'B' and s[i+1] == 'b') or (s[i] == 'b' and s[i+1] == 'B')):
                result.append('B')
                i += 2
            else:
                if s[i] in ['A', 'B', 'a', 'b']:
                    result.append('X')
                else:
                    result.append(s[i])
                i += 1

        return ''.join(result)
    return process_string(clade_model)


def get_model_dic(interspecies_node_list,genetree,sptree):
    model_dic={}
    for i in interspecies_node_list:
        clade=genetree&i
        s=get_model(clade,sptree)
        model_dic.setdefault(s, []).append(i)
    return model_dic

def get_empty_count_dict(sptree: object) -> dict:
    empty_count_dic = {node.name: 0 for node in sptree.traverse()}
    return empty_count_dic

def judge_support(support,support_value):
    if  support <=1 and 0.5 <=support_value <=1:
        if support>=support_value:
            return True
        else:
            return False
        
    elif support <=1 and 50 <= support_value <=100:
        support_value=support_value/100
        if support>=support_value:
            return True
        else:
            return False
    elif support > 1 and 0.5 <=support_value <=1:
        support_value=support_value*100
        if support>=support_value:
            return True
        else:
            return False
    elif support > 1 and 50 <=support_value <=100:
        if support>=support_value:
            return True
        else:
            return False

def gene_pair(clade):
    result_pairs = set()

    children = clade.get_children()
    child1, child2 = sorted(children, key=lambda c: len(c.get_leaf_names()), reverse=True)
    
    leaves1, leaves2 = child1.get_leaf_names(), child2.get_leaf_names()

    for tip1 in leaves1:
        matching_tips = [tip2 for tip2 in leaves2 if tip1.split('_')[0] == tip2.split('_')[0]]
        
        if matching_tips:
            result_pairs.update(f"{tip1}-{tip2}" for tip2 in matching_tips)
        else:
            result_pairs.add(f"{tip1}-null")

    for tip2 in leaves2:
        if all(tip2.split('_')[0] != tip1.split('_')[0] for tip1 in leaves1):
            result_pairs.add(f"null-{tip2}")

    return result_pairs

if __name__ == "__main__":
    support=50
    dup_species_percent = 0.5
    dup_species_num = 2
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap.txt")
    sptree=PhyloTree('30sptree.nwk')
    sptree=rename_species_tree(sptree, voucher2taxa_dic)
    num_tre_node(sptree)
    tre_dic=read_and_return_dict('GF.txt')
    filename = 'result.txt'
    sp_dic=[]
    write_gd_result(sp_dic,filename, tre_dic, support,dup_species_percent, dup_species_num,sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic)
   
