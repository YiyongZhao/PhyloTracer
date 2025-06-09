from __init__ import *
from collections import defaultdict
from Bio import SeqIO
import phyde as hd
import subprocess
from collections import Counter
import time

def format_time(seconds):
    days = seconds // (24 * 3600)
    hours = (seconds % (24 * 3600)) // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return f"{int(days)} d, {int(hours)} h, {int(minutes)} m, {seconds:.2f} s"

def calculate_depth(node1, node2):
    if node1 in node2.iter_ancestors() or node2 in node1.iter_ancestors():
        distance = node1.get_distance(node2, topology_only=True)+2
        return distance
    common_ancestor = node1.get_common_ancestor(node2)
    return abs(node1.get_distance(common_ancestor, topology_only=True) - \
           node2.get_distance(common_ancestor, topology_only=True))

def calculate_distance(node1,node2):
    distance = node1.get_distance(node2)
    return distance

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

def create_fasta_dict(fasta_file,gene2new_named_gene_dic):
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[gene2new_named_gene_dic[record.id]] = str(record.seq)
    return fasta_dict

def get_seq(genes,seq_dic):
    gene_seqs=[]
    for gene in genes :
        seq=seq_dic[gene]
        gene_seqs.append(gene+'\t'+seq)
    return gene_seqs

def get_outgroup_gene(gd_clade, sptree):
    def get_outsps_sort_lst(map_t, sptree):
        sps_depth = {}
        for gene in sptree:
            if gene.name in get_species_list(map_t):
                continue
            sps_depth[gene.name] = gene.get_distance(map_t, topology_only=True)
        sorted_out = sorted(sps_depth.items(), key=lambda item: item[1])
        return [sps for sps, _ in sorted_out]
        
    def get_outsps_min_distance_gene(node, sps):
        dic = {i.name: node.get_distance(i) for i in node if sps in i.name}
        if not dic:
            return None
        return min(dic, key=dic.get)
    
    def get_sister_nodes(node):
        sister_nodes = []
        current_node = node
        while current_node and not current_node.is_root():
            sisters = current_node.get_sisters()
            if sisters:
                sister_nodes.append(sisters[0])
            current_node = current_node.up
        return sister_nodes
    
    def get_target_gene(sisters, sps):
        for sis_node in sisters:
            if sps in get_species_set(sis_node):
                return get_outsps_min_distance_gene(sis_node, sps)
        return None
                
    map_t = sptree.get_common_ancestor(get_species_set(gd_clade))
    sisters = get_sister_nodes(gd_clade)
    sps_lst = get_outsps_sort_lst(map_t, sptree)

    for sps in sps_lst:
        target = get_target_gene(sisters, sps)
        if target:
            return target
    return None
                
    map_t = sptree.get_common_ancestor(get_species_set(node))
    sis = get_sister_nodes(node)

    sps_lst = get_outsps_sort_lst(map_t, sptree)

    for sps in sps_lst:
        target=get_target_gene(sis,sps)
        if target:
            return target


def get_max_deepth(root:object)->int:
    if not root:
        return 0
    
    max_child_depth = 0
    for child in root.children:
        child_depth = get_max_deepth(child)
        max_child_depth = max(max_child_depth, child_depth)
    
    return max_child_depth + 1

def run_hyde_for_subclades(outgroup_gene,node, seq_dic,voucher2taxa_dic):
    def run_hyde_for_clade(leafs, seq_dic, voucher2taxa_dic, outgroup_gene):
        leafs.append(outgroup_gene)  # 仅加一次 outgroup gene
        sps_seq = {}

        with open("temp.imap", 'w') as imap_file:
            for leaf in leafs:
                sps = leaf.split('_')[0]
                sps_0 = voucher2taxa_dic[sps]
                if sps not in sps_seq:
                    imap_file.write(f'{sps_0}\t{sps_0}\n')
                    sps_seq[sps] = seq_dic[leaf]
                else:
                    sps_seq[sps] += seq_dic[leaf]

        max_length = max(len(seq) for seq in sps_seq.values())
        padded_sps_seq = {k: v + '-' * (max_length - len(v)) for k, v in sps_seq.items()}

        with open('temp.phy', 'w') as phy_file:
            for species, sequence in padded_sps_seq.items():
                phy_file.write(f'{voucher2taxa_dic[species]}\t{sequence}\n')

        outgroup_taxa = voucher2taxa_dic[outgroup_gene.split('_')[0]]
        hyde_result_lst = []
        sps_num = len(sps_seq)

        dat = hd.HydeData("temp.phy", "temp.imap", outgroup_taxa, sps_num, sps_num, max_length)
        res = dat.list_triples()
        for t in res:
            p1, hyb, p2 = t
            result = dat.test_triple(p1, hyb, p2)
            combined_element = (t, result, len(res))
            hyde_result_lst.append(combined_element)

        os.remove('temp.phy')
        os.remove('temp.imap')

        return hyde_result_lst

    children = node.children
    if len(children) != 2:
        raise ValueError("The input node must have exactly two children.")

    subclade_results = []
    for child in children:
        leafs = child.get_leaf_names()
        result = run_hyde_for_clade(leafs, seq_dic, voucher2taxa_dic, outgroup_gene)
        subclade_results.extend(result)

    return subclade_results

def run_hyde_for_gdclade(outgroup_gene,node, seq_dic,voucher2taxa_dic):
    leafs = node.get_leaf_names()
    leafs.extend([outgroup_gene, outgroup_gene])
    sps_seq = {}

    with open("temp.imap", 'w') as imap_file:
        for leaf in leafs:
            sps = leaf.split('_')[0]
            sps_0 = voucher2taxa_dic[sps]
            if sps not in sps_seq:
                imap_file.write(f'{sps_0}\t{sps_0}\n')
                sps_seq[sps] = seq_dic[leaf]
            else:
                sps_seq[sps] += seq_dic[leaf]

    max_length = max(len(seq) for seq in sps_seq.values())
    padded_sps_seq = {k: v + '-' * (max_length - len(v)) for k, v in sps_seq.items()}


    with open('temp.phy', 'w') as phy_file:
        for species, sequence in padded_sps_seq.items():
            phy_file.write(f'{voucher2taxa_dic[species]}\t{sequence}\n')
    
    outgroup_taxa = voucher2taxa_dic[outgroup_gene.split('_')[0]]
    

    hyde_result_lst=[]
    sps_num=len(sps_seq)
    dat = hd.HydeData("temp.phy", "temp.imap", outgroup_taxa, sps_num, sps_num, max_length)
    res = dat.list_triples()
    for t in res:
        p1,hyb,p2=t
        result=dat.test_triple(p1,hyb,p2)
        combined_element = (t, result,len(res))
        hyde_result_lst.append(combined_element)
    os.remove('temp.phy')
    os.remove('temp.imap')
    return hyde_result_lst

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

def get_gd_count_dic_and_gd_type_dic(tre_dic,gene2new_named_gene_dic,rename_sptree):
    gd_num = 0
    gd_count_dic={}
    gd_type_dic={}
    for tre_id, tre_path in tre_dic.items():
        t = read_phylo_tree(tre_path)
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        
        gds = find_dup_node(t1, rename_sptree, 50)

        for gd in gds:
            type_str=get_model(gd,rename_sptree)
            if gd.map in gd_count_dic:
                gd_count_dic[gd.map].append((tre_id+'-'+str(gd_num),gd))
                gd_type_dic[gd.map].append(type_str)
            else:
                gd_count_dic[gd.map]=[(tre_id+'-'+str(gd_num),gd)]
                gd_type_dic[gd.map]=[type_str]
            gd_num+=1
    return gd_count_dic,gd_type_dic

def get_process_gd_clade(gd_type_dic,gd_count_dic):
    gd_clade=[]
    for k,v in gd_type_dic.items():
        if (v['ABB']+v['AAB']+v['ABAB'])!=0:
            type_ratio=(v['ABB']+v['AAB'])/(v['ABB']+v['AAB']+v['ABAB'])
            if (v['ABB']+v['AAB']+v['ABAB'])>=10 and type_ratio>=0.1 :
                gd_clade.append((k,gd_count_dic[k]))
    return gd_clade

def write_out(out, triple, outfile):
    """
    Take the output from test_triple() and write it to file.
    """
    print(triple[0],"\t",triple[1],"\t",triple[2],"\t",sep="",end="",file=outfile,)
    print(out["Zscore"], "\t", sep="", end="", file=outfile)
    print(out["Pvalue"], "\t", sep="", end="", file=outfile)
    print(out["Gamma"], "\t", sep="", end="", file=outfile)
    print(out["AAAA"], "\t", sep="", end="", file=outfile)
    print(out["AAAB"], "\t", sep="", end="", file=outfile)
    print(out["AABA"], "\t", sep="", end="", file=outfile)
    print(out["AABB"], "\t", sep="", end="", file=outfile)
    print(out["AABC"], "\t", sep="", end="", file=outfile)
    print(out["ABAA"], "\t", sep="", end="", file=outfile)
    print(out["ABAB"], "\t", sep="", end="", file=outfile)
    print(out["ABAC"], "\t", sep="", end="", file=outfile)
    print(out["ABBA"], "\t", sep="", end="", file=outfile)
    print(out["BAAA"], "\t", sep="", end="", file=outfile)
    print(out["ABBC"], "\t", sep="", end="", file=outfile)
    print(out["CABC"], "\t", sep="", end="", file=outfile)
    print(out["BACA"], "\t", sep="", end="", file=outfile)
    print(out["BCAA"], "\t", sep="", end="", file=outfile)
    print(out["ABCD"], "\n", sep="", end="", file=outfile)



def build_concatenated_matrix(all_genes):
    # 统计所有基因的顺序和每个基因的长度
    gene_order = list(all_genes.keys())
    gene_len = {gene: len(seq) for gene, seq in all_genes.items()}
    # 获取所有物种
    species_set = set([gene.split('_')[0] for gene in all_genes.keys()])
    # 初始化物种到各基因序列的映射
    species2seq = {sp: [] for sp in species_set}
    for gene in gene_order:
        sp = gene.split('_')[0]
        seq = all_genes[gene]
        for s in species_set:
            if s == sp:
                species2seq[s].append(seq)
            else:
                species2seq[s].append('-' * gene_len[gene])
    # 拼接每个物种的所有基因序列
    concat_matrix = {sp: ''.join(seq_list) for sp, seq_list in species2seq.items()}
    seq_length = max(len(seq) for seq in concat_matrix.values())
    return concat_matrix, seq_length

def run_hyde_from_matrix(all_genes:dict,voucher2taxa_dic:dict,outgroup:str):
    concat_matrix,max_length = build_concatenated_matrix(all_genes)
    with open("temp.imap", 'w') as imap_file:
        for spsecies, seq in concat_matrix.items():
            sp = spsecies.split('_')[0]
            sp_0 = voucher2taxa_dic[sp]
            imap_file.write(f'{sp_0}\t{sp_0}\n')

    with open('temp.phy', 'w') as phy_file:
        for species, sequence in concat_matrix.items():
            phy_file.write(f'{voucher2taxa_dic[species]}\t{sequence}\n')

    hyde_result_lst=[]
    sps_num=len(concat_matrix)
    dat = hd.HydeData("temp.phy", "temp.imap", outgroup, sps_num, sps_num, max_length)
    res = dat.list_triples()
    for t in res:
        p1,hyb,p2=t
        result=dat.test_triple(p1,hyb,p2)
        combined_element = (t, result,len(res))
        hyde_result_lst.append(combined_element)

    header = ("P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\t""ABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n")
    with open("hyde_out.txt", 'w') as out_file, open('hyde_filtered_out.txt', 'w') as filtered_out_file:
        print(header, end="", file=out_file)
        print(header, end="", file=filtered_out_file)

        for t, res, t_num in hyde_result_lst:
            write_out(res, t, out_file)
            if is_filtered(res, t_num):
                write_out(res, t, filtered_out_file)

    # os.remove('temp.phy')
    # os.remove('temp.imap')
    return hyde_result_lst

def hyde_main(tre_dic, seq_path_dic, rename_sptree, gene2new_named_gene_dic, voucher2taxa_dic, taxa2voucher_dic,new_named_gene2gene_dic,target_node=None,split_gd:bool=False):
    hyde_tuple_lst=[]
    all_genes = {}
    all_seq_dic = {}
    gd_lst=[]
    gd_group=5
    
    for tre_id, tre_path in tre_dic.items():
        t = read_phylo_tree(tre_path)
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        seq_dic = create_fasta_dict(seq_path_dic[tre_id], gene2new_named_gene_dic)
        all_seq_dic.update(seq_dic)
        gds = find_dup_node(t1, rename_sptree, 50)
        gd_lst.extend(gds)

    for gd in gd_lst:

        for gene in gd.get_leaf_names():
            all_genes[gene] = all_seq_dic[gene]
  
    vo='CEC'
    hyde_result_lst=run_hyde_from_matrix(all_genes,voucher2taxa_dic,vo)
    hyde_tuple_lst.extend(hyde_result_lst)
    
    
    # group_hyde_result_lst=[]
    # for i in range(0, len(gd_lst), gd_group):
    #     group = gd_lst[i:i+gd_group]
        
    #     group_genes = {}
    #     sps=set()
    #     out=[]
    #     for gd in group:
    #         outgroup_gene = get_outgroup_gene(gd, rename_sptree)
    #         out.append(outgroup_gene.split('_')[0] if outgroup_gene else None)
    #         for gene in gd.get_leaf_names():
    #             group_genes[gene] = all_seq_dic[gene]
    #             sps.add(gene.split('_')[0] )
                
    #     out_filtered = [voucher2taxa_dic.get(x,'') for x in out if x is not None]
    #     if 'Cercis_chinensis' in out_filtered:
    #         vo='Cercis_chinensis'
    #         if vo in group_genes:
               
    #             hyde_result_lst=run_hyde_from_matrix(group_genes,voucher2taxa_dic,vo)
    #             group_hyde_result_lst.extend(hyde_result_lst)
    #     else:
    #         print('no outgroup')
    #         continue
            
            

        
        








def hyde_main1(tre_dic, seq_path_dic, rename_sptree, gene2new_named_gene_dic, voucher2taxa_dic, taxa2voucher_dic,new_named_gene2gene_dic,target_node=None,split_gd:bool=False):
    hyde_tuple_lst=[]
    gd_count_dic,gd_type_dic=get_gd_count_dic_and_gd_type_dic(tre_dic,gene2new_named_gene_dic,rename_sptree)
    data=count_elements_in_lists(gd_type_dic)
    gd_clades=get_process_gd_clade(data,gd_count_dic)

    for gd in gd_clades:
        gd_name,gds=gd


        if target_node and gd_name != target_node:
            continue

        start_time = time.time()
        print(f'{gd_name} is processing')

        clade=rename_sptree&gd_name

        for gd_clade_set in gds:
            gdid=gd_clade_set[0]
            gd_clade=gd_clade_set[1]
            outfile =gdid
            tre_id1=outfile.split('-')[0]
            seq_dic = create_fasta_dict(seq_path_dic[tre_id1], gene2new_named_gene_dic)
            outgroup_gene = get_outgroup_gene(gd_clade, rename_sptree)
            if outgroup_gene!=None:
                if split_gd:
                    tup = run_hyde_for_subclades(outgroup_gene,gd_clade, seq_dic,voucher2taxa_dic)
                else:
                    tup = run_hyde_for_gdclade(outgroup_gene,gd_clade, seq_dic,voucher2taxa_dic)
                hyde_tuple_lst.extend(tup)
            else:
                print(f'{gdid} not finde outgroup_gene')

        output_files = {"normal": ("hyde_out.txt", "hyde_filtered_out.txt"),"split": ("split_hyde_out.txt", "split_hyde_filtered_out.txt")}
        file_key = "split" if split_gd else "normal"
        out_file_name, filtered_out_file_name = output_files[file_key]
        header = ("P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\t""ABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n")
        with open(out_file_name, 'w') as out_file, open(filtered_out_file_name, 'w') as filtered_out_file:
            print(header, end="", file=out_file)
            print(header, end="", file=filtered_out_file)

            for t, res, t_num in hyde_tuple_lst:
                write_out(res, t, out_file)
                if is_filtered(res, t_num):
                    write_out(res, t, filtered_out_file)

def is_filtered(res, t_num):
    """Check if the result passes filtering criteria."""
    return (
        res["Pvalue"] < (0.05 / t_num)
        and abs(res["Zscore"]) != 99999.9
        and 0.0 < res["Gamma"] < 1.0
    )


















