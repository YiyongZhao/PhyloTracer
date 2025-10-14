from signal import SIGUSR2
from __init__ import *
from collections import defaultdict
from Bio import SeqIO
import phyde as hd
import subprocess
from collections import Counter
import time


def process_start_node(file_path: str, sptree: object) -> str:
    """
    Reads species names from a file and returns the name of their common ancestor node in the species tree.

    Args:
        file_path (str): Path to a file containing species names, one per line.
        sptree (object): The species tree object (an ete3.Tree or PhyloTree object).

    Returns:
        str: The name of the common ancestor node.
    """
    species_list = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                species_list.append(line.strip())
    except FileNotFoundError:
        print(f"Error: Species list file not found at {file_path}")
        return None
    except Exception as e:
        print(f"Error reading species list file: {e}")
        return None

    if not species_list:
        print("Warning: Species list file is empty.")
        return None

    try:
        # Find the common ancestor of the species in the list
        common_ancestor = sptree.get_common_ancestor(species_list)
        return common_ancestor.name
    except Exception as e:
        print(f"Error finding common ancestor in species tree: {e}")
        return None

def create_fasta_dict(fasta_file,gene2new_named_gene_dic):
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[gene2new_named_gene_dic[record.id]] = str(record.seq)
    return fasta_dict

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


def run_hyde_from_matrix(all_genes:dict,voucher2taxa_dic:dict,clade:object):
    target_sps=clade.get_leaf_names()
    concat_matrix,max_length = build_concatenated_matrix(all_genes)
    
    imap={}
    with open("temp.imap", 'w') as imap_file:
        for spsecies, seq in concat_matrix.items():
            sp = spsecies.split('_')[0]
            sp_0 = voucher2taxa_dic[sp]
            if sp not in target_sps:
                imap[sp_0]='out'
                imap_file.write(f'{sp_0}\tout\n')
            else:
                imap[sp_0]=sp_0
                imap_file.write(f'{sp_0}\t{sp_0}\n')

    with open('temp.phy', 'w') as phy_file:
        for species, sequence in concat_matrix.items():
            phy_file.write(f'{voucher2taxa_dic[species]}\t{sequence}\n')

    hyde_result_lst=[]
    sps_num=len(imap.keys())
    taxa_num=len(set(imap.values()))
    dat = hd.HydeData("temp.phy", "temp.imap", 'out', sps_num, taxa_num, max_length,quiet=True)
    res = dat.list_triples()
    for t in res:
        p1,hyb,p2=t
        result=dat.test_triple(p1,hyb,p2)
        combined_element = (t, result,len(res))
        hyde_result_lst.append(combined_element)

    #os.remove('temp.phy')
    #os.remove('temp.imap')
    return hyde_result_lst


            
def hyde_main(tre_dic, seq_path_dic, rename_sptree, gene2new_named_gene_dic, voucher2taxa_dic, taxa2voucher_dic,new_named_gene2gene_dic,target_node=None,gd_group:int=1):

    hyde_tuple_lst=[]
    gd_count_dic,gd_type_dic=get_gd_count_dic_and_gd_type_dic(tre_dic,gene2new_named_gene_dic,rename_sptree)
    data=count_elements_in_lists(gd_type_dic)
    gd_clades=get_process_gd_clade(data,gd_count_dic)

    sps_dic={}
    for gd in gd_clades:
        gd_name,gds=gd
        if target_node and gd_name != target_node:
            continue

        start_time = time.time()
        print(f'{gd_name} is processing')
        # Display gene family processing statistics
        print(f"Total gene duplication events in dataset: {len(gds)}")
        print(f"Number of parallel processing groups: {gd_group}")
        
        split_gds = np.array_split(gds, gd_group)
        
        # Display parallelization partitioning results
        print(f"Gene duplication events partitioned into {len(split_gds)} processing batches")
        for i, sub_group in enumerate(split_gds):
            print(f"Batch {i+1}: {len(sub_group)} gene duplication events")
        for sub_group in split_gds:
            hyde_result_lst=process_gd_group(sub_group, rename_sptree, gd_name, seq_path_dic, gene2new_named_gene_dic,new_named_gene2gene_dic, voucher2taxa_dic,sps_dic,target_node)
            hyde_tuple_lst.extend(hyde_result_lst)
    header = ("P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\t""ABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n")
    with open('hyde_out.txt', 'w') as out_file, open('hyde_filtered_out.txt', 'w') as filtered_out_file:
        print(header, end="", file=out_file)
        print(header, end="", file=filtered_out_file)

        for t, res, t_num in hyde_tuple_lst:
            write_out(res, t, out_file)
            if is_filtered(res, t_num):
                write_out(res, t, filtered_out_file)

def process_gd_group(gds, rename_sptree, gd_name, seq_path_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic,sps_dic,target_node):
    """
    处理一组GD数据
    :param gds: GD数据列表
    :param rename_sptree: 重命名的物种树
    :param gd_name: GD名称
    :param seq_path_dic: 序列路径字典
    :param gene2new_named_gene_dic: 基因重命名字典
    :param voucher2taxa_dic: 凭证到分类的映射
    :return: 返回hyde处理结果
    """
    target_clade=rename_sptree&target_node
    all_matrices = []  # 收集所有的矩阵
    col_counter = 1
    seq_dic_all={}
    for gd_clade_set in gds:
        gdid = gd_clade_set[0]
        gd_clade = gd_clade_set[1]

        outfile = gdid
        tre_id1 = outfile.split('-')[0]
        
        try:
            seq_dic = create_fasta_dict(seq_path_dic[tre_id1], gene2new_named_gene_dic)
            seq_dic_all.update(seq_dic)
        except KeyError as e:
            print(f"KeyError encountered for GD {gdid}: {e}. Skipping this GD.")
            continue 

        outgroup_gene = get_outgroup_gene(gd_clade, rename_sptree)
        if outgroup_gene is None:
            continue
        if voucher2taxa_dic[outgroup_gene.split('_')[0]]=='MED':

            one_gd_matrix=process_one_gd(gd_clade, outgroup_gene)
            num_cols = one_gd_matrix.shape[1]
            new_col_names = [f"col{col_counter + i}" for i in range(num_cols)]
            one_gd_matrix.columns = new_col_names

            # 更新列计数器
            col_counter += num_cols

            all_matrices.append(one_gd_matrix)
        # print(one_gd_matrix)
    
    # Check if all_matrices is empty before concatenating
    if not all_matrices:
        print(f"Warning: No valid matrices found for GD group {gd_name}. Returning empty result.")
        return []
        
    merged_matrix = pd.concat(all_matrices, axis=1)
    merged_matrix = merged_matrix.fillna('-')

    clean_matrix=clean_matrix_by_dash_count(merged_matrix)
    hyde_result_lst=run_hyde_from_matrix_integrated(clean_matrix, seq_dic_all, voucher2taxa_dic, target_clade)
    return hyde_result_lst

def run_hyde_from_matrix_integrated(clean_matrix: pd.DataFrame, seq_dic: dict, voucher2taxa_dic: dict, clade: object, trim: bool = True):
    """
    整合 matrix_to_phy 和 HyDe 分析的完整流程
    
    参数:
        clean_matrix: DataFrame，清洗后的基因矩阵
        seq_dic: dict，基因ID到序列的映射
        voucher2taxa_dic: dict，voucher到taxa的映射
        clade: object，目标分支
        trim: bool，是否裁剪序列
    """
    target_sps = clade.get_leaf_names()
    
    # 1. 使用 matrix_to_phy 生成 PHYLIP 文件
    matrix_to_phy(clean_matrix, seq_dic, voucher2taxa_dic, "temp.phy", trim=trim)
    
    # 2. 创建 temp.imap 文件
    imap = {}
    with open("temp.imap", 'w') as imap_file:
        # 从 clean_matrix 的索引（物种名）创建 imap
        for species in clean_matrix.index:
            sp_0 = voucher2taxa_dic[species]
            if species not in target_sps:
                imap[sp_0] = 'out'
                imap_file.write(f'{sp_0}\tout\n')
            else:
                imap[sp_0] = sp_0
                imap_file.write(f'{sp_0}\t{sp_0}\n')
    
    # 3. 读取生成的 PHYLIP 文件获取序列长度
    with open("temp.phy", 'r') as f:
        first_line = f.readline().strip()
        num_species, max_length = map(int, first_line.split())
    
    # 4. 运行 HyDe 分析
    hyde_result_lst = []
    sps_num = len(imap.keys())
    taxa_num = len(set(imap.values()))
    
    dat = hd.HydeData("temp.phy", "temp.imap", 'out', sps_num, taxa_num, max_length, quiet=True)
    res = dat.list_triples()
    
    for t in res:
        p1, hyb, p2 = t
        result = dat.test_triple(p1, hyb, p2)
        combined_element = (t, result, len(res))
        hyde_result_lst.append(combined_element)
    
    # 5. 清理临时文件
    os.remove('temp.phy')
    os.remove('temp.imap')
    
    return hyde_result_lst

def clean_matrix_by_dash_count(matrix: pd.DataFrame, max_allowed_dashes: int = 1):
    """
    清洗矩阵：
    1. 删除重复列；
    2. 删除含有 `'-'` 数量大于 max_allowed_dashes 的列；
    3. 重命名列为 col1, col2, ...；
    4. 保存为 CSV。

    参数:
        matrix: 输入 DataFrame
        max_allowed_dashes: 允许列中最多出现几次 '-'（默认 0 表示完全不允许）
        output_file: 输出文件路径
    """


    # 1. 删除重复列
    matrix = matrix.T.drop_duplicates().T
    # 2. 删除含有超过 max_allowed_dashes 个 '-' 的列
    dash_counts = (matrix == '-').sum(axis=0)
    matrix = matrix.loc[:, dash_counts <= max_allowed_dashes]

    # 3. 重命名列
    matrix.columns = [f"col{i+1}" for i in range(matrix.shape[1])]

    return matrix

def matrix_to_phy(clean_matrix: pd.DataFrame, seq_dic: dict, voucher2taxa_dic,output_file: str = "output.phy", trim: bool = True):
    """
    将清洗后的基因矩阵转换为 PHYLIP 格式，并可选择是否裁剪全为 '-' 的列。

    参数:
        clean_matrix: DataFrame，行是物种，列是不同 gene 编号（如 col1, col2...）
        seq_dic: dict，gene_id -> 蛋白序列
        output_file: str，输出的 .phy 文件路径
        trim: bool，是否裁剪拼接序列中全是 '-' 的列
    """
    # 1. 构建每个物种的完整拼接序列
    species_seqs = defaultdict(str)

    for col in clean_matrix.columns:
        for species, gene_id in clean_matrix[col].items():
            sp=voucher2taxa_dic[species]
            if gene_id != '-' and gene_id in seq_dic:
                seq = seq_dic[gene_id]
            else:
                # 缺失或找不到就填充一个合适长度的 '-'
                seq = '-' * max(len(s) for s in seq_dic.values())
            species_seqs[sp] += seq

    # 2. 可选裁剪：去除全是 '-' 的列
    if trim:
        def trimal_matrix(concat_matrix: dict):
            """
            去除全为 '-' 的列。自动修复拼接长度不一致的问题。

            参数:
                concat_matrix: dict，物种名 -> 拼接序列

            返回:
                filtered_matrix: dict，去除空列后的新矩阵
                new_len: int，保留的列数
            """
            # 获取最长序列长度
            max_len = max(len(seq) for seq in concat_matrix.values())

            # 填补短序列为相同长度（右侧补 '-')
            for sp in concat_matrix:
                seq = concat_matrix[sp]
                if len(seq) < max_len:
                    concat_matrix[sp] = seq + '-' * (max_len - len(seq))

            # 找到不是全 '-' 的列索引
            valid_cols = [i for i in range(max_len) if any(concat_matrix[sp][i] != '-' for sp in concat_matrix)]

            # 构建新矩阵
            filtered_matrix = {
                sp: ''.join(concat_matrix[sp][i] for i in valid_cols)
                for sp in concat_matrix
            }

            return filtered_matrix, len(valid_cols)


        species_seqs, trimmed_len = trimal_matrix(species_seqs)
    else:
        trimmed_len = len(next(iter(species_seqs.values())))

    # 3. 写入 PHYLIP 文件
    num_species = len(species_seqs)
    with open(output_file, "w") as f:
        f.write(f"{num_species} {trimmed_len}\n")
        for species, seq in species_seqs.items():
            name = species[:10].ljust(10)
            f.write(f"{name} {seq}\n")

def process_one_gd(gd_clade:object, outgroup_gene:str):
    result = []
    def is_dup_node(node): #改判断重复有问题 应该是上下分枝有重复的物种才能算作重复 不能简单的gene和species比对
        genes = get_species_list(node)
        species = get_species_set(node)
        return len(genes) != len(species)

    def traverse(node):
        if is_dup_node(node):
            for child in getattr(node, 'children', []):
                traverse(child)
        else:
            genes = node.get_leaf_names() if hasattr(node, 'get_leaf_names') else []
            tuple_genes = tuple(sorted(genes + [outgroup_gene], key=lambda x: x.split('_')[0]))
            result.append(tuple_genes)
            
    traverse(gd_clade)
    
    # 将结果转换为矩阵格式
    if not result:
        return pd.DataFrame()
    
    # 获取所有唯一的物种
    all_species = set()
    for triple in result:
        for gene in triple:
            species = gene.split('_')[0]
            all_species.add(species)
    
    all_species = sorted(list(all_species))
    
    # 创建矩阵，行为物种，列为每个triple
    matrix_data = []
    for i, triple in enumerate(result):
        column_data = []
        for species in all_species:
            # 找到该物种在当前triple中的基因
            genes_in_species = [gene for gene in triple if gene.split('_')[0] == species]
            if genes_in_species:
                column_data.append(';'.join(genes_in_species))  # 如果有多个基因，用分号分隔
            else:
                column_data.append('-')  # 如果该物种不在当前triple中
        matrix_data.append(column_data)
    
    # 转置矩阵，使得行为物种，列为triple
    matrix_df = pd.DataFrame(matrix_data).T
    matrix_df.index = all_species
    matrix_df.columns = [f'col_{i+1}' for i in range(len(result))]
    
    return matrix_df

def is_filtered(res, t_num):
    """Check if the result passes filtering criteria."""
    return (
        res["Pvalue"] < (0.05 / t_num)
        and abs(res["Zscore"]) != 99999.9
        and 0.0 < res["Gamma"] < 1.0
    )
    
















