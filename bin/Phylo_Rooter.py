from Ortho_Retriever import *
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick
import math
import pandas as pd
import os
import shutil
from tqdm import tqdm

# --------------------------
# 1. 辅助函数区
# --------------------------

def get_species_map_and_depth(species_tree):
    """预计算物种树节点的深度，用于快速查询。"""
    species_depth = {}
    for node in species_tree.traverse("preorder"):
        if node.is_root():
            depth = 0
        else:
            depth = species_depth[node.up] + 1
        species_depth[node] = depth
    return species_depth

def annotate_mapped_depths(gene_tree: object, species_tree: object) -> object:
    """
    将基因树节点映射到物种树，并标记深度。
    使用缓存机制加速。
    """
    if not gene_tree or not species_tree:
        return gene_tree
    
    # 获取预计算的物种树深度字典
    sp_depths = get_species_map_and_depth(species_tree)
    
    cache = {}
    for node in gene_tree.traverse():
        species_set = get_species_set(node)
        if not species_set:
            node.add_feature('mapped_depth', 0)
            continue
            
        key = frozenset(species_set)
        if key in cache:
            node.add_feature('mapped_depth', cache[key])
            continue
            
        try:
            if len(species_set) == 1:
                species_name = list(species_set)[0]
                mapped_node = species_tree & species_name
            else:
                mapped_node = species_tree.get_common_ancestor(species_set)
            
            # 直接查字典获取深度，比 while 循环快
            depth = sp_depths.get(mapped_node, 0)
            
            cache[key] = depth
            node.add_feature('mapped_depth', depth)
        except Exception:
            node.add_feature('mapped_depth', 0)
    return gene_tree

def get_species_tree_basal_set(species_tree_obj: object) -> set:
    """获取物种树的基部类群集合。"""
    if len(species_tree_obj.get_children()) < 2:
        return set(species_tree_obj.get_leaf_names())
    children = species_tree_obj.get_children()
    child_0_leaves = set(children[0].get_leaf_names())
    child_1_leaves = set(children[1].get_leaf_names())
    return child_0_leaves if len(child_0_leaves) < len(child_1_leaves) else child_1_leaves

def get_dynamic_basal_set(gene_tree_species_set: set, species_tree_root: object) -> set:
    """
    根据基因树中存在的物种，在物种树上动态寻找“相对基部类群”。
    
    Args:
        gene_tree_species_set: 当前基因树包含的所有物种集合 (set of strings)
        species_tree_root: 物种树的根节点 (ETE3 TreeNode)
        
    Returns:
        set: 适用于当前基因树的有效基部物种集合。
    """
    # 递归终止条件：到达叶子节点
    if species_tree_root.is_leaf():
        return set([species_tree_root.name])
    
    children = species_tree_root.get_children()
    
    # 假设物种树是二叉的（如果有三叉，逻辑需微调，这里按标准二叉树处理）
    # 如果有多个子节点，我们需要判断基因树跨越了哪次分化
    
    # 获取物种树当前节点两个子分支下的所有物种名
    # 注意：这里需要缓存优化，否则大树会慢。但物种树通常不大，直接遍历也可。
    clade_sets = []
    for child in children:
        clade_sets.append(set(child.get_leaf_names()))
        
    # 检查基因树的物种落在哪些分支里
    present_flags = [] # 记录每个分支里是否有基因树的物种
    for c_set in clade_sets:
        # 如果基因树物种与该分支有交集，记为 True
        if not gene_tree_species_set.isdisjoint(c_set):
            present_flags.append(True)
        else:
            present_flags.append(False)
            
    # --- 核心判断逻辑 ---
    
    # 情况 A：基因树的物种跨越了当前的分歧点（即至少在两个子分支里都有分布）
    # 例如：既有单子叶，又有双子叶。
    if sum(present_flags) >= 2:
        # 此时，这就是我们要找的最早分歧点。
        # 我们选择叶子数量较少的那个分支作为“相对外群”。
        # (通常 Outgroup 分支物种数 < Ingroup 分支)
        
        # 找到包含基因树物种的所有分支，并按物种树上的叶子总数排序
        valid_children = []
        for i, has_gene in enumerate(present_flags):
            if has_gene:
                valid_children.append((children[i], len(clade_sets[i])))
        
        # 按该分支在物种树上的大小排序（从小到大）
        valid_children.sort(key=lambda x: x[1])
        
        # 返回最小分支的所有物种作为 Basal Set
        # 这意味着：如果 Os(单子叶) 和 AT(双子叶) 都在，且单子叶分支小，则 Os 是基部。
        best_outgroup_node = valid_children[0][0]
        return set(best_outgroup_node.get_leaf_names())

    # 情况 B：基因树的物种完全落在一个子分支里（例如 Amborella 丢了，都在另一个分支里）
    # 此时我们需要“深入”那个分支继续找
    elif sum(present_flags) == 1:
        for i, has_gene in enumerate(present_flags):
            if has_gene:
                return get_dynamic_basal_set(gene_tree_species_set, children[i])
                
    # 情况 C：基因树物种都不在（理论不应发生），返回空
    return set()

def get_all_rerooted_trees_filtered(tree: object, basal_species_set: set) -> list:
    """基部类群过滤 + 候选树生成。"""
    rerooted_trees = []
    # 确保节点有唯一名字以便查找
    num_tre_node(tree)

    for node in tree.traverse("preorder"):
        if node.is_root(): continue
        # if node.is_leaf(): continue # 通常我们不在叶子上定根，除非是单系
        
        # 获取当前分支下的物种（去除基因名后缀）
        node_species = set([leaf.split('_')[0] for leaf in node.get_leaf_names()])
        
        # 过滤：如果当前分支不包含任何基部物种，跳过
        if not node_species.issubset(basal_species_set):
            continue
            
        tree_copy = tree.copy()
        try:
            target_node = tree_copy.search_nodes(name=node.name)[0]
            tree_copy.set_outgroup(target_node)
            tree_copy.name = node.name
            rerooted_trees.append(tree_copy)
        except Exception:
            continue
    return rerooted_trees

def calculate_tree_statistics(tree: object, renamed_species_tree: object) -> tuple:
    """
    Stage 1 快速统计 (不含 RF)。
    """
    up_clade = tree.children[1]
    down_clade = tree.children[0]
    
    var_up = compute_tip_to_root_branch_length_variance(up_clade)
    var_down = compute_tip_to_root_branch_length_variance(down_clade)
    var = abs(var_up - var_down)
    
    # 使用预计算的 mapped_depth
    if len(up_clade.get_leaf_names()) > len(down_clade.get_leaf_names()):
        deep = getattr(down_clade, 'mapped_depth', 0)
    else:
        deep = getattr(up_clade, 'mapped_depth', 0)
        
    species_overlap, GD = calculate_species_overlap_gd_num(tree)
    return deep, var, GD, species_overlap

def calculate_species_overlap_gd_num(gene_tree: object) -> float:
    dup_nodes = find_dup_node_simple(gene_tree)
    if not dup_nodes:
        return 0.0, 0
    largest_tree = max(dup_nodes, key=lambda node: len(node.get_leaves()))
    up_clade = largest_tree.children[1]
    down_clade = largest_tree.children[0]
    species_list_a = get_species_list(up_clade)
    species_list_b = get_species_list(down_clade)
    
    union_len = len(set(species_list_a) | set(species_list_b))
    if union_len == 0: return 0.0, len(dup_nodes)
    
    overlap_ratio = len(set(species_list_a) & set(species_list_b)) / union_len
    return overlap_ratio, len(dup_nodes)

# --------------------------
# 2. 核心 RF 计算逻辑 (集成用户代码)
# --------------------------

def calculate_rf_strategy(
    tree: object, 
    renamed_species_tree: object,
    renamed_length_dict: dict,
    gene_to_new_name: dict, 
    new_name_to_gene: dict,
    tree_path: str, 
    tree_id: str
) -> float:
    """
    根据单拷贝/多拷贝选择不同的 RF 计算策略。
    """
    species_list = get_species_list(tree)
    species_set = get_species_set(tree)
    
    # --- 情况 A: 单拷贝 ---
    if len(species_list) == len(species_set):
        return calculate_RF_distance(tree, renamed_species_tree)
    
    # --- 情况 B: 多拷贝 (集成用户逻辑) ---
    else:
        # 1. 拆分 Principal Gene Set 和 Offcuts
        principal_gene_set, filtered_offcut_ev_seqs = offcut_tre(tree, renamed_length_dict)
        minor_orthologs = []
        minor_orthologs = iterator(filtered_offcut_ev_seqs, tree, gene_to_new_name, minor_orthologs, tree_path, renamed_length_dict)
        ordered_name_OG_list = rename_OGs_tre_name(principal_gene_set, minor_orthologs, tree_id)
        RF = 0
        for tre_name, OG_set in ordered_name_OG_list:
            OG_list = [new_name_to_gene[OG] for OG in OG_set]
            phylo_tree_0 = read_phylo_tree(tree_path)
            phylo_tree = root_tre_with_midpoint_outgroup(phylo_tree_0)
            phylo_tree_OG_list = extract_tree(OG_list, phylo_tree)
            phylo_tree_OG_list = rename_input_tre(phylo_tree_OG_list, gene_to_new_name)
            RF += calculate_RF_distance(phylo_tree_OG_list, renamed_species_tree)
            
        return RF

def calculate_RF_distance(Phylo_t_OG_L: object, sptree: object) -> int:
    tcopy = Phylo_t_OG_L.copy()
    for leaf in tcopy:
        if "_" in leaf.name:
            leaf.name = leaf.name.split('_')[0]
    try:
        rf, _ = tcopy.robinson_foulds(sptree)[:2]
        return rf
    except:
        return 100# 错误惩罚

# --------------------------
# 3. 评分与归一化
# --------------------------

def normalize_and_score(df, weights, include_rf=True):
    # 辅助 lambda: 防止分母为 0
    norm = lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0
    
    # deep/var/GD 越小越好；overlap 越大越好
    norm_deep = norm(df["deep"])
    norm_var = norm(df["var"])
    norm_GD = norm(df["GD"])
    norm_species_overlap = norm(df["species_overlap"]) # Overlap 越大越好

    weighted_norm_RF = 0
    if include_rf and "RF" in df.columns:
        norm_RF = norm(df["RF"])
        weighted_norm_RF = norm_RF * weights.get("RF", 0)

    weighted_norm_overlap = norm_species_overlap * weights.get("species_overlap", 0)
    score = (
        norm_deep * weights.get("deep", 0) +
        norm_var * weights.get("var", 0) +
        norm_GD * weights.get("GD", 0) -
        weighted_norm_overlap +
        weighted_norm_RF
    )
    return score

# --------------------------
# 4. 主流程
# --------------------------

def root_main(
    tree_dict: dict,
    gene_to_new_name: dict,
    renamed_length_dict: dict,
    new_name_to_gene: dict,
    renamed_species_tree: object,
    voucher_to_taxa: dict
) -> None:
    dir_path = os.path.join(os.getcwd(), "rooted_trees/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    
    # 预处理：获取基部类群
    basal_species_set = get_species_tree_basal_set(renamed_species_tree)
    print(f"Basal Filter: {len(basal_species_set)} species")

    pbar = tqdm(total=len(tree_dict), desc="Processing trees", unit="tree")
    stat_matrix = []
    
    stage1_weights_multi = {"deep": 0.3, "var": 0.1, "RF": 0.0, "GD": 0.5, "species_overlap": 0.1}
    stage1_weights_single = {"deep": 0.7, "var": 0.3, "RF": 0.0}

    try:
        for tree_id, tree_path in tree_dict.items():
            pbar.set_description(f"Processing {tree_id}")
            
            # 1. 读取与初步处理
            Phylo_t0 = read_phylo_tree(tree_path)
            Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
            Phylo_t2 = rename_input_tre(Phylo_t1, gene_to_new_name)

            # 2. 简单情况处理
            if len(get_species_set(Phylo_t2)) == 1 or len(get_species_list(Phylo_t2)) <= 3:
                rename_output_tre(Phylo_t2, new_name_to_gene, tree_id, dir_path)
                pbar.update(1)
                continue
            
            # 3. 深度注释 & 基部过滤
            Phylo_t2 = annotate_mapped_depths(Phylo_t2, renamed_species_tree) 

            current_gene_species = set(get_species_list(Phylo_t2))
            dynamic_basal_set = get_dynamic_basal_set(current_gene_species, renamed_species_tree)
            root_list = get_all_rerooted_trees_filtered(Phylo_t2, dynamic_basal_set)
            
            if not root_list:
                print("No valid root candidates after basal filter.")
                rename_output_tre(Phylo_t2, new_name_to_gene, tree_id, dir_path)
                pbar.update(1)
                continue
            
            # ==========================================
            # Stage 1: 快速筛选 (无 RF)
            # ==========================================
            tree_objects = {} 
            temp_stats = []
            
            for n, tree in enumerate(root_list):
                # tree 已经是 rename 过的 Taxa ID 格式
                tree_key = f"{tree_id}_{n+1}"
                tree_objects[tree_key] = tree 
                # print(rename_input_tre(tree, new_name_to_gene))
                deep, var, GD, species_overlap = calculate_tree_statistics(tree, renamed_species_tree)
                # print(deep, var, GD, species_overlap)
                temp_stats.append({
                    "Tree": tree_key,
                    "deep": deep, "var": var, "GD": GD, 
                    "species_overlap": species_overlap, "RF": 0,
                    "tree_obj_ref": tree # 暂存对象引用，Stage 2 用
                })

            current_df = pd.DataFrame(temp_stats)
            
            # 权重选择
            is_multi_copy = len(get_species_list(Phylo_t2)) != len(get_species_set(Phylo_t2))
            used_weights = stage1_weights_multi if is_multi_copy else stage1_weights_single
            
            # 初步评分
            current_df["score"] = normalize_and_score(current_df, used_weights, include_rf=False)
            
            # ==========================================
            # Stage 2: 精细筛选 (Top N 算复杂 RF)
            # ==========================================
            total_candidates = len(current_df)
            # 策略：取前 40% 或前 20 名 (保持之前的宽松漏斗，防止漏掉 RF 好的树)
            top_n = max(20, math.ceil(total_candidates * 0.8)) 
            top_n = min(top_n, total_candidates)
            
            # 1. 选出初试成绩最好的 Top N
            top_candidates_df = current_df.nsmallest(top_n, "score").copy()
            # print('-'*50)
            # for _, row in top_candidates_df.iterrows():
            #     tk = row["Tree"]
            #     t = tree_objects[tk]
            #     t_print = t.copy('newick')
            #     for nd in t_print.traverse():
            #         if nd.name in new_name_to_gene:
            #             nd.name = new_name_to_gene[nd.name]
            #     print(t_print)
            
            # 2. 对 Top N 计算 RF (面试)
            for idx, row in top_candidates_df.iterrows():
                target_tree = row["tree_obj_ref"]
                
                # 计算 RF
                rf_val = calculate_rf_strategy(
                    tree=target_tree,
                    renamed_species_tree=renamed_species_tree,
                    renamed_length_dict=renamed_length_dict,
                    gene_to_new_name=gene_to_new_name,
                    new_name_to_gene=new_name_to_gene,
                    tree_path=tree_path,
                    tree_id=tree_id
                )
                
                top_candidates_df.at[idx, "RF"] = rf_val
            
            # 3. 【核心修改】直接排序选最优，不再重新加权
            # 逻辑：优先看 RF (越小越好)，RF 如果一样，看 Stage 1 的 Score (越小越好)
            best_row = top_candidates_df.sort_values(
                by=["RF", "score"], 
                ascending=[True, True]
            ).iloc[0]
            
            # ==========================================
            # 输出
            # ==========================================
            best_tree_key = best_row["Tree"]
            best_tree = tree_objects[best_tree_key]
            
            rename_output_tre(best_tree, new_name_to_gene, tree_id, dir_path)
            
           
            record = best_row.to_dict()
            del record["tree_obj_ref"]
            stat_matrix.append(record)
            
            pbar.update(1)

        # 保存报表
        stat_df = pd.DataFrame(stat_matrix)
        if not stat_df.empty:
            stat_df["tree_id"] = stat_df["Tree"].apply(lambda x: "_".join(x.split("_")[:-1]))
            cols = ["Tree", "score", "deep", "var", "GD", "species_overlap", "RF"]
            stat_df = stat_df.sort_values(by=["tree_id", "score"])[cols]
            stat_df.to_csv("stat_matrix.csv", index=False)
            
    finally:
        pbar.close()

def rename_output_tre(tree: object, name_mapping: dict, tree_id: str, output_dir: str) -> None:
    # 还原名字
    for node in tree.traverse():
        if node.name in name_mapping:
            node.name = name_mapping[node.name]
    tree_str = tree.write(format=0)
    write_tree_to_newick(tree_str, tree_id, output_dir)
