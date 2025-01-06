#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#Usage: python dotplot.py AHC.gff AAS.gff AHC.lens AAS.lens GD_pairs.txt
#Usage: 20240310, Yiyogn Zhao modified

#AHC.gff
#for gff, the last column is the gene count for each chromosome.
#chr1   AHC_0004600 12837   26777   +   1
#chr1   AHC_0004073 33171   35791   +   2
#chr1   AHC_0004946 46794   47258   -   3

#AHC.lens
#chr1   23021215
#chr2   18778401
#chr3   19321588

#AAS.gff
#for gff, the last column is the gene count for each chromosome.
#Chr1   AAS_0012495 3631    5899    +   1
#Chr1   AAS_0011322 5928    8737    -   2
#Chr1   AAS_0007791 11649   13714   -   3

#AAS.lens
#Chr1   30425192
#Chr2   19696821
#Chr3   23459804

#GD_pairs.txt
#if first iput is AHC.gff, so the GD pairs first column is AHC genes; 
#AHC_0000001    AAS_0011214 NONE
#AHC_0000001    AAS_0011214 NONE
#AHC_0000001    AAS_0011214 NONE
from __init__ import *
import matplotlib 
matplotlib.use('Agg') 
import sys 
import re 
import numpy as np 
import matplotlib.pyplot as plt; plt.rcdefaults() 
import matplotlib.patches as mpatches 
import time 
import gc
from collections import defaultdict,Counter


def read_gff(fn):
    f=open(fn) 
    data,dict=[],{} 
    for line in f.readlines(): 
        a=line.strip().split("\t") 
        dict[a[1]]=[a[0],a[2],a[3],a[4],a[5]]
        data.append(a) 
    return data,dict

def read_lens(fn, chrs=None): 
    data = []
    
    chrs_lst = None
    if chrs:
        with open(chrs, 'r') as f:
            chrs_lst = {i.strip() for i in f.readlines()} 
    
    with open(fn, 'r') as fp:
        for row in fp.readlines():
            r1, r2 = row.strip().split('\t')           
            if chrs_lst is None or r1 in chrs_lst:
                data.append([r1, r2])
    
    return data

def plot_chr1(lens,gl,gl2,mark,name): 
    total_lens=sum([float(k[1]) for k in lens]) 
    step=gl/float(total_lens) 
    gl_start,n,start_x=0.95,0,0.05 
    mark_y=0.04 
    align = dict(family='Times New Roman',style='normal',horizontalalignment ="center", verticalalignment="center") 
    for k in lens: 
        n+=float(k[1]) 
        mark_new=str(mark)+str(k[0]) 
        x=gl_start-float(n)*step 
        mark_x=x+0.5*float(k[1])*step 
        plt.plot([start_x,start_x+gl2],[x,x],linestyle = '-',color='black', linewidth=0.5) 
        plt.text(mark_y,mark_x,mark_new,color='black',fontsize = 12,rotation = 90,weight='semibold',**align) 
    plt.plot([start_x,start_x+gl2],[gl_start,gl_start],linestyle='-', color='black',linewidth=1) 
    plt.text(mark_y-0.02,0.5*(2*gl_start-gl),name,color='black',fontsize = 18,rotation = 90,weight='semibold',**align) 
    return step

def plot_chr2(lens,gl,gl2,mark,name): 
    total_lens=sum([float(k[1]) for k in lens]) 
    step=gl/float(total_lens) 
    gl_start,n,start_x=0.05,0,0.95 
    mark_y=0.96 
    align = dict(family='Times New Roman',style='normal',horizontalalignment="center", verticalalignment="center") 
    for k in lens: 
        n+=float(k[1]) 
        mark_new=str(mark)+str(k[0]) 
        x=gl_start+float(n)*step 
        mark_x=x-0.5*float(k[1])*step 
        plt.plot([x,x],[start_x,start_x-gl2],linestyle='-',color='black', linewidth=0.5) 
        plt.text(mark_x,mark_y,mark_new,color='black',fontsize = 12,rotation = 0,weight='semibold',**align) 
    plt.plot([gl_start,gl_start],[start_x,start_x-gl2],linestyle='-', color='black',linewidth=1) 
    plt.text(0.5*(2*gl_start+gl),mark_y+0.02,name,color='black',fontsize = 18,rotation = 0,weight='semibold',**align) 
    return step

def gene_location(gff, lens, step): 
    loc_gene,dict_chr,n={},{},0 
    for i in lens: 
        dict_chr[i[0]]=n 
        n+=float(i[1]) 
    for k in gff: 
        if k[0] not in dict_chr.keys(): 
            continue 
        loc=(float(dict_chr[k[0]])+float(k[5]))*step 
        loc_gene[k[1]]=loc 
    return loc_gene

def read_gd_pairs(gd_lst):
    dict_gd,dict_gd1={},{}
    for line in gd_lst:
        a=line.strip().split("\t")
        dict_gd1[str(a[0])+":"+str(a[1])]=str(a[2])
    for pair, wgd_notch in dict_gd1.items():
        fd=pair.strip().split(":")
        gene1=fd[0]
        gene2=fd[1]
        index= gene1+":"+gene2
        reverse_index=gene2+":"+ gene1
        dict_gd[index]=dict_gd1[index]
        if reverse_index not in dict_gd1.keys():
            dict_gd[reverse_index]=dict_gd1[index]
    return dict_gd
    
def plot_dot(root, data, loc1, loc2, gl,dict_gd,size=0.001):
    gl_start1, gl_start2 = 0.95, 0.05
    
    for pair, wgd_notch in dict_gd.items():
        fd = pair.strip().split(":")
        gene1 = fd[0]
        gene2 = fd[1]
        index = gene1 + ":" + gene2
        if gene1 in loc1.keys() and gene2 in loc2.keys():
            x, y = loc1[gene1], loc2[gene2]
            x, y = gl_start1 - x, gl_start2 + y
            if dict_gd[index] == "NONE":
                DrawCircle(root, [y, x], size, 'gray', 0.6)

    for pair, wgd_notch in dict_gd.items():
        fd = pair.strip().split(":")
        gene1 = fd[0]
        gene2 = fd[1]
        index = gene1 + ":" + gene2
        if gene1 in loc1.keys() and gene2 in loc2.keys():
            x, y = loc1[gene1], loc2[gene2]
            x, y = gl_start1 - x, gl_start2 + y
            if dict_gd[index] == "red":
                DrawCircle(root, [y, x], size, 'red', 0.6)
            elif dict_gd[index] == "green":
                DrawCircle(root, [y, x], size, 'green', 0.6)
            elif dict_gd[index] == "purple":
                DrawCircle(root, [y, x], size, 'purple', 0.6)
            elif dict_gd[index] == "blue":
                DrawCircle(root, [y, x], size, 'blue', 0.6)   

def DrawCircle(ax, loc, radius,color, alpha):
    circle = mpatches.Circle(loc, radius, edgecolor="none", facecolor=color, alpha=alpha)
    ax.add_patch(circle)


def generate_dotplot(gff1,gff2,lens1,lens2,gd_pairs,spe1,spe2,file_name,target_chr1=None, target_chr2=None, size=None):
    plt.figure(figsize=(10, 10)) 
    root = plt.axes([0, 0, 1, 1]) 
    align = dict(family='Arial',style='normal',horizontalalignment="center", verticalalignment="center") 
    t1=time.time() 
    print(f"Dotplot of {file_name} are ready to begin") 
    gff_1,dict_gff1=read_gff(gff1) 
    gff_2,dict_gff2=read_gff(gff2) 
    t2=time.time() 
    print("Reading gff took "+str(t2-t1)+" second") 
    if target_chr1 and target_chr2:
        lens_1=read_lens(lens1,chrs1)
        lens_2=read_lens(lens2,chrs2)
    else:
        lens_1=read_lens(lens1)
        lens_2=read_lens(lens2)
    t3=time.time() 
    print("Reading lens took "+str(t3-t2)+" second") 
    gl1,gl2=0.92,0.92 
    step_1=plot_chr1(lens_1,gl1,gl2,'',spe1) 
    step_2=plot_chr2(lens_2,gl2,gl1,'',spe2)
    dict_gd=read_gd_pairs(gd_pairs)
    t4=time.time() 
    print("Reading lebel_pairs took "+str(t4-t3)+" second")  
    gene_loc_1=gene_location(gff_1,lens_1,step_1)
    gene_loc_2=gene_location(gff_2,lens_2,step_2) 


    gene_conversion_list=find_gene_conversion(dict_gd, dict_gff1, dict_gff2, lens_1, lens_2,file_name)
    find_gene_pair_info(gene_conversion_list, dict_gd, dict_gff1, dict_gff2,file_name)
    t5=time.time() 
    print("Dealing lebel_pairs took "+str(t5-t4)+" second") 
    gc.collect() 

    if size:
        plot_dot(root,dict_gd,gene_loc_1,gene_loc_2,gl1,dict_gd,size) 
    else:
        plot_dot(root,dict_gd,gene_loc_1,gene_loc_2,gl1,dict_gd) 
    t6=time.time() 
    print("Ploting dot took "+str(t6-t5)+" second") 
    root.set_xlim(0, 1) 
    root.set_ylim(0, 1) 
    root.set_axis_off() 
    plt.savefig(file_name+"_dotplot.pdf",dpi=500) 
    plt.savefig(file_name+"_dotplot.png",dpi=500) 
    t7=time.time() 
    print(f"{file_name} Dotplot totaly took "+str(t7-t1)+" second")
#####################################
def assign_colors_by_alignment(alignments, alignment_scores):
    gene_to_alignments = defaultdict(list)  

   
    for alignment, gene_pairs in alignments.items():
        for gene_a, gene_b in gene_pairs:
            gene_to_alignments[gene_a].append((alignment, gene_b))

    result = set()
    gene_pair_set=set()
    
    for gene_a, alignments_list in gene_to_alignments.items():
        
        alignments_list.sort(key=lambda x: alignment_scores[x[0]], reverse=True)

       
        for i, (alignment, gene_b) in enumerate(alignments_list):
            if i == 0:
                color = "red"
            else:
                color = "blue"
            gene_pair=f"{gene_a}\t{gene_b}"
            if gene_pair not in gene_pair_set:
                result.add(f"{gene_pair}\t{color}")
                gene_pair_set.add(gene_pair)
            else:
                continue

    sorted_list = sorted(result)
    return sorted_list

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

#####################################
def find_dup_node(Phylo_t:object,sp1,sp2,new_named_gene2gene_dic,processed_lines,written_results,support_value=50)->list:
    """
    Find all duplication nodes in the phylogenetic tree.
    
    Args:
        Phylo_t (object): The phylogenetic tree object.
    
    Returns:
        list: A list of duplication node names.
    """
    events = Phylo_t.get_descendant_evol_events()
    dup_node_list = []
    for event in events:
        if event.etype == "D":
            event_nodes = set(event.in_seqs) | set(event.out_seqs)
            dup_node = Phylo_t.get_common_ancestor(list(event_nodes))
            if judge_support(dup_node.support,support_value):
                child1,child2=dup_node.get_children()

                if judge_support(child1.support,support_value) and judge_support(child2.support,support_value):
                    dup_node_list.append(dup_node)
        else:
            event_nodes = set(event.in_seqs) | set(event.out_seqs)
            dup_node = Phylo_t.get_common_ancestor(list(event_nodes))
            orthology=get_gene_pairs(dup_node,sp1,sp2)
            for i in orthology:
                gene1, gene2 = i
                gene_a=new_named_gene2gene_dic[gene1]
                gene_b=new_named_gene2gene_dic[gene2]
                result = f"{gene_a}\t{gene_b}\tred\n"
                s = f"{gene2}-r"
                if s not in written_results:
                    processed_lines.append(result)
                    written_results.add(s)
    return dup_node_list


def generate_combinations(list1, list2):
    gene_pairs=[]
    for elem1 in list1:
        for elem2 in list2:
            gene_pairs.append((elem1,elem2))
    return gene_pairs

def get_gene_pairs(t,sp1,sp2):
    sp_count = defaultdict(list)
    for i in t.get_leaf_names():
        sps = i.split('_')[0]
        sp_count[sps].append(i)

    sp1_gene_list = list(sp_count[sp1])
    sp2_gene_list = list(sp_count[sp2])
    gene_pairs = generate_combinations(sp1_gene_list, sp2_gene_list)
    return gene_pairs

def find_independent_dup_nodes(dup_node_list, Phylo_t: object) -> list:
    dup_node_list.sort(key=lambda x: len(x.get_leaf_names()))
    node_dic={}

    for node in dup_node_list:
        if len(get_species_set(node))==1:
            continue
        else:
            tips=len(node)
            if tips in node_dic:
                node_dic[tips].append(node)
            else:
                node_dic[tips]=[node]
    if node_dic:
        min_key_value = min(node_dic.items(), key=lambda x: x[0])

        return min_key_value[1]       
    else:
        return []

def find_gene_conversion(dict_gd, dict_gff1, dict_gff2, lens_1, lens_2,gd_pairs):
    chrs_combinations = [(chr_a, chr_b) for chr_a in [i[0] for i in lens_1] 
                         for chr_b in [i[0] for i in lens_2]]

    block_list = []

    for chr_a, chr_b in chrs_combinations:
        
        if chr_a==chr_b: 
            red_count = 0
            blue_count = 0

            for pair, color in dict_gd.items():
                gene_a, gene_b = pair.split(":")
                chr_gene_a = dict_gff1.get(gene_a, [None])[0]
                chr_gene_b = dict_gff2.get(gene_b, [None])[0]

                if chr_gene_a == chr_a and chr_gene_b == chr_b:
                    if color == "red":
                        red_count += 1
                    elif color == "blue":
                        blue_count += 1

            total = red_count + blue_count
            red_ratio = red_count / total if total > 0 else 0
            blue_ratio = blue_count / total if total > 0 else 0

            block_list.append({
                "chr_a": chr_a,
                "chr_b": chr_b,
                "red_count": red_count,
                "blue_count": blue_count,
                "red_ratio": red_ratio,
                "blue_ratio": blue_ratio
            })

    filter_block_list = [
        block for block in block_list
        if block["blue_ratio"] > 0 and block["red_ratio"] / block["blue_ratio"] < 10
    ]
    with open(f'chr_conversion_list_{gd_pairs}.txt', 'w') as f:
        for i in filter_block_list:
            f.write(f'{i}\n')
    return filter_block_list


def find_gene_pair_info(gene_conversion_list, dict_gd, dict_gff1, dict_gff2,gd_pairs):
    with open(f'gene_conversion_{gd_pairs}.txt', 'w') as f:
        for block in gene_conversion_list:
            chr_a = block["chr_a"]
            chr_b = block["chr_b"]

            for pair, color in dict_gd.items():
                gene_a, gene_b = pair.split(":")
                chr_gene_a = dict_gff1.get(gene_a, [None])[0]
                chr_gene_b = dict_gff2.get(gene_b, [None])[0]

                if chr_gene_a == chr_a and chr_gene_b == chr_b:
                    f.write(f'{gene_a}\t{gene_b}\t{color}\n')

def process_gd_result(gf,imap,sp1,sp2,support):
    tre_dic=read_and_return_dict(gf)
    gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(imap)
    rename_sp1=taxa2voucher_dic[sp1]
    rename_sp2=taxa2voucher_dic[sp2]
    processed_lines=[]
    written_results = set()
    for tre_ID, tre_path in tre_dic.items():
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t0.resolve_polytomy(recursive=True)
        Phylo_t0.sort_descendants()
        Phylo_t1 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)

        if len(get_species_set(Phylo_t1)) == 1:
            continue
        sps_tol=get_species_set(Phylo_t1)

        if len(sps_tol)==2 and {sp1,sp2} not in sps_tol:
            continue

        dup_node_list = find_dup_node(Phylo_t1,rename_sp1, rename_sp2,new_named_gene2gene_dic,processed_lines,written_results,support)
        

        
        for i in dup_node_list:
            if len(i) ==3:
                gene_pairs = get_gene_pairs(i, rename_sp1, rename_sp2)
                
                gd_clade1, gd_clade2 = i.get_children()
                gd_tips1 = set(gd_clade1.get_leaf_names())
                gd_tips2 = set(gd_clade2.get_leaf_names())

                for item in gene_pairs:
                    gene1, gene2 = item
                    item_set = set(item)
                    gene_a=new_named_gene2gene_dic[gene1]
                    gene_b=new_named_gene2gene_dic[gene2]
                    
                    if item_set <= gd_tips1 or item_set <= gd_tips2:
                        result = f"{gene_a}\t{gene_b}\tred\n"
                        s = f"{gene2}-r"

                    else:
                        result = f"{gene_a}\t{gene_b}\tblue\n"
                        s = f"{gene2}-b"


                    if s not in written_results:
                        processed_lines.append(result)
                        written_results.add(s)

            elif len(i) ==4:
                gd_clade1, gd_clade2 = i.get_children()
                tips1 = set(get_species_set(gd_clade1))
                tips2 = set(get_species_set(gd_clade2))
                if len(tips1)==len(tips2)==2:
                    gene_pairs = get_gene_pairs(i, rename_sp1, rename_sp2)
                    gd_tips1 = set(gd_clade1.get_leaf_names())
                    gd_tips2 = set(gd_clade2.get_leaf_names())
                    for item in gene_pairs:
                        gene1, gene2 = item
                        item_set = set(item)
                        gene_a=new_named_gene2gene_dic[gene1]
                        gene_b=new_named_gene2gene_dic[gene2]
                        
                        if item_set <= gd_tips1 or item_set <= gd_tips2:
                            result = f"{gene_a}\t{gene_b}\tred\n"
                            s = f"{gene2}-r"

                        else:
                            result = f"{gene_a}\t{gene_b}\tblue\n"
                            s = f"{gene2}-b"


                        if s not in written_results:
                            processed_lines.append(result)
                            written_results.add(s)

    sorted_lines = sorted(processed_lines, key=lambda x: x.split('\t')[0])
    return sorted_lines

if __name__=="__main__":
    gff1=sys.argv[1]
    gff2=sys.argv[2]
    lens1=sys.argv[3]
    lens2=sys.argv[4]
    blastp_pairs=sys.argv[5]
    synteny_pairs=sys.argv[6]
    spe1=sys.argv[7]
    spe2=sys.argv[8]
    num=int(sys.argv[9])
    gf=sys.argv[10]
    imap=sys.argv[11]
    target_chr1=sys.argv[12]
    target_chr2=sys.argv[13]
    size=sys.argv[14]
    

    process_blastp_pairs=process_blastp_result(blastp_pairs, num)
    alignments, alignment_scores = parse_synteny_file(synteny_pairs)
    process_synteny_pairs = assign_colors_by_alignment(alignments, alignment_scores)
    process_gd_pairs=process_gd_result(gf,imap,spe1,spe2)
    generate_dotplot(gff1,gff2,lens1,lens2,process_blastp_pairs,spe1,spe2,'blastp_pairs',target_chr1,target_chr2,size)
    print('-'*30)
    generate_dotplot(gff1,gff2,lens1,lens2,process_synteny_pairs,spe1,spe2,'synteny_pairs',target_chr1,target_chr2,size)
    print('-'*30)
    generate_dotplot(gff1,gff2,lens1,lens2,process_gd_pairs,spe1,spe2,'gd_pairs',target_chr1,target_chr2,size)
    total_pairs=[process_blastp_pairs,process_synteny_pairs,process_gd_pairs]
    total_lst=process_total_color_list(total_pairs)
    generate_dotplot(gff1,gff2,lens1,lens2,total_lst,spe1,spe2,'total_pairs',target_chr1,target_chr2,size)
