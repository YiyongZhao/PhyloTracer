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
from collections import defaultdict


def read_gff(fn):
    f=open(fn) 
    data,dict=[],{} 
    for line in f.readlines(): 
        a=line.strip().split("\t") 
        dict[a[1]]=[a[0],a[2],a[3],a[4],a[5]]
        data.append(a) 
    return data,dict
def read_lens(fn): 
    fp=open(fn) 
    data=[] 
    for row in fp.readlines(): 
        r1,r2=row.split() 
        data.append([r1,r2]) 
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
        loc=(float(dict_chr[k[0]])+float(k[3]))*step 
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
    
def plot_dot(root,data,loc1,loc2,gl,dict_gd):
    gl_start1,gl_start2=0.95,0.05
    for pair, wgd_notch in dict_gd.items():
        fd=pair.strip().split(":")
        gene1=fd[0]
        gene2=fd[1]
        index= gene1+":"+gene2
        if gene1 in loc1.keys() and gene2 in loc2.keys():
            x,y=loc1[gene1],loc2[gene2]
            x,y=gl_start1-x,gl_start2+y
            if index in dict_gd.keys():
                if dict_gd[index] == "red":
                    DrawCircle(root,[y,x],0.001,'red',0.6)
                if dict_gd[index] == "green":
                    DrawCircle(root,[y,x],0.001,'green',0.6)
                if dict_gd[index] == "orange":
                    DrawCircle(root,[y,x],0.001,'orange',0.6)
                if dict_gd[index] == "blue":
                    DrawCircle(root,[y,x],0.001,'blue',0.6)
                if dict_gd[index] == "NONE":
                    DrawCircle(root,[y,x],0.001,'gray',0.6)   

def DrawCircle(ax, loc, radius,color, alpha):
    circle = mpatches.Circle(loc, radius, edgecolor="none", facecolor=color, alpha=alpha)
    ax.add_patch(circle)


def generate_dotplot(gff1,gff2,lens1,lens2,gd_pairs,spe1,spe2,file_name):
    plt.figure(figsize=(10, 10)) 
    root = plt.axes([0, 0, 1, 1]) 
    align = dict(family='Arial',style='normal',horizontalalignment="center", verticalalignment="center") 
    t1=time.time() 
    print(f"Dotplot of {file_name} are ready to begin") 
    gff_1,dict_gff1=read_gff(gff1) 
    gff_2,dict_gff2=read_gff(gff2) 
    t2=time.time() 
    print("Reading gff took "+str(t2-t1)+" second") 
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
    t5=time.time() 
    print("Dealing lebel_pairs took "+str(t5-t4)+" second") 
    gc.collect() 
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
def process_blastp_result(input_file, num: int = 5):
    col_count = {}
    processed_lines = []
    
    with open(input_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            col1 = cols[0]
            
            if col1 in col_count:
                if col_count[col1] < num:
                    processed_lines.append(f"{cols[0]}\t{cols[1]}\tNONE")
                    col_count[col1] += 1
            else:
                processed_lines.append(f"{cols[0]}\t{cols[1]}\torange")
                col_count[col1] = 1

    sorted_lines = sorted(processed_lines, key=lambda x: x.split('\t')[2])
    return sorted_lines

#####################################
def process_synteny_result(input_file):
    with open(input_file, 'r') as file:
        processed_lines = [f"{cols[0]}\t{cols[2]}\tblue\n" for line in file if not line.startswith('#') for cols in [line.strip().split('\t')]]
    sorted_lines = sorted(processed_lines, key=lambda x: x.split('\t')[2])
    
    return sorted_lines

#####################################
def find_dup_node(Phylo_t:object)->list:#After searching for duplication events, the list of node names where duplication events occurred is as follows:
    events = Phylo_t.get_descendant_evol_events()
    dup_node_name_list = []
    for ev in events:
        if ev.etype == "D":
            i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
            events_node_name_list = i.split(',')
            common_ancestor_node_name = Phylo_t.get_common_ancestor(events_node_name_list)
            dup_node_name_list.append(common_ancestor_node_name)
    return dup_node_name_list


def filter_gd(gene_pair, paral_clade_list, new_named_gene2gene_dic):
    item1_set = set(gene_pair)

    for gd in paral_clade_list:
        gd_clade1, gd_clade2 = gd.get_children()
        gd_tips1 = [leaf for leaf in gd_clade1.get_leaf_names()]
        gd_tips2 = [leaf for leaf in gd_clade2.get_leaf_names()]
        set1 = set(gd_tips1)
        set2 = set(gd_tips2)
        total = set(gd_tips1 + gd_tips2)

        if item1_set <= total:
            if item1_set <= set1 or item1_set <= set2:
                return new_named_gene2gene_dic[gene_pair[0]] + '\t' + new_named_gene2gene_dic[gene_pair[1]] + '\t' + 'green'+'\n'
            else:
                return new_named_gene2gene_dic[gene_pair[0]] + '\t' + new_named_gene2gene_dic[gene_pair[1]] + '\t' + 'red'+'\n'
        else:
            return new_named_gene2gene_dic[gene_pair[0]] + '\t' + new_named_gene2gene_dic[gene_pair[1]] + '\t' + 'NONE'+'\n'  # Return empty string if no match found

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

def process_gd_result(gf,imap,sp1,sp2):
    tre_dic=read_and_return_dict(gf)
    gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic= gene_id_transfer(imap)
    rename_sp1=taxa2voucher_dic[sp1]
    rename_sp2=taxa2voucher_dic[sp2]
    processed_lines=[]
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

        paral_clade_list = find_dup_node(Phylo_t1)
        gene_pairs = get_gene_pairs(Phylo_t1,rename_sp1,rename_sp2)

        
        for item in gene_pairs:
            result = filter_gd(item, paral_clade_list, new_named_gene2gene_dic)
            if result:
                processed_lines.append(result)
    sorted_lines = sorted(processed_lines, key=lambda x: x.split('\t')[2])
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
 

    process_blastp_pairs=process_blastp_result(blastp_pairs,num)
    process_synteny_pairs=process_synteny_result(synteny_pairs)
    process_gd_pairs=process_gd_result(gf,imap,spe1,spe2)
    generate_dotplot(gff1,gff2,lens1,lens2,process_blastp_pairs,spe1,spe2,'blastp_pairs')
    print('-'*30)
    generate_dotplot(gff1,gff2,lens1,lens2,process_synteny_pairs,spe1,spe2,'synteny_pairs')
    print('-'*30)
    generate_dotplot(gff1,gff2,lens1,lens2,process_gd_pairs,spe1,spe2,'gd_pairs')
    total_pairs=process_blastp_pairs+process_synteny_pairs+process_gd_pairs
    generate_dotplot(gff1,gff2,lens1,lens2,total_pairs,spe1,spe2,'total_pairs')
