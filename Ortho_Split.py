from Base import *

def rename_len_dic(len_dic:dict, gene2new_named_gene_dic:dict) -> dict:
    return {gene2new_named_gene_dic[key]: value for key,value in len_dic.items() if key in gene2new_named_gene_dic.keys()}

def count_sps_num(ev_seqs:set)->int:
    sps_set={gene[0:3] for gene in ev_seqs }
    return len(sps_set) 
  
def find_tre_dup(Phylo_t:object) -> list:   #seperator either "@" or "_"
    tre_ParaL=[]
    GF_leaves_S = set(Phylo_t.get_leaf_names())
    events = Phylo_t.get_descendant_evol_events()
    for ev in events:
        if ev.etype == "D":
            tre_ParaL.append(",".join(ev.in_seqs)+"<=>"+",".join(ev.out_seqs))
    return tre_ParaL,GF_leaves_S
  
def extract_tree(ev_seqs:set,Phylo_t:object)-> object:
    Phylo_t.prune(ev_seqs)
    return Phylo_t

def offcut_tre(Phylo_t:object,renamed_len_dic:dict)-> list:
    tre_ParaL,GF_leaves_S = find_tre_dup(Phylo_t)
    offcut_ev_seqs_L0=[]
    for ev_seqs in tre_ParaL:
        fd=ev_seqs.strip().split("<=>")
        ev_seqs1,ev_seqs2=set(fd[0].split(",")),set(fd[1].split(","))
        sps_num1,seq_num1,sps_num2,seq_num2 = count_sps_num(ev_seqs1),len(ev_seqs1),count_sps_num(ev_seqs2),len(ev_seqs2)
        if sps_num1 == sps_num2:
            if seq_num1 > seq_num2:
                offcut_ev_seqs_L0.append(ev_seqs2)
            elif seq_num1 < seq_num2:
                offcut_ev_seqs_L0.append(ev_seqs1)
            elif seq_num1 == seq_num2:
                ev_seqs1_avg_len,ev_seqs2_avg_len = sum([renamed_len_dic[Gene] for Gene in ev_seqs1 ])/sps_num1,sum([renamed_len_dic[Gene] for Gene in ev_seqs2 ])/sps_num2
                if ev_seqs1_avg_len >= ev_seqs2_avg_len:
                    offcut_ev_seqs_L0.append(ev_seqs2)
                else:
                    offcut_ev_seqs_L0.append(ev_seqs1)
        else:
            if sps_num1 > sps_num2:
                offcut_ev_seqs_L0.append(ev_seqs2)
            else:
                offcut_ev_seqs_L0.append(ev_seqs1)
    offcut_Gene_S=set()
    for ev_seqs in offcut_ev_seqs_L0:
        offcut_Gene_S = offcut_Gene_S | ev_seqs
    principal_gene_S=GF_leaves_S - offcut_Gene_S
    filtered_offcut_ev_seqs_L0=[]
    for ev_seqs in offcut_ev_seqs_L0:
        sps_num = count_sps_num(ev_seqs)
        if sps_num >= 2:
            filtered_offcut_ev_seqs_L0.append(ev_seqs)
    filtered_offcut_ev_seqs_L0=rm_dup(filtered_offcut_ev_seqs_L0)
    return principal_gene_S,filtered_offcut_ev_seqs_L0# filtered_offcut_ev_seqs_L0=offcut_ev_seqs_L0

def rm_dup(paralogs_L:list)->list: #[{1,2},{2,3}]
    for ev_seqs1 in paralogs_L:
        for ev_seqs2 in paralogs_L:
            if ev_seqs1.issubset(ev_seqs2):
                if ev_seqs1 in paralogs_L:
                    paralogs_L.remove(ev_seqs1)
            elif ev_seqs2.issubset(ev_seqs1):
                if ev_seqs2 in paralogs_L:
                    paralogs_L.remove(ev_seqs2)
    return paralogs_L

def split_offcut_ev_seqs(offcut_ev_seqs_L0:list) ->list:  #[{1,2},{2,3}]
    othologs_L=[]
    paralogs_L=[]
    for ev_seqs in offcut_ev_seqs_L0:
        sps_num,seq_num = count_sps_num(ev_seqs),len(ev_seqs)
        if sps_num >=2:
            if seq_num > sps_num:
                paralogs_L.append(ev_seqs)
            else:
                othologs_L.append(ev_seqs)
    paralogs_L=rm_dup(paralogs_L)
    return othologs_L,paralogs_L

def iterator(offcut_ev_seqs_L0:list,Phylo_t:object,new_named_gene2gene_dic:dict,minor_othologs_L:list,tre_path:str,renamed_len_dic:dict) ->list:
    
    othologs_L,paralogs_L=split_offcut_ev_seqs(offcut_ev_seqs_L0)
    minor_othologs_L+=(othologs_L)
    if paralogs_L !=[]:    
        for i,ev_seqs in enumerate(paralogs_L):
            Phylo_t=read_tree(tre_path)
            Phylo_t=rename_input_tre(Phylo_t,new_named_gene2gene_dic)
            Phylo_t = extract_tree(ev_seqs,Phylo_t)
            principal_gene_S,offcut_ev_seqs_L0=offcut_tre(Phylo_t,renamed_len_dic)
            minor_othologs_L.append(principal_gene_S)
            iterator(offcut_ev_seqs_L0,Phylo_t,new_named_gene2gene_dic,minor_othologs_L,tre_path,renamed_len_dic)
    return minor_othologs_L
      
def rename_OGs_tre_name(principal_gene_S:list,minor_othologs_L:list,tre_ID:str)->list:
    OG_count = len(minor_othologs_L)+1
    tre_name_L =["T"+str(i+1) for i in range(OG_count)]
    sps_num_L=[len(OG) for OG in minor_othologs_L]
        
    tre_name2sps_num_dict=dict(zip(tre_name_L,sps_num_L))
    tre_name2OG_dict=dict(zip(tre_name_L,minor_othologs_L))
    orderd_OG_sps_num_L = sorted(tre_name2sps_num_dict.items(), key=lambda x: x[1],reverse =True)
        
    ordered_name_OG_L=[]
    ordered_name_OG_L.append((str(tre_ID)+"_"+tre_name_L[0]+"_"+str(len(principal_gene_S)),principal_gene_S))
    index=0
    for name,count in orderd_OG_sps_num_L:
        index+=1
        ordered_name_OG_L.append((str(tre_ID)+"_"+tre_name_L[index]+"_"+str(count),tre_name2OG_dict[name]))
    return ordered_name_OG_L
#ordered_name_OG_L=rename_OGs_tre_name(principal_gene_S,minor_othologs_L,tre_ID)
##########################################################################
def split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic,out_file_name):
    o=open(out_file_name,'w')
    for tre_ID,tre_path in tre_dic.items():
        Phylo_t0=read_tree(tre_path)
        Phylo_t1=root_tre_with_midpoint_outgroup(Phylo_t0)
        Phylo_t1=rename_input_tre(Phylo_t1,gene2new_named_gene_dic)
        principal_gene_S,filtered_offcut_ev_seqs_L0=offcut_tre(Phylo_t1,renamed_len_dic)
        minor_othologs_L=[]
        minor_othologs_L=iterator(filtered_offcut_ev_seqs_L0,Phylo_t1,gene2new_named_gene_dic,minor_othologs_L,tre_path,renamed_len_dic)
        ordered_name_OG_L=rename_OGs_tre_name(principal_gene_S,minor_othologs_L,tre_ID)
        for tre_name,OG_S in ordered_name_OG_L: 
            OG_L = [new_named_gene2gene_dic[OG] for OG in OG_S]
            Phylo_t0=read_tree(tre_path)
            Phylo_t=root_tre_with_midpoint_outgroup(Phylo_t0)
            Phylo_t_OG_L=extract_tree(OG_L,Phylo_t)
            o.write(tre_name+"\t"+Phylo_t_OG_L.write()+"\n")

if __name__ == "__main__":
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    len_dic=read_and_return_dict('length_list')
    renamed_len_dic=rename_len_dic(len_dic,gene2new_named_gene_dic)
    tre_dic=read_and_return_dict('GF_list')   
    split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic,out_file_name)