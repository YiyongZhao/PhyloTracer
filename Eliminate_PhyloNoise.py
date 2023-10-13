import pandas as pd
from ete3 import TextFace,TreeStyle,Tree,NodeStyle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from PyPDF4 import PdfFileReader, PdfFileWriter
import os

def read_and_return_dict(filename, separator="\t") -> dict:
    df=pd.read_csv(filename,sep=separator,header=None)
    return df.set_index([0])[1].to_dict()

def read_tree(tre_path:str) -> object:
    return PhyloTree(tre_path)

def rename_input_tre(Phylo_t:object, gene2new_named_gene_dic:dict) -> object:
    for node in Phylo_t :
        sps=node.name[0:3]
        if sps in gene2new_named_gene_dic:
            node.name=gene2new_named_gene_dic[sps]+'_'+node.name
    return Phylo_t

def num_tre_node(Phylo_t:object)->object:#Numbering the nodes in the tree
    i = 1
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            node.name = "node" + str(i)
            i += 1
    return Phylo_t

def calculate_species_num(node:object)->int:# Obtain the number of species under a node
    leaf_list=node.get_leaf_names()
    species_num=len(set(i.split('_')[0] for i in leaf_list)) 
    return species_num

def get_species_list(Phylo_t:object)->list:
    leaves = Phylo_t.get_leaf_names()
    species_lst = [leaf.split("_")[0] for leaf in leaves]
    return species_lst

def get_pdf(t,c_dic):
    c_color_dict=get_color_dict(c_dic)
    for i in t.traverse():
        nstyle=NodeStyle()
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        i.set_style(nstyle)
        if  i.is_leaf():
            v=i.name.split('_')[1]
            if v in c_color_dict:
                color1 = c_color_dict[v].split('-')[-1]
                face = TextFace(''+i.name, fgcolor=color1,fstyle='italic')
                i.add_face(face, column=0)
    return t

def get_color_dict(dictory:dict)->dict:
    colormap = plt.get_cmap("gist_rainbow")
    # color dictionary\n",
    unique_values=set(dictory.values())
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_list = {k: v + '-' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}

    return sps_color_list 

def get_single_clades(Phylo_t,empty_set):
    if calculate_species_num(Phylo_t)==1:
        empty_set.add(Phylo_t)
        return
    for i in Phylo_t.get_children():
        get_single_clades(i,empty_set)

def get_node_single_taxa_dict(Phylo_t):     
    def get_single_taxa_caldes(empty_set):
        single_taxa_dict = {}
        for clade in empty_set:
            single_taxa = get_species_list(clade)[0]
            single_taxa_dict.setdefault(single_taxa, []).append(clade)
        return single_taxa_dict   
    empty_set=set()
    get_single_clades(Phylo_t,empty_set)
    single_taxa_dict =get_single_taxa_caldes(empty_set)      
    
    return single_taxa_dict                         

def judge(t):
    single_taxa_dict = get_node_single_taxa_dict(t)
    single_taxa_dict.pop('basal angiosperms')
    def check_single(dictionary):
        for value in dictionary.values():
            if len(value) != 1:
                return False
        return True

    return check_single(single_taxa_dict)

def get_single_inform_between_clade(clade1,clade2):
    ancestor=t.get_common_ancestor(clade1.name,clade2.name)
    empty_set=set()
    get_single_clades(ancestor,empty_set)
    empty_set -= {clade1, clade2}
    return empty_set

def paraphyletic_process(clades,rm_group,single_taxa)->list:
    clade1=clades[0]
    clade2=clades[1]
    clades_between_two_clade=get_single_inform_between_clade(clade1,clade2)
    if len(clades_between_two_clade) ==1:
        in_clade=list(clades_between_two_clade)[0]
        sps=get_species_list(in_clade)[0]
        if sps in single_taxa:
            if len(clade1)>len(clade2):
                rm_group.add(clade2.name) 
            else:
                rm_group.add(clade1.name)
        else:
                rm_group.add(in_clade.name)
    else:
        if len(clade1)>len(clade2):
            rm_group.add(clade2.name)
        elif len(clade1)==len(clade2):
            rm_group.add(clade1.name)

def polyphyletic_process(clades,rm_group,single_taxa):
    paraphyletic_combinations=set()
    for i in range(len(clades)):
        for j in range(i+1, len(clades)):
            paraphyletic_combinations.add((clades[i], clades[j]))

    for combination in  paraphyletic_combinations:     
        paraphyletic_process(combination,rm_group,single_taxa)
 
def merge_pdfs_side_by_side(file1, file2, output_file):
    pdf_writer = PdfFileWriter()

    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        pdf1 = PdfFileReader(f1)
        pdf2 = PdfFileReader(f2)
        page1 = pdf1.getPage(0)
        page2 = pdf2.getPage(0)

        new_page_width = page1.mediaBox.getWidth() + page2.mediaBox.getWidth()
        new_page_height = max(page1.mediaBox.getHeight(), page2.mediaBox.getHeight())

        new_page = pdf_writer.addBlankPage(new_page_width, new_page_height)
        new_page.mergeTranslatedPage(page1, 0, 0)
        new_page.mergeTranslatedPage(page2, page1.mediaBox.getWidth(), 0)

        with open(output_file, 'wb') as output:
            pdf_writer.write(output)

    f1.close()
    f2.close()

              

dic=read_and_return_dict('taxa.txt')
c_dic={k:v.split('_')[0] for k,v in dic.items()}
t=Tree('polytypic_tree/OGs_100000_T1_41_single_copy.treefile')
t.ladderize()
t.resolve_polytomy(recursive=True)
t.sort_descendants("support")
num_tre_node(t)
rename_input_tre(t,c_dic)
k='OGs_100000_T1_40_copy.treefile'
get_pdf(t,c_dic)
ts=TreeStyle()
ts.show_leaf_name=False
ts.title.add_face(TextFace(k+'_before', fsize=10), column=0)
t.render(file_name=k+'_before.pdf',tree_style=ts)
prune_single(t)
ts1=TreeStyle()
ts1.show_leaf_name=False
ts1.title.add_face(TextFace(k+'_after', fsize=10), column=0)
t.render(file_name=k+'_after.pdf',tree_style=ts1)
merge_pdfs_side_by_side(k+'_before.pdf', k+'_after.pdf', k+'.pdf')
os.remove(k+'_before.pdf')
os.remove(k+'_after.pdf')
 
