from __init__ import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from PyPDF4 import PdfFileReader, PdfFileWriter


def rename_input_single_tre(Phylo_t:object, gene2new_named_gene_dic:dict) -> object:
    for node in Phylo_t :
        sps=node.name.split('_')[0]
        if sps in gene2new_named_gene_dic:
            node.name=gene2new_named_gene_dic[sps]+'_'+node.name
    return Phylo_t

def rename_output_tre(Phylo_t:object) -> object:
    for node in Phylo_t :
        parts = node.name.split("_")  
        new_str = "_".join(parts[1:]) 
        node.name=new_str
    return Phylo_t


def is_multi_copy_tree(Phylo_t:object)->bool:
    leafs=Phylo_t.get_leaf_names()
    uniq_species=get_species_set(Phylo_t)
    if len(leafs) !=len(uniq_species):
        return True

def set_style(Phylo_t:object,color_dict:dict)->object:
    for node in Phylo_t.traverse():
        nstyle=NodeStyle()
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        node.set_style(nstyle)
        if  node.is_leaf():
            species_name=node.name.split('_')[1]
            if species_name in color_dict:
                color = color_dict[species_name].split('-')[-1]
                face = TextFace(''+node.name, fgcolor=color,fstyle='italic')
                node.add_face(face, column=0)
    return Phylo_t

def get_color_dict(dictory:dict)->dict:
    colormap = plt.get_cmap("gist_rainbow")
    # color dictionary\n",
    unique_values=set(dictory.values())
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_dict = {k: v + '-' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}

    return sps_color_dict

def get_single_clades(Phylo_t:object,empty_set:set):
    if calculate_species_num(Phylo_t)==1:
        empty_set.add(Phylo_t)
        return
    for i in Phylo_t.get_children():
        get_single_clades(i,empty_set)

def get_node_single_taxa_dict(Phylo_t:object)->dict:
    def get_single_taxa_caldes(empty_set):
        single_taxa_dict = {}
        for clade in empty_set:
            single_taxa = get_species_list(clade)[0]
            single_taxa_dict.setdefault(single_taxa, []).append(clade)
        return single_taxa_dict   
    empty_set=set()
    get_single_clades(Phylo_t,empty_set)
    single_taxa_dict =get_single_taxa_caldes(empty_set)      
    sorted_dict = dict(sorted(single_taxa_dict.items(), key=lambda item: len(item[1]), reverse=True))
    return sorted_dict                         

def is_single_tree(Phylo_t:object)->bool:
    single_taxa_dict = get_node_single_taxa_dict(Phylo_t)
    if 'basal angiosperms'  in single_taxa_dict.keys():
        single_taxa_dict.pop('basal angiosperms')
    if 'Outgroup' in single_taxa_dict.keys():
        single_taxa_dict.pop('Outgroup')
    def check_single(dictionary):
        for value in dictionary.values():
            if len(value) != 1:
                return False
        return True

    return check_single(single_taxa_dict)

   
def prune_single(Phylo_t:object):
    rm_list=[]
    single_taxa_dict=get_node_single_taxa_dict(Phylo_t)
    total_leafs=Phylo_t.get_leaf_names()
    for taxa,clades in single_taxa_dict.items():
        if len(clades) ==1 :
            pass
        elif len(clades)==2:
            clade1=clades[0]
            clade2=clades[1]
            if len(clade1)>len(clade2):
                leafs=clade2.get_leaf_names()
                diff = [leaf for leaf in total_leafs if leaf not in leafs]
                Phylo_t.prune(diff)

            else:
                leafs=clade1.get_leaf_names()
                diff = [leaf for leaf in total_leafs if leaf not in leafs]
                Phylo_t.prune(diff)
        else:
            result=[]
            for clade in clades :
                insertion_index=calculate_insertion_index(i)
                result.append(insertion_index)
            min_score=min(result)
            d=[]
            for j in v :
                insertion_index1=calculate_insertion_index(j)
                if insertion_index1 ==min_score:
                    branch_length_relative_score=calculate_branch_length_relative_score(Phylo_t,j)
                    d.append(branch_length_relative_score)
            min_score1=min(d)
            for h in v :
                branch_length_relative_score1=calculate_branch_length_relative_score(Phylo_t,h)
                if branch_length_relative_score1 ==min_score1:
                    rm_list.append(h.name)
                    
    for i in rm_list :
        clade=Phylo_t&i
        if clade.is_leaf():
            clade.delete()
        else:
            leafs=clade.get_leaf_names()
            diff = [a for a in total_leafs if a not in set(leafs)]
            Phylo_t.prune(diff)
    
                    
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

def calculate_avg_length(node):
    total_length = 0.0
    leaf_count = 0

    for leaf in node.iter_leaves():
        total_length += node.get_distance(leaf)
        leaf_count += 1

    if leaf_count > 0:
        avg_length = total_length / leaf_count
    else:
        avg_length = 0.0

    return avg_length+(node.dist/leaf_count)
    
def calculate_branch_length_relative_score(node):
    if not node.is_leaf():
        avg_length=calculate_avg_length(node)
        sister = node.get_sisters()[0] if not node.is_root() else None
        if not sister.is_leaf():
            sister_length = calculate_avg_length(sister) if sister else 0.0
        else:
            sister_length=sister.dist
        if sister_length != 0:
            return avg_length / sister_length
        
        else:
            return 0.0
       
    else:
        branch_length=node.dist
        sister = node.get_sisters()[0] if not node.is_root() else None
        if not sister.is_leaf():
            sister_length = calculate_avg_length(sister) if sister else 0.0
        else:
            sister_length=sister.dist
        if sister_length != 0:
            return branch_length / sister_length
        else:
            return 0.0

            
        
def get_root_relative_branch_ratio(leaf,avg_length):
    branch_length=leaf.dist
    return (branch_length-avg_length)/avg_length

def get_sister_relative_branch_ratio(leaf,sister):
    branch_length=leaf.dist
    if not sister.is_leaf():
        sister_length = calculate_avg_length(sister) if sister else 0.0
    else:
        sister_length=sister.dist
    if sister_length != 0:
        return (branch_length-sister_length)/sister_length 
    else:
        return 0.0
    
    


def calculate_insertion_index(node):
    insertion_index = 1.0 
    current_node = node
    while current_node.up:
        current_tips = len(current_node.get_leaves())
        parent_node = current_node.up
        parent_tips = len(parent_node.get_leaves())
        
        if parent_tips == 0:
            insertion_index = 0.0  
            break
        current_insertion_index = len(node.get_leaves())/ current_tips
        insertion_index = min(insertion_index, current_insertion_index)
        
        current_node = parent_node  
    
    return insertion_index


def prune_mc_main(tre_dic,taxa_dic,long_brancch_index):
    color_dic=get_color_dict(taxa_dic)
    dir_path1 = os.path.join(os.getcwd(), "pruned_tree")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    dir_path2 = os.path.join(os.getcwd(), "pruned_tree_pdf")
    if os.path.exists(dir_path2):
         shutil.rmtree(dir_path2)
    os.makedirs(dir_path2)
    dir_path3 = os.path.join(os.getcwd(), "long_branch_gene")
    if os.path.exists(dir_path3):
        shutil.rmtree(dir_path3)
    os.makedirs(dir_path3)
    for k,v in tre_dic.items():
        t=Tree(v)
        t.ladderize()
        t.resolve_polytomy(recursive=True)
        t.sort_descendants("support")
        num_tre_node(t)
        if is_multi_copy_tree(t):
            
            o = open(os.path.join(dir_path3, k + '_delete_gene.txt'), 'w')
            o.write('tre_ID'+'\t'+'delete_label'+'\t'+'gene'+'\t'+'root_relative_branch_ratio'+'\t'+'sister_relative_branch_ratio'+'\n')
            rename_input_single_tre(t,taxa_dic)

            set_style(t,color_dic)
            ts=TreeStyle()
            ts.show_leaf_name=False
            ts.title.add_face(TextFace(k+'_before', fsize=10), column=0)
            t.render(file_name=k+'_before.pdf',tree_style=ts)

            avg_length=sum([leaf.dist for leaf in t])/len(t)
            for leaf in t :
                sps_gene='_'.join(leaf.name.split('_')[1:])
                sister=leaf.get_sisters()[0] if not leaf.is_root() else None
                a=get_root_relative_branch_ratio(leaf,avg_length)
                b=get_sister_relative_branch_ratio(leaf,sister)
               
                if a  >long_brancch_index or b >long_brancch_index :
                     o.write(k+'\t'+'*'+'\t'+sps_gene+'\t'+str(a)+'\t'+str(b)+'\t'+'\n')
                     leaf.delete()
                else:
                     o.write(k+'\t'+'\t'+'\t'+sps_gene+'\t'+str(a)+'\t'+str(b)+'\t'+'\n')
                     
            
            ts1=TreeStyle()
            ts1.show_leaf_name=False
            ts1.title.add_face(TextFace(k+'_after', fsize=10), column=0)
            t.render(file_name=k+'_after.pdf',tree_style=ts1)
            merge_pdfs_side_by_side(k+'_before.pdf', k+'_after.pdf', os.path.join(dir_path2, k + '.pdf'))
            os.remove(k+'_before.pdf')
            os.remove(k+'_after.pdf')

            t1=rename_output_tre(t)
            t1.write(outfile=os.path.join(dir_path1, k + '.nwk'),format=0,dist_formatter='%.10f')
            o.close()

        else:
            print(k+' is not multi copy tree')
            continue
    
if __name__ == "__main__":
    os.makedirs(os.path.join(os.getcwd(), "pruned_tree"))
    taxa_dic=read_and_return_dict('taxa.txt')
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    long_brancch_index=5
    prune_main(tre_dic,taxa_dic,long_brancch_index)
   

