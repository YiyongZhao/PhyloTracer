from __init__ import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from PyPDF4 import PdfFileReader, PdfFileWriter
from BranchLength_NumericConverter import write_tree_to_newick

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
        if  node.is_leaf() and node.name !='Outgroups':
            species_name=node.name.split('_')[1]
            if species_name in color_dict:
                color = color_dict[species_name].split('*')[-1]
                face = TextFace(''+node.name, fgcolor=color,fstyle='italic')
                node.add_face(face, column=0)
    return Phylo_t

def get_color_dict(dictory:dict)->dict:
    colormap = plt.get_cmap("gist_rainbow")
    unique_values=set(dictory.values())
    colors_lst = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict=dict(zip(unique_values,colors_lst)) 
    sps_color_dict = {k: v + '*' + color_dict.get(v) for k, v in dictory.items() if v in color_dict}
    return sps_color_dict

def get_single_clades(Phylo_t:object,empty_set:set)->None:
    if calculate_species_num(Phylo_t)==1:
        empty_set.add(Phylo_t)
        return
    for child in Phylo_t.get_children():
        get_single_clades(child,empty_set)                        
                
def merge_pdfs_side_by_side(file1:str, file2:str, output_file:str)->None:
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

def calculate_avg_length(node:object)->int:
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
    
def calculate_branch_length_relative_score(node:object)->int:
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

def get_root_relative_branch_ratio(leaf:object,avg_length:int)->int:
    branch_length=leaf.dist
    return (branch_length-avg_length)/avg_length

def get_sister_relative_branch_ratio(leaf:object,sister:object)->int:
    branch_length=leaf.dist
    if not sister.is_leaf():
        sister_length = calculate_avg_length(sister) if sister else 0.0
    else:
        sister_length=sister.dist
    if sister_length != 0:
        return (branch_length-sister_length)/sister_length 
    else:
        return 0.0
    
def calculate_insertion_depth(clade: object, node: object) -> int:
    if clade is None:
        return 0
    return clade.get_distance(node, topology_only=True)

def calculate_insertion_coverage(clade: object, node: object) -> float:
    if clade is None or node is None:
        return 0.0  
    return len(node) / len(clade)

def calculate_insertion(clade:object, node:object)->int:
    depth = calculate_insertion_depth(clade, node)
    coverage = calculate_insertion_coverage(clade, node)
    return depth / coverage if coverage != 0 else 0

def get_target_clade(clade:object)->object:
    node2root = clade.get_ancestors()
    target_clade = None
    for ancestor in node2root:
        if len(get_species_set(ancestor)) > 2:
            break
        elif len(get_species_set(ancestor)) == 2:
            target_clade = ancestor      
    return target_clade

def is_ancestor_sister_same(node:object)->bool:
    sps = get_species_list(node)[0]
    sis=node.get_sisters()[0]
    sis_name=get_species_list(sis)[0]
    if sps==sis_name:
        return False
    else:
        an=node.up
        if an.is_root():
            return False
        else:
            an_si=an.get_sisters()[0]
            an_name=get_species_list(an_si)[0]
            return sis_name==an_name

def remove_long_gene(Phylo_t:object,long_branch_index:int,outfile:str,tre_ID:str)->object:
    Phylo_t1=Phylo_t.copy()
    remove_gene_set=set()
    avg_length=sum([leaf.dist for leaf in Phylo_t1])/len(Phylo_t1)
    for leaf in Phylo_t1 :
        sps_gene='_'.join(leaf.name.split('_')[1:])
        sister=leaf.get_sisters()[0] if not leaf.is_root() else None
        root2clade_radio=get_root_relative_branch_ratio(leaf,avg_length)
        sister2clade_radio=get_sister_relative_branch_ratio(leaf,sister)
        if  root2clade_radio  >long_branch_index or sister2clade_radio >long_branch_index :
            outfile.write(tre_ID+'\t'+'*'+'\t'+sps_gene+'\t'+str(root2clade_radio)+'\t'+str(sister2clade_radio)+'\t'+'\n')
            remove_gene_set.add(leaf.name)
        else:
            outfile.write(tre_ID+'\t'+'\t'+'\t'+sps_gene+'\t'+str(root2clade_radio)+'\t'+str(sister2clade_radio)+'\t'+'\n')
    total_leafs_set=set(Phylo_t1.get_leaf_names())
    diff=total_leafs_set - remove_gene_set
    Phylo_t1.prune(diff,preserve_branch_length=True)
    return Phylo_t1

def remove_insert_gene(Phylo_t:object,long_branch_index:int,outfile:str,tre_ID:str)->object:
    Phylo_t1 = Phylo_t.copy()    
    taxa_clade = set()
    get_single_clades(Phylo_t1, taxa_clade)
    remove_gene_set = set()

    for clade in taxa_clade:
        if is_ancestor_sister_same(clade):
            target_clade = get_target_clade(clade)
            if target_clade==None:
                continue
            else:
                index_num = calculate_insertion_depth(target_clade, clade)
                index_over = calculate_insertion_coverage(target_clade, clade)
                insert=calculate_insertion(target_clade, clade)
                outfile.write(tre_ID+'\t'+'@'+'\t'+clade.name+'\t'+str(index_num)+'\t'+str(index_over)+'\t'+str(insert)+'\n')     
                remove_gene_set.add(clade.name)              
    total_leafs_set=set(Phylo_t1.get_leaf_names())
    diff=total_leafs_set - remove_gene_set
    Phylo_t1.prune(diff,preserve_branch_length=True)
    return Phylo_t1

def generate_pdf_before(tre_ID,Phylo_t)->None:
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(tre_ID + '_before', fsize=10), column=0)
    Phylo_t.render(file_name=tre_ID + '_before.pdf', tree_style=ts)

def generate_pdf_after(tre_ID,Phylo_t)->None:
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(tre_ID + '_after', fsize=10), column=0)
    Phylo_t.render(file_name=tre_ID + '_after.pdf', tree_style=ts)

def prune_mc_main(tre_dic:dict, taxa_dic:dict, long_branch_index:int, insert_branch_index:int)->None:
    color_dic = get_color_dict(taxa_dic)
    dir_path1 = os.path.join(os.getcwd(), "output/pruned_tree/")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    dir_path2 = os.path.join(os.getcwd(), "output/pruned_tree_pdf/")
    if os.path.exists(dir_path2):
        shutil.rmtree(dir_path2)
    os.makedirs(dir_path2)
    dir_path3 = os.path.join(os.getcwd(), "output/long_branch_gene/")
    if os.path.exists(dir_path3):
        shutil.rmtree(dir_path3)
    os.makedirs(dir_path3)

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t = Tree(tre_path)
        t.ladderize()
        t.resolve_polytomy(recursive=True)
        t.sort_descendants("support")
        num_tre_node(t)
        if is_multi_copy_tree(t):
            o = open(os.path.join(dir_path3, tre_ID + '_delete_gene.txt'), 'w')
            o.write('tre_ID' + '\t' + 'long_branch_label' + '\t' + 'gene' + '\t' +'root_relative_branch_ratio' + '\t' + 'sister_relative_branch_ratio' + '\n')
            rename_input_single_tre(t, taxa_dic)
            set_style(t, color_dic)
            generate_pdf_before(tre_ID,t)
            t1 = remove_long_gene(t, long_branch_index, o, tre_ID)
            o.write('\n')
            if len(get_species_set(t1)) !=1:
                o.write('tre_ID' + '\t' + 'insert_branch_label' + '\t' + 'gene' + '\t' + 'insertion_depth' + '\t' + 'insertion_coverage' + '\t' + 'calculate_insertion' + '\n')
                t2 = remove_insert_gene(t1, insert_branch_index, o, tre_ID)
                o.close()
                generate_pdf_after(tre_ID,t2)
                merge_pdfs_side_by_side(tre_ID + '_before.pdf', tre_ID + '_after.pdf', os.path.join(dir_path2, tre_ID + '.pdf'))
                os.remove(tre_ID + '_before.pdf')
                os.remove(tre_ID + '_after.pdf')
                t3 = rename_output_tre(t2)
                tree_str = t3.write(format=0)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)
            else:
                generate_pdf_after(tre_ID,t1)
                merge_pdfs_side_by_side(tre_ID + '_before.pdf', tre_ID + '_after.pdf', os.path.join(dir_path2, tre_ID + '.pdf'))
                os.remove(tre_ID + '_before.pdf')
                os.remove(tre_ID + '_after.pdf')
                t3 = rename_output_tre(t1)
                tree_str = t3.write(format=0)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)


        else:
            print(tre_ID + ' is not multi copy tree')
            continue
        pbar.update(1)
    pbar.close()
    
if __name__ == "__main__":
    os.makedirs(os.path.join(os.getcwd(), "pruned_tree"))
    taxa_dic=read_and_return_dict('taxa.txt')
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    long_brancch_index=5
    prune_mc_main(tre_dic,taxa_dic,long_brancch_index,insert_branch_index)
