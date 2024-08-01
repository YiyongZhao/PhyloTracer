from __init__ import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from PyPDF4 import PdfFileReader, PdfFileWriter
from BranchLength_NumericConverter import trans_branch_length, write_tree_to_newick
import os
import shutil
from tqdm import tqdm

def has_multiple_copies(phylo_t: object) -> bool:
    """Check if the tree has multiple copies based on unique species."""
    leaf_names = phylo_t.get_leaf_names()
    unique_species = get_species_set(phylo_t)
    return len(leaf_names) != len(unique_species)

def style_tree(phylo_t: object, color_dict: dict, new_named_gene2gene_dic: dict) -> object:
    """Set styles for tree nodes based on gene colors."""
    for node in phylo_t.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        node.set_style(nstyle)
        if node.is_leaf():
            parts = node.name.split("_")
            species_name = parts[0]
            gene = new_named_gene2gene_dic[node.name]
            if species_name in color_dict:
                color = color_dict[species_name].split('*')[-1]
                face = TextFace(f'{gene}', fgcolor=color, fstyle='italic')
                node.add_face(face, column=0)
    return phylo_t

def create_color_mapping(dictionary: dict) -> dict:
    """Create a color mapping for unique values in the dictionary."""
    colormap = plt.get_cmap("gist_rainbow")
    unique_values = set(dictionary.values())
    colors_list = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 1, len(unique_values))]
    color_dict = dict(zip(unique_values, colors_list)) 
    return {k: v + '*' + color_dict.get(v) for k, v in dictionary.items() if v in color_dict}

def merge_pdfs(file1: str, file2: str, output_file: str) -> None:
    """Merge two PDF files side by side."""
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

def get_average_tip_length(phylo_t: object) -> float:
    """Calculate the average length of tips in the phylogenetic tree."""
    return sum(leaf.dist for leaf in phylo_t)/len(phylo_t)

def get_average_node_length(phylo_t: object) -> float:
    """Calculate the average length of nodes in the phylogenetic tree."""
    return sum(phylo_t.get_distance(i)+phylo_t.dist for i in phylo_t)/ len(phylo_t)

def remove_long_branches(phylo_t: object, absolute_branch_length: int, relative_branch_length: int,outfile: str, tre_ID: str, new_named_gene2gene_dic: dict) -> object:
    """Remove branches longer than the specified index from the phylogenetic tree."""
    phylo_t1 = phylo_t.copy()
    remove_gene_set = set()
    tips_avg_length = get_average_tip_length(phylo_t1)
  
    for leaf in phylo_t1:
        if leaf.dist!=0:
            sps_gene = leaf.name
            distance = leaf.dist
            distance_to_root_ratio = (distance-tips_avg_length) / tips_avg_length
            sister = leaf.get_sisters()[0]
            if sister.is_leaf():
                sister_avg_length=sister.dist
                if sister_avg_length==0:
                    leaf_to_sister_ratio=4
                else:
                    leaf_to_sister_ratio = (distance - sister_avg_length) / sister_avg_length
            else:
                sister_avg_length=get_average_node_length(sister)
                if sister_avg_length==0:
                    leaf_to_sister_ratio=4
                else:
                    leaf_to_sister_ratio = (distance - sister_avg_length) / sister_avg_length
            
            if distance_to_root_ratio >= absolute_branch_length:
                if leaf_to_sister_ratio >= relative_branch_length:
                    outfile.write(f"{tre_ID}\t*\t{new_named_gene2gene_dic[sps_gene]}\t{distance_to_root_ratio}\t{leaf_to_sister_ratio}\n")
                    remove_gene_set.add(leaf.name)
                else:
                    outfile.write(f"{tre_ID}\t\t{new_named_gene2gene_dic[sps_gene]}\t{distance_to_root_ratio}\t{leaf_to_sister_ratio}\n")

                #is_modified = True  
            else:
                outfile.write(f"{tre_ID}\t\t{new_named_gene2gene_dic[sps_gene]}\t{distance_to_root_ratio}\t{leaf_to_sister_ratio}\n")
        else:
            outfile.write(f"{tre_ID}\t\t{new_named_gene2gene_dic[sps_gene]}\t{distance_to_root_ratio}\t{leaf_to_sister_ratio}\n")
            continue

    # if not is_modified:
    #     break
        
    total_leafs_set = set(phylo_t1.get_leaf_names())
    diff = total_leafs_set - remove_gene_set
    phylo_t1.prune(diff, preserve_branch_length=True)
    
    return phylo_t1

def generate_pdf_before(tre_ID: str, phylo_t: object) -> None:
    """Generate a PDF from the phylogenetic tree before pruning."""
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(f'{tre_ID}_before', fsize=10), column=0)
    phylo_t.render(file_name=f'{tre_ID}_before.pdf', tree_style=ts)

def generate_pdf_after(tre_ID: str, phylo_t: object) -> None:
    """Generate a PDF from the phylogenetic tree after pruning."""
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(f'{tre_ID}_after', fsize=10), column=0)
    phylo_t.render(file_name=f'{tre_ID}_after.pdf', tree_style=ts)

def prune_main_LB(tre_dic: dict, voucher2taxa_dic: dict, gene2new_named_gene_dic: dict, new_named_gene2gene_dic: dict, absolute_branch_length:int=5, relative_branch_length: int=5, visual: bool = False) -> None:
    """Main function to prune long branches from phylogenetic trees and visualize the results."""
    color_dic = create_color_mapping(voucher2taxa_dic)

    dir_path_pruned = os.path.join(os.getcwd(), "orthofilter_lb/pruned_tree/")
    shutil.rmtree(dir_path_pruned, ignore_errors=True)
    os.makedirs(dir_path_pruned)

    if visual:
        dir_path_pdf = os.path.join(os.getcwd(), "orthofilter_lb/pruned_tree_pdf/")
        shutil.rmtree(dir_path_pdf, ignore_errors=True)
        os.makedirs(dir_path_pdf)

    dir_path_long_branch = os.path.join(os.getcwd(), "orthofilter_lb/long_branch_gene/")
    shutil.rmtree(dir_path_long_branch, ignore_errors=True)
    os.makedirs(dir_path_long_branch)

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t0 = Tree(tre_path)
        t0.ladderize()
        t0.resolve_polytomy(recursive=True)
        t0.sort_descendants("support")
        t=rename_input_tre(t0,gene2new_named_gene_dic)
        num_tre_node(t)

        output_file = open(os.path.join(dir_path_long_branch, f'{tre_ID}_delete_gene.txt'), 'w')
        output_file.write('tre_ID\tlong_branch_label\tgene\troot_relative_branch_ratio\tsister_relative_branch_ratio\n')        
        if has_multiple_copies(t):
            if visual:
                style_tree(t, color_dic, new_named_gene2gene_dic)
                generate_pdf_before(tre_ID, t)

            t1 = remove_long_branches(t, absolute_branch_length, relative_branch_length,output_file, tre_ID, new_named_gene2gene_dic)
            output_file.close()

            if visual:
                generate_pdf_after(tre_ID, t1)
                merge_pdfs(f"{tre_ID}_before.pdf", f"{tre_ID}_after.pdf", os.path.join(dir_path_pdf, f"{tre_ID}.pdf"))
                os.remove(f"{tre_ID}_before.pdf")
                os.remove(f"{tre_ID}_after.pdf")

            t2 = rename_input_tre(t1, new_named_gene2gene_dic)
            tree_str = trans_branch_length(t2)
            write_tree_to_newick(tree_str, tre_ID, dir_path_pruned)
        else:
            if visual:
                style_tree(t, color_dic, new_named_gene2gene_dic)
                generate_pdf_before(tre_ID, t)

            t1 = remove_long_branches(t, absolute_branch_length, relative_branch_length,output_file, tre_ID, new_named_gene2gene_dic)
            output_file.close()

            if visual:
                generate_pdf_after(tre_ID, t1)
                merge_pdfs(f"{tre_ID}_before.pdf", f"{tre_ID}_after.pdf", os.path.join(dir_path_pdf, f"{tre_ID}.pdf"))
                os.remove(f"{tre_ID}_before.pdf")
                os.remove(f"{tre_ID}_after.pdf")

            t2 = rename_input_tre(t1, new_named_gene2gene_dic)
            tree_str = trans_branch_length(t2)
            write_tree_to_newick(tree_str, tre_ID, dir_path_pruned)

        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    tre_dic = read_and_return_dict('gf')
    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic,taxa2voucher_dic = gene_id_transfer('imap.txt')
    absolute_branch_length=4
    relative_branch_length=1.5
    prune_main_LB(tre_dic, voucher2taxa_dic,gene2new_named_gene_dic, new_named_gene2gene_dic, absolute_branch_length,relative_branch_length, visual=True)
