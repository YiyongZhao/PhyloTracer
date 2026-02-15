# PhyloTracer

PhyloTracer is a command-line toolkit for phylogenomics and comparative genomics.

It is built for practical workflows on gene-family trees: rooting, noise filtering, ortholog retrieval, gene-duplication (GD) and loss analysis, topology summaries, hybridization screening, and visualization.

---

## What

Given gene trees, a species tree, and mapping files, PhyloTracer helps you:

- Root gene trees with species-aware scoring.
- Remove long-branch and non-monophyletic noise.
- Detect GD events and classify GD types.
- Track post-GD loss patterns at event/path level.
- Summarize topology frequencies across gene families.
- Build GD-aware hybridization tests (HyDe workflow).
- Generate publication-ready tree and species-tree visualizations.

All tools are exposed through one executable:

```bash
phylotracer <SUBCOMMAND> [options]
```

---

## Key Features

- **Modular**: each subcommand can run independently.
- **Pipeline-friendly**: simple TSV/Newick inputs and plain-text outputs.
- **Tree-focused**: consistent handling of gene and species trees.
- **GD-centric**: detection, typing, visualization, and loss tracking.
- **Reproducible**: deterministic CLI behavior and stable output files.

---

## Installation

### 1) Conda (recommended)

```bash
git clone https://github.com/YiyongZhao/PhyloTracer.git
cd PhyloTracer
conda env create -f environment.yml
conda activate phylotracer
```

### 2) PyPI

```bash
pip install PhyloTracer
```

### 3) Optional headless rendering fix

```bash
export QT_QPA_PLATFORM=offscreen
```

---

## Basic Usage

### Show help

```bash
phylotracer -h
```

### Example: GD detection

```bash
phylotracer GD_Detector \
  --input_GF_list example_data/GD_Detector/GF_ID2path.imap \
  --input_imap example_data/GD_Detector/gene2sps.imap \
  --input_sps_tree example_data/GD_Detector/sptree.nwk \
  --gd_support 50 \
  --subclade_support 50 \
  --dup_species_proportion 0 \
  --dup_species_num 2 \
  --deepvar 1
```

### Example: Rooting

```bash
phylotracer Phylo_Rooter \
  --input_GF_list example_data/Phylo_Rooter/GF_ID2path.imap \
  --input_imap example_data/Phylo_Rooter/gene2sps.imap \
  --input_gene_length example_data/Phylo_Rooter/gene2length.imap \
  --input_sps_tree example_data/Phylo_Rooter/sptree.nwk
```

---

## Typical Workflow

Recommended analysis order:

1. **Preprocess trees**
   - `PhyloSupport_Scaler`, `BranchLength_NumericConverter`, `PhyloTree_CollapseExpand`
2. **Root and denoise**
   - `Phylo_Rooter`, then `OrthoFilter_LB` / `OrthoFilter_Mono`
3. **Detect and summarize GD**
   - `GD_Detector`, then `GD_Visualizer`
4. **Track GD-associated loss**
   - `GD_Loss_Tracker`, then `GD_Loss_Visualizer`
5. **Downstream summaries**
   - `TreeTopology_Summarizer`, `Hybrid_Tracer`, `Hybrid_Visualizer`, `Ortho_Retriever`, `HaploFinder`

---

## Core Input Files

Most modules use two-column TSV files:

- `GF_ID2path.imap`: `<gene_family_id><TAB><gene_tree_path>`
- `gene2sps.imap`: `<gene_id><TAB><species_name>`
- `gene2length.imap`: `<gene_id><TAB><gene_length>`

Optional module-specific files include:

- taxonomy/category maps (`gene_id<TAB>label`)
- alignment maps (`GF_ID<TAB>alignment_path`)
- GFF/lens files
- expression matrix (`.csv/.xls/.xlsx`)

---

## Command Line Options

PhyloTracer follows the OrthoFinder-style CLI pattern:

```bash
phylotracer <SUBCOMMAND> [OPTIONS]
```

Get options for one module:

```bash
phylotracer <SUBCOMMAND> -h
```

Below is the full command-line option reference for all 16 modules.

### 1) `PhyloTree_CollapseExpand`

**Purpose**
- Collapse internal nodes by support threshold, optionally revert to binary form.

**Command line options (required)**
- `--input_GF_list`: gene-tree list (`GF_ID<TAB>tree_path`)
- `--support_value`: collapse cutoff (0–100)

**Command line options (optional)**
- `--revert`: expand comb-like structures back to binary trees

**Input**
- Newick gene trees from `--input_GF_list`

**Output**
- Directory: `collapse_expand_tree/`
- Files: one `.nwk` per input family

---

### 2) `PhyloSupport_Scaler`

**Purpose**
- Rescale branch support values between `[0,1]` and `[0,100]`.

**Command line options (required)**
- `--input_GF_list`
- `--scale_to` (`1` or `100`)

**Input**
- Gene trees with support values

**Output**
- Directory: `support_scaler_tree/`
- Files: one `.nwk` per input family

---

### 3) `BranchLength_NumericConverter`

**Purpose**
- Standardize branch-length decimal precision.

**Command line options (required)**
- `--input_GF_list`

**Command line options (optional)**
- `--decimal_place` (default `10`)

**Input**
- Gene trees

**Output**
- Directory: `converter_tree/`
- Files: one `.nwk` per input family

---

### 4) `Phylo_Rooter`

**Purpose**
- Root gene trees using species-tree context and gene lengths.

**Command line options (required)**
- `--input_GF_list`
- `--input_imap` (`gene_id<TAB>species_name`)
- `--input_gene_length` (`gene_id<TAB>gene_length`)
- `--input_sps_tree` (Newick)

**Input**
- Gene trees + species map + species tree + gene lengths

**Output**
- Directory: `rooted_trees/` (rooted gene trees)
- File: `stat_matrix.csv` (rooting statistics)

---

### 5) `OrthoFilter_LB`

**Purpose**
- Remove long-branch outlier tips.

**Command line options (required)**
- `--input_GF_list`
- `--input_imap`
- `--absolute_branch_length`
- `--relative_branch_length`

**Command line options (optional)**
- `--visual`

**Input**
- Gene trees + species map

**Output**
- Directory: `orthofilter_lb/pruned_tree/` (pruned trees)
- Directory: `orthofilter_lb/long_branch_gene/` (`*_delete_gene.txt` logs)
- Optional directory: `orthofilter_lb/pruned_tree_pdf/` (before/after PDFs)

---

### 6) `OrthoFilter_Mono`

**Purpose**
- Remove alien/non-monophyletic tips inside dominant lineages.

**Command line options (required)**
- `--input_GF_list`
- `--input_taxa` (`gene_id<TAB>clade_label`)
- `--input_imap`
- `--input_sps_tree`

**Command line options (optional)**
- `--purity_cutoff` (default `0.95`)
- `--max_remove_fraction` (default `0.5`)
- `--visual`

**Input**
- Gene trees + taxon labels + species map + species tree

**Output**
- Directory: `orthofilter_mono/pruned_tree/`
- Directory: `orthofilter_mono/insert_gene/` (`*_insert_gene.txt` scoring logs)
- Optional directory: `orthofilter_mono/visual/`

---

### 7) `TreeTopology_Summarizer`

**Purpose**
- Count and rank absolute/relative topologies.

**Command line options (required)**
- `--input_GF_list`
- `--input_imap`

**Command line options (optional)**
- `--visual_top` (default `10`)

**Input**
- Gene trees + species map

**Output**
- Files: `absolute_topology.txt`, `relative_topology.txt` style summaries
- Files: `absolute_*.txt`, `relative_*.txt`
- Figures: `merge_absolutely_top*.png`, `merge_relative_top*.png`

---

### 8) `Tree_Visualizer`

**Purpose**
- Render gene trees with optional annotations (categories, GD, family, expression).

**Command line options (required)**
- `--input_GF_list`
- `--input_imap`

**Command line options (optional)**
- `--gene_categories`
- `--keep_branch` (`0`/`1`)
- `--tree_style` (`r`/`c`)
- `--gene_family`
- `--input_sps_tree` (required with `--gene_family`)
- `--gene_expression`
- `--visual_gd`

**Input**
- Gene trees + species map + optional annotation files

**Output**
- Directory: `tree_visualizer/` (per-tree PDFs)
- File: `genefamily_map2_sptree.pdf` (when species-tree family mapping is used)

---

### 9) `GD_Detector`

**Purpose**
- Detect GD events and classify GD types.

**Command line options (required)**
- `--input_GF_list`
- `--input_imap`
- `--gd_support`
- `--subclade_support`
- `--dup_species_proportion`
- `--dup_species_num`
- `--input_sps_tree`
- `--deepvar`

**Command line options (optional)**
- `--gdtype_mode` (`relaxed` default, `strict`)

**Input**
- Gene trees + species map + species tree

**Output**
- File: `gd_result_relaxed.txt` or `gd_result_strict.txt`
- File: `gd_type_relaxed.tsv` or `gd_type_strict.tsv`
- File: `numed_sptree.nwk` (written by CLI wrapper)

---

### 10) `GD_Visualizer`

**Purpose**
- Plot GD counts/types on a species tree.

**Command line options (required)**
- `--input_sps_tree`
- `--gd_result`
- `--input_imap`

**Input**
- Numbered species tree + GD result table + species map

**Output**
- File: `phylotracer_gd_visualizer.pdf`

---

### 11) `GD_Loss_Tracker`

**Purpose**
- Track post-GD loss at event and path levels.

**Command line options (required)**
- `--input_GF_list`
- `--input_sps_tree`
- `--input_imap`

**Command line options (optional)**
- `--target_species` (repeatable)
- `--mrca_node` (`SP1,SP2`, repeatable)
- `--include_unobserved_species`

**Input**
- Gene trees + species tree + species map

**Output**
- File: `gd_loss_summary.txt`
- File: `gd_loss_count_summary.txt`
- File: `gd_loss.xlsx`

---

### 12) `GD_Loss_Visualizer`

**Purpose**
- Visualize GD-loss summaries on species-tree topology.

**Command line options (required)**
- `--gd_loss_result` (expects `gd_loss_summary.txt` format)

**Command line options (optional)**
- `--input_sps_tree`

**Input**
- GD-loss summary table (+ species tree)

**Output**
- File: `gd_loss_pie_visualizer.PDF`

---

### 13) `Ortho_Retriever`

**Purpose**
- Split paralogous trees and retrieve ortholog groups.

**Command line options (required)**
- `--input_GF_list`
- `--input_imap`
- `--input_gene_length`

**Input**
- Rooted gene trees + species map + gene lengths

**Output**
- File: `ortho_retriever_summary.txt`
- File: `ortholog_trees.tsv`

---

### 14) `Hybrid_Tracer`

**Purpose**
- Build GD-aware matrices and run HyDe-based tests.

**Command line options (required)**
- `--input_GF_list`
- `--input_Seq_GF_list`
- `--input_sps_tree`
- `--input_imap`

**Command line options (optional)**
- `--mrca_node`
- `--split_groups` (default `1`)

**Input**
- Gene trees + alignments + species tree + species map

**Output**
- File: `hyde_out.txt`
- File: `hyde_filtered_out.txt`

---

### 15) `Hybrid_Visualizer`

**Purpose**
- Visualize HyDe results in leaf-mode or node-mode heatmaps.

**Command line options (required)**
- `--hyde_out`
- `--input_sps_tree`

**Command line options (optional)**
- `--node`

**Input**
- HyDe output + species tree

**Output**
- Files with patterns like `*_img_faces.png` and `*_hotmap.png`

---

### 16) `HaploFinder`

**Purpose**
- Detect haplotype-level GD signals (`mode=haplofinder`) or split FASTA (`mode=split`).

**Required (global)**
- `--mode` (`haplofinder` or `split`)

**Required in `haplofinder` mode**
- `--input_GF_list`, `--input_imap`, `--input_sps_tree`
- `--species_a`, `--species_b`
- `--species_a_gff`, `--species_b_gff`
- `--species_a_lens`, `--species_b_lens`

**Optional in `haplofinder` mode**
- `--visual_chr_a`, `--visual_chr_b`
- `--gd_support` (default `50`)
- `--pair_support` (default `50`)
- `--size` (default `0.0005`)

**Required in `split` mode (CLI checks)**
- `--input_GF_list`, `--input_imap`, `--input_fasta`
- `--cluster_file`, `--hyb_sps`, `--parental_sps`, `--species_b_gff`

**Input**
- Mode-specific tree/mapping/annotation/FASTA files

**Output**
- Dotplots: `*_dotplot.pdf`, `*_dotplot.png`
- Conversion summaries: `gene_conversion_*.txt`, `gene_conversion.txt`
- Color labels / split outputs: e.g., `color_label.txt`, split FASTA files

---

## Example Data

Example inputs are available in:

- `example_data/`

Use them to test module behavior and output formats quickly.

---

## FAQ

**Q1: Which file should I run first?**
- Start with `phylotracer -h`, then run one module on `example_data/` to confirm environment setup.

**Q2: How do I see module-specific parameters?**
- Use `phylotracer <SUBCOMMAND> -h`.

**Q3: My plotting step fails on server/headless mode.**
- Set `export QT_QPA_PLATFORM=offscreen` before running visualization modules.

---

## Support

- GitHub: https://github.com/YiyongZhao/PhyloTracer
- PyPI: https://pypi.org/project/PhyloTracer
- Issue tracker: https://github.com/YiyongZhao/PhyloTracer/issues

---

## Citation

If you use PhyloTracer in published work, please cite the software and corresponding manuscript.

---

## License

PhyloTracer is licensed under MIT (`LICENSE`).
