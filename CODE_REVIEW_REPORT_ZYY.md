# PhyloTracer Code Review Report

**Date:** 2026-02-21
**Scope:** All Python source code in `bin/` and all documentation (`README.md`, `docs/`)

---

## Table of Contents

- [Part 1: Critical Errors (Runtime Crashes)](#part-1-critical-errors-runtime-crashes)
- [Part 2: Logic Bugs (Wrong Results)](#part-2-logic-bugs-wrong-results)
- [Part 3: Resource Leaks & Code Quality](#part-3-resource-leaks--code-quality)
- [Part 4: Documentation Errors](#part-4-documentation-errors)
- [Part 5: Priority Summary](#part-5-priority-summary)

---

## Part 1: Critical Errors (Runtime Crashes)

These issues **will cause the program to crash** when the affected code path is executed.

### C1. `GD_Detector.py:9` — Invalid import

```python
from calendar import c
```

`calendar` module has no export named `c`. This raises `ImportError` at import time. **The entire module fails to load.** Should be deleted.

---

### C2. `GD_Detector.py:349` — Unpacking mismatch

```python
gene_to_new_name, new_name_to_gene, voucher_to_taxa = gene_id_transfer("imap.txt")
```

`gene_id_transfer` returns a 4-tuple `(gene2new, new2gene, voucher2taxa, taxa2voucher)`, but only 3 values are unpacked. Raises `ValueError: too many values to unpack`.

**Fix:** Add a 4th variable: `gene_to_new_name, new_name_to_gene, voucher_to_taxa, taxa_to_voucher = gene_id_transfer("imap.txt")`

---

### C3. `Tree_Visualizer.py:1275` — Same unpacking error

```python
gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer("imap")
```

Same issue as C2. `gene_id_transfer` returns 4 values but only 3 are unpacked.

---

### C4. `Tree_Visualizer.py:645,713` — `RectFace` not imported

`RectFace` is used in `add_heat_map_to_node` and `add_color_bar` but is never imported (not in this file, not in `__init__.py`). Raises `NameError` whenever a heatmap DataFrame is provided.

**Fix:** Add `from ete3 import RectFace` to the imports.

---

### C5. `Tree_Visualizer.py:1287-1295` — `view_main` called with wrong arguments

The `__main__` block passes arguments in the wrong order to `view_main`. The `sptree` parameter is missing entirely, and `gene2new_named_gene_dic` is passed where `sptree` should be.

---

### C6. `TreeTopology_Summarizer.py:295-296,330` — Missing `math` and `PIL.Image` imports

```python
math.ceil(...)   # line 295 — math is never imported
math.sqrt(...)   # line 296
Image.open(...)  # line 330 — PIL.Image is never imported
```

Both will raise `NameError` at runtime.

**Fix:** Add `import math` and `from PIL import Image` at the top of the file.

---

### C7. `TreeTopology_Summarizer.py:424` — Wrong number of arguments to `statistical_main`

```python
statistical_main(tre_dic, outfile, gene2new_named_gene_dic, new_named_gene2gene_dic)
```

Function signature requires 5 parameters: `(tre_dic, outfile, gene2new_named_gene_dic, voucher2taxa_dic, top_n)`. Missing `top_n` and `new_named_gene2gene_dic` is passed where `voucher2taxa_dic` is expected.

---

### C8. `__init__.py:829` — `calculate_gd_num` missing required argument

```python
find_dup_node(tree)
```

`find_dup_node` requires both `gene_tree` and `species_tree` as positional arguments. Only one is provided. Raises `TypeError`.

---

### C9. `__init__.py:770` — `species_tree & "unknown"` raises TreeError

When `node.map` is set to `"unknown"` (line 451), the expression `species_tree & node.map` will raise an ete3 `TreeError` because no node named `"unknown"` exists in the species tree. No try/except guard.

---

### C10. `__init__.py:808` — `dup_map` can be `None`

```python
abs(clade_map.depth - dup_map.depth)
```

`map_species_set_to_node()` can return `None` (lines 703, 713). Accessing `.depth` on `None` raises `AttributeError`.

---

### C11. `Phylo_Rooter.py:380-381,415-416` — Assumes root has exactly 2 children

```python
tree.children[1]
tree.children[0]
```

If the root has 1 or 3+ children (e.g., unrooted trifurcating tree), this raises `IndexError`.

---

### C12. `OrthoFilter_Mono.py:850-852` — `os.makedirs(None)` when `visual=False`

When `visual=False`, `out_visual_dir` is `None`. The code calls `os.makedirs(None, exist_ok=True)` which raises `TypeError: expected str, bytes or os.PathLike object, not NoneType`.

**Fix:** Add a guard: `if dir_path is not None:` before `shutil.rmtree` and `os.makedirs`.

---

### C13. `OrthoFilter_LB.py:369-371` — Same `os.makedirs(None)` crash

Same pattern as C12. When `visual=False`, `pdf_dir` is `None`, and `os.makedirs(None)` raises `TypeError`.

---

### C14. `Ortho_Retriever.py:201,204` — `rm_dup` logic inversion

```python
if ev_seqs1 not in paralogs_L:
    paralogs_L.remove(ev_seqs1)
```

Condition is inverted: tries to remove an element only when it is **not** in the list, which always raises `ValueError`.

**Fix:** Change `not in` to `in`:
```python
if ev_seqs1 in paralogs_L:
    paralogs_L.remove(ev_seqs1)
```

---

### C15. `HaploFinder.py:804` — `set` is unhashable, raises `TypeError`

```python
if len(sps_tol) == 2 and {sp1, sp2} not in sps_tol:
```

`sps_tol` is a `set` of strings. `{sp1, sp2}` is a `set`. Checking `set in set` raises `TypeError: unhashable type: 'set'`.

**Fix:**
```python
if len(sps_tol) == 2 and not {sp1, sp2}.issubset(sps_tol):
```

---

### C16. `HaploFinder.py:1123,1134` — Undefined functions in `__main__`

`process_blastp_result`, `parse_synteny_file`, and `process_total_color_list` are called in `__main__` but are never defined anywhere. Raises `NameError`.

---

### C17. `HaploFinder.py:649` — Chromosome name without digits crashes

```python
chr1 = int("".join(filter(str.isdigit, chr_a)))
```

If chromosome name has no digits (e.g., `"chrX"`, `"scaffold_abc"`), the result is an empty string, and `int("")` raises `ValueError`.

---

### C18. `Hybrid_Visualizer.py:900` — `hyde_visual_leaf_main` missing argument

```python
hyde_visual_leaf_main(sptree)
```

Function requires `(out_file_name, sptree)` — two parameters. Only one is passed. Raises `TypeError`.

---

### C19. `GD_Loss_Tracker.py:823-828` — `__main__` uses undefined variables

`sptree_path` and `gf` are never defined. `get_path_str_num_dic` is called with wrong argument count. The entire `__main__` block crashes immediately with `NameError`.

---

## Part 2: Logic Bugs (Wrong Results)

These issues **will not crash** but **produce incorrect results**.

### B1. `__init__.py:582` — `judge_support` misclassifies integer support=1

```python
if support <= 1:
    support *= 100
```

A support value of `1` (meaning 1%) on a 0-100 scale is misinterpreted as `1.0` on a 0-1 scale and becomes `100`. This is indistinguishable from true 100% support.

---

### B2. `__init__.py:767` — Child support compared without normalization

```python
if any(c.support < clade_support for c in children)
```

Main node support is normalized via `judge_support()`, but child support values are compared raw. If scales differ (0-1 vs 0-100), the comparison is wrong.

---

### B3. `Phylo_Tracer.py:333-334` — Branch length thresholds decremented by 1

```python
absolute_branch_length = cli_args.absolute_branch_length - 1
relative_branch_length = cli_args.relative_branch_length - 1
```

- User input of minimum value `1` becomes `0` (disables filter)
- `relative_branch_length` input of `0.0` becomes `-1.0` (logically wrong)

---

### B4. `Phylo_Tracer.py:655` — `all(required_args)` treats 0 as missing

```python
all(required_args)
```

If user passes `--gd_support 0` (a valid value), `all()` treats `0` as falsy, rejecting the valid input as "missing argument".

---

### B5. `Phylo_Tracer.py:666` — `size=0` silently overridden

```python
size = cli_args.size if cli_args.size else 0.001
```

If user explicitly passes `--size 0`, the `0` is falsy so it gets overridden to `0.001`.

---

### B6. `Phylo_Rooter.py:538` — `norm()` returns scalar 0 instead of Series

```python
return 0
```

When all values are the same, should return `pd.Series(0.0, index=values.index)` to maintain consistent types in downstream arithmetic.

---

### B7. `GD_Loss_Tracker.py:354-357` — Operator precedence bug

```python
parts[-1] = taxa_name + "(" + count_part if count_part else ""
```

Python parses this as `(taxa_name + "(" + count_part) if count_part else ""`. When `count_part` is empty, `parts[-1]` becomes `""` instead of `taxa_name`.

**Fix:**
```python
parts[-1] = (taxa_name + "(" + count_part) if count_part else taxa_name
```

---

### B8. `Hybrid_Visualizer.py:135,165` — Mean divided by wrong denominator

`calculate_gamma` and `calculate_pvalue` skip `"nan"` and `"-inf"` entries but still divide by `len(lst)` (full length including skipped entries). This deflates the mean.

**Fix:** Divide by the count of valid (non-skipped) entries.

---

### B9. `Hybrid_Visualizer.py:381-391` — Colormap off-by-one

First loop generates indices 0-4999 (5000 entries). Second loop uses `range(5001, 9999)` generating indices 5001-9998 (4998 entries). Index 5000 is never populated, creating a gap in the colormap.

---

### B10. `SyntenyRate.py:113-114` — Results returned in wrong order

```python
as_completed(futures)
```

`as_completed` returns results in completion order, not submission order. The rates are then misattributed to wrong trees when printing `tree_1`, `tree_2`, etc.

**Fix:** Use `futures` list directly or map results back by key.

---

### B11. `OrthoFilter_LB.py:288` — `get_sisters()[0]` can IndexError

```python
sister = leaf.get_sisters()[0]
```

If a leaf has no sisters, `get_sisters()` returns empty list and `[0]` raises `IndexError`.

---

### B12. `OrthoFilter_LB.py:373,419` — Progress bar double-counted

`tqdm` auto-increments on iteration, then `pbar.update(1)` manually increments again, causing progress to count double.

---

### B13. `Tree_Visualizer.py:1180,1286` — `keep_branch` type mismatch

Function checks `if keep_branch != "1":` but `__main__` passes integer `1`. Since `1 != "1"` is always `True`, `realign_branch_length` is always called.

---

### B14. `HaploFinder.py:1121` — `size` passed as string to matplotlib

`sys.argv[14]` is a string. `mpatches.Circle` expects a float for `radius`, causing `TypeError`.

---

### B15. `Hybrid_Visualizer.py:575` — `gamma_df` references wrong variable

`gamma_df` refers to the last iteration's value, not the averaged `gamma_result_df`. The colorbar is based on a single species' data instead of the aggregate.

---

### B16. `PhyloTree_CollapseExpand.py:44-45` — May delete root node

If root node support < threshold, `node.delete()` is called on the root, which is undefined behavior in ETE3.

---

### B17. `GD_Visualizer.py:171-173` — Species tree leaf names mutated in-place

`mark_sptree` permanently renames species tree leaf names. If the tree object is reused downstream, lookups by original voucher codes will fail.

---

### B18. `Hybrid_Tracer.py:905` — PHYLIP name truncation risk

```python
name = species[:10].ljust(10)
```

Two species sharing the first 10 characters get identical names in the PHYLIP file, corrupting HyDe input.

---

## Part 3: Resource Leaks & Code Quality

### R1. File handles not using `with` statement

| File | Line | Variable |
|------|------|----------|
| `OrthoFilter_Mono.py` | 866 | `log = open(...)` — not closed on exception |
| `Ortho_Retriever.py` | 446 | `o = open(...)` — no `with` |
| `GD_Loss_Tracker.py` | 602 | `out = open(...)` — no `with` |
| `HaploFinder.py` | 49 | `f = open(fn)` — never closed |

---

### R2. `GD_Loss_Tracker.py:80` — `set | None` syntax requires Python 3.10+

```python
family_species_set: set | None = None
```

Should use `Optional[set]` from `typing` for Python 3.8/3.9 compatibility.

---

### R3. `HaploFinder.py:51` — `dict` built-in shadowed

```python
data, dict = [], {}
```

Shadows Python's built-in `dict` type within the function scope.

---

### R4. `Hybrid_Tracer.py:751-782` — Race condition with temp files

`run_hyde_from_matrix_integrated` writes to hardcoded `temp.phy` and `temp.imap`. Concurrent calls overwrite each other.

---

### R5. All files — Fragile `from __init__ import *`

All modules use `from __init__ import *` which depends on `bin/` being on `sys.path`. This is fragile and pollutes the namespace.

---

### R6. `Phylo_Tracer.py:33-55,264` — Module-level side effects

The ASCII banner prints at import time, and `parser.parse_args()` runs at module load. Importing this module from elsewhere will crash.

---

### R7. `Phylo_Rooter.py:589` — Silent deletion of previous output

`shutil.rmtree(dir_path)` deletes the entire `rooted_trees/` directory without confirmation. Re-running destroys all previous results.

---

### R8. `__init__.py:261` — Hardcoded output path in utility function

`sptree.write(outfile="numed_sptree.nwk")` writes to current working directory with a hardcoded filename, and `Phylo_Tracer.py:462-463` writes the same file again (redundant).

---

## Part 4: Documentation Errors

### Fatal Documentation Errors (Users will fail)

| # | File | Line | Issue |
|---|------|------|-------|
| D1 | `docs/installation.rst` | 46 | **`pip install pPhyloTracer`** — extra `p`. Should be `pip install PhyloTracer` |
| D2 | `docs/installation.rst` | 60 | **`bash install_package.sh`** — missing `s`. Should be `install_packages.sh` |
| D3 | `docs/installation.rst` | 33 | **conda install uses commas** (invalid syntax). Should use spaces. Also lists `time` (stdlib, not installable) and `HyDe` (not on conda) |
| D4 | `docs/analyze.rst` | 109-122 | **OrthoFilter_Mono documents nonexistent parameters**: `--branch_length_multiples`, `--insert_branch_index` |
| D5 | `docs/analyze.rst` | 215-223 | **GD_Loss_Tracker documents nonexistent parameters**: `--all`, `--start_node`, `--end_species` |
| D6 | `docs/analyze.rst` | 271 | **Hybrid_Tracer documents `--target_node`** — does not exist. Actual parameter is `--mrca_node` |
| D7 | `docs/analyze.rst` | all examples | **Command invocation uses `Phylo_Tracer.py`** — should be `PhyloTracer` |
| D8 | `docs/input_files.rst` | 8,29,45,59,75,91,107,123 | **All GitHub links point to `/examples/`** — actual directory is `/example_data/`. All links 404 |
| D9 | `docs/analyze.rst` | 195-201 | **GD_Visualizer missing required `--input_imap`** parameter |
| D10 | `docs/analyze.rst` | 232-238 | **GD_Loss_Visualizer missing required `--gd_loss_result`** parameter |

### Undocumented Parameters (code has, docs don't)

| Parameter | Module | Notes |
|-----------|--------|-------|
| `--target_species`, `--mrca_node`, `--include_unobserved_species` | GD_Loss_Tracker | Missing from both README and analyze.rst |
| `--visual_top` | TreeTopology_Summarizer | Missing from both README and analyze.rst |
| `--gdtype_mode` (choices: relaxed/strict) | GD_Detector | Missing from both README and analyze.rst |
| `--mode` (haplofinder/split), `--pair_support`, `--input_sps_tree` | HaploFinder | Missing from analyze.rst |

### Other Documentation Issues

| # | File | Issue |
|---|------|-------|
| D11 | `README.md:207` | Python version listed as `3.0+` — actual minimum is 3.8+ given dependencies |
| D12 | `README.md:28` | Travis CI badge uses `master` branch (current is `main`); `travis-ci.org` is deprecated |
| D13 | `README.md:26` | Documentation badge links to `hybridization-detection.readthedocs.io` (wrong project) |
| D14 | `README.md:19` | `Licence` should be `License` (for US English consistency) |
| D15 | `docs/installation.rst:20` | Recommends Python 3.6 (EOL); Miniconda link `repo.continuum.io` is outdated |
| D16 | `docs/installation.rst:65` vs `README.md:99` | Contradictory QT platform: `linuxfb` vs `offscreen` |
| D17 | `docs/analyze.rst:129` | Typo: `visualiz` should be `visualize` |
| D18 | `docs/analyze.rst:150-151` | `Tree_Visualizer`: `--keep_branch` and `--tree_style` listed as Required, but code has them as optional |
| D19 | `docs/analyze.rst:177` | `--dup_species_proportion` shown with default 0.2, but code has it as `required=True` (no default) |
| D20 | `docs/analyze.rst:313` | `--gd_support` listed as Required, but code has `required=False` with default=50 |
| D21 | `docs/index.rst:17` | Outdated terminology: "inserted genes" — current code uses "purity", "alien lineage" |
| D22 | `README.md:21` | Contact email inconsistent with `Phylo_Tracer.py` banner |

---

## Part 5: Priority Summary

### Highest Priority (Blocks Users)

1. **C1** — Delete `from calendar import c` in `GD_Detector.py`
2. **C2, C3** — Fix all `gene_id_transfer` unpacking to 4 values
3. **C6** — Add missing `math` and `PIL.Image` imports in `TreeTopology_Summarizer.py`
4. **C12, C13** — Guard against `None` path in `OrthoFilter_Mono.py` / `OrthoFilter_LB.py`
5. **C14** — Fix `rm_dup` condition inversion in `Ortho_Retriever.py`
6. **C15** — Fix set unhashable check in `HaploFinder.py`
7. **D1** — Fix `pPhyloTracer` typo in `docs/installation.rst`
8. **D2** — Fix `install_package.sh` filename in `docs/installation.rst`

### High Priority (Wrong Results)

9. **B10** — Fix `as_completed` ordering in `SyntenyRate.py`
10. **B1, B2** — Fix `judge_support` normalization logic
11. **B3** — Remove branch length threshold decrement in `Phylo_Tracer.py`
12. **B8** — Fix gamma/pvalue mean denominator in `Hybrid_Visualizer.py`
13. **B7** — Fix operator precedence in `GD_Loss_Tracker.py`

### Medium Priority (Documentation)

14. **D4-D6** — Update stale parameters in `docs/analyze.rst`
15. **D7** — Change all `Phylo_Tracer.py` to `PhyloTracer` in examples
16. **D8** — Fix `/examples/` to `/example_data/` in `docs/input_files.rst`
17. **D9, D10** — Add missing required parameters to docs

---

*Report generated by Claude Code on 2026-02-21*
