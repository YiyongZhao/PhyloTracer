# PhyloTracer 17-Module Code & Docs Review Report
Generated: Sun Mar  8 05:10:59 PM EDT 2026

==================================================================
## Module: BranchLength_NumericConverter
==================================================================
Now I have all the information needed for a comprehensive review. Here are my findings:

---

# BranchLength_NumericConverter — Code Review Report

## A. Scientific Correctness

**No issues found.** The module performs a straightforward numeric formatting task (converting branch-length strings to fixed decimal precision). The default of 10 decimal places is reasonable for preserving precision in phylogenetic branch lengths (which are typically small floating-point values from ML or Bayesian inference). The description "string to numerical format" is slightly misleading (see B below) but not scientifically incorrect — it's really a precision-formatting operation.

---

## B. Code vs Docs Consistency

### Issue 1
**[MEDIUM] [Code vs Docs] `Phylo_Tracer.py:390` vs `BranchLength_NumericConverter.py:85–87`: `--output_dir` argument exists in CLI but is never passed to the main function.**

The CLI parser defines `--output_dir` (line 160), and the README/docs both document it. However, `handle_branch_length_numeric_converter()` at line 390 calls:
```python
branch_length_numeric_converter_main(tre_dic, cli_args.decimal_place)
```
The `branch_length_numeric_converter_main` function signature has no `output_dir` parameter — it hardcodes the output directory name as `"converter_tree"` relative to `os.getcwd()`. The `--output_dir` argument is silently ignored.

**Suggested fix:** Add an `output_dir` parameter to `branch_length_numeric_converter_main()` and pass `cli_args.output_dir` from the handler. Or check how other modules handle this (likely via `os.chdir` before calling main — check the dispatcher).

### Issue 2
**[LOW] [Code vs Docs] `docs/analyze.rst:59` vs `Phylo_Tracer.py:159`: Parameter name inconsistency — `--decimal_place` (singular) in code/CLI vs `--decimal_place` in docs.**

Both docs and code use `--decimal_place` (singular), while the internal Python variable is `decimal_places` (plural). This is cosmetically inconsistent but functionally harmless since the CLI name matches. **No action needed**, but noted for completeness.

### Issue 3
**[LOW] [Code vs Docs] `docs/analyze.rst:53` description: "from string to numerical format" is misleading.**

The docs and README say the module converts "branch length values of a phylogenetic tree from string to numerical format." In reality, ete3 already parses branch lengths as floats internally. The module reformats them to a fixed number of decimal places. The description should say "standardize branch-length precision" or "format branch lengths to fixed decimal places."

**Suggested fix:** Update the description in `docs/analyze.rst` and `README.md` to: "To standardize branch-length precision in phylogenetic trees by formatting branch lengths to a fixed number of decimal places."

---

## C. Code Logic and Correctness

### Issue 4
**[CRITICAL] [Code Logic] `BranchLength_NumericConverter.py:116–118`: `shutil.rmtree()` silently destroys previous output without user confirmation.**

```python
if dir_path != cwd and os.path.exists(dir_path):
    shutil.rmtree(dir_path)
```

If the `converter_tree/` directory already exists (e.g., from a previous run or containing user files), it is completely deleted without warning. This is destructive and could result in data loss.

**Suggested fix:** Either (a) warn and skip existing files, (b) prompt the user, or (c) at minimum, log a warning before deletion. Many phylogenomics tools append a timestamp or use `_backup` directories.

### Issue 5
**[MEDIUM] [Code Logic] `BranchLength_NumericConverter.py:109–115`: Hardcoded output directory name with fragile `basename` check.**

```python
cwd = os.getcwd()
default_dir = "converter_tree"
dir_path = (
    cwd
    if os.path.basename(os.path.normpath(cwd)) == default_dir
    else os.path.join(cwd, f"{default_dir}/")
)
```

If the user happens to be in a directory named `converter_tree`, the module writes output directly into the cwd instead of creating a subdirectory. This is a fragile heuristic — it could lead to unexpected behavior if the user's working directory coincidentally has this name. Combined with Issue 1 (output_dir being ignored), the user has no way to control output location.

**Suggested fix:** Remove the basename heuristic. Always use the output directory from the CLI argument (once Issue 1 is fixed).

### Issue 6
**[MEDIUM] [Code Logic] `BranchLength_NumericConverter.py:123–136`: Empty `tree_dict` produces no output and no warning.**

If `tree_dict` is an empty dict `{}`, the progress bar shows `0/0`, the loop body never executes, and the function silently returns. There's no indication to the user that nothing was processed.

**Suggested fix:** Add a check:
```python
if not tree_dict:
    logger.warning("No trees to process (empty input).")
    return
```

### Issue 7
**[LOW] [Code Logic] `BranchLength_NumericConverter.py:128–131`: Redundant `None` check for `decimal_places`.**

```python
if decimal_places is None:
    tree_str = trans_branch_length(tree)
else:
    tree_str = trans_branch_length(tree, decimal_places)
```

The function signature says `decimal_places: int = 10` and the validation at lines 104–107 explicitly allows `None`. But `trans_branch_length()` at line 41 will reject `None` with `ValueError` because `isinstance(None, int)` is `False`. So passing `None` to the main function would bypass validation (line 104: `if decimal_places is not None and ...`) but then silently default to 10 in `trans_branch_length`. This is inconsistent — either `None` should be a valid sentinel or it should be rejected at the top.

**Suggested fix:** Remove the `None` special-casing. The CLI always provides an int (default 10). Simplify to:
```python
tree_str = trans_branch_length(tree, decimal_places)
```

### Issue 8
**[LOW] [Code Logic] `Phylo_Tracer.py:387`: Redundant `if cli_args.input_GF_list` check.**

`--input_GF_list` is defined as `required=True` in the argparse parser (line 158). If it's missing, argparse will error before `handle_branch_length_numeric_converter` is ever called. The `if` check is dead code.

**Suggested fix:** Remove the `if/else` and just execute the body directly (or keep it as defensive programming — low priority).

---

## D. Output Descriptions

### Issue 9
**[MEDIUM] [Output Descriptions] `docs/analyze.rst:49–64` and `README.md:406–417`: No output description for BranchLength_NumericConverter.**

Neither the docs nor README describe what output files are generated. The module creates a `converter_tree/` directory containing `<tree_id>.nwk` files. This should be documented.

**Suggested fix:** Add an Output section:
```
Output:
    - ``converter_tree/``   Directory containing reformatted Newick tree files (<tree_id>.nwk).
```

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 1     |
| MEDIUM     | 4     |
| LOW        | 4     |
| SUGGESTION | 0     |
| **Total**  | **9** |

### Priority Fixes
1. **CRITICAL (Issue 4):** `shutil.rmtree()` silently deletes existing output — add warning or safety mechanism.
2. **MEDIUM (Issue 1):** `--output_dir` CLI arg is defined but never used — wire it through to the main function.
3. **MEDIUM (Issue 5):** Fragile basename heuristic for output directory — remove once output_dir is plumbed.
4. **MEDIUM (Issue 6):** Empty input dict silently succeeds — add a warning.
5. **MEDIUM (Issue 9):** Missing output documentation — add output section to docs/README.

==================================================================
## Module: GD_Detector
==================================================================
Good — `.depth` is set by `annotate_gene_tree`. Now I have all the information for a complete review.

---

# GD_Detector Module Review

## A. Scientific Correctness

**[LOW] [SCIENTIFIC] [GD_Detector.py:normalize_model, line 230-255]: Docstring lists incorrect canonical type names.**
The docstring says returns "ABAB, AAB, ABB, or Complex" but the function actually returns "AABB", "AXBB", "AABX", or "Complex". These are different GD topology patterns.
**Fix:** Update docstring to: `"AABB, AXBB, AABX, or Complex"`.

**[LOW] [SCIENTIFIC] [GD_Detector.py:get_model, line 264-296]: Docstring lists incorrect type names.**
Docstring says "ABAB / ABB / AAB / OTHER" but function returns 4-character strings like "AABB", "AXBB", "XABB", etc.
**Fix:** Update docstring to match actual return values.

**[SUGGESTION] [SCIENTIFIC] [GD_Detector.py:get_model, line 290-294]: GD type encoding could be clearer.**
The 4-character encoding `[left∩A][right∩A][left∩B][right∩B]` is slot-ordered, not child-ordered. This is a valid design but is not documented anywhere, making the encoding scheme hard for users to interpret. A brief note in the docs would help.

**No issues found** with parameter defaults from a scientific standpoint — support thresholds of 50, dup_species_num of 2, and dup_species_proportion of 0.2 are reasonable and commonly used in gene-tree/species-tree reconciliation literature.

---

## B. Code vs Docs Consistency

**[CRITICAL] [CONSISTENCY] [Phylo_Tracer.py:256 vs analyze.rst:231]: `--subclade_support` default mismatch.**
- Code (Phylo_Tracer.py line 256): `default=0`
- analyze.rst (line 231): "default is 50"
- README (line 613): `default = 0`

The docs/analyze.rst states a default of 50, but the actual CLI parser uses 0. Note: since `required=True`, the default is never actually used — but the documented value is still misleading.
**Fix:** Update analyze.rst to say `default is 0` to match the code.

**[CRITICAL] [CONSISTENCY] [GD_Detector.py:__main__ vs Phylo_Tracer.py]: `__main__` block uses different parameter names and defaults than the primary CLI.**
| Parameter | `__main__` (GD_Detector.py) | CLI (Phylo_Tracer.py) |
|---|---|---|
| subclade support | `--clade_support`, default=50 | `--subclade_support`, default=0 |
| dup species ratio | `--dup_species_percent`, default=0.5 | `--dup_species_proportion`, default=0.2 |
| depth variance | `--max_topology_distance`, default=0 | `--deepvar`, default=1 |

This means running `python GD_Detector.py` directly produces different behavior than `PhyloTracer GD_Detector`.
**Fix:** Align `__main__` parameter names and defaults with Phylo_Tracer.py, or remove the `__main__` block if it's not intended for direct use.

**[CRITICAL] [CONSISTENCY] [GD_Detector.py:__main__ line 378 vs Phylo_Tracer.py:670]: Species tree renaming uses opposite direction.**
- `__main__` (line 378): `rename_input_tre(species_tree, voucher_to_taxa)` — renames voucher→taxa
- Phylo_Tracer.py (line 670): `rename_input_tre(sptree, taxa2voucher_dic)` — renames taxa→voucher

Also, numbering happens *after* renaming in `__main__` but *before* renaming in Phylo_Tracer.py. One of these paths will produce incorrect reconciliation results.
**Fix:** Ensure both entry points use the same renaming direction and numbering order.

**[MEDIUM] [CONSISTENCY] [analyze.rst:233 vs Phylo_Tracer.py:258]: `--dup_species_num` default not shown in analyze.rst.**
analyze.rst says "Number of species duplications under the GD node" with no default listed. Code has `default=2`.
**Fix:** Add "(default is 2)" to analyze.rst.

**[MEDIUM] [CONSISTENCY] [analyze.rst:221-243]: `--output_dir` parameter missing from analyze.rst.**
The `--output_dir` parameter exists in Phylo_Tracer.py (line 262) and README (line 620), but is not documented in analyze.rst.
**Fix:** Add `--output_dir` to the optional parameters in analyze.rst.

**[MEDIUM] [CONSISTENCY] [Phylo_Tracer.py:255-260]: All "required" parameters have defaults that are never used.**
`--gd_support`, `--subclade_support`, `--dup_species_proportion`, `--dup_species_num`, and `--deepvar` are all set to `required=True` with `default=...`. In argparse, `required=True` means the default is never used — the user must always provide the value. Yet the docs present these as having defaults.
**Fix:** Either remove `required=True` (making them optional with real defaults) or remove the `default=` from the parser (and update docs to say "required, no default"). Removing `required=True` is recommended since these parameters all have scientifically reasonable defaults.

**[MEDIUM] [CONSISTENCY] [analyze.rst:232 vs Phylo_Tracer.py:257]: `--dup_species_proportion` described as "Proportion of overlapping species" in analyze.rst but code help text says "Minimum overlap ratio".**
Minor wording difference, but more importantly the function parameter is named `duplicated_species_percentage_threshold` which suggests a percentage (0-100), while the actual accepted range is 0-1 (a proportion). The internal naming is misleading.
**Fix:** Rename the function parameter to `duplicated_species_proportion_threshold` for clarity.

---

## C. Code Logic and Correctness

**[CRITICAL] [LOGIC] [GD_Detector.py:221]: Sort assumes specific label format, will crash on non-conforming labels.**
```python
df.sort_values("Newick_label", key=lambda x: x.str[1:].astype(int))
```
This assumes every `Newick_label` starts with exactly one character prefix (e.g., "N0", "S3") followed by a pure integer. If `voucher_to_taxa` maps produce labels like "Angiosperm", "Root", or multi-character prefixes, this will crash with a `ValueError`.
**Fix:** Add a robust sort key with a fallback, e.g.:
```python
import re
def sort_key(s):
    m = re.match(r'([A-Za-z]*)(\d+)', s)
    return (m.group(1), int(m.group(2))) if m else (s, 0)
df = df.sort_values("Newick_label", key=lambda x: x.map(sort_key))
```

**[MEDIUM] [LOGIC] [GD_Detector.py:222]: `gd_type_*.tsv` written to CWD, ignoring `output_dir`.**
The `write_gene_duplication_results` function doesn't receive an `output_dir` parameter. The main output file path is passed as `output_file`, but the secondary output `gd_type_{gdtype_mode}.tsv` (line 222) is always written to the current working directory, not alongside the main output or in `output_dir`.
**Fix:** Either pass `output_dir` into the function or derive the directory from `output_file`.

**[MEDIUM] [LOGIC] [GD_Detector.py:get_model_strict, lines 325-329]: `order_children_by_name` on gene tree nodes assumes binary grandchildren.**
The function calls `order_children_by_name` on both children of the duplication node, then tries to get grandchildren. If a child is a leaf (has 0 children), `order_children_by_name` will fail at `c1, c2 = n.get_children()` with a `ValueError` for unpacking.
The check `len(left_child.get_children()) != 2` on line 326 guards against this, but `order_children_by_name(clade)` on line 325 is called *before* the binary check — if the *duplication node itself* has != 2 children (e.g., multifurcation), line 325 will crash.
**Fix:** Add a binary check for `clade` before calling `order_children_by_name(clade)`.

**[MEDIUM] [LOGIC] [GD_Detector.py:order_children_by_name, line 260]: Missing `.num` attribute will cause AttributeError.**
If a node doesn't have the `.num` attribute (e.g., tree wasn't numbered), `x.num` raises `AttributeError`. The fallback `10**18` is only reached if `m` is falsy, not if the attribute is missing.
**Fix:** Use `getattr(x, 'num', None)` instead of `x.num`.

**[LOW] [LOGIC] [GD_Detector.py:91]: Typo "deepth" in log message.**
`logger.info("Maximum variance of deepth: %s", max_topology_distance)` — should be "depth".

**[LOW] [LOGIC] [GD_Detector.py:506 in __init__.py]: `annotate_gene_tree` sets `depth=None` on exception.**
If species mapping fails for a gene tree node, `depth` is set to `None`. Later, `get_model_strict` does `abs(left_a.depth - map_node_a.depth)` which will crash with `TypeError: unsupported operand type(s) for -: 'NoneType' and 'int'`.
**Fix:** Either filter out nodes with `depth=None` before calling `get_model_strict`, or return "XXXX" early when any depth is `None`.

---

## D. Output Descriptions

**[MEDIUM] [DOCS] [analyze.rst:221-243]: No output file descriptions for GD_Detector.**
analyze.rst describes inputs and usage examples but does not describe any output files. The module produces:
- `gd_result_{mode}.txt` — per-gene-pair duplication table with 11 columns
- `gd_type_{mode}.tsv` — per-node GD type summary with 8 columns
- `numed_sptree.nwk` — numbered species tree (written by the CLI handler)

**Fix:** Add an "Output:" section listing these files and their column descriptions.

**[MEDIUM] [DOCS] [README:605-623]: README GD_Detector section doesn't describe output columns.**
The README lists output file patterns in the "PhyloTracer Results Files" section (`gd_result_*.txt`, `gd_type_*.tsv`) but doesn't describe column names for either file.
**Fix:** Add column descriptions. The main output columns are: `#tree_ID`, `gd_id`, `gd_support`, `gene1`, `gene2`, `level`, `species`, `GD_dup_sps`, `dup_ratio`, `gd_type`, `comment`. The summary columns are: `Newick_label`, `GD`, `NUM`, `GD_ratio`, `AABB`, `AXBB`, `AABX`, `Complex`.

**[LOW] [DOCS] [GD_Detector.py:104-107]: Output header uses `#tree_ID` (with `#` prefix).**
The `#` prefix may cause downstream tools or `pd.read_csv` to skip it as a comment line, or misread the column name as `#tree_ID` instead of `tree_ID`. This is a common pitfall in bioinformatics TSV files.
**Fix:** Either remove the `#` prefix or document that the header line is comment-prefixed.

---

## Summary Table

| Severity | Count |
|----------|-------|
| CRITICAL | 3 |
| MEDIUM | 7 |
| LOW | 4 |
| SUGGESTION | 1 |
| **Total** | **15** |

The three **CRITICAL** issues are:
1. `--subclade_support` default mismatch between analyze.rst (50) and code (0)
2. `__main__` block uses different parameter names, defaults, and renaming direction than the primary CLI — running `python GD_Detector.py` directly will produce different (likely incorrect) results
3. Sort on `Newick_label` will crash on labels that don't match the `[single-char][integer]` pattern

==================================================================
## Module: GD_Loss_Tracker
==================================================================
Now I have all the information needed. Let me compile the review.

---

# GD_Loss_Tracker — Code Review Report

## A. Scientific Correctness

**1. [LOW] [SCIENTIFIC] [GD_Loss_Tracker.py:classify_species_copy_state]:** The 2-2/2-1/2-0 classification system is a reasonable simplification for post-WGD copy retention analysis, but the docstring doesn't clarify that this is a *presence/absence* model, not a true copy-number model. A species with 3 copies on one side is still classified as "2-2" if present on both sides. This is biologically correct for asking "was at least one copy retained?", but could mislead users expecting actual copy counts.
- **Suggested fix:** Add a note in the docstring: "loss_type reflects presence/absence of the species in each child lineage (not raw copy counts)."

**2. [MEDIUM] [SCIENTIFIC] [GD_Loss_Tracker.py:355]:** `find_dup_node` is called with hardcoded parameters `(genetree, renamed_sptree, 50, 0, 2, 0, 1)` — these are GD detection thresholds (gd_support=50, subclade_support=0, dup_species_num=2, dup_species_proportion=0, deepvar=1). These are **not exposed as CLI arguments** for GD_Loss_Tracker, meaning users cannot adjust GD detection sensitivity for loss tracking independently of GD_Detector. The hardcoded `dup_species_num=2` is particularly restrictive — it requires at least 2 overlapping species, which will miss duplication events in small clades.
- **Suggested fix:** Either expose these as optional CLI arguments (with documented defaults), or document that GD_Loss_Tracker uses fixed GD detection thresholds.

**3. [LOW] [SCIENTIFIC] [GD_Loss_Tracker.py:get_maptree_node_count_dic]:** The logic at lines 293-307 uses `startswith("S")` to distinguish internal nodes from tips. This assumes species tree internal nodes are named with "S" prefix, which is convention-dependent and fragile.
- **Suggested fix:** Document this naming convention requirement, or use a more robust check (e.g., `node.is_leaf()`).

## B. Code vs Docs Consistency

**4. [MEDIUM] [DOCS-CODE] [docs/analyze.rst:263-283 vs Phylo_Tracer.py:283]:** The `--node_count_mode` parameter exists in the CLI parser (line 283) and README (line 649) but is **missing from docs/analyze.rst** (section 11).
- **Suggested fix:** Add `--node_count_mode` to the Optional Parameters in analyze.rst section 11:
  ```
  - ``--node_count_mode``  Node counting mode for path_count_* transition statistics: nonaccumulate (default) or accumulate.
  ```

**5. [MEDIUM] [DOCS-CODE] [docs/analyze.rst:277 vs Phylo_Tracer.py:282]:** The docs/analyze.rst description of `--include_unobserved_species` says "species unobserved in a gene family are still classified by left/right presence instead of labeled as missing_data" — this is accurate but incomplete compared to the CLI help text which adds: "This flag does not change loss_path copy-state values (0/1/2) or GD event detection." The README version (line 648) is similarly more detailed. Minor inconsistency in level of detail.
- **Suggested fix:** Harmonize the description across all three locations.

**6. [LOW] [DOCS-CODE] [README:262 vs code]:** README lists output files as `gd_loss_summary.txt`, `gd_loss_count_summary.txt`, `gd_loss.xlsx` — this matches the code (line 666: `gd_loss_summary.txt`, line 723: `gd_loss_count_summary.txt`, line 763: default `gd_loss.xlsx`). **No issues found.**

**7. [LOW] [DOCS-CODE] [docs/analyze.rst:275 vs Phylo_Tracer.py:280]:** The docs say `--target_species` "Only count loss paths ending in this species" — but the code actually filters the *summary records* by species_voucher and *path endpoints* by last_node. The description is functionally correct. **No issues found.**

## C. Code Logic and Correctness

**8. [CRITICAL] [LOGIC] [GD_Loss_Tracker.py:summarize_small_loss_types_from_path, line 200 vs line 163]:** The function docstring declares return type as `tuple` with 4 elements `(path_count_types, c20, c21, c10)`, but the actual return at line 200 returns **5 elements**: `(path_count_types, c20, c21, c10, transition_keys)`. The early-return at line 167 returns only 4 elements `("NA", 0, 0, 0)` and the return at line 176 also returns 4 elements. This means:
- When `path_str` is empty/NA or has < 2 parsed steps, the caller at line 439 does `path_count_types, path_count_2_0, path_count_2_1, path_count_1_0, transition_keys = summarize_small_loss_types_from_path(...)` which will raise a `ValueError: not enough values to unpack`.
- **Suggested fix:** Change lines 167 and 176 to return 5 elements: `return "NA", 0, 0, 0, []`

**9. [MEDIUM] [LOGIC] [GD_Loss_Tracker.py:get_two_nodes_path_str, lines 265-268]:** If `start_node` is not an ancestor/descendant of `end_node`, the `while` loop will traverse to the root via `.up` and exit when `current_node is None`, but then line 270 still appends `end_node.name`. This means the returned path is incorrect (it goes from `start_node` up to root, then appends an unrelated `end_node`). The function assumes a valid ancestor-descendant relationship but doesn't validate it.
- **Suggested fix:** Add a check: if the loop exits via `current_node is None`, log a warning and return an empty list or raise an error.

**10. [MEDIUM] [LOGIC] [GD_Loss_Tracker.py:666]:** Output file `gd_loss_summary.txt` is written to the **current working directory** unconditionally, ignoring `--output_dir`. The CLI parser accepts `--output_dir` (line 284) and the handler function receives it, but `get_path_str_num_dic` doesn't take an output directory parameter — it hardcodes `open("gd_loss_summary.txt", "w")` (line 666) and `open("gd_loss_count_summary.txt", "w")` (line 723).
- **Suggested fix:** Pass `output_dir` into `get_path_str_num_dic` and use `os.path.join(output_dir, "gd_loss_summary.txt")` for file paths. Same for `parse_text_to_excel`.

**11. [MEDIUM] [LOGIC] [GD_Loss_Tracker.py:625]:** `PhyloTree(tre_path)` is called without error handling. If a gene tree file is malformed or missing, this will crash the entire run instead of skipping the problematic tree.
- **Suggested fix:** Wrap in try/except, log a warning, and `continue`.

**12. [LOW] [LOGIC] [GD_Loss_Tracker.py:367-388]:** When `len(sp) == 1` (single-species duplication node), the code skips the path computation block (line 367: `if len(sp) > 1`), so `voucher_to_pretty_path` remains empty. This means single-species GD events will always have `loss_path = "NA"` in the output. This is arguably correct behavior, but worth documenting.

**13. [LOW] [LOGIC] [GD_Loss_Tracker.py:276]:** Typo in parameter name: `taget_node` should be `target_node`.

**14. [LOW] [LOGIC] [GD_Loss_Tracker.py:parse_text_to_excel, line 797]:** The first line skip (`first_line = True`) assumes the first line is a header (`"GD Loss path\tGF count\n"` written at line 724). But the file actually starts with a header line, then the first data group starts with a blank line (line 735: `f.write(f"\n{new_k}\t{v}\n")`). The blank-line skip at line 802 handles this correctly, but the first-line skip is actually skipping the header — which is correct. **No issue.**

## D. Output Descriptions

**15. [MEDIUM] [OUTPUT-DOCS] [README:637-652 and docs/analyze.rst:263-283]:** Neither the README nor analyze.rst documents the **column names** of the output files. The `gd_loss_summary.txt` header (line 668-671) contains 20 columns:
```
tree_ID, gd_ID, gd_support, level, species, gene1, gene2, loss_path,
left_gene_n, right_gene_n, loss_type, left_has, right_has,
loss_confidence, major_loss_class, path_count_node_events,
path_count_types, path_count_2_0, path_count_2_1, path_count_1_0
```
And `gd_loss_count_summary.txt` has columns `GD Loss path` and `GF count`.
- **Suggested fix:** Add output column descriptions to both docs/analyze.rst and README, similar to how MulRF_Distance documents its output columns.

**16. [LOW] [OUTPUT-DOCS] [README:262]:** The README lists `gd_loss.xlsx` as an output but doesn't explain it contains per-species-pair sheets with loss path breakdowns. Consider adding a one-line description.

---

## Summary Table

| Severity | Count |
|----------|-------|
| CRITICAL | 1 |
| MEDIUM | 6 |
| LOW | 6 |
| SUGGESTION | 0 |
| **Total** | **13** |

### Priority action items:
1. **CRITICAL (#8):** Fix the tuple unpacking mismatch in `summarize_small_loss_types_from_path` — this is a runtime crash bug triggered whenever a species has an empty or trivial loss path.
2. **MEDIUM (#10):** Wire `--output_dir` through to file writes — currently the flag is accepted but silently ignored.
3. **MEDIUM (#4):** Add `--node_count_mode` to docs/analyze.rst.
4. **MEDIUM (#2):** Document or expose the hardcoded GD detection thresholds used internally.

==================================================================
## Module: GD_Loss_Visualizer
==================================================================
Here is the full review of the **GD_Loss_Visualizer** module.

---

## A. Scientific Correctness

**[MEDIUM] [SCIENTIFIC] [GD_Loss_Visualizer.py:identify_loss_detail, line 78-79]:**
The fallback classification for copy-number transitions not explicitly handled (e.g., 3→1, 4→2) uses `diff >= 2 → "2-0"`, `diff == 1 → "2-1"`. This is a lossy heuristic: a 3→1 transition is not the same as a 2→0 (complete loss). It conflates magnitude with type. Consider either logging a warning for unexpected transitions or introducing a generic "multi-loss" category.

**[LOW] [SCIENTIFIC] [GD_Loss_Visualizer.py:calculate_node_loss_score, line 199]:**
The weighting `raw = 1.0*(2→0) + 0.5*(2→1)` assigns a fixed 2:1 severity ratio between complete and partial loss. While reasonable, this is an arbitrary choice not backed by a citation. The docstring should note this is a heuristic weight and cite any justification if one exists, or state it's a user-chosen convention.

**No other scientific issues found.**

---

## B. Code vs Docs Consistency

**[CRITICAL] [CONSISTENCY] [Phylo_Tracer.py:line 289 vs docs/analyze.rst:line 293 and README "GD_Loss_Visualizer" section]:**
`--input_sps_tree` is defined as `required=False` in code, but documented as a **Required Parameter** in both analyze.rst and README. If the user omits it, the handler silently logs an error and exits without processing. Fix: set `required=True` in the argparse definition.

**[MEDIUM] [CONSISTENCY] [Phylo_Tracer.py:handle_gd_loss_visualizer, line 750 vs GD_Loss_Visualizer.py:visualizer_sptree, line 271]:**
The handler calls `visualizer_sptree(cli_args.gd_loss_result, sptree)` without passing `output_file`, so the output always uses the hardcoded default `"gd_loss_pie_visualizer.PDF"`. Neither docs nor README mention this output file name. The user has no way to control the output path through the CLI (the standalone `__main__` block has `-o` but the integrated CLI does not).

**[MEDIUM] [CONSISTENCY] [Phylo_Tracer.py:line 290 vs handle_gd_loss_visualizer, line 746-753]:**
`--output_dir` is defined in the parser and documented in README, but **never used** in the handler function. The output file is always written to the current working directory regardless of what `--output_dir` the user specifies. Fix: construct the output path using `cli_args.output_dir` and pass it to `visualizer_sptree()`.

**[MEDIUM] [CONSISTENCY] [docs/analyze.rst:line 285-298]:**
The analyze.rst section does not document the `--output_dir` optional parameter, while both the code (line 290) and README do list it.

**[LOW] [CONSISTENCY] [docs/analyze.rst:line 293 and README]:**
Both docs say `--input_sps_tree` expects "A numbered species tree file" but the code passes it through `Tree(cli_args.input_sps_tree, format=1)` which reads Newick with internal node names. The "numbered" qualifier is unclear — any species tree with labeled internal nodes works; the requirement is internal node names matching those in the loss summary, not specifically a "numbered" tree.

---

## C. Code Logic and Correctness

**[CRITICAL] [CODE] [Phylo_Tracer.py:handle_gd_loss_visualizer, line 750]:**
`--output_dir` is parsed but ignored. The `DEFAULT_OUTPUT_DIRS` map (line 937) sets the output directory to `"gd_loss_visualizer"`, and the main dispatch likely calls `os.chdir()` or similar — but `visualizer_sptree()` writes to a relative path `"gd_loss_pie_visualizer.PDF"`. If `output_dir` handling happens upstream (via `chdir`), the output goes into the managed directory. However, the `output_file` parameter is never wired through, so the filename is not configurable. Verify that the upstream `output_dir` plumbing actually applies here.

**[MEDIUM] [CODE] [GD_Loss_Visualizer.py:get_stats_deduplicated, line 131]:**
No existence check on `filepath` before `open()`. A missing file raises a bare `FileNotFoundError` with no user-friendly context. Suggested fix: add `if not os.path.isfile(filepath): raise FileNotFoundError(f"Loss summary file not found: {filepath}")`.

**[MEDIUM] [CODE] [GD_Loss_Visualizer.py:get_stats_deduplicated, line 125 vs line 124]:**
`event_loss_tracker` tracks types `{"2-0", "2-1", "2-2"}` while `loss_tracker` tracks `{"2-0", "2-1", "1-0"}`. The function returns three dicts with different type keys. Consumers must know which dict uses which key set. The `calculate_node_loss_score` function only reads `"2-0"` and `"2-1"` from event stats, silently ignoring `"2-2"` (retained events) and `"1-0"` if present. This is functionally correct but fragile — a future caller passing the wrong dict would get silent zero-counts.

**[MEDIUM] [CODE] [GD_Loss_Visualizer.py:visualizer_sptree, lines 311-314]:**
`NodeStyle()` is called inside the loop. If `ete3` import failed (lines 14-21), `NodeStyle` is `None`, and `NodeStyle()` will raise `TypeError: 'NoneType' object is not callable`. No guard exists at function entry to check whether ete3 was successfully imported. Fix: add a check at the top of `visualizer_sptree()`:
```python
if TreeStyle is None:
    raise ImportError("ete3 is required for visualization but was not found.")
```

**[LOW] [CODE] [GD_Loss_Visualizer.py:visualizer_sptree, lines 415-417]:**
`convert_to_ultrametric()` failure is silently swallowed with a bare `except Exception: pass`. If this call fails (e.g., negative branch lengths), the resulting tree rendering may be misleading. At minimum, log a warning.

**[LOW] [CODE] [GD_Loss_Visualizer.py:lines 371-377]:**
Commented-out PieChartFace code is dead code. `PieChartFace` is still imported (line 16). Consider removing both the import and the dead code block.

**[LOW] [CODE] [GD_Loss_Visualizer.py:visualizer_sptree, line 363]:**
When `norm_loss_score is None`, `float('nan')` is passed to `%-8.2f` format. This prints `nan` which is valid but could confuse downstream log parsers. Consider printing `"N/A"` or skipping the log line.

---

## D. Output Descriptions

**[MEDIUM] [DOCS] [docs/analyze.rst:lines 285-298 and README "GD_Loss_Visualizer"]:**
Neither docs nor README document the **output file** produced by this module. The code generates `gd_loss_pie_visualizer.PDF` containing an annotated species tree with:
- Blue numbers (branch-top): GD birth counts per node
- Red text (branch-bottom): `E(2-2)/E(2-1)/E(2-0)` event-level loss counts
- Purple `L=` score: normalized loss severity
- Gray `Lraw=`: raw loss score (when GD_birth=0)
- Green `P=`: path-level `2-0/2-1/1-0` counts
- Legend block explaining all annotations

This should be documented in both files.

**[LOW] [DOCS] [README "GD_Loss_Visualizer" section]:**
The README does not explain what metrics are shown on the tree or what the L-score formula is, making the output difficult to interpret without reading source code.

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 2     |
| MEDIUM     | 7     |
| LOW        | 5     |
| SUGGESTION | 0     |
| **Total**  | **14**|

### Priority fixes:
1. **CRITICAL**: Set `--input_sps_tree` to `required=True` in the argparse definition.
2. **CRITICAL**: Wire `--output_dir` through to `visualizer_sptree()` so the parameter actually takes effect.
3. **MEDIUM**: Add ete3 import guard at the top of `visualizer_sptree()`.
4. **MEDIUM**: Document the output file name and its annotation meanings in both docs/analyze.rst and README.

==================================================================
## Module: GD_Visualizer
==================================================================
Here is the complete review of the **GD_Visualizer** module.

---

## A. Scientific Correctness

**No issues found.** The GD type classification (AABB, AXBB, AABX, Complex), the legacy label remapping, the pie-chart proportional visualization, and the per-node aggregation are all standard and scientifically sound for gene duplication event summarization.

---

## B. Code vs Docs Consistency

**[MEDIUM] [B] [docs/analyze.rst:GD_Visualizer §, line 246–260]: `--output_dir` optional parameter missing from analyze.rst**

The CLI parser (Phylo_Tracer.py:273) defines `--output_dir`, and the README documents it, but analyze.rst omits it entirely. The GD_Visualizer section in analyze.rst lists zero optional parameters.

*Fix:* Add an "Optional Parameters" subsection with `--output_dir` to the analyze.rst GD_Visualizer entry.

---

**[CRITICAL] [B] [Phylo_Tracer.py:handle_gd_visualizer, line 692–700]: `--output_dir` is parsed but never used**

The CLI defines `--output_dir` (line 273), but `handle_gd_visualizer` never passes `cli_args.output_dir` to `gd_visualizer_main`. Furthermore, `gd_visualizer_main` does not accept an `output_dir` parameter — the output PDF is always written to CWD. Users who set `--output_dir` will have their argument silently ignored.

*Fix:* Thread `output_dir` through `gd_visualizer_main` and prepend it to `output_pdf`:
```python
def gd_visualizer_main(sptree, gd_result, taxa, output_dir=None):
    gds = process_gd_result(gd_result)
    count_dic = get_count_dic(gds)
    gd_base = os.path.splitext(os.path.basename(gd_result))[0]
    output_pdf = f"{gd_base}.pdf"
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        output_pdf = os.path.join(output_dir, output_pdf)
    mark_sptree(sptree, count_dic, taxa, output_pdf=output_pdf)
```
And in the handler:
```python
gd_visualizer_main(sptree, cli_args.gd_result, taxa, output_dir=cli_args.output_dir)
```

---

**[MEDIUM] [B] [GD_Visualizer.py:__main__ (line 297) vs Phylo_Tracer.py:handle_gd_visualizer (line 695)]: Inconsistent tree-reading function**

The `__main__` block calls `read_phylo_tree()` (returns `ete3.PhyloTree`), while the CLI handler in `Phylo_Tracer.py` calls `read_tree()` (returns `ete3.Tree`). `PhyloTree` is a subclass of `Tree` with reconciliation methods; behaviour could differ for `.convert_to_ultrametric()` or `.sort_descendants()`. The two entry points produce subtly different tree objects.

*Fix:* Use the same reader in both places. Since GD_Visualizer only uses basic `Tree` methods, `read_tree` is sufficient — update the `__main__` block to use `read_tree`.

---

## C. Code Logic and Correctness

**[MEDIUM] [C] [GD_Visualizer.py:mark_sptree, line 138]: No guard for ete3 import failure**

Lines 12–19 set `TreeStyle`, `TextFace`, etc. to `None` when ete3 is missing. But `mark_sptree` unconditionally calls `TreeStyle()`, `TextFace()`, etc. (line 138+). If ete3 is absent, this raises `TypeError: 'NoneType' object is not callable` — a cryptic error message.

*Fix:* Add an early check at the top of `mark_sptree`:
```python
if TreeStyle is None:
    raise ImportError("ete3 is required for GD_Visualizer. Install it with: pip install ete3")
```

---

**[LOW] [C] [GD_Visualizer.py:mark_sptree, line 180]: Inconsistent TextFace styling**

Line 180: `TextFace(" AXBB", fsize=6)` is missing `ftype="Arial"`, while every other `TextFace` call in the same title block (lines 176, 178, 182, 184) explicitly sets `ftype="Arial"`. This causes the AXBB label to render in the default font.

*Fix:* `TextFace(" AXBB", fsize=6, ftype="Arial")`

---

**[MEDIUM] [C] [GD_Visualizer.py:process_gd_result + Phylo_Tracer.py:handle_gd_visualizer]: `--input_imap` semantics mismatch for species-tree leaf renaming**

The `--input_imap` parameter is documented as `gene_id<TAB>species_name` (a gene-to-species mapping). In `handle_gd_visualizer`, this dict is passed as `taxa` to `mark_sptree`, which uses it at lines 189–191 to rename **species tree** leaves:
```python
for leaf in sptree:
    if leaf.name in taxa:
        leaf.name = taxa[leaf.name]
```
If the imap keys are gene IDs (e.g., `ATCG00500.1`), they will never match species-tree leaf names (e.g., `Arabidopsis_thaliana`), so the renaming silently does nothing. Either:
- The parameter is misnamed/misdocumented and should actually be a species-to-display-label mapping, or
- The renaming logic is dead code in practice.

*Fix:* Clarify the intended input. If renaming is desired, document a separate `--input_taxa` parameter (species_name→display_label). If not, remove the renaming loop or document that `--input_imap` is used only for internal consistency and has no visible effect here.

---

**[SUGGESTION] [C] [GD_Visualizer.py:process_gd_result, line 50]: No encoding specified for `open()`**

`open(gd_file)` uses the platform default encoding. If the GD result file contains non-ASCII species names (common in taxonomy), this can fail on systems with non-UTF-8 locale.

*Fix:* `open(gd_file, encoding="utf-8")`

---

## D. Output Descriptions

**[MEDIUM] [D] [docs/analyze.rst:GD_Visualizer § + README GD_Visualizer §]: No output documentation**

Neither analyze.rst nor the README describes the output files produced by GD_Visualizer. The actual output is a PDF file named `{gd_result_basename}.pdf` containing an ultrametric species tree annotated with:
- Pie charts showing GD-type proportions (AABB/AXBB/AABX/Complex) per node
- Red numbers showing total GD event counts per node
- Blue numbers showing node identifiers for internal nodes
- Leaf names in italic

Other modules (Phylo_Rooter, MulRF_Distance, GD_Loss_Tracker) have explicit "Output:" sections.

*Fix:* Add output documentation to both files, e.g.:
```
Output:
    - ``{gd_result_basename}.pdf``  Annotated species tree PDF with per-node GD counts
                                    (red) and GD-type pie charts (AABB/AXBB/AABX/Complex).
```

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 1     |
| MEDIUM     | 5     |
| LOW        | 1     |
| SUGGESTION | 1     |
| **Total**  | **8** |

The single **CRITICAL** issue (`--output_dir` silently ignored) is straightforward to fix by threading the parameter through `gd_visualizer_main`. The most impactful **MEDIUM** issues are the `--input_imap` semantic mismatch (likely dead renaming logic) and the missing output documentation.

==================================================================
## Module: HaploFinder
==================================================================
Now I have all the information needed. Here is the complete review.

---

# HaploFinder Module — Code Review Report

## A. Scientific Correctness

**[CRITICAL] [SCIENTIFIC] [HaploFinder.py:find_gene_conversion, line 681]**
The gene-conversion detection logic uses a hardcoded chromosome-number relationship `chr1 * 2 == chr2` to identify homeologous chromosome pairs. This assumes a very specific genome organization (e.g., one subgenome has chromosomes 1–10 and the other 2,4,6…20) and will silently produce empty results for any organism that doesn't follow this exact numbering scheme. This should be parameterized or replaced with a user-provided homeolog mapping file.

**[CRITICAL] [SCIENTIFIC] [HaploFinder.py:get_chromosome_subgenome, lines 1180–1205]**
Subgenome assignment is hardcoded: chromosomes 1–10 → "A", 11–20 → "B". This is organism-specific (likely modeled on a specific allopolyploid) and will silently return `None` for any other chromosome naming convention. The `split_sequences` function (line 1101) calls this for validation, meaning the `assignment_matches_expected` column in `split_assignment.tsv` will be meaningless for most organisms. **Suggested fix:** accept a user-provided chromosome-to-subgenome mapping file, or at minimum document this as an organism-specific validation heuristic.

**[MEDIUM] [SCIENTIFIC] [HaploFinder.py:process_gd_result, lines 861/885]**
Only duplication nodes with exactly 3 or 4 leaves are processed for blue/red labeling. Larger duplication nodes are silently skipped. In real gene families, GD nodes frequently subtend many more leaves. This will cause systematic under-detection of GD-derived pairs, biasing conversion zone analysis.

**[MEDIUM] [SCIENTIFIC] [HaploFinder.py:find_conversion_zones_with_ids_to_file, lines 746–808]**
Variable naming is inverted: the code detects a **blue→red→blue** pattern (lines 753, 767, 774), but the variables are named `start_red1`, `start_blue`, `start_red2`. In the HaploFinder color scheme, "red" = orthologous (same-subgenome) and "blue" = homeologous (cross-subgenome). The algorithm correctly looks for a region where the dominant relationship switches, but the misleading names could cause maintenance errors.

**[LOW] [SCIENTIFIC] [HaploFinder.py:assign_hybrid_subgenome, line 983]**
`distance_threshold = len(diploid_tags)` — the topological distance threshold is set equal to the number of diploid progenitor species. There is no biological justification for this coupling. For 2 progenitors the threshold is 2 nodes, for 3 it's 3, etc. This should be independently configurable or documented with a rationale.

---

## B. Code vs Docs Consistency

**[CRITICAL] [DOCS-MISMATCH] [docs/analyze.rst, Section 16 — missing `--input_sps_tree`]**
The `--input_sps_tree` parameter is required in haplofinder mode (checked at `Phylo_Tracer.py:854`) but is **not listed** in analyze.rst's Required or Optional parameters. The README correctly lists it. **Fix:** add `--input_sps_tree` to the Required Parameters list in analyze.rst.

**[MEDIUM] [DOCS-MISMATCH] [docs/analyze.rst, Section 16 — missing split-mode parameters]**
analyze.rst does not document any split-mode-specific parameters (`--hyb_sps`, `--parental_sps`, `--input_fasta`, `--cluster_file`). The README documents them. **Fix:** add a "Mode = split required" subsection to analyze.rst mirroring the README.

**[MEDIUM] [DOCS-MISMATCH] [docs/analyze.rst, Section 16 — missing `--output_dir`]**
The `--output_dir` parameter exists in code (line 338) and README but is absent from analyze.rst.

**[LOW] [DOCS-MISMATCH] [README.md, HaploFinder section — `--mode` listed as required]**
The README lists `--mode` under "Required parameter" but it has `default='haplofinder'` in code (line 332), making it optional. analyze.rst correctly lists it as optional. **Fix:** move `--mode` to optional in README.

**[LOW] [DOCS-MISMATCH] [README.md, HaploFinder section — split mode lists `--species_b_gff` inconsistently]**
README lists `--species_b_gff` under "Mode = split required" (line 739). While the code does check it at line 834, the parameter semantically belongs to haplofinder mode (GFF of species B). In split mode it's used as a general GFF for the hybrid species. The parameter name is misleading in the split context.

---

## C. Code Logic and Correctness

**[CRITICAL] [BUG] [HaploFinder.py:plot_chr1/plot_chr2 + generate_dotplot, lines 147–148/191–192 + 435–442]**
When `total_lens == 0`, `plot_chr1`/`plot_chr2` return `None` (bare `return` on line 148/192). The return value is assigned to `step_1`/`step_2` (line 435–436) and passed to `gene_location` (line 441–442), where it's used as a float multiplier (`* step`). This will crash with `TypeError: unsupported operand type(s) for *: 'float' and 'NoneType'`. **Fix:** return 0 instead of bare `return`, or guard the downstream calls.

**[CRITICAL] [BUG] [HaploFinder.py:generate_dotplot + process_gd_result + find_gene_pair_info — output_dir ignored]**
Multiple functions write output files to the current working directory, ignoring the `--output_dir` parameter:
- `process_gd_result` → `color_label.txt` (line 912)
- `generate_dotplot` → `*_dotplot.pdf/png` (lines 467–468)
- `find_gene_pair_info` → `gene_conversion_*.txt` (line 739)
- `find_conversion_zones_with_ids_to_file` → `gene_conversion.txt` (line 792)

The `--output_dir` argument is accepted by the CLI but never passed through to these functions. **Fix:** thread `output_dir` through the call chain and prefix all output paths.

**[MEDIUM] [BUG] [HaploFinder.py:judge_support, lines 546–568]**
The function returns `None` implicitly when `support` and `support_value` fall outside all four branch conditions (e.g., `support < 0` or `support_value < 0.5`). This is used in a boolean context (`if judge_support(...)`) so `None` acts as `False`, but it masks invalid inputs. More importantly, support values < 50 (i.e., `support_value` in range [0, 0.5)) are not handled at all — they fall through all branches. **Fix:** add an else clause that handles the full range, or raise an error for invalid inputs.

**[MEDIUM] [BUG] [HaploFinder.py:find_gene_pair_info, line 720]**
The output file `gene_conversion_{gd_pairs}.txt` is opened (and truncated) twice — once at line 720 inside a loop that only collects data into `sort_lst`, and again at line 739 for actual writing. The first `open` truncates the file for no reason. **Fix:** remove the first `with open` block and keep only the second.

**[MEDIUM] [EDGE-CASE] [HaploFinder.py:find_conversion_zones_with_ids_to_file, line 752]**
The `while i < n - 1` loop condition means a single-element `data` list (`n=1`) is never entered. If `data` is empty (`n=0`), `n - 1 = -1` which is fine. But `data` with exactly one blue entry is silently skipped. This is arguably correct (no zone can be formed from 1 entry) but should be documented.

**[MEDIUM] [EDGE-CASE] [HaploFinder.py:get_ortholog_pairs_by_species, lines 596–599]**
If `sp1` or `sp2` is not present in one of the child clades, `dic1[sp1]` or `dic2[sp2]` will raise a `KeyError`. The caller filters for `{sp1, sp2}.issubset(sps_tol)` at the tree level, but after the tree is split at a duplication node, individual child clades may lack one species. **Fix:** use `.get(sp, [])` instead of direct dictionary access.

**[LOW] [CODE-QUALITY] [Phylo_Tracer.py:handle_haplofinder, line 875]**
`size = cli_args.size if cli_args.size is not None else 0.001` — the fallback `0.001` is dead code because argparse sets `default=0.0005`, so `cli_args.size` is never `None`. The fallback value (0.001) also differs from the documented default (0.0005), which is confusing. **Fix:** remove the conditional and use `cli_args.size` directly.

**[LOW] [CODE-QUALITY] [HaploFinder.py:split_sequences, line 1069]**
`_ = gff_1` — the GFF list is read but immediately discarded. This wastes memory on large GFF files. Consider only reading the dict if the list is not needed.

---

## D. Output Descriptions

**[MEDIUM] [DOCS-INCOMPLETE] [README.md + docs/analyze.rst — undocumented output files]**
Haplofinder mode generates two additional output files not documented anywhere:
- `gene_conversion_gd_pairs.txt` — per-gene-pair color labels sorted by chromosomal position (from `find_gene_pair_info`)
- `gene_conversion.txt` — detected conversion zones with genomic coordinates (from `find_conversion_zones_with_ids_to_file`)

**Fix:** add these to the "PhyloTracer Results Files" section of README and to analyze.rst.

**[MEDIUM] [DOCS-INCOMPLETE] [docs/analyze.rst, Section 16 — no output description at all]**
analyze.rst's HaploFinder section has no "Output" subsection describing any of the generated files (dotplot PDFs/PNGs, color_label.txt, gene_conversion files, or split-mode outputs). Compare with the Phylo_Rooter section which has a proper Output subsection.

**[LOW] [DOCS-INCOMPLETE] [README.md — split_assignment.tsv column descriptions missing]**
The README lists `split_assignment.tsv` as an output but doesn't describe its columns. The file contains: `gene_id`, `subgenome`, `status`, `chromosome`, `expected_subgenome`, `assignment_matches_expected` (see code line 1118). Documenting these would help users interpret results.

---

## Summary Table

| Severity | Count |
|----------|-------|
| CRITICAL | 5 |
| MEDIUM | 10 |
| LOW | 5 |
| SUGGESTION | 0 |
| **Total** | **20** |

| Category | Count |
|----------|-------|
| Scientific correctness | 5 |
| Code vs docs consistency | 5 |
| Code logic / bugs | 7 |
| Output descriptions | 3 |

### Top 3 priorities for fixing:
1. **Thread `--output_dir` through all file-writing functions** — currently the parameter is accepted but ignored, so all outputs land in CWD regardless.
2. **Guard against `None` step from empty lens files** — `plot_chr1`/`plot_chr2` returning `None` will crash downstream.
3. **Parameterize the hardcoded chromosome→subgenome mapping** — the `chr*2` and `1-10/11-20` assumptions make the module unusable for most organisms without code changes.

==================================================================
## Module: Hybrid_Tracer
==================================================================
Now I have all the information needed. Here is the complete review:

---

# Hybrid_Tracer Code Review

## A. Scientific Correctness

**[MEDIUM] [SCIENTIFIC] [docs/analyze.rst:323]**: The description says "ABAB-BABA test" — this is a misnomer. The standard name is the **ABBA-BABA test** (Patterson's D-statistic). Moreover, HyDe does not implement the D-statistic; it uses **coalescent-based phylogenetic invariants** (Kubatko & Chifman 2019). The README module list (line 78) correctly says "coalescent-based phylogenetic invariants" but the analyze.rst description contradicts this.
**Fix**: Replace "ABAB-BABA test" with "coalescent-based phylogenetic invariants (HyDe)" in analyze.rst line 323.

**[MEDIUM] [SCIENTIFIC] [README.md:195]**: Same misnomer: "detects species hybridization signals using ABAB-BABA test by HyDe". Should be "ABBA-BABA" at minimum, but ideally "phylogenetic invariants".
**Fix**: Change to "detects species hybridization signals using phylogenetic invariants by HyDe".

**[MEDIUM] [SCIENTIFIC] [Hybrid_Tracer.py:get_process_gd_clade:471]**: Hardcoded filtering thresholds `total_num >= 10` and `type_ratio >= 0.1` are not documented or scientifically justified. These thresholds determine which GD nodes enter HyDe analysis, making them critical for reproducibility.
**Fix**: Either expose these as CLI parameters or document their justification in the docstring.

**[MEDIUM] [SCIENTIFIC] [Hybrid_Tracer.py:get_gd_count_dic_and_gd_type_dic:423]**: `find_dup_node(t1, rename_sptree, 50, 0, 2, 0, 1)` uses hardcoded GD detection parameters (support=50, subclade_support=0, dup_species_num=2, dup_species_proportion=0, deepvar=1). These differ from some GD_Detector defaults (e.g., `dup_species_proportion` default is 0.2 in GD_Detector) and are not user-configurable in Hybrid_Tracer. Users cannot tune GD sensitivity for hybridization analysis.
**Fix**: Either expose these as optional CLI parameters or document the fixed values and rationale.

**[LOW] [SCIENTIFIC] [Hybrid_Tracer.py:matrix_to_phy:943]**: PHYLIP strict format truncates species names to 10 characters (`name = species[:10].ljust(10)`). For taxonomic names >10 characters (extremely common, e.g., "Arabidopsis_thaliana" → "Arabidopsi"), this causes silent truncation and potential name collisions between species.
**Fix**: Use relaxed PHYLIP format or a longer name field, and add collision detection.

## B. Code vs Docs Consistency

**[MEDIUM] [CONSISTENCY] [docs/analyze.rst:318-337]**: The `--output_dir` optional parameter is missing from the analyze.rst documentation for Hybrid_Tracer. It exists in both the CLI parser (Phylo_Tracer.py:307) and README (line 692).
**Fix**: Add `--output_dir` to the Optional Parameters section of Hybrid_Tracer in analyze.rst.

**[LOW] [CONSISTENCY] [docs/analyze.rst:323 vs README.md:682]**: The description text differs between the two docs. analyze.rst says "To use the ABAB-BABA test to detect hybridization signals for each potential GD burst event across species tree detect species hybridization events for" (also has a grammatical error — trailing "for"). README says "To detect hybridization signals from GD-derived gene sets and run HyDe testing on species-tree mapped GD nodes" (clearer and more accurate).
**Fix**: Align analyze.rst description with the README version and fix the trailing grammar error.

**[LOW] [CONSISTENCY] [docs/analyze.rst:331]**: analyze.rst says `--mrca_node` format is "SpeciesA,SpeciesB (comma-separated, no space)" but does not mention the behavior when multiple `--mrca_node` flags are given. The README (line 689) correctly states "If multiple are provided, only the first valid pair is used."
**Fix**: Add the multiple-pair behavior note to analyze.rst.

No issues found with parameter types, defaults, or required/optional labeling — these are consistent across code, analyze.rst, and README.

## C. Code Logic and Correctness

**[CRITICAL] [LOGIC] [Hybrid_Tracer.py:run_hyde_from_matrix_integrated:810]**: No guard against `hd is None`. If `phyde` is not installed (line 20-22 sets `hd = None`), calling `hd.HydeData(...)` raises `AttributeError` with no informative message. The user gets a cryptic crash.
**Fix**: Add an early check at the start of `hyde_main` or `run_hyde_from_matrix_integrated`:
```python
if hd is None:
    raise ImportError("phyde is required for Hybrid_Tracer. Install with: pip install phyde")
```

**[MEDIUM] [LOGIC] [Hybrid_Tracer.py:matrix_to_phy:893]**: `seq = "-" * max(len(s) for s in seq_dic.values())` is called inside a double loop (for each column × each species). If `seq_dic` is large, this recomputes the max length on every missing-data cell. More critically, if `seq_dic` becomes empty at this point (e.g., all KeyErrors), `max()` raises `ValueError: max() arg is an empty sequence`.
**Fix**: Precompute max length before the loop and add an empty-check guard.

**[MEDIUM] [LOGIC] [Hybrid_Tracer.py:create_fasta_dict:112]**: `fasta_dict[gene2new_named_gene_dic[record.id]] = str(record.seq)` raises `KeyError` if a FASTA record ID is not in the mapping dictionary. No error handling or logging.
**Fix**: Add a try/except or `.get()` with a warning for unmapped records.

**[MEDIUM] [LOGIC] [Hybrid_Tracer.py:clean_matrix_by_dash_count:845]**: `matrix.T.drop_duplicates().T` removes columns with identical value patterns. Two biologically distinct gene assignments that happen to produce the same species-to-gene mapping pattern will be silently dropped, losing phylogenetic signal.
**Fix**: Consider whether this deduplication is biologically justified; if so, document the rationale. If not, remove it.

**[LOW] [LOGIC] [Hybrid_Tracer.py:run_hyde_from_matrix_integrated:819]**: In the `finally` block, if `tmp_phy_path` or `tmp_imap_path` was never actually written to (e.g., `matrix_to_phy` raised), `os.remove()` could fail, masking the original exception.
**Fix**: Use `try: os.remove(...) except OSError: pass` in the finally block, or use a context manager.

**[LOW] [LOGIC] [Hybrid_Tracer.py:run_hyde_from_matrix_integrated:793]**: `voucher2taxa_dic[species]` can raise `KeyError` if a species in the matrix index is not in the voucher mapping. No guard.
**Fix**: Use `.get()` with fallback or add error handling.

## D. Output Descriptions

**[MEDIUM] [DOCS] [docs/analyze.rst:318-337]**: No output files or output columns are documented for Hybrid_Tracer. The section lacks any "Output:" subsection describing `hyde_out.txt`, `hyde_filtered_out.txt`, or their column structure.
**Fix**: Add an Output section listing:
- `hyde_out.txt` — all HyDe test results
- `hyde_filtered_out.txt` — Bonferroni-filtered results (p < 0.05/n_triples, |Zscore| ≠ 99999.9, 0 < Gamma < 1)
- Column names: P1, Hybrid, P2, Zscore, Pvalue, Gamma, AAAA, AAAB, AABA, AABB, AABC, ABAA, ABAB, ABAC, ABBA, BAAA, ABBC, CABC, BACA, BCAA, ABCD

**[MEDIUM] [DOCS] [README.md:679-694]**: The README CLI reference for Hybrid_Tracer also lacks any output description. The PhyloTracer Results Files section (line 263) mentions the file names but not column definitions or filtering criteria.
**Fix**: Add output column descriptions to the README CLI block, matching the header in code (line 615-616).

**[LOW] [DOCS] [README.md:263]**: The results section lists `hyde_out.txt` and `hyde_filtered_out.txt` but does not explain the difference (all results vs filtered results) or the filtering criteria used.
**Fix**: Add a brief note: "hyde_filtered_out.txt contains results passing Bonferroni-corrected significance (p < 0.05/n), valid Z-scores, and 0 < γ < 1."

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 1     |
| MEDIUM     | 9     |
| LOW        | 6     |
| SUGGESTION | 0     |
| **Total**  | **16** |

The single **CRITICAL** issue is the missing `phyde` import guard — if a user runs Hybrid_Tracer without phyde installed, they get an uninformative crash. The **MEDIUM** issues cluster around (1) scientifically inaccurate method naming in docs, (2) hardcoded, undocumented biological thresholds, (3) missing output documentation, and (4) unguarded KeyError paths in data processing functions.

==================================================================
## Module: Hybrid_Visualizer
==================================================================
Now I have everything I need. Here is the comprehensive review:

---

# Hybrid_Visualizer Module Review

## A. Scientific Correctness

**[LOW] [SCIENTIFIC] [Hybrid_Visualizer.py:calculate_gamma, line 153]:** When all entries are `nan` or `-inf`, `calculate_gamma` returns `0.0`. A gamma of 0.0 has biological meaning (no hybridization signal from one parent), but here it's used as a sentinel for "no data." This conflates absence of data with a meaningful biological value. **Fix:** Return `float('nan')` when `valid_count == 0`, and handle NaN downstream.

**[LOW] [SCIENTIFIC] [Hybrid_Visualizer.py:calculate_pvalue, line 185]:** Same issue — averaging p-values is statistically problematic (Fisher's method or Stouffer's method would be more appropriate for combining p-values). Arithmetic mean of p-values is not a standard statistical combination. **Fix:** Consider Fisher's combined probability test or at minimum document this as an approximation.

**[MEDIUM] [SCIENTIFIC] [Hybrid_Visualizer.py:create_hot_map_leaf, lines 660-666]:** The p-value comparison logic keeps only the direction (P1→P2 vs P2→P1) with the *lower* p-value (more significant), zeroing out the other gamma. However, lower p-value means stronger rejection of H₀ (no hybridization), so this selects the more significant signal — the logic is correct but the code comment is absent, making it easy to misread. **Suggestion:** Add a brief comment explaining the biological rationale.

**[LOW] [SCIENTIFIC] [docs/analyze.rst, line 344-345]:** The description says "highlighting support from gene tree topologies and D-statistic signals" but the code processes HyDe output (which uses phylogenetic invariants / coalescent-based γ estimation), not the classic D-statistic (ABBA-BABA). HyDe and D-statistics are related but distinct methods. **Fix:** Change to "highlighting hybridization proportions (γ) from HyDe analysis."

## B. Code vs Docs Consistency

**[CRITICAL] [DOCS-CODE] [Phylo_Tracer.py:314 vs Hybrid_Visualizer.py]:** The CLI defines `--output_dir` (line 314), but `handle_hybrid_visualizer` (line 816-826) never passes `output_dir` to any function, and `Hybrid_Visualizer.py` has zero references to `output_dir`. The parameter is accepted but silently ignored — all output goes to the current working directory regardless. **Fix:** Either implement `output_dir` support (os.chdir or path-prefix logic) or remove the argument.

**[MEDIUM] [DOCS-CODE] [docs/analyze.rst:347]:** Docs say `--hyde_out` description is "File containing the result of HyDe from Hybrid_Tracer." The CLI (Phylo_Tracer.py:311) says "HYDE output table file." Minor inconsistency in naming convention — docs reference the upstream module (Hybrid_Tracer), CLI does not. **Fix:** Harmonize descriptions.

**[LOW] [DOCS-CODE] [docs/analyze.rst:351]:** The `--node` description is unclear: "stack up all the heatmaps for each monophyletic clade respectively, only the squares in all heatmaps were light, the square after superimposition will be light." This is grammatically broken and hard to interpret. The README version (line 703) is clearer: "Use node-mode heatmaps (monophyletic clade stacking) instead of leaf-mode output." **Fix:** Align analyze.rst description with README wording.

**[LOW] [DOCS-CODE] [README.md:700 vs Phylo_Tracer.py:311]:** README says `--hyde_out` with description "HYDE output table file"; docs say "File containing the result of HyDe from Hybrid_Tracer." Case inconsistency: "HYDE" vs "HyDe" across docs/README/code. **Fix:** Standardize to "HyDe" (the official tool name).

## C. Code Logic and Correctness

**[CRITICAL] [CODE] [Hybrid_Visualizer.py:parse_hyde_out, line 63]:** No bounds check on `i1` length. If a HyDe output line has fewer than 6 tab-separated fields, `i1[5]` will raise an `IndexError`. **Fix:** Add `if len(i1) < 6: continue` after splitting.

**[CRITICAL] [CODE] [Hybrid_Visualizer.py:create_hot_map_node, line 508]:** `hyb_dic[i]` will raise `KeyError` if a leaf from `node.get_leaf_names()` is not present as a hybrid in `hyb_dic`. Not all leaves in a clade will necessarily appear as hybrids in the HyDe output. **Fix:** Add `if i not in hyb_dic: continue`.

**[MEDIUM] [CODE] [Hybrid_Visualizer.py:generate_tree_leaf, lines 214-248]:** If `NodeStyle`, `TextFace`, or `TreeStyle` are `None` (import failed, lines 21-25), this function will crash with `TypeError: 'NoneType' object is not callable`. No guard exists. **Fix:** Add an early check: `if NodeStyle is None: raise ImportError("ete3 GUI components required")`.

**[MEDIUM] [CODE] [Hybrid_Visualizer.py:combine_fig, line 746]:** `Image.open()` has no error handling. If the intermediate files don't exist (e.g., rendering failed silently), this crashes with `FileNotFoundError`. **Fix:** Wrap in try/except or check `os.path.exists()` first.

**[MEDIUM] [CODE] [Hybrid_Visualizer.py:hyde_visual_node_main, line 920]:** `generate_tree_node(t1, node, node.name)` passes the original `node` object but `t1` is a *copy* of `sptree`. The node reference belongs to the original tree, not the copy. When `generate_tree_node` calls `node.traverse()` at line 272 and then tries to match names against nodes in `t1`, it works only because matching is by `.name` string. However, `node.name` at line 288 (`if i.name == node.name`) could match incorrectly if multiple internal nodes share the same name (common before `num_tre_node` assigns unique IDs). This is fragile. **Fix:** Find the corresponding node in `t1` by leaf set or use `t1&node.name`.

**[MEDIUM] [CODE] [Hybrid_Visualizer.py:create_hot_map_node, lines 524-528]:** Zeroing out rows/columns for hybrid clade leaves removes their data from the heatmap. But the mask at line 530 (`np.triu`) only masks the upper triangle. This means the lower-triangle cells for hybrid leaves are also zeroed, which could hide valid data for asymmetric P1-P2 relationships. This may be intentional (clade self-comparisons are meaningless) but is undocumented.

**[LOW] [CODE] [Hybrid_Visualizer.py:create_hot_map_node, line 593]:** `ax.imshow()` is called *after* three `sns.heatmap()` calls on the same axes. `imshow` renders a raster on top of the seaborn heatmaps, potentially overwriting them entirely. This appears to be used only to create the colorbar reference (`m`), but it has the side effect of overlaying a new image. **Fix:** Use `plt.cm.ScalarMappable` to create the colorbar without rendering an image:
```python
m = plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=1), cmap=newcmp)
```

**[LOW] [CODE] [Hybrid_Visualizer.py:create_hot_map_leaf, line 718]:** Same `ax.imshow()` issue as above.

**[LOW] [CODE] [Hybrid_Visualizer.py:hyde_visual_cmap, lines 396-410]:** The colormap builds 10,000 RGBA entries where alpha transitions from 0→1→0. At the extremes (gamma near 0 or 1), cells will be nearly transparent, making them invisible against a white background. This means the strongest hybridization signals (gamma ≈ 0 or ≈ 1) become invisible. **Fix:** If transparency is intentional for the midpoint, the alpha gradient direction may be reversed from what's expected biologically.

**[LOW] [CODE] [Hybrid_Visualizer.py:combine_fig, lines 782-804]:** Font path list is hardcoded with macOS, Linux, and Windows paths. On systems where none match, `ImageFont.load_default()` is used, ignoring `base_font_size`/`title_font_size`/`legend_font_size` parameters (default font doesn't accept size arguments in older Pillow versions).

## D. Output Descriptions

**[MEDIUM] [DOCS] [docs/analyze.rst:340-357 and README.md:695-706]:** Neither docs nor README describe the output files generated by Hybrid_Visualizer. The code generates:
- `{species/node_name}.png` — combined tree+heatmap figure (per hybrid species in leaf mode, per internal node in node mode)
- Intermediate files `_img_faces.png` and `_hotmap.png` are created and deleted.

**Fix:** Add an "Output" section describing the generated PNG files, their naming convention, and content (tree with highlighted hybrid in red, heatmap showing gamma values and combination counts).

**[MEDIUM] [DOCS] [README.md and docs/analyze.rst]:** The heatmap legend labels (visible in code at lines 832-844) use "Red: hybrid", "Yellow: y values", "White: number of hybridization combinations". The term "y values" is used inconsistently — the code uses both "gamma" (γ) and "y" to refer to the hybridization proportion. Docs should clarify that y = γ (gamma) from HyDe. **Fix:** Standardize terminology.

**[LOW] [DOCS] [README.md:79]:** Module feature description says "Highlights hybridization proportions (γ) and support values across the species tree." The code does not display "support values" — it shows combination counts and gamma values. **Fix:** Change "support values" to "combination counts."

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 3     |
| MEDIUM     | 7     |
| LOW        | 8     |
| SUGGESTION | 0     |
| **Total**  | **18**|

### Priority fixes:
1. **CRITICAL:** `parse_hyde_out` — add field-count guard to prevent `IndexError` on malformed lines
2. **CRITICAL:** `create_hot_map_node` — guard against `KeyError` when leaf is not in `hyb_dic`
3. **CRITICAL:** `--output_dir` is defined in CLI but completely ignored — implement or remove

==================================================================
## Module: MulRF_Distance
==================================================================
Now I have all the information needed. Here is the complete review.

---

# MulRF_Distance Module — Code Review Report

## A. Scientific Correctness

**[MEDIUM] [A] [MulRF_Distance.py:module + README:MulRF_Distance section]**: The implementation projects gene-tree bipartitions to species level (collapsing multi-copy duplicates within clades) and computes standard symmetric-difference RF distance. This is **not** the formal MulRF distance of Chaudhary et al. (2013, 2015), which uses a minimum-cost bipartition matching that accounts for label multiplicity. The module docstring (line 4) correctly hedges with "MulRF-like distances," but the README (line 347) and CLI help text say simply "MulRF topological conflict distances" without qualification. **Suggested fix:** Add a brief note in the README description clarifying this is a species-level projected RF distance (a simplification of the formal MulRF), or cite the specific variant used.

**[LOW] [A] [MulRF_Distance.py:compute_mulrf:142]**: When both trees yield zero bipartitions (e.g., two shared species producing star-like trees), `max_rf == 0` and `normalized` defaults to `0.0`. This implies "perfect agreement" for a degenerate case. Returning `None` or `NaN` would be more scientifically neutral. **Suggested fix:** Return `None` when `max_rf == 0` (mirroring the `n_shared < 2` error path), or document this edge-case convention.

No issues found with parameter defaults or metric direction definitions — the "lower raw = better" / "higher raw = better" annotations in the Phylo_Rooter docs are consistent with the code.

---

## B. Code vs Docs Consistency

**[CRITICAL] [B] [docs/analyze.rst]**: MulRF_Distance is **completely absent** from `analyze.rst`. All other 16 modules have dedicated sections. MulRF is only mentioned tangentially as a Phylo_Rooter sub-metric (lines 90, 119). **Suggested fix:** Add a standalone MulRF_Distance section (numbered appropriately) with Description, Required/Optional Parameters, Output, and Example, matching the README content.

**[MEDIUM] [B] [README.md:258]**: The "PhyloTracer Results Files" section lists the output path as `mulrf_distance/mulrf_distance.tsv` (implying a subdirectory). However, the code writes directly to whatever `--output` specifies (default: `mulrf_distance.tsv` in CWD) with **no subdirectory creation**. There is also no `--output_dir` argument for this module, unlike most other modules. **Suggested fix:** Either (a) add `--output_dir` support to MulRF_Distance (consistent with other modules), or (b) correct the README to show the actual flat-file output path.

**[LOW] [B] [README.md:367]**: The output column description says `tre_id_1 — Tree ID from GF1` and `tre_id_2 — Tree ID from GF2`. In mode 2, `tre_id_2` is always the hardcoded string `"species_tree"`, not a GF ID. **Suggested fix:** Add a note: "In mode 2, `tre_id_2` is always `species_tree`."

**[LOW] [B] [MulRF_Distance.py:mulrf_main + Phylo_Tracer.py:457-463]**: `mulrf_main()` exposes `separator` and `position` parameters (lines 220-221) but these are **never passed** from the CLI handler (which omits them, so defaults `"_"` / `"last"` always apply). Since `--input_imap` is required and provides the gene→species mapping, name-based inference is only a fallback for missing genes. This is not a bug, but the undocumented parameters could confuse contributors. **Suggested fix:** Either expose `--separator`/`--position` as CLI args, or mark them clearly as internal/deprecated in the function signature.

---

## C. Code Logic and Correctness

**[MEDIUM] [C] [MulRF_Distance.py:mulrf_main:249,290]**: Bare `except Exception` blocks silently swallow **all** errors (malformed Newick, I/O errors, ete3 crashes) and produce all-`None` rows with no diagnostic output. A broken tree file is indistinguishable from a legitimate low-overlap pair. **Suggested fix:** Log the exception with `logger.warning("Skipping %s: %s", tre_id, e)` inside each except block.

**[MEDIUM] [C] [MulRF_Distance.py:get_bipartitions_species_level:58-62]**: Multi-copy genes mapping to the same species are collapsed within each clade (correct for species-level projection), but **different internal nodes can produce identical species frozensets**, which are deduplicated by the Python `set`. This is intentional and correct for RF-like comparison. However, for very large multi-copy families (e.g., transposable element families), the species-level collapse can yield very few distinct bipartitions, producing artificially low distances. No code fix needed, but a **documentation caveat** about multi-copy saturation would help users interpret results on highly duplicated families.

**[LOW] [C] [Phylo_Tracer.py:handle_mulrf_distance:449]**: The guard `if cli_args.input_GF_list and cli_args.input_imap and cli_args.mode` is redundant because all three arguments are `required=True` in argparse — argparse will error before this function is reached if any are missing. Not harmful, just dead code. **Suggested fix:** Remove the guard or replace with a simple assertion.

**[SUGGESTION] [C] [MulRF_Distance.py:mulrf_main:269]**: `sp_tree = Tree(species_tree_path)` is called once and then reused across all gene trees. Since `prune()` is called on a `copy()` (line 115-131), the original `sp_tree` is safe. Correct — no issue, just confirming.

---

## D. Output Descriptions

**[MEDIUM] [D] [README.md:258 + docs/analyze.rst]**: As noted in B above, the documented output path (`mulrf_distance/mulrf_distance.tsv`) does not match actual code behavior (`mulrf_distance.tsv` in CWD or `--output` value). The column list in the README (lines 366-379) **does** accurately match the code's `columns` list (lines 311-325).

No other output description issues found — column names, types, and ordering are consistent between code and README.

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 1     |
| MEDIUM     | 5     |
| LOW        | 4     |
| SUGGESTION | 1     |
| **Total**  | **11**|

### Top priorities:
1. **CRITICAL**: Add MulRF_Distance section to `docs/analyze.rst`
2. **MEDIUM**: Fix output path mismatch in README (or add `--output_dir` support)
3. **MEDIUM**: Log exceptions instead of silently swallowing them
4. **MEDIUM**: Clarify that this is a species-level projected RF, not formal MulRF

==================================================================
## Module: OrthoFilter_LB
==================================================================
Good — `output_dir` is handled generically by `chdir` before calling the handler. Now I have everything needed for the review.

---

# OrthoFilter_LB — Code Review Report

## A. Scientific Correctness

### [CRITICAL] [Scientific] [README.md:OrthoFilter_LB §2 SRBR formula]
**SRBR formula in docs does not match code implementation.**

The README states:
> SRBR = (Branch Length − Sister Branch Length) / Sister Branch Length

This implies SRBR uses terminal branch lengths only (`leaf.dist`). However, the code (`OrthoFilter_LB.py:297-317`) uses **root-to-tip distances**:

```python
tip_root_distance = pruned_tree.get_distance(leaf)          # root-to-tip
srbr_value = (tip_root_distance - sister_avg) / sister_avg   # sister_avg is also root-based
```

The code comment on line 316 even confirms: `"SRBR is also defined on root-to-tip distances."` The README formula should be updated to match the code, using "Root-to-tip distance" and "Sister average root-to-tip distance" instead of "Branch Length" and "Sister Branch Length".

**Fix:** Update the README SRBR formula and the "Where" definitions to use root-to-tip distance terminology consistently with the code.

---

### [CRITICAL] [Scientific] [OrthoFilter_LB.py:105-108 `get_average_node_length`]
**`sister.dist` is double-counted when computing average root-to-tip distance for an internal sister subtree.**

`get_average_node_length(sister)` computes:
```python
mean(sister.get_distance(leaf) + sister.dist for leaf in sister)
```
This equals `mean(sister_to_leaf) + sister.dist`, i.e. the mean distance from **sister's parent** to each leaf.

Then on line 314:
```python
sister_avg = get_average_node_length(sister) + pruned_tree.get_distance(sister)
```
`pruned_tree.get_distance(sister)` = root → sister (already includes `sister.dist`).

Total = `mean(sister_to_leaf) + sister.dist + root_to_sister_parent + sister.dist` = correct value **+ one extra `sister.dist`**.

This inflates `sister_avg`, making SRBR smaller than it should be, causing the filter to be slightly **less aggressive** than intended for tips with internal sisters.

**Fix:** Either remove `+ subtree.dist` from `get_average_node_length`, or replace line 314 with:
```python
sister_avg = mean(pruned_tree.get_distance(leaf) for leaf in sister)
```

---

## B. Code vs Docs Consistency

### [MEDIUM] [Consistency] [docs/analyze.rst:134-152]
**`--lb_mode` parameter missing from analyze.rst.**

The CLI parser (Phylo_Tracer.py:197) defines `--lb_mode` with choices `[or, and]` and default `or`. The README documents it. But analyze.rst lists only `--visual` as optional and omits `--lb_mode` entirely.

**Fix:** Add `--lb_mode` to the Optional Parameters in analyze.rst section 5.

---

### [MEDIUM] [Consistency] [docs/analyze.rst:134-152]
**`--output_dir` parameter missing from analyze.rst.**

The CLI parser (Phylo_Tracer.py:199) defines `--output_dir`. The README documents it. But analyze.rst omits it.

**Fix:** Add `--output_dir` to the Optional Parameters in analyze.rst section 5.

---

### [LOW] [Consistency] [docs/analyze.rst:143-144]
**`--rrbr_cutoff` and `--srbr_cutoff` listed as "Required Parameters" but have defaults.**

In the CLI parser, both are `required=True` with `default=5` / `default=2.5`. Since they are `required=True`, listing them as required is technically correct. However, this is confusing — a required parameter with a default is unusual. Consider either making them optional (remove `required=True` since defaults exist) or clarifying in the docs that they are required but have sensible defaults.

**Fix:** Clarify in the docs or remove `required=True` from the parser since defaults are provided.

---

## C. Code Logic and Correctness

### [CRITICAL] [Logic] [OrthoFilter_LB.py:314]
**Python operator precedence bug — `or 1e-6` guard may not work as intended.**

```python
sister_avg = get_average_node_length(sister) + pruned_tree.get_distance(sister) or 1e-6
```

Due to Python precedence, this evaluates as:
```python
sister_avg = (get_average_node_length(sister) + pruned_tree.get_distance(sister)) or 1e-6
```

This only triggers the `1e-6` fallback when the sum is exactly `0` (falsy). A negative sum (possible if tree data is malformed) would pass through and could cause unexpected SRBR values. More importantly:

**Fix:** Use explicit guarding:
```python
sister_avg = get_average_node_length(sister) + pruned_tree.get_distance(sister)
if sister_avg == 0:
    sister_avg = 1e-6
```

---

### [MEDIUM] [Logic] [OrthoFilter_LB.py:311-317]
**Division by zero when sister is a leaf with root-to-tip distance of 0.**

When `sister.is_leaf()` is `True`, `sister_avg = pruned_tree.get_distance(sister)`. If this is `0`, line 317 divides by zero. The `or 1e-6` guard only exists for the internal-sister branch (line 314), not for the leaf-sister branch (line 312).

**Fix:** Add a zero guard after line 312:
```python
if sister.is_leaf():
    sister_avg = pruned_tree.get_distance(sister) or 1e-6
```

---

### [LOW] [Logic] [OrthoFilter_LB.py:440-448]
**PDF cleanup assumes CWD — fragile when `--output_dir` changes CWD.**

`generate_tree_pdf` writes `{tree_id}_before.pdf` and `{tree_id}_after.pdf` to CWD (line 206), and cleanup on lines 447-448 also references CWD. This works because `--output_dir` changes CWD before the handler runs. However, if the function is ever called from a different context (e.g., tests, library use), the temp PDFs will litter an unexpected directory.

**Fix (suggestion):** Use `tempfile` or write temp PDFs to `pdf_dir` directly.

---

### [LOW] [Logic] [OrthoFilter_LB.py:44-63]
**`has_multiple_copies` is defined but never called anywhere in the module.**

Dead code. It's not used in `remove_long_branches` or `prune_main_LB`.

**Fix:** Remove it, or document it as a utility for external callers.

---

## D. Output Descriptions

### [MEDIUM] [Output] [docs/analyze.rst:134-152]
**No output files described for OrthoFilter_LB in analyze.rst.**

The README mentions `orthofilter_lb/pruned_tree/` and `orthofilter_lb/delete_gene/` in the results summary, but analyze.rst has no "Output" section for this module. The per-gene `.delete_gene.tsv` file columns (`tre_ID`, `gene`, `root_relative_branch_ratio`, `sister_relative_branch_ratio`) and the pruned Newick files are not documented in the manual.

**Fix:** Add an Output section in analyze.rst describing the two output directories and the TSV columns.

---

### [LOW] [Output] [README.md:259]
**`visual/` output directory not mentioned in README results summary.**

When `--visual` is used, a `orthofilter_lb/visual/` directory is created with merged before/after PDFs. The README results section only lists `pruned_tree/` and `delete_gene/`.

**Fix:** Add `orthofilter_lb/visual/` to the README results listing.

---

## Summary

| Severity   | Count |
|------------|-------|
| CRITICAL   | 3     |
| MEDIUM     | 4     |
| LOW        | 4     |
| **Total**  | **11**|

The three critical issues are:
1. SRBR formula in README contradicts the code (root-to-tip vs terminal branch)
2. `sister.dist` double-counted in `get_average_node_length` inflating sister averages
3. Operator precedence bug with `or 1e-6` on line 314

==================================================================
## Module: OrthoFilter_Mono
==================================================================
Now I have everything needed for a comprehensive review. Here is the full report:

---

# OrthoFilter_Mono — Code Review Report

## A. Scientific Correctness

### [CRITICAL] [A] README.md OrthoFilter_Mono §2 "Phylogenetic Distance Score" formula: **Formula does not match code**

README states:
$$\text{PhyloDist}=\text{Depth}(\text{MRCA}(\textbf{alien}))-\text{Depth}(\text{MRCA}(\text{target}\cup\text{alien}))$$

Code (`OrthoFilter_Mono.py:506`) computes:
```python
phylo_distance = brass_depth - union_depth
```
where `brass_depth = Depth(MRCA(target))`, not `Depth(MRCA(alien))`.

**Fix:** Change the README formula's first operand from `MRCA(alien)` to `MRCA(target)`:
$$\text{PhyloDist}=\text{Depth}(\text{MRCA}(\textbf{target}))-\text{Depth}(\text{MRCA}(\text{target}\cup\text{alien}))$$

---

### [MEDIUM] [A] OrthoFilter_Mono.py:2–5 module docstring: **Brassicaceae-specific language in generic module**

The docstring says *"detects dominant Brassicaceae lineages"*. The module is fully generic — it works on any user-supplied clade labels.

**Fix:** Replace "Brassicaceae" with "target-clade" or remove the family-specific reference.

---

### [LOW] [A] OrthoFilter_Mono.py:483,488 variable names: **Brassicaceae-specific naming in generic code**

`brass_species_set` and `brass_depth` refer to the target clade generically but are named after Brassicaceae.

**Fix:** Rename to `target_species_set` and `target_depth`.

---

## B. Code vs Docs Consistency

### [MEDIUM] [B] docs/analyze.rst:155–175: **`--output_dir` not documented**

The CLI defines `--output_dir` (`Phylo_Tracer.py:210`) and it is documented in README, but `analyze.rst` section 6 omits it entirely.

**Fix:** Add `--output_dir` to the Optional Parameters list in the analyze.rst OrthoFilter_Mono section.

---

### [MEDIUM] [B] docs/analyze.rst:163 vs Phylo_Tracer.py:204: **`--input_taxa` help text inconsistent**

- analyze.rst: *"File with taxonomic information for species."*
- Phylo_Tracer.py: *"Two-column mapping file (gene_id<TAB>clade_or_lineage_label)"*
- README: *"Two-column mapping file (gene_id<TAB>clade_or_lineage_label)"*

The analyze.rst description is vague and doesn't match the more precise CLI/README wording.

**Fix:** Update analyze.rst to: `Two-column mapping file (gene_id<TAB>clade_or_lineage_label)`.

---

### [MEDIUM] [B] OrthoFilter_Mono.py:975–976 / README / analyze.rst: **Hardcoded `dominant_purity=0.9` and `min_tips_per_clade=2` not documented**

`prune_main_Mono` calls `prune_all_clades` with `dominant_purity=0.9` and `min_tips_per_clade=2`. These are not exposed as CLI parameters and not mentioned anywhere in docs. Users have no way to know:
- That a subtree needs ≥90% purity to qualify as a "dominant lineage"
- That clades with fewer than 2 tips in a tree are skipped entirely

**Fix:** Either expose these as optional CLI parameters (recommended), or at minimum document them in README/analyze.rst as fixed internal thresholds.

---

### [LOW] [B] README.md "PhyloTracer Results Files" section: **`visual/` output directory not documented**

The code creates `orthofilter_mono/visual/` when `--visual` is set (`OrthoFilter_Mono.py:928`), but the results section only lists `pruned_tree/` and `delete_gene/`.

**Fix:** Add `orthofilter_mono/visual/` to the results list.

---

## C. Code Logic and Correctness

### [CRITICAL] [C] OrthoFilter_Mono.py:1018 `__main__` block: **`gene_id_transfer` not imported**

The `__main__` block calls `gene_id_transfer(args.input_imap)` but `gene_id_transfer` is not imported in this module. This will raise a `NameError` at runtime.

**Fix:** Add `gene_id_transfer` to the import from `phylotracer`:
```python
from phylotracer import (
    gene_id_transfer,
    num_tre_node,
    read_and_return_dict,
    rename_input_tre,
    stable_color_for_label,
)
```

---

### [MEDIUM] [C] OrthoFilter_Mono.py:671–673 `prune_all_clades`: **No guard against pruning all leaves**

If `global_remove` covers every leaf in the tree, `keep` becomes an empty set and `t.prune(keep)` will crash (ETE raises a `ValueError` or similar).

**Fix:** Add a guard:
```python
if global_remove:
    keep = set(t.get_leaf_names()) - global_remove
    if not keep:
        return t  # or log a warning
    t.prune(keep, preserve_branch_length=True)
```

---

### [MEDIUM] [C] OrthoFilter_Mono.py:526–531 `score_alien_candidates`: **Single-candidate scoring always yields 0.0**

When there is exactly one candidate, `normalize_values` returns `[0.0]` for both `phylo_norm` and `deep_norm`, making `combined_score = 0.0 * 0.0 * coverage_term = 0.0`. While the candidate can still be selected in `select_alien_tips_for_removal` (the score just affects sort order), a score of 0.0 is misleading in log output — it suggests no signal rather than "only candidate."

**Fix:** Consider returning `[1.0]` when there is exactly one candidate (or `[0.5]`), or document this behavior.

---

### [LOW] [C] OrthoFilter_Mono.py:757–773 log writing: **Log header/data semantic mismatch**

The TSV header defines `candidate_id`, `candidate_type`, `candidate_leaf_count` columns, but the data rows always write `removed_gene`, `"tip"`, and `1` respectively. When a candidate is a multi-leaf subtree MRCA, the per-tip logging loses the candidate-node identity.

**Fix:** Either (a) change the header to match per-tip semantics (`removed_gene_id`, `removed_type`, ...), or (b) also log the candidate node ID in a separate column.

---

## D. Output Descriptions

### [MEDIUM] [D] README.md OrthoFilter_Mono section / docs/analyze.rst: **Delete-gene TSV columns not documented**

The code writes a detailed per-tree TSV (`{tre_ID}.delete_gene.tsv`) with 16 columns (lines 954–959), but neither README nor analyze.rst documents the column names or their meanings.

**Fix:** Add an "Output" subsection documenting the TSV columns, similar to how `Phylo_Rooter` documents `stat_matrix.csv`.

---

## Summary

| Severity | Count |
|----------|-------|
| CRITICAL | 2 |
| MEDIUM | 7 |
| LOW | 3 |
| SUGGESTION | 0 |
| **Total** | **12** |

**Most urgent fixes:**
1. README PhyloDist formula uses `MRCA(alien)` but code uses `MRCA(target)` — scientifically misleading
2. `gene_id_transfer` missing import in `__main__` — runtime crash when module is run standalone

==================================================================
## Module: Ortho_Retriever
==================================================================
Good - `output_dir` is handled globally via `os.chdir`. Now I have all the information needed.

---

# Ortho_Retriever Code Review

## A. Scientific Correctness

**[MEDIUM] [SCIENTIFIC] [Ortho_Retriever.py:offcut_tre, lines 143-166]:** Average length computed by dividing by `sps_num` (species count) instead of `seq_num` (sequence count). The variable is named `ev_seqs1_avg_len` suggesting a per-sequence average, but it computes sum/species_count. In the specific context where this is reached (`sps_num1 == sps_num2` AND `seq_num1 == seq_num2`), the comparison result is unaffected because both sides divide by the same denominator. However, if the code is ever refactored to use these averages independently, the semantics are wrong. **Suggested fix:** divide by `seq_num` (i.e., `len(ev_seqs1)`) and rename variable for clarity.

**[LOW] [SCIENTIFIC] [Ortho_Retriever.py:count_sps_num, line 72]:** Species identity is inferred from the first 3 characters of renamed gene IDs (`gene[0:3]`). This is consistent with PhyloTracer's internal renaming scheme but is undocumented in user-facing docs. If the renaming scheme ever changes, this silently breaks. **Suggested fix:** add an inline comment or docstring noting this depends on the 3-character species prefix from `gene_id_transfer`.

**[LOW] [SCIENTIFIC] [Ortho_Retriever.py:offcut_tre, line 186]:** Offcut groups with fewer than 2 species are silently discarded. This is a reasonable filter but is not documented anywhere. Users cannot control this threshold. **Suggested fix:** document this minimum-species-2 filter in the docs, or expose as a parameter if flexibility is desired.

**No issues found** with the overall algorithmic approach: recursive paralog splitting by comparing species count → sequence count → average length is a standard heuristic for ortholog retrieval from gene family trees.

## B. Code vs Docs Consistency

**[MEDIUM] [DOCS-CONSISTENCY] [docs/analyze.rst:301-315]:** The `--output_dir` parameter is missing from the analyze.rst Ortho_Retriever section. It is defined in the CLI parser (Phylo_Tracer.py:297) and documented in the README (line 675). **Suggested fix:** add `--output_dir` to the Optional Parameters section in analyze.rst.

**[LOW] [DOCS-CONSISTENCY] [docs/analyze.rst:301-315]:** No "Output" section for Ortho_Retriever. Other modules in analyze.rst (e.g., Phylo_Rooter at line 116) describe their output files. The module produces `ortho_retriever_summary.txt` and `ortholog_trees.tsv`. **Suggested fix:** add an Output subsection listing both files and their column descriptions.

**[LOW] [DOCS-CONSISTENCY] [README.md:667-678 vs docs/analyze.rst:301-315]:** README documents the default output directory as `ortho_retriever` (via `DEFAULT_OUTPUT_DIRS` in Phylo_Tracer.py:938), but neither README nor analyze.rst explicitly mentions this default. **Suggested fix:** note that output defaults to `ortho_retriever/` subdirectory when `--output_dir` is not specified.

All three CLI arguments (`--input_GF_list`, `--input_imap`, `--input_gene_length`) are consistent across code, README, and analyze.rst in terms of names, types, and required status.

## C. Code Logic and Correctness

**[CRITICAL] [CODE-LOGIC] [Ortho_Retriever.py:rm_dup, lines 210-218]:** The function modifies the list it is iterating over. When `paralogs_L.remove(x)` is called during `for ev_seqs1 in paralogs_L`, Python's list iterator shifts indices, causing elements to be skipped. This can result in incomplete deduplication — subset relationships may survive.

```python
# Current (buggy):
for ev_seqs1 in paralogs_L:
    for ev_seqs2 in paralogs_L:
        if ev_seqs1.issubset(ev_seqs2):
            if ev_seqs1 in paralogs_L:
                paralogs_L.remove(ev_seqs1)
```

**Suggested fix:** build a removal set, then filter:
```python
def rm_dup(paralogs_L):
    to_remove = set()
    for i, s1 in enumerate(paralogs_L):
        for j, s2 in enumerate(paralogs_L):
            if i != j:
                if s1.issubset(s2):
                    to_remove.add(i)
                elif s2.issubset(s1):
                    to_remove.add(j)
    return [s for i, s in enumerate(paralogs_L) if i not in to_remove]
```

**[CRITICAL] [CODE-LOGIC] [Ortho_Retriever.py:rm_dup, line 212]:** When `ev_seqs1 == ev_seqs2` (same set object, i.e., self-comparison), `ev_seqs1.issubset(ev_seqs2)` is `True`, so the function removes the set from the list on self-comparison. This means every unique set that appears exactly once will be removed during the self-comparison iteration. Combined with the iteration-mutation bug above, the behavior is unpredictable.

**Suggested fix:** add `if ev_seqs1 is not ev_seqs2:` guard (or use index-based iteration with `i != j`).

**[MEDIUM] [CODE-LOGIC] [Ortho_Retriever.py:iterator, line 255]:** Parameter named `new_named_gene2gene_dic` (suggesting renamed→original mapping) but actually receives `gene2new_named_gene_dic` (original→renamed) from `get_single_copy_trees` line 409. The function works correctly because `rename_input_tre` on line 296 expects original→renamed, but the misleading name is a maintenance hazard. **Suggested fix:** rename parameter to `gene2new_named_gene_dic` to match what is actually passed.

**[MEDIUM] [CODE-LOGIC] [Ortho_Retriever.py:rename_OGs_tre_name, lines 342-363]:** When `principal_gene_S` is empty (all leaves are offcuts), it is still added to the output list with a name like `TREEID_T1_0`. Downstream, `extract_tree` would call `Phylo_t.prune([])` on an empty set, which raises an ETE3 error. **Suggested fix:** guard with `if not principal_gene_S: return ordered_name_OG_L` early, or skip appending empty sets.

**[MEDIUM] [CODE-LOGIC] [Ortho_Retriever.py:rename_OGs_tre_name, lines 344-346]:** `sps_num_L` computes `len(OG)` for each minor ortholog group — this counts total genes, not species. The variable name `sps_num_L` is misleading. This propagates into the sorting (line 348-352) and naming (line 361-363). The sort and output naming use this as a "count" suffix, which is gene count, not species count. **Suggested fix:** either rename to `gene_count_L` or actually compute species counts via `count_sps_num(OG)`.

**[LOW] [CODE-LOGIC] [Ortho_Retriever.py:offcut_tre, line 130]:** The check `if None in ev_seqs1 or None in ev_seqs2` guards against None values from the split, but `"".split(",")` returns `['']` not `[None]`, so this guard never triggers. The None check is dead code. **Suggested fix:** check for empty strings instead: `if '' in ev_seqs1 or '' in ev_seqs2`.

**[LOW] [CODE-LOGIC] [Ortho_Retriever.py:split_main, line 466]:** `len(Phylo_t0)` on an ETE tree returns the number of leaves. `len(Phylo_t0) == 0` would only occur for an entirely empty tree object, which `read_phylo_tree` would likely fail on before reaching here. The `len(Phylo_t0.children) == 0` check is more meaningful (single-node tree). Not harmful, but the `len == 0` check is effectively unreachable. **Suggested fix:** simplify to `if len(Phylo_t0.children) == 0`.

## D. Output Descriptions

**[MEDIUM] [OUTPUT-DOCS] [docs/analyze.rst:301-315]:** No output files are described for Ortho_Retriever. The module produces two files:
- `ortho_retriever_summary.txt` — columns: `tre_name`, `single_copy_tree` (Newick string)
- `ortholog_trees.tsv` — columns: `tre_ID`, `OG_1`, `OG_2`, ... (comma-separated leaf names per ortholog group)

**Suggested fix:** add an Output section to analyze.rst matching the format used by Phylo_Rooter (line 116-119).

**[LOW] [OUTPUT-DOCS] [README.md:264]:** README lists `ortho_retriever_summary.txt` and `ortholog_trees.tsv` in the Results Files section but does not describe column contents. **Suggested fix:** add column descriptions inline or reference the analyze.rst documentation.

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 2     |
| MEDIUM     | 5     |
| LOW        | 6     |
| SUGGESTION | 0     |
| **Total**  | **13** |

The two CRITICAL issues are both in `rm_dup` — the list-mutation-during-iteration bug and the self-comparison bug. Together they can cause incorrect ortholog group deduplication, potentially producing redundant or missing ortholog sets in downstream output. These should be fixed before any production use.

==================================================================
## Module: Phylo_Rooter
==================================================================
Now I have everything. Here is the complete review:

---

# Phylo_Rooter Code Review

## A. Scientific Correctness

### No issues found.

The six metrics (OD, BLV, GD, SO, GDC, MulRF), their direction classifications (cost vs. benefit), the direction-aware min-max normalization, the EWM algorithm, and the composite scoring formula are all scientifically accurate and consistent with the literature. Default weights giving higher emphasis to OD and GD (0.30 each) are biologically reasonable (parsimony-driven duplication minimization and outgroup depth are the most informative rooting signals). The neutral value of 0.5 for degenerate metrics (max == min) is a sound choice. The dynamic basal-set logic correctly walks up the species tree to find the deepest divergence point covered by the gene tree. No claims in the docs appear overstated.

---

## B. Code vs Docs Consistency

**[LOW] [CODE-VS-DOCS] [docs/analyze.rst:line 283, Phylo_Tracer.py:line 283]:** The `--node_count_mode` parameter for `GD_Loss_Tracker` is defined in `Phylo_Tracer.py` (line 283) but is **not documented** in `docs/analyze.rst` (section 11). The docs only mention `--target_species`, `--mrca_node`, and `--include_unobserved_species` as optional parameters.

> **Suggested fix:** Add `--node_count_mode` to the GD_Loss_Tracker optional parameters in `docs/analyze.rst`.

*Note: This is outside Phylo_Rooter proper, but was found during the cross-check.*

### Phylo_Rooter-specific: No issues found.

All three required parameters (`--input_GF_list`, `--input_imap`, `--input_sps_tree`) and all three optional parameters (`--weights`, `--weight-strategy`, `--output_dir`) are consistent across code, `docs/analyze.rst`, and `README.md`. Default values match exactly: weights default `[0.30, 0.10, 0.30, 0.10, 0.10, 0.10]`, weight-strategy default `empirical`. Parameter types (6 floats for `--weights`, choice for `--weight-strategy`) are correctly described. Required vs optional labeling is accurate in all three locations.

---

## C. Code Logic and Correctness

**[MEDIUM] [CODE-LOGIC] [Phylo_Rooter.py:root_main, line 834]:** The `tree_id` column is computed via `lambda x: "_".join(x.split("_")[:-1])` to strip the candidate suffix from the `Tree` column (e.g., `"OG0001_3"` → `"OG0001"`). However, if the original tree ID itself contains underscores (e.g., `"OG_0001"`), this splits incorrectly (`"OG_0001_3"` → `"OG_0001"` is actually fine, but `"OG_00_01_3"` → `"OG_00_01"` is also fine). The logic is correct for stripping the last `_N` suffix. However, the `tree_id` column is computed but **never written to the output CSV** — it is used only for sorting. The `cols` list on line 835 does not include `tree_id`, so the sort uses a column that is silently discarded. This is not a bug per se, but it makes the sort by `tree_id` cosmetic-only (ordering rows in `stat_matrix.csv` by gene family). **No real issue** — just slightly dead logic since only the best candidate per tree is kept anyway (one row per tree), making the sort trivially stable.

**[LOW] [CODE-LOGIC] [Phylo_Rooter.py:calculate_species_overlap_gd_num, line 482]:** If the largest duplication node has fewer than 2 children (a polytomy or degenerate node), `children[1]` / `children[0]` could still succeed for a 3+ child case but the function only examines children 0 and 1, ignoring additional children. This is consistent with the bifurcating assumption stated in the docstring and the `resolve_polytomy` call at line 795, so it is safe in practice.

**[LOW] [CODE-LOGIC] [Phylo_Rooter.py:root_main, line 701-702]:** When `dir_path != cwd`, the code does `shutil.rmtree(dir_path)` on an existing `rooted_trees/` directory without warning. If a user re-runs the command, previously rooted trees are silently deleted. This is intentional overwrite-on-rerun behavior but could surprise users.

> **Suggested fix:** Consider logging a warning before `rmtree`, e.g., `logger.warning("Removing existing output directory: %s", dir_path)`.

**[LOW] [CODE-LOGIC] [Phylo_Rooter.py:_tip_distance_mad, line 369]:** If all branch lengths are zero (e.g., topology-only tree), `get_distance` returns 0 for all leaves, MAD = 0.0, and all candidates tie. The function handles this gracefully (returns 0.0), but users should be aware that MAD/MinVar candidates are uninformative for trees without branch lengths. This is an inherent limitation, not a bug.

---

## D. Output Descriptions

### No issues found.

Both `docs/analyze.rst` (lines 116-119) and `README.md` correctly describe the two outputs:
- `rooted_trees/` — directory of rooted gene trees in Newick format.
- `stat_matrix.csv` — per-tree scoring table.

The column names documented (`Tree, score, deep (OD), var (BLV), GD, species_overlap (SO), gd_consistency (GDC), RF (MulRF)`) match exactly the `cols` list at `Phylo_Rooter.py:835`:
```python
cols = ["Tree", "score", "deep", "var", "GD", "species_overlap", "gd_consistency", "RF"]
```

The parenthetical aliases in the docs (e.g., "deep (OD)") correctly map to column names.

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 0     |
| MEDIUM     | 0*    |
| LOW        | 3     |
| SUGGESTION | 0     |

\* The one MEDIUM-level finding (dead `tree_id` sort column) is functionally harmless since only one row per tree is ever stored in `stat_matrix`.

**Overall assessment:** The Phylo_Rooter module is well-implemented and thoroughly documented. Code, CLI definitions, docs, and README are highly consistent. The scientific claims are accurate and appropriately scoped. The one cross-module documentation gap (`--node_count_mode` missing from `analyze.rst` for GD_Loss_Tracker) is outside the Phylo_Rooter scope. No critical or high-severity issues were found.

==================================================================
## Module: PhyloSupport_Scaler
==================================================================
Here is the complete review of the **PhyloSupport_Scaler** module.

---

## A. Scientific Correctness

**[CRITICAL] [Scientific] [PhyloSupport_Scaler.py:55 `scale_support`]: Boundary condition at `max_support == 1.0` causes incorrect scale detection**

The heuristic `if max_support < 1` to determine whether support values are in [0,1] vs [0,100] fails when the maximum support is exactly `1.0`. A tree with posterior probabilities in [0,1] where the max is 1.0 will be treated as if it's already in [0,100] scale. Requesting `--scale_to 1` would then divide all values by 100, producing `0.01` instead of leaving them unchanged. Conversely, `--scale_to 100` would log "already in range" and do nothing.

**Suggested fix:** Change the threshold to `max_support <= 1` (i.e., `if max_support <= 1:`), or better, use a threshold like `max_support <= 1.0 + epsilon` to account for float imprecision. Alternatively, document that the user must know their current scale and the function will not auto-detect — but then the auto-detection logic should be removed entirely.

---

**[LOW] [Scientific] [PhyloSupport_Scaler.py:48 `scale_support`]: Auto-detection heuristic cannot distinguish all-zero support from [0,1] scale**

If all support values are 0 (uninformative/missing), `max_support == 0 < 1`, so the code classifies them as [0,1] scale. Scaling to 100 multiplies all values by 100 but `0 * 100 = 0`, so no data corruption occurs. However, the silent assumption is misleading and a warning would be appropriate.

**Suggested fix:** Add a warning when `max_support == 0`.

---

**[SUGGESTION] [Scientific] [PhyloSupport_Scaler.py:25 `scale_support`]: Docstring assumption is imprecise**

The docstring states "Support values are either in [0, 1] or [0, 100]" but the detection logic can only distinguish "max < 1" from "max >= 1" — it cannot handle intermediate cases (e.g., a Bayesian tree where all PP values happen to be < 1 but the tree is intended to be on a 0-100 scale with low support). This is inherent to the heuristic approach and should be explicitly documented as a limitation.

---

## B. Code vs Docs Consistency

**[MEDIUM] [Docs-vs-Code] [docs/analyze.rst section 2 "PhyloSupport_Scaler"]: `--output_dir` parameter missing from analyze.rst**

The CLI parser (Phylo_Tracer.py:154) defines `--output_dir` for PhyloSupport_Scaler, and the README documents it, but docs/analyze.rst does not list it under Optional Parameters.

**Suggested fix:** Add to analyze.rst:
```
Optional Parameters:
    - ``--output_dir``      Output directory (default: current working directory).
```

---

**[MEDIUM] [Docs-vs-Code] [docs/analyze.rst section 2 + README "PhyloSupport_Scaler"]: Output directory/files not documented**

Neither docs/analyze.rst nor the README PhyloTracer Results Files section mentions the output directory `support_scaler_tree/` or describes the output files (rescaled Newick tree files). Users have no way to know where their results go or what format they are in.

**Suggested fix:** Add an Output section to both docs/analyze.rst and README:
```
Output:
    - ``support_scaler_tree/``  Directory containing rescaled gene trees in Newick format.
```

---

**[LOW] [Docs-vs-Code] [docs/analyze.rst vs README]: `--scale_to` description phrasing inconsistency**

- docs/analyze.rst uses directional language: *"Specify '1' to scale support values from 1-100 to 0-1, or '100' to scale from 0-1 to 1-100"*
- README uses target language: *"Target support scale: '1' for [0,1], '100' for [0,100]"*

Both are correct but could confuse users comparing the two sources. The README phrasing is cleaner and matches the CLI help text.

**Suggested fix:** Align analyze.rst to match the README/CLI phrasing.

---

**[LOW] [Docs-vs-Code] [docs/analyze.rst section 2]: `--scale_to` is documented as Required but not explicitly flagged as `required=True`-style in the RST layout**

Both `--input_GF_list` and `--scale_to` appear under "Required Parameters" in analyze.rst, which matches the code (`required=True` on both). No actual inconsistency — just noting this is correct.

**No issues found** regarding parameter types, choices, or default values — `choices=['1', '100']` in code matches both docs.

---

## C. Code Logic and Correctness

**[MEDIUM] [Logic] [PhyloSupport_Scaler.py:106-107 `support_scaler_main`]: `shutil.rmtree` silently destroys previous output**

If the output directory `support_scaler_tree/` already exists, it is unconditionally deleted via `shutil.rmtree(dir_path)` before writing new results. This can destroy previous results without any warning or confirmation.

**Suggested fix:** Either warn the user before deleting, or use a timestamped/versioned directory, or skip deletion and only overwrite individual files.

---

**[MEDIUM] [Logic] [PhyloSupport_Scaler.py:48 `scale_support`]: Empty tree traversal edge case**

If `phylo_tree.traverse()` yields zero nodes (degenerate/empty tree), `max()` on an empty sequence raises `ValueError`, which is caught by the generic `except Exception` and re-raised as `RuntimeError("Failed to get support values: ...")`. While it doesn't crash, the error message is misleading — the real issue is an empty tree, not a failure to read support.

**Suggested fix:** Add an explicit check: `nodes = list(phylo_tree.traverse()); if not nodes: raise ValueError("Tree has no nodes")`.

---

**[LOW] [Logic] [PhyloSupport_Scaler.py:48 vs 58]: Inconsistent support attribute access pattern**

Line 48 uses `getattr(node, "support", 0)` (with a fallback default), while line 58-68 uses `if hasattr(node, "support")` before accessing `node.support`. In ete3, all `TreeNode` objects have a `support` attribute (default 1.0), so `hasattr` is always `True` and the guard is unnecessary. The inconsistency is cosmetic but confusing.

**Suggested fix:** Use a consistent pattern throughout — either always use `getattr` with default, or drop the `hasattr` guards since ete3 guarantees the attribute.

---

**[LOW] [Logic] [PhyloSupport_Scaler.py:120 `support_scaler_main`]: Error in individual tree processing is silently logged, not raised**

If `read_tree(tree_path)` fails (e.g., corrupt Newick file), the error is caught and logged at `logger.error` level, but execution continues to the next tree. This is reasonable for batch processing, but the user gets no summary of how many trees failed at the end.

**Suggested fix:** Track failed tree count and log a summary after the loop.

---

## D. Output Descriptions

**[MEDIUM] [Output] [docs/analyze.rst section 2 + README "PhyloTracer Results Files"]: No output description for PhyloSupport_Scaler**

The module writes rescaled trees to `support_scaler_tree/<tree_id>.nwk` (via `write_tree_to_newick`), but neither docs/analyze.rst nor the README Results Files section documents this output. Users cannot discover the output location or format from the documentation.

**Suggested fix:** Add to both locations:
- Output directory: `support_scaler_tree/`
- Output files: one Newick file per input tree, with rescaled support values
- Newick format: ete3 `format=0` (all attributes)

---

## Summary Table

| Severity | Count |
|----------|-------|
| CRITICAL | 1 |
| MEDIUM | 4 |
| LOW | 4 |
| SUGGESTION | 1 |
| **Total** | **10** |

The most urgent fix is the **CRITICAL boundary bug** at `max_support == 1.0` — this produces silently wrong results for any tree with posterior probabilities where the maximum is exactly 1.0, which is an extremely common case in Bayesian phylogenetics.

==================================================================
## Module: PhyloTree_CollapseExpand
==================================================================
Now I have all the information needed for a thorough review. Let me compile the findings.

---

## PhyloTree_CollapseExpand — Code Review Report

---

### A. Scientific Correctness

**[MEDIUM] [SCIENTIFIC] [docs/analyze.rst:22, README.md:387]**: **Collapse threshold semantics mismatch with code.**
- Docs (analyze.rst line 22) state: *"Nodes whose support is **less than or equal to** this value will be collapsed"*
- README (line 387) states: *"Node support cutoff used for collapsing internal branches"* (ambiguous)
- Code (`PhyloTree_CollapseExpand.py:52`) uses **strict less-than**: `if getattr(node, "support", 0) < node_support_value`
- **Impact**: A user setting `--support_value 50` expects nodes with support=50 to be collapsed (per docs), but the code retains them.
- **Fix**: Either change docs to say "less than" or change code to `<=`. The `<=` semantics is more conventional (bootstrap < threshold means poorly supported), so updating docs to "less than" is recommended since the code behavior (retaining nodes at exactly the threshold) is the more conservative and standard choice.

**[LOW] [SCIENTIFIC] [PhyloTree_CollapseExpand.py:24, Phylo_Tracer.py:146]**: **`node_support_value` typed as `int`, but support values can be floats.**
- Posterior probabilities (e.g., Bayesian PP in [0,1]) are floats. The function signature and argparse `type=bounded_int(0, 100)` restrict to integers.
- **Fix**: Accept `float` in both the function signature and the CLI parser (`bounded_float(0, 100)`) to support PP-scaled trees, or document that trees must first be scaled to integer range via `PhyloSupport_Scaler`.

No other scientific issues found.

---

### B. Code vs Docs Consistency

**[CRITICAL] [CONSISTENCY] [Phylo_Tracer.py:146 vs Phylo_Tracer.py:148 vs PhyloTree_CollapseExpand.py:62-66]**: **`--output_dir` exists in CLI but is not passed to `collapse_expand_main`.**
- The CLI defines `--output_dir` (line 148) and the dispatch infrastructure resolves it and `os.chdir`s into it. However, `collapse_expand_main` hardcodes its own output directory name `"collapse_expand_tree"` (line 90-95) and creates/deletes it internally. This means:
  - If a user passes `--output_dir /some/path`, the process `chdir`s there, but then `collapse_expand_main` creates a **nested** `collapse_expand_tree/` subdirectory inside it.
  - The hardcoded `shutil.rmtree` (line 98) silently deletes prior results without warning.
- **Fix**: Remove the hardcoded directory logic from `collapse_expand_main` and accept an `output_dir` parameter, or document the nested directory behavior.

**[MEDIUM] [CONSISTENCY] [Phylo_Tracer.py:146 vs docs/analyze.rst:22]**: **`--support_value` is `required=True` in CLI but docs say "default is 50".**
- The CLI parser sets `required=True` with no `default=` value (line 146), so the user **must** provide it.
- Docs (analyze.rst line 22) and the standalone `__main__` block (line 132) both say "default is 50".
- README (line 387) says "default = 50".
- **Fix**: Either add `default=50` and remove `required=True` in Phylo_Tracer.py, or remove the "default is 50" claim from docs. Given that the standalone script uses default=50, adding the default is the right approach.

**[MEDIUM] [CONSISTENCY] [docs/analyze.rst:21-22]**: **`--support_value` listed as "Required" in analyze.rst but has a stated default.**
- Parameters with defaults are conventionally optional. Listing it as required while stating "default is 50" is contradictory.
- **Fix**: Move to Optional if a default is added, or remove the default claim.

**[LOW] [CONSISTENCY] [docs/analyze.rst:25 vs README.md:389]**: **`--output_dir` is documented in README but missing from analyze.rst.**
- README line 392 shows `[--output_dir DIR]` in usage.
- analyze.rst has no mention of `--output_dir` for PhyloTree_CollapseExpand.
- **Fix**: Add `--output_dir` to the Optional Parameters section of analyze.rst.

No other consistency issues found.

---

### C. Code Logic and Correctness

**[MEDIUM] [LOGIC] [PhyloTree_CollapseExpand.py:97-98]**: **Unconditional `shutil.rmtree` on existing output directory is destructive.**
- If a user has an existing `collapse_expand_tree/` directory with other files, it is silently deleted before processing.
- **Fix**: Either warn before deletion, or use a timestamped subdirectory, or only overwrite individual output files rather than nuking the whole directory.

**[LOW] [LOGIC] [PhyloTree_CollapseExpand.py:62-66]**: **`collapse_expand_main` does not use `output_dir` parameter — function signature mismatch with CLI.**
- The function signature is `(tree_dict, support_value, revert)` but the handler at Phylo_Tracer.py:370 only passes these three. If the infrastructure `chdir`s into `output_dir`, the hardcoded directory creation still creates a subdirectory.
- Already covered in B above; noting the code-logic side here.

**[LOW] [LOGIC] [PhyloTree_CollapseExpand.py:103]**: **Empty `tree_dict` proceeds without warning.**
- If the input file is empty (no trees), the progress bar shows 0 total and the function silently creates an empty output directory.
- **Fix**: Add a warning log or early return for empty input.

**[SUGGESTION] [LOGIC] [PhyloTree_CollapseExpand.py:52]**: **`getattr(node, "support", 0)` treats missing support as 0 (collapse).**
- If a tree has no support values at all (e.g., ML tree without bootstraps), every internal node gets support=0 and is collapsed.
- This is defensible but could surprise users. Consider logging a warning if no nodes have support values.

No off-by-one errors, division-by-zero, or encoding issues found.

---

### D. Output Descriptions

**[LOW] [OUTPUT] [docs/analyze.rst, README.md]**: **No output section for PhyloTree_CollapseExpand.**
- Neither analyze.rst nor README describes what output files are generated. The code writes Newick files to `collapse_expand_tree/` using `write_tree_to_newick`, but this is not documented anywhere.
- **Fix**: Add an Output section stating that collapsed (and optionally re-expanded) trees are written as individual Newick files in the `collapse_expand_tree/` directory, one per input tree, named by tree ID.

**[LOW] [OUTPUT] [README.md "PhyloTracer Results Files" section]**: **PhyloTree_CollapseExpand is not listed in the results files overview.**
- The "PhyloTracer Results Files" section (README ~line 253-267) lists outputs for most modules but omits PhyloTree_CollapseExpand.
- **Fix**: Add `- Collapse/expand outputs: collapse_expand_tree/` to the results listing.

---

### Summary Table

| Severity | Count |
|----------|-------|
| CRITICAL | 1 |
| MEDIUM | 4 |
| LOW | 5 |
| SUGGESTION | 1 |
| **Total** | **11** |

**Key takeaways:**
1. The **most impactful issue** is the `--output_dir` / hardcoded directory mismatch (CRITICAL) — users passing `--output_dir` get unexpected nested directories and silent deletion of prior results.
2. The **collapse threshold semantics** (`<` vs `<=`) discrepancy between code and docs could cause reproducibility issues.
3. The `--support_value` required/default contradiction should be resolved.
4. Output documentation for this module is entirely absent.

==================================================================
## Module: TreeTopology_Summarizer
==================================================================
Now I have all the information needed. Here is the complete review.

---

# TreeTopology_Summarizer — Module Review

## A. Scientific Correctness

### Issues Found

**[MEDIUM] [SCIENTIFIC] [TreeTopology_Summarizer.py:group_trees_by_topology_with_ids, line 141]:**
`robinson_foulds()` on unrooted trees with different leaf sets will raise an error. The code groups species-only simplified trees, but trees with different species compositions (different leaf sets) will fail at the RF comparison. The code implicitly assumes all trees share the same leaf set, but `statistical_main` only filters out single-species trees — trees with different species subsets will still be compared.
**Suggested fix:** Add `unrooted_trees=True` parameter and handle the `TreeError` when leaf sets differ (either skip or use `correct_by_polytomy_size`), or explicitly filter to trees with identical leaf sets before grouping.

**[LOW] [SCIENTIFIC] [TreeTopology_Summarizer.py:group_trees_by_topology_with_ids, line 131–148]:**
The grouping algorithm compares all trees against the *first* tree in the list. This is correct for partitioning by RF=0, but the order-dependence means the representative topology for each group depends on input order. This is algorithmically fine but worth documenting — the first tree encountered becomes the "pivot" for each group.

**[SUGGESTION] [SCIENTIFIC] [TreeTopology_Summarizer.py:visualize_top_trees, line 456]:**
`convert_to_ultrametric(tree_length=1)` is applied before rendering. This silently distorts branch-length proportions. Since the code already sets `force_topology = False` (line 445, comment says "keep branch-length geometry"), there is a contradiction: the ultrametric conversion overrides the original branch lengths. Either set `force_topology = True` (topology-only view) or remove `convert_to_ultrametric` to preserve real branch lengths.

## B. Code vs Docs Consistency

### Issues Found

**[MEDIUM] [DOCS-CODE] [README.md:line 565 vs TreeTopology_Summarizer.py:visualize_top_trees, line 339]:**
README states: *"each topology panel is rendered at A4 width and merged in a single-column, top-to-bottom layout"*. The code actually uses a **2-column grid** layout (`cols = 2` on line 345), not single-column.
**Suggested fix:** Change README description to "2-column grid layout" or update code to match docs.

**[MEDIUM] [DOCS-CODE] [README.md:line 265 vs code output]:**
README "PhyloTracer Results Files" section says topology summaries produce *"merged PNG summaries"*. The code generates **merged PDF files** (`merge_relative_top{N}.pdf` and `merge_absolutely_top{N}.pdf`), not PNGs.
**Suggested fix:** Change "merged PNG summaries" to "merged PDF summaries" in README line 265.

**[LOW] [DOCS-CODE] [docs/analyze.rst:line 185 vs Phylo_Tracer.py:line 214]:**
analyze.rst describes `--input_GF_list` as *"File containing paths to gene tree files, one per line"*. The actual CLI help and README describe it as *"Tab-delimited mapping file (GF_ID<TAB>gene_tree_path)"*. The analyze.rst description omits the tab-delimited two-column format.
**Suggested fix:** Update analyze.rst to say "Tab-delimited mapping file (GF_ID<TAB>gene_tree_path); one gene tree per line."

**[LOW] [DOCS-CODE] [docs/analyze.rst:line 189 vs code]:**
analyze.rst does not document `--output_dir`, which exists in the CLI parser (Phylo_Tracer.py:217).
**Suggested fix:** Add `--output_dir` to the analyze.rst Optional Parameters section.

**[LOW] [DOCS-CODE] [README.md TreeTopology_Summarizer section vs code output files]:**
Neither README nor analyze.rst document the actual output file names: `absolute_topology.txt`, `relative_topology.txt`, `merge_absolutely_top{N}.pdf`, `merge_relative_top{N}.pdf`. The README just says "absolute_*.txt, relative_*.txt" generically.
**Suggested fix:** Document the exact output file names and their column headers (Topology_ID, Topology_num, Topology, Tree_IDs).

## C. Code Logic and Correctness

### Issues Found

**[CRITICAL] [LOGIC] [TreeTopology_Summarizer.py:visualize_top_trees, line 442]:**
`len(count)` is used where `count` is actually the list of `(tree_id, tree)` tuples (the variable name is misleading from the unpacking `for i, (tree_str, count) in enumerate(sorted_trees)`). The variable `count` is the list of trees, not a count integer. `len(count)` happens to give the correct number, but the `TextFace` label says `Count: {len(count)}` which is correct by accident. The real bug is the variable naming — `count` should be `trees_list` or similar. Not a functional bug but highly confusing and error-prone for maintenance.

**[MEDIUM] [LOGIC] [TreeTopology_Summarizer.py:statistical_main, line 527]:**
`exit(1)` is called on non-single-copy detection instead of `sys.exit(1)`. The bare `exit()` builtin is intended for interactive use and may behave differently in some environments. More importantly, the entire program terminates on the *first* non-single-copy tree rather than skipping it with a warning, which is very harsh for large-scale analyses.
**Suggested fix:** Use `sys.exit(1)` or, better, skip the offending tree with a warning and continue processing.

**[MEDIUM] [LOGIC] [TreeTopology_Summarizer.py:write_relative_summary_with_ids, line 214]:**
`if top_n:` evaluates to `False` when `top_n=0`, which is correct, but also when `top_n` is `None`. The CLI default is `10`, so `top_n` will always be truthy from the CLI path. However, when called programmatically with `top_n=0` intending "no visualization", this works correctly. The issue is that there is no way to disable visualization from the CLI since `--visual_top` defaults to 10 and `bounded_int(1)` requires >= 1.
**Suggested fix:** Consider adding a `--no_visual` flag or allowing `--visual_top 0` to disable visualization.

**[MEDIUM] [LOGIC] [TreeTopology_Summarizer.py:group_trees_by_topology_with_ids, line 140–141]:**
The first tree is compared against itself (RF distance to self = 0), so it always goes into `same_rf_trees`. This is functionally correct but wastes one RF computation per group. Not a bug, just unnecessary work.

**[LOW] [LOGIC] [TreeTopology_Summarizer.py:statistical_main, line 508–528]:**
If all trees are single-species (skipped at line 515–516), `trees_with_ids` will be empty. The downstream functions `write_absolute_summary` and `process_tree_with_ids` will run on an empty list. `write_absolute_summary` will write a header-only file (benign). `process_tree_with_ids` will not enter the while loop (benign). No crash, but no warning to the user that zero trees were processed.
**Suggested fix:** Add a warning/early return if `trees_with_ids` is empty after filtering.

**[LOW] [LOGIC] [TreeTopology_Summarizer.py:visualize_top_trees, line 428]:**
If `NodeStyle`, `TextFace`, or `TreeStyle` are `None` (import failed), the function will crash with `TypeError` at line 428 (`NodeStyle()`). There's no guard for missing ete3 visualization dependencies.
**Suggested fix:** Add an early check: `if NodeStyle is None: logger.warning(...); return`.

**[LOW] [LOGIC] [TreeTopology_Summarizer.py:write_relative_summary_with_ids, line 220]:**
`Topology_ID` starts at 0. This is consistent with `write_absolute_summary` (line 290) which also starts at 0. Fine, but docs should clarify that IDs are 0-indexed.

## D. Output Descriptions

### Issues Found

**[MEDIUM] [OUTPUT] [docs/analyze.rst + README.md]:**
Neither document describes the output file columns. The output files have columns: `Topology_ID`, `Topology_num`, `Topology`, `Tree_IDs`. These are not documented anywhere.
**Suggested fix:** Add an "Output" section to the TreeTopology_Summarizer docs describing each output file and its columns.

**[LOW] [OUTPUT] [README.md:line 265]:**
The output description doesn't mention the PDF visualization outputs (`merge_relative_top{N}.pdf`, `merge_absolutely_top{N}.pdf`).
**Suggested fix:** Add these to the results files section.

---

## Summary Table

| Severity | Count |
|----------|-------|
| CRITICAL | 0* |
| MEDIUM | 6 |
| LOW | 7 |
| SUGGESTION | 1 |
| **Total** | **14** |

*\*The `len(count)` issue at line 442 is functionally correct by accident, so I classified it LOW rather than CRITICAL.*

### Priority Fixes

1. **README "single-column" → "2-column grid"** and **"PNG" → "PDF"** — straightforward doc corrections
2. **RF distance on mismatched leaf sets** — add error handling in `group_trees_by_topology_with_ids` for trees with different species compositions
3. **Missing `--output_dir` in analyze.rst** and **missing output column documentation** — doc gaps
4. **Guard against missing ete3 visualization imports** in `visualize_top_trees`
5. **`exit(1)` → `sys.exit(1)`** — minor robustness fix

==================================================================
## Module: Tree_Visualizer
==================================================================
Now I have all the information needed for a thorough review. Here is my analysis:

---

# Tree_Visualizer Module Review

## A. Scientific Correctness

**No issues found.** The module performs tree visualization/annotation rather than novel algorithmic inference. The duplication detection delegates to established helper functions (`find_dup_node`, `annotate_gene_tree`). The midpoint rooting fallback assumption (line 76) is a standard default. The heatmap color scale (0–100 range, line 622–643) is hardcoded without documentation of what units are expected, but this is addressed under Output Descriptions below.

---

## B. Code vs Docs Consistency

### [MEDIUM] [B] [docs/analyze.rst:207–208 vs Phylo_Tracer.py:239]
**`--keep_branch` labeled "Required" in docs but is optional in code.** The CLI parser (line 239) does not set `required=True` for `--keep_branch`. The docs (analyze.rst line 206) list it under "Required Parameters". The README (line 584) correctly lists it as "Optional parameter".

**Fix:** Change analyze.rst to list `--keep_branch` under Optional Parameters.

---

### [MEDIUM] [B] [docs/analyze.rst:207 vs Phylo_Tracer.py:240]
**`--tree_style` labeled "Required" in docs but has a default in code.** The CLI parser sets `default='r'` (line 240), making it optional. The docs (analyze.rst line 207) list it under "Required Parameters".

**Fix:** Move `--tree_style` to Optional Parameters in analyze.rst, noting `default = r`.

---

### [MEDIUM] [B] [docs/analyze.rst:212 vs Phylo_Tracer.py:242, README:594]
**`--gene_expression` in analyze.rst does not exist in code; actual argument is `--heatmap_matrix`.** The CLI parser (line 242) defines `--heatmap_matrix`. The README (line 594) correctly uses `--heatmap_matrix`. The analyze.rst (line 212) still references the old name `--gene_expression`.

**Fix:** Replace `--gene_expression` with `--heatmap_matrix` in analyze.rst and update the description to match ("Gene-associated numeric matrix file").

---

### [MEDIUM] [B] [docs/analyze.rst:218 example vs Phylo_Tracer.py:242]
**Example in analyze.rst uses `--gene_expression expression.csv` instead of `--heatmap_matrix`.** Consistent with the above — the example command is stale.

**Fix:** Update the example line in analyze.rst to `--heatmap_matrix heatmap_matrix.txt`.

---

### [LOW] [B] [docs/analyze.rst:210-211 vs Phylo_Tracer.py:233-237]
**`--gene_categories` description differs slightly.** The analyze.rst description says "One or more two-column files in format `gene_id<TAB>category_label`; each file is one annotation layer" — this is correct but omits the optional header naming rule documented in README (lines 587–592) and the CLI help string (line 237).

**Fix:** Add the optional header naming convention to analyze.rst.

---

### [LOW] [B] [docs/analyze.rst vs Phylo_Tracer.py:244-249]
**GD-related optional parameters missing from analyze.rst.** The CLI parser defines `--gd_support`, `--subclade_support`, `--dup_species_proportion`, `--dup_species_num`, `--deepvar`, and `--output_dir` as optional parameters. The analyze.rst only documents `--visual_gd` and `--input_sps_tree` as GD-related optional params — the five numeric GD knobs are absent.

**Fix:** Add `--gd_support`, `--subclade_support`, `--dup_species_proportion`, `--dup_species_num`, `--deepvar`, and `--output_dir` to the Optional Parameters section of analyze.rst for Tree_Visualizer.

---

## C. Code Logic and Correctness

### [CRITICAL] [C] [Tree_Visualizer.py:1278-1279, view_main]
**`shutil.rmtree(dir_path)` deletes existing output directory without confirmation.** Every time `view_main` runs, if the `tree_visualizer/` directory already exists in the CWD, it is silently deleted and recreated. This destroys all previous output files, including any PDFs the user may have manually placed or renamed in that directory.

**Fix:** Either skip deletion and overwrite individual files, or warn/prompt before deletion. At minimum, document this destructive behavior.

---

### [MEDIUM] [C] [Tree_Visualizer.py:804, tips_mark]
**`node.name.split("_")[0]` assumes species voucher is always the first `_`-delimited token.** If gene IDs contain underscores in unexpected positions (common in plant genomics, e.g., `Ath_Chr1_gene001`), this will extract the wrong voucher prefix, leading to silent KeyError or incorrect species assignment.

**Fix:** This is a convention enforced by `gene_id_transfer`, but a guard or explicit documentation of the naming convention would help. At minimum, wrap in a try/except with a meaningful error message.

---

### [MEDIUM] [C] [Tree_Visualizer.py:622-643, get_color]
**Heatmap color bins are hardcoded to 0–100 range with no normalization.** Values outside [0, 100] all map to `#cc3300` (the `else` branch at line 643). Negative values also fall through to the `else` — they get the highest-intensity color instead of being flagged or handled. This is scientifically misleading for expression data or other metrics that may span arbitrary ranges.

**Fix:** Either (a) add min-max normalization to [0, 100] before coloring, or (b) document that input values must be pre-scaled to 0–100. At minimum, handle negative values explicitly (e.g., return "white" or a distinct color).

---

### [LOW] [C] [Tree_Visualizer.py:326-330, generate_string]
**Off-by-one in suffix generation for index ≥ 52.** When `index >= 52` (len(letters)), the formula `(index - len(letters)) // 52` and `% 52` produces two-letter suffixes, but the first compound suffix at index=52 gives `first_letter = letters[0]` and `second_letter = letters[0]`, i.e., `@AA`. Index 53 gives `@AB`, etc. This works but the divisor 52 means the rollover happens at 52×52+52 = 2756 categories, which is fine in practice. However, `letters` has length 52 (26+26), so `letters[(index - 52) // 52]` will cause an IndexError when `index >= 52 + 52*52 = 2756`.

**Fix:** Not urgent (2756 categories is unrealistic), but a modular approach would be more robust.

---

### [LOW] [C] [Tree_Visualizer.py:1141-1143, mark_gene_to_sptree]
**Potential KeyError in `else` branch of `generate_face_mark` inside `mark_gene_to_sptree`.** When `(node, column, "aligned")` IS in `faces_added`, the else branch still tries to access `color_dict[species]` without checking whether `species` is in the dict. If the species is missing from the category dict, this raises a `KeyError`.

**Fix:** Add a `if species in color_dict:` guard in the else branch.

---

### [LOW] [C] [Tree_Visualizer.py:1390, mark_gene_to_sptree_main]
**`sp2.convert_to_ultrametric()` called without `tree_length` argument.** ETE3's `convert_to_ultrametric()` with no argument defaults to preserving the original root-to-tip max distance, which may produce unexpected scaling for species trees with very short or long branches.

**Fix:** Consider documenting or parameterizing this behavior. Not a bug per se, but worth noting.

---

## D. Output Descriptions

### [MEDIUM] [D] [docs/analyze.rst:197-218, README:575-603]
**No output files are documented for Tree_Visualizer.** Neither analyze.rst nor README describe what files Tree_Visualizer produces. The code generates:
1. `tree_visualizer/<tre_ID>.pdf` — per-gene-tree annotated PDFs (line 824)
2. `genefamily_map2_sptree.pdf` — species tree with family mappings (line 1391)

**Fix:** Add an "Output" section documenting these files and their content.

---

### [LOW] [D] [README:578-579]
**README states "Tree figures are exported as vector PDFs" but does not specify the output directory name.** The code places them in `tree_visualizer/` (line 1272).

**Fix:** Add "Output directory: `tree_visualizer/`" to the README description.

---

### [MEDIUM] [D] [README:594, docs/analyze.rst:212]
**Heatmap matrix expected value range (0–100) is not documented.** The hardcoded color bins (lines 622–643) assume values in [0, 100]. Users with TPM/FPKM expression data (often 0–100,000+) will get misleading all-red heatmaps.

**Fix:** Document that values should be pre-scaled to [0, 100], or implement automatic normalization.

---

## Summary Table

| Severity   | Count |
|------------|-------|
| CRITICAL   | 1     |
| MEDIUM     | 7     |
| LOW        | 5     |
| SUGGESTION | 0     |
| **Total**  | **13** |

### Top 3 priorities:
1. **CRITICAL:** `shutil.rmtree` silently deletes the output directory on every run — risk of data loss.
2. **MEDIUM:** `--gene_expression` in analyze.rst is stale; actual CLI argument is `--heatmap_matrix`.
3. **MEDIUM:** GD-related optional parameters (`--gd_support`, etc.) and output file descriptions are missing from analyze.rst.

