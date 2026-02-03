import argparse
import os
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

from ete3 import Tree


def parse_tree_line(line: str) -> Tree:
    s = line.strip()
    try:
        return Tree(s, format=0)
    except Exception:
        return Tree(s, format=1)


def split_species_gene(name: str) -> tuple:
    if '_' in name:
        parts = name.split('_')
        species = parts[0]
        gene = '_'.join(parts[1:]) if len(parts) > 1 else parts[0]
        return species, gene
    if '|' in name:
        p = name.split('|')
        return p[0], p[1]
    m = re.match(r"([^\.\s]+)[\.|\s]+(.+)", name)
    if m:
        return m.group(1), m.group(2)
    return name, name


def build_synteny_index(dir_path: str) -> dict:
    files = os.listdir(dir_path)
    idx = defaultdict(list)
    for f in files:
        lf = f.lower()
        tokens = re.split(r"[^a-z0-9]+", lf)
        uniq = [t for t in tokens if t]
        for i in range(len(uniq)):
            for j in range(i + 1, len(uniq)):
                k = frozenset({uniq[i], uniq[j]})
                idx[k].append(os.path.join(dir_path, f))
    return idx


def choose_file(candidates: list, sp_a: str, sp_b: str) -> str:
    if not candidates:
        return None
    prefer = [p for p in candidates if re.search(r"syn|anchor|pair", os.path.basename(p), re.I)]
    if prefer:
        return max(prefer, key=lambda x: os.path.getsize(x))
    return max(candidates, key=lambda x: os.path.getsize(x))


@lru_cache(maxsize=128)
def load_synteny_pairs(fp: str) -> set:
    s = set()
    if fp is None:
        return s
    with open(fp, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            a = re.split(r"[\t,\s]+", line.strip())
            if len(a) < 2:
                continue
            g1 = a[0]
            g2 = a[1]
            s.add((g1, g2))
            s.add((g2, g1))
    return s


def synteny_rate_for_tree(t: Tree, syn_idx: dict) -> float:
    leaves = t.get_leaf_names()
    sg = [split_species_gene(x) for x in leaves]
    by_sp = defaultdict(list)
    for sp, g in sg:
        by_sp[sp].append(g)
    species = list(by_sp.keys())
    total = 0
    hit = 0
    for i in range(len(species)):
        for j in range(i + 1, len(species)):
            sp_a = species[i]
            sp_b = species[j]
            key = frozenset({sp_a.lower(), sp_b.lower()})
            fps = choose_file(syn_idx.get(key, []), sp_a, sp_b)
            syn_pairs = load_synteny_pairs(fps)
            a_list = by_sp[sp_a]
            b_list = by_sp[sp_b]
            for ga in a_list:
                for gb in b_list:
                    total += 1
                    if (ga, gb) in syn_pairs:
                        hit += 1
    return (hit / total) if total > 0 else 0.0


def compute_all_rates(trees_path: str, syn_dir: str, threads: int) -> tuple:
    syn_idx = build_synteny_index(syn_dir)
    rates = []
    with open(trees_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    def worker(line):
        t = parse_tree_line(line)
        return synteny_rate_for_tree(t, syn_idx)
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [ex.submit(worker, line) for line in lines]
        for fut in as_completed(futs):
            rates.append(fut.result())
    avg = sum(rates) / len(rates) if rates else 0.0
    return rates, avg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trees', required=True)
    ap.add_argument('--syn_dir', required=True)
    ap.add_argument('--threads', type=int, default=os.cpu_count())
    args = ap.parse_args()
    rates, avg = compute_all_rates(args.trees, args.syn_dir, args.threads)
    for i, r in enumerate(rates, start=1):
        print(f"tree_{i}\tsynteny_rate\t{r:.6f}")
    print(f"global_avg\t{avg:.6f}")


if __name__ == '__main__':
    main()
