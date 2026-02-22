# PhyloTracer Benchmarks

Standalone benchmarking script for measuring performance of key PhyloTracer operations.

## Requirements

- Python 3.8+
- `ete3`
- `numpy`

PhyloTracer itself is optional. When installed, the benchmark uses the actual
PhyloTracer functions. Otherwise it falls back to equivalent built-in
implementations and prints a notice.

## Usage

Run with default settings (sizes 50, 100, 500, 1000 taxa; 3 repetitions):

```bash
python benchmark.py
```

Specify custom tree sizes:

```bash
python benchmark.py --sizes 50 100 500
```

Increase repetitions for more stable averages:

```bash
python benchmark.py --sizes 50 100 500 1000 --repeat 5
```

## Benchmarked Operations

| Operation                  | Description                                                       |
|----------------------------|-------------------------------------------------------------------|
| Tree read/parse            | Parse a Newick string into an ete3 Tree object                    |
| Tree rooting (midpoint)    | Root a tree using midpoint outgroup selection                     |
| Species set extraction     | Extract unique species labels from every node in the tree         |
| Gene duplication detection | Annotate a gene tree against a species tree and detect GD nodes   |
| Tree topology comparison   | Compute Robinson-Foulds distance between two trees                |

## Output

The script prints a formatted table with four columns:

- **Operation** -- name of the benchmarked function
- **Input Size** -- number of taxa in the synthetic tree
- **Time (s)** -- average wall-clock time across repetitions
- **Peak Mem (MB)** -- peak memory allocated during the operation (via `tracemalloc`)
