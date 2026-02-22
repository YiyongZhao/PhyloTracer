"""
Shared pytest fixtures for the PhyloTracer test suite.
"""

import os
import sys
import pytest

# Ensure the project root is on the import path.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from ete3 import Tree, PhyloTree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False


@pytest.fixture
def simple_tree():
    """A simple rooted binary Newick tree with species_gene leaf names."""
    if not HAS_ETE3:
        pytest.skip("ete3 is not installed")
    return Tree("((A_1:0.1,B_1:0.2):0.3,(C_1:0.4,D_1:0.5):0.6):0.0;")


@pytest.fixture
def unrooted_tree():
    """An unrooted tree (trifurcating root) with species_gene leaf names."""
    if not HAS_ETE3:
        pytest.skip("ete3 is not installed")
    return Tree("(A_1:0.1,B_1:0.2,(C_1:0.3,D_1:0.4):0.5);")


@pytest.fixture
def species_tree():
    """A simple rooted species tree with plain species names."""
    if not HAS_ETE3:
        pytest.skip("ete3 is not installed")
    return Tree("((A:1,B:1):1,(C:1,D:1):1);")


@pytest.fixture
def phylo_species_tree():
    """A PhyloTree version of the species tree."""
    if not HAS_ETE3:
        pytest.skip("ete3 is not installed")
    return PhyloTree("((A:1,B:1):1,(C:1,D:1):1);")


@pytest.fixture
def gene_tree_with_duplication():
    """A gene tree where species A has two copies (duplication signal)."""
    if not HAS_ETE3:
        pytest.skip("ete3 is not installed")
    return Tree("((A_1:0.1,B_1:0.2):0.3,(A_2:0.1,C_1:0.4):0.6);")


@pytest.fixture
def mapping_file(tmp_path):
    """Create a temporary two-column mapping file."""
    content = "gene1\tspeciesA\ngene2\tspeciesB\ngene3\tspeciesA\n"
    fpath = tmp_path / "mapping.tsv"
    fpath.write_text(content)
    return str(fpath)
