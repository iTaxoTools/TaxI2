from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from itaxotools.taxi2.trees import Tree, Trees

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


class ReadTest(NamedTuple):
    fixture: Callable[[], Trees]
    input: str

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> Trees:
        return self.fixture()

    def validate(self):
        generated = Trees.fromPath(self.input_path)
        assert generated == self.fixed


def trees_single() -> Trees:
    return Trees(
        [
            Tree("(A,B,(C,D));"),
        ]
    )


def trees_multiple() -> Trees:
    return Trees(
        [
            Tree("(A,B,(C,D));"),
            Tree("(A,B,(D,C));"),
            Tree("(B,A,(C,D));"),
        ]
    )


def trees_complex() -> Trees:
    return Trees(
        [
            Tree("[&U](A:1,B:2,(C,D));"),
        ]
    )


read_tests = [
    ReadTest(trees_single, "single.tree"),
    ReadTest(trees_multiple, "multiple.tree"),
    ReadTest(trees_complex, "complex.tree"),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_trees(test: ReadTest) -> None:
    test.validate()


def test_read_bad_tree() -> None:
    test = ReadTest(None, "bad.tree")
    with pytest.raises(ValueError):
        test.validate()


def test_newick_semicolon_add():
    tree = Tree.from_newick_string("(A,B)")
    string = tree.get_newick_string(lengths=False, semicolon=True, comments=False)
    assert string == "(A,B);"


def test_newick_semicolon_remove():
    tree = Tree.from_newick_string("(A,B);")
    string = tree.get_newick_string(lengths=False, semicolon=False, comments=False)
    assert string == "(A,B)"


def test_newick_preserve_comments():
    tree = Tree.from_newick_string("[&U](A,B[123])")
    string = tree.get_newick_string(lengths=False, semicolon=False, comments=True)
    assert string == "[&U](A,B[123])"


def test_newick_remove_comments():
    tree = Tree.from_newick_string("[&U](A,B[123])")
    string = tree.get_newick_string(lengths=False, semicolon=False, comments=False)
    assert string == "(A,B)"


def test_newick_preserve_lengths():
    tree = Tree.from_newick_string("(A:1,B:2):3")
    string = tree.get_newick_string(lengths=True, semicolon=False, comments=False)
    assert string == "(A:1,B:2):3"


def test_newick_remove_lengths():
    tree = Tree.from_newick_string("(A:1,B:2):3")
    string = tree.get_newick_string(lengths=False, semicolon=False, comments=False)
    assert string == "(A,B)"


def test_newick_names_rooted():
    tree = Tree.from_newick_string("((A,B),(C,D))")
    string = tree.get_newick_string(lengths=False, semicolon=False, comments=False)
    assert string == "((A,B),(C,D))"
    names = tree.get_node_names()
    assert names == ["A", "B", "C", "D"]


def test_newick_names_unrooted():
    tree = Tree.from_newick_string("(A,B,(C,D))")
    string = tree.get_newick_string(lengths=False, semicolon=False, comments=False)
    assert string == "(A,B,(C,D))"
    names = tree.get_node_names()
    assert names == ["A", "B", "C", "D"]


def test_newick_scientific():
    Tree.from_newick_string("(A:1E1,B:0.0):-1e42")


def test_newick_underscores():
    Tree.from_newick_string("(a_sample,b_sample_42)")


def test_newick_single():
    Tree.from_newick_string("a")


def test_newick_extra_parentheses():
    Tree.from_newick_string("((A,B))")


def test_newick_extra_parentheses_with_distances():
    Tree.from_newick_string("(((A,B),C):1)")


def test_newick_strip():
    Tree.from_newick_string("\n(A,B);\r\n")


def test_newick_empty():
    with pytest.raises(ValueError):
        Tree.from_newick_string("")


def test_newick_bad_parentesis():
    with pytest.raises(ValueError):
        Tree.from_newick_string("(A,B")


def test_newick_bad_poly_root():
    with pytest.raises(ValueError):
        Tree.from_newick_string("(A,B,C,D)")


def test_newick_bad_characters():
    with pytest.raises(ValueError):
        Tree.from_newick_string("(A!,B)")
