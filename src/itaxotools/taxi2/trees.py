from __future__ import annotations

import re
from pathlib import Path
from typing import Literal, NamedTuple

from .handlers import FileHandler, ReadHandle, WriteHandle


class Tree(NamedTuple):
    """Uses Newick strings as a base"""

    newick: str

    @classmethod
    def from_newick_string(cls, newick: str) -> Tree:
        """Validate the string before construction"""
        valid, _ = cls._process_newick_string(newick)
        if not valid:
            raise ValueError("Cannot parse Newick string!")
        return cls(newick)

    @classmethod
    def _process_newick_string(cls, newick: str) -> tuple[bool, list[str]]:
        newick = cls._format_newick_string(newick, False, False, False)
        names = set()

        # recursively collapse leaves into branches
        pattern = re.compile(r"\(([\w\.\-]+?),([\w\.\-]+?)\)")
        while True:
            hit = pattern.search(newick)
            if hit is None:
                break
            names.add(hit.group(1))
            names.add(hit.group(2))
            newick = newick.replace(hit.group(0), hit.group(1))

        # remove remaining parentheses around root
        while newick.startswith("(") and newick.endswith(")"):
            newick = newick[1:-1]

        # final unrooted tree may contain at most three nodes
        pattern_unrooted = re.compile(r"^([\w\.\-]+?),([\w\.\-]+?),([\w\.\-]+?)$")
        hit = pattern_unrooted.search(newick)
        if hit is not None:
            names.add(hit.group(1))
            names.add(hit.group(2))
            names.add(hit.group(3))

        names = list(sorted(names))

        pattern_single = re.compile(r"^[\w\.\-]+?$")

        # check if the collapsed tree is valid
        if re.fullmatch(pattern_single, newick) or re.fullmatch(
            pattern_unrooted, newick
        ):
            return True, names

        return False, names

    @staticmethod
    def _format_newick_string(
        newick: str,
        lengths: bool,
        semicolon: bool,
        comments: bool,
    ) -> str:
        newick = newick.strip()

        if semicolon and not newick.endswith(";"):
            newick += ";"
        if not semicolon and newick.endswith(";"):
            newick = newick[:-1]

        if not comments:
            newick = re.sub(r"\[[^\]]*\]", "", newick)

        if not lengths:
            newick = re.sub(r":-?\d*\.?\d+(-?[Ee]\d+)?", "", newick)

        return newick

    def get_newick_string(
        self,
        lengths=True,
        semicolon=False,
        comments=False,
    ) -> str:
        """Get a formatted newick string"""
        return self._format_newick_string(self.newick, lengths, semicolon, comments)

    def get_node_names(self) -> list[str]:
        """Get the names of all nodes"""
        _, names = self._process_newick_string(self.newick)
        return names


class Trees(list[Tree]):
    @classmethod
    def fromPath(cls, path: Path, *args, **kwargs) -> list[Tree]:
        return cls(tree for tree in NewickTreeHandler(path, *args, **kwargs))


class NewickTreeHandler(FileHandler[Tree]):
    """Strict parser for Newick files, one tree per line"""

    def _open(self, path: Path, mode: Literal["r", "w"] = "r", *args, **kwargs):
        self.filter = filter
        super()._open(path, mode, *args, **kwargs)

    def _iter_write(self) -> WriteHandle[Tree]:
        raise NotImplementedError()

    def _iter_read(self, *args, **kwargs) -> ReadHandle[Tree]:
        with open(self.path) as file:
            yield self
            for line in file:
                line = line.strip()
                yield Tree.from_newick_string(line)
