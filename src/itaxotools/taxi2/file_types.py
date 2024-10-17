from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from .types import Type


class FileFormat(Enum):
    Ali = "Ali", ".ali"
    Fasta = "Fasta", ".fas"
    FastQ = "FastQ", ".fq"
    Tabfile = "Tabfile", ".tsv"
    Spart = "Spart", ".spart"
    Excel = "Excel", ".xlsx"
    Newick = "Newick", ".tree"
    Unknown = "Unknown", None

    def __init__(self, label, extension):
        self.label = label
        self.extension = extension

    def __repr__(self):
        return f"<{type(self).__name__}.{self._name_}>"


@dataclass
class FileInfo(Type):
    path: Path
    format: FileFormat
    size: int


@dataclass
class Tabular(FileInfo):
    headers: list[str]
    header_individuals: str | None
    header_sequences: str | None
    header_organism: str | None
    header_species: str | None
    header_genus: str | None


@dataclass
class Tabfile(Tabular, FileInfo):
    pass


@dataclass
class Excel(Tabular, FileInfo):
    pass


@dataclass
class Fasta(FileInfo):
    has_subsets: bool
    subset_separator: str


@dataclass
class Spart(FileInfo):
    spartitions: list[str]
    is_matricial: bool
    is_xml: bool


@dataclass
class Newick(FileInfo):
    count: int
    names: set[str]
