from __future__ import annotations

from typing import Callable
from enum import Enum, auto
from pathlib import Path
from re import fullmatch
from dataclasses import dataclass, asdict

from .types import Type
from .handlers import FileHandler
from .partitions import PartitionHandler
from .encoding import sanitize


class Identifier:
    def isFasta(path: Path) -> bool:
        with path.open() as file:
            return bool(file.read(1) == '>')

    def isTabFile(path: Path) -> bool:
        with path.open() as file:
            line = file.readline()
            return bool(fullmatch(r'([^\t]+\t)+[^\t]+', line))


@dataclass
class FileInfo(Type):
    path: Path
    format: FileFormat
    size: int

    @classmethod
    def from_path(cls, path: Path) -> FileInfo:
        format = FileFormat.identify(path)
        return format.info.get(path)

    @classmethod
    def get(cls, path: Path):
        return cls(
            path = path,
            format = FileFormat.identify(path),
            size = path.stat().st_size,
        )


@dataclass
class Tabfile(FileInfo):
    headers: list[str]
    header_individuals: str | None
    header_sequences: str | None
    header_organism: str | None
    header_species: str | None
    header_genus: str | None

    @classmethod
    def get(cls, path: Path):
        headers = FileHandler.Tabfile(path, has_headers=True).headers
        headers = [sanitize(header) for header in headers]

        return cls(
            headers = list(headers),
            header_individuals = 'seqid' if 'seqid' in headers else None,
            header_sequences = 'sequence' if 'sequence' in headers else None,
            header_organism = 'organism' if 'organism' in headers else None,
            header_species = 'species' if 'species' in headers else None,
            header_genus = 'genus' if 'genus' in headers else None,
            **asdict(FileInfo.get(path))
        )


@dataclass
class Fasta(FileInfo):
    has_subsets: bool

    @classmethod
    def get(cls, path: Path):
        has_subsets = PartitionHandler.Fasta.has_subsets(path,)
        return cls(
            has_subsets = has_subsets,
            **asdict(FileInfo.get(path))
        )


class FileFormat(Enum):
    Fasta = 'Fasta', '.fas', Identifier.isFasta, FileInfo.Fasta
    Tabfile = 'Tabfile', '.tsv', Identifier.isTabFile, FileInfo.Tabfile
    Excel = 'Excel', '.xlsx', None, None
    Unknown = 'Unknown', None, None, None

    def __init__(self, label, extension, identifier, info):
        self.label = label
        self.extension = extension
        self.identifier = identifier or (lambda _: False)
        self.info = info or FileInfo

    def __repr__(self):
        return f'<{type(self).__name__}.{self._name_}>'

    @classmethod
    def identify(cls, path: Path) -> FileFormat:
        for format in cls:
            if format.identifier(path):
                return format
        return cls.Unknown
