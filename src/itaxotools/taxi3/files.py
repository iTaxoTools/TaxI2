from __future__ import annotations

from typing import Callable
from enum import Enum, auto
from pathlib import Path
from re import fullmatch
from dataclasses import dataclass, asdict

from .types import Type
from .handlers import FileHandler


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
    header_id: str | None
    header_seq: str | None
    header_species: str | None
    header_genus: str | None

    @classmethod
    def get(cls, path: Path):
        headers = FileHandler.Tabfile(path, has_headers=True).headers
        header_id = 'id' if 'id' in headers else None
        header_seq = 'seq' if 'seq' in headers else None
        header_species = 'organism' if 'organism' in headers else None
        header_genus = 'organism' if 'organism' in headers else None

        return cls(
            headers = list(headers),
            header_id = header_id,
            header_seq = header_seq,
            header_species = header_species,
            header_genus = header_genus,
            **asdict(FileInfo.get(path))
        )


class FileFormat(Enum):
    Fasta = 'Fasta', Identifier.isFasta, FileInfo
    Tabfile = 'Tabfile', Identifier.isTabFile, FileInfo.Tabfile
    Excel = 'Excel', None, None
    Unknown = 'Unknown', None, None

    def __init__(self, label, identifier, info):
        self.label = label
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
