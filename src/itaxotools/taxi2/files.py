from __future__ import annotations

from typing import Callable
from enum import Enum, auto
from pathlib import Path
from re import fullmatch
from dataclasses import dataclass, asdict

from itaxotools.spart_parser.main import Spart as SpartParserSpart, is_path_xml

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

    def isSpart(path: Path) -> bool:
        try:
            SpartParserSpart.fromPath(path)
        except Exception:
            return False
        return True


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

        header_individuals = 'seqid' if 'seqid' in headers else None
        header_sequences = 'sequence' if 'sequence' in headers else None
        header_organism = 'organism' if 'organism' in headers else None
        header_species = 'species' if 'species' in headers else None
        header_genus = 'genus' if 'genus' in headers else None

        species_is_binomen = False
        if 'species' in headers:
            index = headers.index('species')
            with FileHandler.Tabfile(path, columns=[index], has_headers=True) as file:
                first = file.read()
                if first is not None:
                    parts = first[0].split(' ')
                    species_is_binomen = bool(len(parts) > 1)

        if species_is_binomen:
            if 'organism' not in headers and 'genus' not in headers:
                header_organism = 'species'
                header_species = None
                header_genus = None

        return cls(
            headers = headers,
            header_individuals = header_individuals,
            header_sequences = header_sequences,
            header_organism = header_organism,
            header_species = header_species,
            header_genus = header_genus,
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


@dataclass
class Spart(FileInfo):
    spartitions: list[str]
    is_matricial: bool
    is_xml: bool

    @classmethod
    def get(cls, path: Path):
        is_xml = is_path_xml(path)

        if is_xml:
            spart = SpartParserSpart.fromXML(path)
        else:
            spart = SpartParserSpart.fromMatricial(path)

        spartitions = spart.getSpartitions()

        return cls(
            spartitions = spartitions,
            is_matricial = not is_xml,
            is_xml = is_xml,
            **asdict(FileInfo.get(path))
        )


class FileFormat(Enum):
    Fasta = 'Fasta', '.fas', Identifier.isFasta, FileInfo.Fasta
    Tabfile = 'Tabfile', '.tsv', Identifier.isTabFile, FileInfo.Tabfile
    Spart = 'Spart', '.spart', Identifier.isSpart, FileInfo.Spart
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
