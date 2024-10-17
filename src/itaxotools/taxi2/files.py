from __future__ import annotations

from pathlib import Path
from re import fullmatch
from typing import Callable

from itaxotools.common.utility import DecoratorDict
from itaxotools.spart_parser.main import Spart as SpartParserSpart
from itaxotools.spart_parser.main import is_path_xml

from .encoding import sanitize
from .file_types import FileFormat, FileInfo
from .handlers import FileHandler
from .partitions import PartitionHandler
from .trees import Trees

FormatIdentifier = Callable[[Path], bool]
InfoGetter = Callable[[Path, FileFormat], bool]

identifier = DecoratorDict[FileFormat, FormatIdentifier]()
info_getter = DecoratorDict[FileFormat, InfoGetter]()


def identify_format(path: Path):
    for format in identifier:
        if identifier[format](path):
            return format
    return FileFormat.Unknown


def get_info(path: Path, format: FileFormat = None):
    if format is None:
        format = identify_format(path)
    if format not in info_getter:
        format = FileFormat.Unknown
    return info_getter[format](path, format)


@identifier(FileFormat.Ali)
def is_ali(path: Path) -> bool:
    with path.open() as file:
        infos = False
        for line in file:
            if not line.strip():
                continue
            if line.startswith("#"):
                infos = True
                continue
            if line.startswith(">"):
                return bool(infos)
    return False


@identifier(FileFormat.Fasta)
def is_fasta(path: Path) -> bool:
    with path.open() as file:
        for line in file:
            if not line.strip():
                continue
            if line.startswith(";"):
                continue
            if line.startswith(">"):
                return True
    return False


@identifier(FileFormat.FastQ)
def is_fastq(path: Path) -> bool:
    with open(path, "r") as file:
        id = False
        for line in file:
            if not line.strip():
                continue
            if line.startswith("@"):
                id = True
            if line.startswith("+"):
                return bool(id)
    return False


@identifier(FileFormat.Tabfile)
def is_tabfile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r"([^\t]+\t)+[^\t]+", line))


@identifier(FileFormat.Spart)
def is_spart(path: Path) -> bool:
    try:
        SpartParserSpart.fromPath(path)
    except Exception:
        return False
    return True


@identifier(FileFormat.Newick)
def is_newick(path: Path) -> bool:
    try:
        trees = Trees.fromPath(path)
    except Exception:
        return False
    if not len(trees):
        return False
    return True


@info_getter(FileFormat.Fasta)
def get_fasta_info(path: Path, format: FileFormat) -> bool:
    subset_separator = PartitionHandler.Fasta.guess_subset_separator(path)
    has_subsets = PartitionHandler.Fasta.has_subsets(path, subset_separator)
    return FileInfo.Fasta(
        path=path,
        format=format,
        size=path.stat().st_size,
        has_subsets=has_subsets,
        subset_separator=subset_separator,
    )


@info_getter(FileFormat.Tabfile)
def get_tabfile_info(path: Path, format: FileFormat) -> bool:
    headers = FileHandler.Tabfile(path, has_headers=True).headers
    headers = [sanitize(header) for header in headers]

    header_individuals = "seqid" if "seqid" in headers else None
    header_sequences = "sequence" if "sequence" in headers else None
    header_organism = "organism" if "organism" in headers else None
    header_species = "species" if "species" in headers else None
    header_genus = "genus" if "genus" in headers else None

    species_is_binomen = False
    if "species" in headers:
        index = headers.index("species")
        with FileHandler.Tabfile(path, columns=[index], has_headers=True) as file:
            first = file.read()
            if first is not None:
                parts = first[0].split(" ")
                species_is_binomen = bool(len(parts) > 1)

    if species_is_binomen:
        if "organism" not in headers and "genus" not in headers:
            header_organism = "species"
            header_species = None
            header_genus = None

    return FileInfo.Tabfile(
        path=path,
        format=format,
        size=path.stat().st_size,
        headers=headers,
        header_individuals=header_individuals,
        header_sequences=header_sequences,
        header_organism=header_organism,
        header_species=header_species,
        header_genus=header_genus,
    )


@info_getter(FileFormat.Spart)
def get_spart_info(path: Path, format: FileFormat) -> bool:
    is_xml = is_path_xml(path)

    if is_xml:
        spart = SpartParserSpart.fromXML(path)
    else:
        spart = SpartParserSpart.fromMatricial(path)

    spartitions = spart.getSpartitions()

    return FileInfo.Spart(
        path=path,
        format=format,
        size=path.stat().st_size,
        spartitions=spartitions,
        is_matricial=not is_xml,
        is_xml=is_xml,
    )


@info_getter(FileFormat.Newick)
def get_newick_info(path: Path, format: FileFormat) -> bool:
    trees = Trees.fromPath(path)
    return FileInfo.Newick(
        path=path,
        format=format,
        size=path.stat().st_size,
        count=len(trees),
        names=set(name for tree in trees for name in tree.get_node_names()),
    )


@info_getter(FileFormat.Unknown)
def get_general_info(path: Path, format: FileFormat) -> bool:
    return FileInfo(
        path=path,
        format=format,
        size=path.stat().st_size,
    )
