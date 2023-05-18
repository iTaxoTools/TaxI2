from __future__ import annotations

from pathlib import Path
from re import fullmatch
from typing import Callable

from itaxotools.spart_parser.main import Spart as SpartParserSpart
from itaxotools.spart_parser.main import is_path_xml

from .encoding import sanitize
from .handlers import FileHandler
from .partitions import PartitionHandler
from .file_types import FileFormat, FileInfo


FormatIdentifier = Callable[[Path], bool]
FormatIdentifierDecorator = Callable[[FormatIdentifier], FormatIdentifier]

InfoGetter = Callable[[Path, FileFormat], bool]
InfoGetterDecorator = Callable[[InfoGetter], InfoGetter]

format_identifiers: Dict[FileFormat, CallableTest] = {}
info_getters: Dict[FileFormat, CallableTest] = {}


def identifier(format: FileFormat) -> FormatIdentifierDecorator:
    def decorator(func: FormatIdentifier) -> FormatIdentifier:
        format_identifiers[format] = func
        return func
    return decorator


def info_getter(format: FileFormat) -> InfoGetterDecorator:
    def decorator(func: InfoGetter) -> InfoGetter:
        info_getters[format] = func
        return func
    return decorator


def identify_format(path: Path):
    for format in format_identifiers:
        if format_identifiers[format](path):
            return format
    return FileFormat.Unknown


def get_info(path: Path, format: FileFormat = None):
    if format is None:
        format = identify_format(path)
    if format not in info_getters:
        format = FileFormat.Unknown
    return info_getters[format](path, format)


@identifier(FileFormat.Fasta)
def is_fasta(path: Path) -> bool:
    with path.open() as file:
        return bool(file.read(1) == '>')


@identifier(FileFormat.Tabfile)
def is_tabfile(path: Path) -> bool:
    with path.open() as file:
        line = file.readline()
        return bool(fullmatch(r'([^\t]+\t)+[^\t]+', line))


@identifier(FileFormat.Spart)
def is_spart(path: Path) -> bool:
    try:
        SpartParserSpart.fromPath(path)
    except Exception:
        return False
    return True


@info_getter(FileFormat.Fasta)
def get_fasta_info(path: Path, format: FileFormat) -> bool:
    subset_separator = PartitionHandler.Fasta.guess_subset_separator(path)
    has_subsets = PartitionHandler.Fasta.has_subsets(path, subset_separator)
    return FileInfo.Fasta(
        path = path,
        format = format,
        size = path.stat().st_size,
        has_subsets = has_subsets,
        subset_separator = subset_separator,
    )


@info_getter(FileFormat.Tabfile)
def get_tabfile_info(path: Path, format: FileFormat) -> bool:
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

    return FileInfo.Tabfile(
        path = path,
        format = format,
        size = path.stat().st_size,
        headers = headers,
        header_individuals = header_individuals,
        header_sequences = header_sequences,
        header_organism = header_organism,
        header_species = header_species,
        header_genus = header_genus,
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
        path = path,
        format = format,
        size = path.stat().st_size,
        spartitions = spartitions,
        is_matricial = not is_xml,
        is_xml = is_xml,
    )

@info_getter(FileFormat.Unknown)
def get_general_info(path: Path, format: FileFormat) -> bool:
    return FileInfo(
        path = path,
        format = format,
        size = path.stat().st_size,
    )
