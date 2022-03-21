#!/usr/bin/env python3
from __future__ import annotations

from typing import List, Set, Any, Optional, Iterator, Type
from abc import ABC, abstractmethod
from pathlib import Path
import os
import re
import itertools
from enum import Enum, auto

import pandas as pd
import openpyxl

from itaxotools.DNAconvert.library.fasta import Fastafile
from itaxotools.DNAconvert.library.genbank import GenbankFile
from itaxotools.DNAconvert.library.record import Record


class InvalidPath(Exception):
    """
    Raised by ValidFilePath
    """

    pass


class ValidFilePath:
    """
    PathLike that contains a path to an existing file
    """

    def __init__(self, *args: Any, **kwargs: Any):
        """
        raise InvalidPath if the arguments are a path that doesn't exists
        or doesn't point to a file
        """
        self._path = Path(*args, **kwargs)
        self.exists = self._path.exists
        self.is_file = self._path.is_file
        self.stat = self._path.stat
        if not self.exists() or not self.is_file():
            raise InvalidPath

    def __fspath__(self) -> str:
        return str(self._path)


os.PathLike.register(ValidFilePath)


def _dataframe_append(dataframe: pd.DataFrame, file: Path) -> None:
    try:
        header = file.stat().st_size == 0
    except FileNotFoundError:
        header = True
    dataframe.to_csv(file, header=header, sep="\t")


class _Header:

    RENAME_DICT = {
        "organism": "species",
        "specimen_identifier": "specimen_voucher",
        "specimen identifier": "specimen_voucher",
        "specimen voucher": "specimen_voucher",
    }

    def __init__(self, header: List[str]):
        self.names = list(map(_Header.rename, map(_Header.to_ascii, header)))

    @staticmethod
    def to_ascii(name: str) -> str:
        regex = r"([a-zA-Z0-9_-]+)"
        return "".join(re.findall(regex, name))

    @staticmethod
    def rename(name: str) -> str:
        if "sequence" in name:
            return "sequence"
        else:
            try:
                return _Header.RENAME_DICT[name]
            except KeyError:
                return name


class DataType(ABC):
    """
    Class of TaxI input and output data
    """

    @classmethod
    @abstractmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> DataType:
        """
        Create data from file at `path` using `protocol` to read it
        """
        pass

    @abstractmethod
    def get_dataframe(self) -> pd.DataFrame:
        pass


class DuplicatedValues(Exception):
    """
    Abstract base class for exception when setting a column with duplicated values as index.
    """

    pass


class DuplicatedSeqids(DuplicatedValues):
    """
    "seqid" column cannot be an index because of duplicated values
    """

    pass


class DuplicatedSpecies(DuplicatedValues):
    """
    "species" column cannot be an index because of duplicated values
    """

    pass


class SequenceData(DataType):
    """
    Represents a collection of sequences.

    Only loads data on demand.
    """

    def __init__(self, path: Optional[ValidFilePath], protocol: FileReader):
        self.path = path
        self.protocol = protocol
        self.normalize = False
        self.dataframe: Optional[pd.DataFrame] = None

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> SequenceData:
        return cls(path, protocol)

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> SequenceData:
        self = cls(None, TabfileReader())
        assert set(dataframe.columns) == {"sequence"}
        assert dataframe.index.name == "seqid"
        self.dataframe = dataframe
        return self

    @staticmethod
    def _normalize_sequences(table: pd.DataFrame) -> None:
        assert "sequence" in table.columns
        table["sequence"] = (
            table["sequence"]
            .str.upper()
            .str.replace("?", "N", regex=False)
            .str.replace("-", "", regex=False)
        )

    def normalize_sequences(self) -> None:
        self.normalize = True
        if self.dataframe is not None:
            self._normalize_sequences(self.dataframe)

    def get_dataframe(self) -> pd.DataFrame:
        """
        Loads the data, if necessary, and returns it.

        Returns pandas DataFrame with index "seqid" and column "sequence"
        """
        if self.dataframe is None:
            assert self.path is not None
            self.dataframe = self.protocol.read(
                self.path, columns=["seqid", "sequence"]
            )
            try:
                self.dataframe.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSeqids()
            if self.normalize:
                self._normalize_sequences(self.dataframe)
        return self.dataframe

    def get_dataframe_chunks(self) -> Iterator[pd.DataFrame]:
        """
        Loads the data in chunks

        Yields pandas DataFrame with index "seqid" and column "sequence"
        """
        if self.dataframe is not None:
            yield self.dataframe
            return
        assert self.path is not None
        for dataframe in self.protocol.read(self.path, columns=["seqid", "sequence"]):
            try:
                dataframe.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSeqids()
            else:
                if self.normalize:
                    self._normalize_sequences(dataframe)
                yield dataframe


class CompleteData(DataType):
    """
    Represents a collection of sequences with additional data.

    Only loads data on demand.
    """

    def __init__(self, path: Optional[ValidFilePath], protocol: FileReader):
        self.path = path
        self.protocol = protocol
        self.normalize = False
        self.dataframe: Optional[pd.DataFrame] = None

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> CompleteData:
        return cls(path, protocol)

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame) -> CompleteData:
        self = cls(None, TabfileReader())
        assert "sequence" in dataframe.columns
        assert dataframe.index.name == "seqid"
        self.dataframe = dataframe
        return self

    @staticmethod
    def _normalize_sequences(table: pd.DataFrame) -> None:
        assert "sequence" in table.columns
        table["sequence"] = (
            table["sequence"]
            .str.upper()
            .str.replace("?", "N", regex=False)
            .str.replace("-", "", regex=False)
        )

    def normalize_sequences(self) -> None:
        self.normalize = True
        if self.dataframe is not None:
            self._normalize_sequences(self.dataframe)

    def get_dataframe(self) -> pd.DataFrame:
        """
        Loads the data, if necessary, and returns it.

        Returns pandas DataFrame with index "seqid" and with a column "sequence"
        """
        if self.dataframe is None:
            assert self.path is not None
            self.dataframe = self.protocol.read(self.path, columns=[])
            if not {"seqid", "sequence"} in set(self.dataframe.columns):
                raise ColumnsNotFound({"seqid", "sequence"})
            try:
                self.dataframe.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSeqids()
            if self.normalize:
                self._normalize_sequences(self.dataframe)
        return self.dataframe

    def get_dataframe_chunks(self) -> Iterator[pd.DataFrame]:
        """
        Loads the data in chunks

        Yields pandas DataFrame with index "seqid" and column "sequence"
        """
        if self.dataframe is not None:
            yield self.dataframe
            return
        assert self.path is not None
        for dataframe in self.protocol.read(self.path, columns=[]):
            if not {"seqid", "sequence"} in set(dataframe.columns):
                raise ColumnsNotFound({"seqid", "sequence"})
            try:
                dataframe.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSeqids()
            else:
                if self.normalize:
                    self._normalize_sequences(dataframe)
                yield dataframe

    def get_chunks(self) -> Iterator[CompleteData]:
        for dataframe in self.get_dataframe_chunks():
            yield CompleteData.from_dataframe(dataframe)

    def get_sequences(self) -> SequenceData:
        sequences = SequenceData(self.path, self.protocol)
        if self.dataframe:
            sequences.dataframe = self.dataframe[["sequence"]].copy()
        return sequences

    def split_sequences(self, sequences: SequenceData) -> SequenceData:
        """
        Removes rows with seqids that are not in `sequences`.

        Returns the removed rows.
        """
        sequences_table = sequences.get_dataframe()
        dataframe = self.get_dataframe()
        removed = SequenceData.from_dataframe(
            dataframe.reindex(dataframe.index.difference(sequences_table.index))
        )
        self.dataframe = dataframe[dataframe.index.intersection(sequences_table.index)]
        return removed

    def append_to_file(self, file: Path) -> None:
        for dataframe in self.get_dataframe_chunks():
            _dataframe_append(dataframe, file)


class Metric(Enum):
    """
    Types of metric for calculating distance between sequences
    """

    Uncorrected = auto()
    Kimura2P = auto()
    JukesCantor = auto()
    UncorrectedWithGaps = auto()

    def __str__(self) -> str:
        return {
            Metric.Uncorrected: "pairwise uncorrected distance",
            Metric.Kimura2P: "Jukes-Cantor distance",
            Metric.JukesCantor: "Kimura-2-Parameter distance",
            Metric.UncorrectedWithGaps: "pairwise uncorrected distance counting gaps",
        }[self]


class SequenceDistanceMatrix(DataType):
    """
    Represents pairwise distances between sequences.
    """

    def __init__(self, distance_matrix: pd.DataFrame):
        assert list(distance_matrix.index.names) == ["seqid_target", "seqid_query"]
        assert set(distance_matrix.columns) <= set(Metric)
        assert distance_matrix.index.is_unique
        self.distance_matrix = distance_matrix

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> SequenceDistanceMatrix:
        raise NotImplementedError

    def get_dataframe(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with MultiIndex ["seqid_target", "seqid_query"]
        and columns, whose names are a subset of `Metric`
        """
        return self.distance_matrix


class SpeciesPartition(DataType):
    """
    Represents partition of sequences into species
    """

    def __init__(self, species_partition: pd.DataFrame):
        assert list(species_partition.index.names) == ["seqid"]
        assert list(species_partition.columns) == ["species"]
        assert species_partition.index.is_unique()
        self.species_partition = species_partition

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> SpeciesPartition:
        species_partition = protocol.read(path, columns=["seqid", "species"])
        try:
            species_partition.set_index("seqid", inplace=True, verify_integrity=True)
        except ValueError:
            raise DuplicatedSeqids
        return cls(species_partition)

    def get_dataframe(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "seqid" and column "species"
        """
        return self.species_partition


class SubsubspeciesPartition(DataType):
    """
    Represents partition of sequences into subspecies
    """

    def __init__(self, subspecies_partition: pd.DataFrame):
        assert list(subspecies_partition.index.names) == ["seqid"]
        assert list(subspecies_partition.columns) == ["subspecies"]
        assert subspecies_partition.index.is_unique()
        self.subspecies_partition = subspecies_partition

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> SubsubspeciesPartition:
        subspecies_partition = protocol.read(path, columns=["seqid", "subspecies"])
        try:
            subspecies_partition.set_index("seqid", inplace=True, verify_integrity=True)
        except ValueError:
            raise DuplicatedSeqids
        return cls(subspecies_partition)

    def get_dataframe(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "seqid" and column "subspecies"
        """
        return self.subspecies_partition


class NonBinomialSpecies(Exception):
    """
    List of non-binomial name encountered
    """

    def __init__(self, args: List[str]):
        super().__init__()
        assert isinstance(args, list)
        self.args = (args,)


class GenusPartition(DataType):
    """
    Represents partition of sequence into species with binomial names
    """

    def __init__(self, genus_partition: pd.DataFrame):
        assert list(genus_partition.index.names) == ["seqid"]
        assert list(genus_partition.columns) == ["genus", "species"]
        assert genus_partition.index.is_unique()
        self.genus_partition = genus_partition

    @staticmethod
    def _separate_species(species: pd.DataFrame) -> pd.DataFrame:
        try:
            species.set_index("seqid", inplace=True, verify_integrity=True)
        except ValueError:
            raise DuplicatedSpecies
        genus_partition = species["species"].str.split(
            r" |_", n=1, expand=True, regex=True
        )
        genus_partition.columns = ["genus", "species"]
        genus_partition["binomial"] = species["species"]

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> GenusPartition:
        try:
            genus_partition = protocol.read(path, columns=["seqid", "species", "genus"])
            try:
                genus_partition.set_index("seqid", inplace=True, verify_integrity=True)
            except ValueError:
                raise DuplicatedSpecies
        except ColumnsNotFound:
            genus_partition = GenusPartition._separate_species(
                protocol.read(path, columns=["seqid", "species"])
            )
        return cls(genus_partition)

    def get_dataframe(self) -> GenusPartition:
        """
        Returns a pandas DataFrame with indes "seqid" and columns "genus" and "species".

        "genus" contains the generic name of the species.
        "species" contains the specific epithet of the species.
        """
        return self.genus_partition


class FileReader(ABC):
    """
    Abstract base class for readers of specific file types
    """

    DEFAULT_CHUNK_SIZE: int = 1000

    @staticmethod
    @abstractmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        """
        Try to read a pandas DataFrame with columns in `columns` from `path`.

        If `columns` is [], all columns will be loaded.

        Can raise `ColumnsNotFound`
        """
        pass

    @staticmethod
    @abstractmethod
    def read_chunks(
        path: ValidFilePath, *, columns: List[str], chunksize: int = DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        """
        Try to read a pandas DataFrame with columns in `columns` from `path` in chunks.

        If `columns` is [], all columns will be loaded.

        Can raise `ColumnsNotFound`
        """
        pass

    @staticmethod
    @abstractmethod
    def read_data(path: ValidFilePath) -> List[DataType]:

        pass


class ColumnsNotFound(Exception):
    """
    Set of requested columns that are missing
    """

    def __init__(self, columns: Set[str]):
        if not isinstance(columns, set):
            raise TypeError("ColumnsNotFound argument should be a set")
        self.args = (columns,)


class TabfileReader(FileReader):
    """
    Reader for tab-separated files
    """

    @staticmethod
    def _verify_columns(path: ValidFilePath, *, columns: List[str]) -> List[str]:
        with open(path, errors="replace") as file:
            header = _Header(file.readline().rstrip("\n").split("\t"))
        if not set(columns) <= set(header.names):
            raise ColumnsNotFound(set(columns) - set(header.names))
        return header.names

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        names = TabfileReader._verify_columns(path, columns=columns)
        return pd.read_table(
            path,
            header=0,
            names=names,
            usecols=columns or names,
            na_filter=False,
            dtype=str,
            encoding_errors="replace",
        )

    @staticmethod
    def read_chunks(
        path: ValidFilePath,
        *,
        columns: List[str],
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        names = TabfileReader._verify_columns(path, columns=columns)
        yield from pd.read_table(
            path,
            header=0,
            names=names,
            usecols=columns or names,
            na_filter=False,
            dtype=str,
            encoding_errors="replace",
            chunksize=chunksize,
        )

    @staticmethod
    def read_data(path: ValidFilePath) -> List[DataType]:
        data: List[DataType] = []

        datatypes: List[Type[DataType]] = [
            SequenceData,
            SpeciesPartition,
            SubsubspeciesPartition,
            GenusPartition,
        ]
        for datatype in datatypes:
            try:
                data.append(datatype.from_path(path, TabfileReader()))
            except Exception:
                pass
        return data


class XlsxReader(FileReader):
    """
    Reader for Excel's xlsx files

    Reads the first sheet
    """

    @staticmethod
    def _verify_columns(path: ValidFilePath, *, columns: List[str]) -> List[str]:
        worksheet = openpyxl.load_workbook(path).worksheets[0]
        header = _Header(list(map(lambda cell: cell.value, next(worksheet.rows))))
        if not set(columns) <= set(header.names):
            raise ColumnsNotFound(set(columns) - set(header.names))
        return header.names

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        names = XlsxReader._verify_columns(path, columns=columns)
        return pd.read_excel(
            path,
            sheet_name=0,
            header=0,
            names=names,
            usecols=columns or names,
            na_filter=False,
            dtype=str,
            engine="openpyxl",
            encoding_errors="replace",
        )

    @staticmethod
    def read_chunks(
        path: ValidFilePath,
        *,
        columns: List[str],
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        names = XlsxReader._verify_columns(path, columns=columns)
        yield from pd.read_table(
            path,
            sheet_name=0,
            header=0,
            names=names,
            usecols=columns or names,
            na_filter=False,
            dtype=str,
            engine="openpyxl",
            encoding_errors="replace",
            chunksize=chunksize,
        )

    @staticmethod
    def read_data(path: ValidFilePath) -> List[DataType]:
        data: List[DataType] = []

        datatypes: List[Type[DataType]] = [
            SequenceData,
            SpeciesPartition,
            SubsubspeciesPartition,
            GenusPartition,
        ]
        for datatype in datatypes:
            try:
                data.append(datatype.from_path(path, XlsxReader()))
            except Exception:
                pass
        return data


class FastaReader(FileReader):
    """
    Reader for Fasta files
    """

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        if not set(columns) <= {"seqid", "sequence"}:
            raise ColumnsNotFound(set(columns) - {"seqid", "sequence"})
        with open(path, errors="replace") as file:
            _, records = Fastafile.read(file)
            return pd.DataFrame(
                ([record["seqid"], record["sequence"]] for record in records()),
                columns=["seqid", "sequence"],
            )

    @staticmethod
    def read_chunks(
        path: ValidFilePath,
        *,
        columns: List[str],
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE
    ) -> pd.DataFrame:
        if not set(columns) <= {"seqid", "sequence"}:
            raise ColumnsNotFound(set(columns) - {"seqid", "sequence"})
        with open(path, errors="replace") as file:
            _, records = Fastafile.read(file)
            records = records()
            while True:
                table = pd.DataFrame(
                    (
                        record["seqid", "sequence"]
                        for record in itertools.islice(records, chunksize)
                    ),
                    columns=["seqid", "sequence"],
                )
                if table.empty:
                    break
                yield table

    @staticmethod
    def read_data(path: ValidFilePath) -> List[DataType]:
        return [SequenceData.from_path(path, FastaReader())]


class GenbankReader(FileReader):
    """
    Reader for Genbank flat files
    """

    @staticmethod
    def select_columns(columns: List[str], record: Record) -> List[str]:
        return list(
            record[column] if column != "species" else record["organism"]
            for column in columns
        )

    @staticmethod
    def _verify_columns(path: ValidFilePath, *, columns: List[str]) -> List[str]:
        with open(path, errors="replace") as file:
            names, _ = GenbankFile.read(file)
        if not set(columns) <= set(names):
            raise ColumnsNotFound(set(columns) - set(names))
        return names

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        names = GenbankReader._verify_columns(path, columns=columns)
        columns = columns or names
        with open(path, errors="replace") as file:
            _, records = GenbankFile.read(file)
            return pd.DataFrame(
                (GenbankReader.select_columns(columns, record) for record in records()),
                columns=columns or names,
            )

    @staticmethod
    def read_chunks(
        path: ValidFilePath,
        *,
        columns: List[str],
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE
    ) -> Iterator[pd.DataFrame]:
        names = GenbankReader._verify_columns(path, columns=columns)
        columns = columns or names
        with open(path, errors="replace") as file:
            _, records = GenbankFile.read(file)
            records = records()
            while True:
                table = pd.DataFrame(
                    (
                        GenbankReader.select_columns(columns, record)
                        for record in itertools.islice(records, chunksize)
                    ),
                    columns=columns,
                )
                if table.empty:
                    break
                yield table

    @staticmethod
    def read_data(path: ValidFilePath) -> List[DataType]:
        data: List[DataType] = []

        datatypes: List[Type[DataType]] = [
            SequenceData,
            SpeciesPartition,
            SubsubspeciesPartition,
            GenusPartition,
        ]
        for datatype in datatypes:
            try:
                data.append(datatype.from_path(path, GenbankReader()))
            except Exception:
                pass
        return data
