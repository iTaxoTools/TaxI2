#!/usr/bin/env python3
from __future__ import annotations

from typing import List, Set, Any, Optional, Iterator, Type, Dict
from abc import ABC, abstractmethod
from pathlib import Path
import os
import re
import itertools
from enum import Enum, auto
from dataclasses import dataclass

import pandas as pd
import openpyxl

from itaxotools.DNAconvert.library.fasta import Fastafile
from itaxotools.DNAconvert.library.genbank import GenbankFile
from itaxotools.DNAconvert.library.record import Record

from .sequence_statistics import sequence_statistics, sequence_statistics_with_gaps


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
        self.read_text = self._path.read_text
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
    dataframe.to_csv(file, header=header, sep="\t", mode="a")


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
        for dataframe in self.protocol.read_chunks(
            self.path, columns=["seqid", "sequence"]
        ):
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
        for dataframe in self.protocol.read_chunks(self.path, columns=[]):
            if not {"seqid", "sequence"} <= set(dataframe.columns):
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
        if self.dataframe is not None:
            sequences.dataframe = self.dataframe[["sequence"]].copy()
        return sequences

    def split_sequences(self, sequences: SequenceData) -> CompleteData:
        """
        Removes rows with seqids that are not in `sequences`.

        Returns the removed rows.
        """
        sequences_table = sequences.get_dataframe()
        dataframe = self.get_dataframe()
        removed = CompleteData.from_dataframe(
            dataframe.reindex(dataframe.index.difference(sequences_table.index))
        )
        self.dataframe = dataframe.loc[
            dataframe.index.intersection(sequences_table.index)
        ]
        return removed

    def append_to_file(self, file: Path) -> None:
        for dataframe in self.get_dataframe_chunks():
            _dataframe_append(dataframe, file)


class DecontaminateSummary(DataType):
    """
    Represent the summary of decontaminate task
    """

    def __init__(self, dataframe: pd.DataFrame):
        assert dataframe.index.name == "seqid_query"
        assert set(dataframe.columns) == {
            "closest possible contaminant",
            "distance",
            "is_contaminant",
        }
        self.dataframe = dataframe

    def get_dataframe(self) -> pd.DataFrame:
        return self.dataframe

    @classmethod
    def from_path(
        self, path: ValidFilePath, protocol: FileReader
    ) -> DecontaminateSummary:
        raise NotImplementedError

    def append_to_file(self, file: Path) -> None:
        _dataframe_append(self.dataframe, file)


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

    @staticmethod
    def short_names() -> Dict[Metric, str]:
        return {
            Metric.Uncorrected: "p-distance",
            Metric.Kimura2P: "JC distance",
            Metric.JukesCantor: "K2P distance",
            Metric.UncorrectedWithGaps: "p-distance with gaps",
        }


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

    def set_self_distance(self, self_distance: float) -> None:
        seqids = self.distance_matrix.index.to_frame()
        self.distance_matrix.loc[
            seqids["seqid_target"] == seqids["seqid_query"]
        ] = self_distance


class SpeciesPartition(DataType):
    """
    Represents partition of sequences into species
    """

    def __init__(self, species_partition: pd.DataFrame):
        assert list(species_partition.index.names) == ["seqid"]
        assert list(species_partition.columns) == ["species"]
        assert species_partition.index.is_unique
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


class VoucherPartition(DataType):
    """
    Represents partition of sequences into specimens
    """

    def __init__(self, voucher_partition: pd.DataFrame):
        assert list(voucher_partition.index.names) == ["seqid"]
        assert list(voucher_partition.columns) == ["specimen_voucher"]
        assert voucher_partition.index.is_unique
        self.voucher_partition = voucher_partition

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> VoucherPartition:
        voucher_partition = protocol.read(path, columns=["seqid", "specimen_voucher"])
        try:
            voucher_partition.set_index("seqid", inplace=True, verify_integrity=True)
        except ValueError:
            raise DuplicatedSeqids
        return cls(voucher_partition)

    def get_dataframe(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with index "seqid" and column "species_voucher"
        """
        return self.voucher_partition


class SubsubspeciesPartition(DataType):
    """
    Represents partition of sequences into subspecies
    """

    def __init__(self, subspecies_partition: pd.DataFrame):
        assert list(subspecies_partition.index.names) == ["seqid"]
        assert list(subspecies_partition.columns) == ["subspecies"]
        assert subspecies_partition.index.is_unique
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
        assert genus_partition.index.is_unique
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
        return genus_partition

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


class RowOrdering(Enum):
    Default = auto()
    Alphabetical = auto()
    Species = auto()

    def description(self) -> str:
        return {
            RowOrdering.Default: "",
            RowOrdering.Alphabetical: "Alphabetical order",
            RowOrdering.Species: "Ordered by species",
        }[self]


class SequenceDistanceOutput(DataType):
    """
    Represents tables of distances between sequence in a form suitable for output.
    """

    def __init__(
        self,
        dataframe: pd.DataFrame,
        *,
        metric: Metric,
        ordering: RowOrdering,
        in_percent: bool,
    ):
        assert dataframe.index.equals(dataframe.columns)
        if ordering is RowOrdering.Species:
            assert list(dataframe.index.names) == ["species", "seqid"]
        else:
            assert dataframe.index.name == "seqid"
        self.dataframe = dataframe
        self.ordering = ordering
        self.in_percent = in_percent
        self.metric = metric

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> SequenceDistanceOutput:
        raise NotImplementedError

    def get_dataframe(self) -> pd.DataFrame:
        return self.dataframe

    def description(self) -> str:
        if self.ordering is RowOrdering.Default:
            ordering_s = ""
        else:
            ordering_s = f" ({self.ordering.description()})"
        return f"{self.metric} between sequences{ordering_s}"

    def make_file_name(self) -> str:
        return (
            self.description().replace(" ", "_")
            + ("_in_percent" if self.in_percent else "")
            + ".txt"
        )

    def append_to_file(self, file: Path) -> None:
        l = len(self.dataframe.index.names)
        self.dataframe.index.names = [None] * l
        self.dataframe.columns.names = [None] * l
        self.dataframe.to_csv(file, float_format="%.4g", sep="\t", mode="a")


class SimpleSequenceStatistic(DataType):
    """
    A column of simple sequence statistics
    """

    def __init__(self, column: pd.Series):
        assert column.name is None
        assert set(column.index) == set.union(
            set(sequence_statistics), set(sequence_statistics_with_gaps)
        )
        self.column = column

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> SimpleSequenceStatistic:
        raise NotImplementedError

    def get_dataframe(self) -> pd.DataFrame:
        return self.column.to_frame()

    def append_to_file(self, path: Path) -> None:
        with open(path, mode="a") as output_file:
            for stat_name, stat_value in self.column.items():
                print(stat_name, f"{stat_value:.4g}", sep=": ", file=output_file)


class SimpleSpeciesStatistic(DataType):
    """
    A table of simple sequence statistics by species
    """

    def __init__(self, table: pd.DataFrame):
        assert table.index.name == "species"
        assert set(table.columns) == set.union(
            set(sequence_statistics), set(sequence_statistics_with_gaps)
        )
        self.dataframe = table

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> SimpleSequenceStatistic:
        raise NotImplementedError

    def get_dataframe(self) -> pd.DataFrame:
        return self.dataframe

    def append_to_file(self, path: Path) -> None:
        self.dataframe.to_csv(path, sep="\t", float_format="%.4g", mode="a")


class TaxonRank(Enum):
    Species = auto()
    Genus = auto()

    def __str__(self) -> str:
        return {
            TaxonRank.Species: "species",
            TaxonRank.Genus: "genus",
        }[self]

    def plural_str(self) -> str:
        return {
            TaxonRank.Species: "species",
            TaxonRank.Genus: "genera",
        }[self]


class Aggregation(Enum):
    Mean = auto()
    Min = auto()
    Max = auto()


@dataclass(frozen=True)
class AggregatedMetric:
    aggregation: Aggregation
    metric: Metric

    @classmethod
    def min(cls, metric: Metric) -> AggregatedMetric:
        return cls(Aggregation.Min, metric)

    @classmethod
    def max(cls, metric: Metric) -> AggregatedMetric:
        return cls(Aggregation.Max, metric)

    @classmethod
    def mean(cls, metric: Metric) -> AggregatedMetric:
        return cls(Aggregation.Mean, metric)


class MeanMinMaxDistances(DataType):
    """
    Three column dataframes of mean, mininum and maximum distances
    """

    def __init__(
        self,
        mean_distances: pd.DataFrame,
        min_distances: pd.DataFrame,
        max_distances: pd.DataFrame,
        *,
        metric: Metric,
        taxon_rank: TaxonRank,
        in_percent: bool,
    ):
        assert mean_distances.index.equals(min_distances.index)
        assert mean_distances.index.equals(max_distances.index)
        assert len(mean_distances.index.names) <= 2
        self.is_square = len(mean_distances.index.names) == 2
        self.metric = metric
        self.taxon_rank = taxon_rank
        self.in_percent = in_percent
        if self.is_square:
            index_names = [
                str(self.taxon_rank) + "_target",
                str(self.taxon_rank) + "_query",
            ]
            assert list(mean_distances.index.names) == index_names
            assert list(min_distances.index.names) == index_names
            assert list(max_distances.index.names) == index_names
        else:
            assert mean_distances.index.name == str(self.taxon_rank)
            assert min_distances.index.name == str(self.taxon_rank)
            assert max_distances.index.name == str(self.taxon_rank)
        assert mean_distances.columns == [self.metric]
        assert min_distances.columns == [self.metric]
        assert max_distances.columns == [self.metric]
        self.mean = mean_distances
        self.min = min_distances
        self.max = max_distances

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> MeanMinMaxDistances:
        raise NotImplementedError

    def get_dataframe(self) -> pd.DataFrame:
        """
        Returns a dataframe with the index:
            str(self.taxon_rank) or
            [str(self.taxon_rank) + "_target", str(self.taxon_rank) + "_query"]
        and columns:
            AggregatedMetric.mean(metric)
            AggregatedMetric.min(metric)
            AggregatedMetric.max(metric)
        """
        mean_distances = self.mean.rename(columns=AggregatedMetric.mean)
        min_distances = self.min.rename(columns=AggregatedMetric.min)
        max_distances = self.max.rename(columns=AggregatedMetric.max)
        return mean_distances.join(min_distances).join(max_distances)


class Source(Enum):
    Query1 = auto()
    Query2 = auto()
    Reference = auto()

    def __str__(self) -> str:
        return {
            Source.Query1: "query 1",
            Source.Query2: "query 2",
            Source.Reference: "reference",
        }[self]


@dataclass(frozen=True)  # frozen=True is required for hashing
class SourcedColumn:
    name: str
    source: Source

    @classmethod
    def query1(cls, name: str) -> SourcedColumn:
        return cls(name=name, source=Source.Query1)

    @classmethod
    def query2(cls, name: str) -> SourcedColumn:
        return cls(name=name, source=Source.Query2)

    @classmethod
    def reference(cls, name: str) -> SourcedColumn:
        return cls(name=name, source=Source.Reference)

    def __str__(self) -> str:
        return f"{self.name} ({self.source})"


class VersusAllSummary(DataType):
    def __init__(self, dataframe: pd.DataFrame):
        assert set(dataframe.columns) <= {
            SourcedColumn("seqid", Source.Query1),
            SourcedColumn("specimen_voucher", Source.Query1),
            SourcedColumn("genus", Source.Query1),
            SourcedColumn("species", Source.Query1),
            SourcedColumn("seqid", Source.Query2),
            SourcedColumn("specimen_voucher", Source.Query2),
            SourcedColumn("genus", Source.Query2),
            SourcedColumn("species", Source.Query2),
            "comparison_type",
            Metric.Uncorrected,
            Metric.JukesCantor,
            Metric.Kimura2P,
            Metric.UncorrectedWithGaps,
        }
        self.dataframe = dataframe

    @classmethod
    def from_path(cls, path: ValidFilePath, protocol: FileReader) -> VersusAllSummary:
        raise NotImplementedError

    def get_dataframe(self) -> pd.DataFrame:
        return self.dataframe

    def to_file(self, path: Path) -> None:
        column_order = [
            column
            for column in [
                SourcedColumn("seqid", Source.Query1),
                SourcedColumn("specimen_voucher", Source.Query1),
                SourcedColumn("genus", Source.Query1),
                SourcedColumn("species", Source.Query1),
                "comparison_type",
                SourcedColumn("seqid", Source.Query2),
                SourcedColumn("specimen_voucher", Source.Query2),
                SourcedColumn("genus", Source.Query2),
                SourcedColumn("species", Source.Query2),
                Metric.Uncorrected,
                Metric.Kimura2P,
                Metric.JukesCantor,
                Metric.UncorrectedWithGaps,
            ]
            if column in self.dataframe.columns
        ]
        with open(path, mode="w") as file:
            print("Summary statistics", file=file)
        different_seqids = (
            self.dataframe[SourcedColumn.query1("seqid")]
            != self.dataframe[SourcedColumn.query2("seqid")]
        )
        self.dataframe.loc[different_seqids, column_order].rename(
            columns=Metric.short_names()
        ).to_csv(path, sep="\t", mode="a", index=False, float_format="%.4g")


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
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE,
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
            VoucherPartition,
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
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE,
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
            VoucherPartition,
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
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE,
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
        chunksize: int = FileReader.DEFAULT_CHUNK_SIZE,
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
            VoucherPartition,
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
