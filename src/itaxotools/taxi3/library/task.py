#!/usr/bin/env python3
from __future__ import annotations

from typing import Callable, Generic, TypeVar, Optional, List, Iterator, Union
from abc import ABC, abstractmethod
from enum import Enum, auto
from dataclasses import dataclass
from collections import defaultdict

import pandas as pd
import networkx as nx

from .datatypes import (
    DataType,
    FileReader,
    SequenceDistanceMatrix,
    SequenceData,
    Metric,
    CompleteData,
    DecontaminateSummary,
    VersusAllSummary,
    SourcedColumn,
    Source,
    VoucherPartition,
    GenusPartition,
    SpeciesPartition,
    SubsubspeciesPartition,
    SequenceDistanceOutput,
    RowOrdering,
    SimpleSequenceStatistic,
    SimpleSpeciesStatistic,
)
from .rust_backend import calc, make_aligner
from .sequence_statistics import (
    sequence_statistics,
    sequence_statistics_with_gaps,
)

from .alfpy_distance import alfpy_distance_array, alfpy_distance_array2

_Result = TypeVar("_Result")


class Alignment(Enum):
    Pairwise = auto()
    AlignmentFree = auto()
    AlreadyAligned = auto()


@dataclass()
class SequencesPair:
    target: SequenceData
    query: SequenceData

    def normalize_sequences(self) -> None:
        self.target.normalize_sequences()
        self.query.normalize_sequences()


class MissingArgument(Exception):
    pass


WarningHandler = Callable[[Warning], None]


class Task(ABC, Generic[_Result]):
    """
    To run:
        * Set the argument attributes (depends on the subclass)
        * call self.start()
        * self.result contains the output

    If self is instance of Task[T], self.result is instance of T
    """

    def __init__(self, warn: WarningHandler):
        self.warn = warn
        self.result: Optional[_Result] = None

    @abstractmethod
    def start(self) -> None:
        pass


class CalculateDistances(Task[SequenceDistanceMatrix]):
    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.sequences: Union[None, SequenceData, SequencesPair] = None
        self.alignment: Optional[Alignment] = None
        self.metrics: List[Metric] = []

    def _alignment_free_start(self) -> None:
        assert self.sequences is not None
        if isinstance(self.sequences, SequenceData):
            sequences = self.sequences.get_dataframe()
            distances_array = alfpy_distance_array(sequences["sequence"])
            seqids1 = sequences.index.copy()
            seqids2 = sequences.index.copy()
        elif isinstance(self.sequences, SequencesPair):
            reference = self.sequences.target.get_dataframe()
            sequences = self.sequences.query.get_dataframe()
            distances_array = alfpy_distance_array2(
                reference["sequence"], sequences["sequence"]
            )
            seqids1 = reference.index.copy()
            seqids2 = sequences.index.copy()
        else:
            assert False
        index = pd.MultiIndex.from_product(
            [seqids1, seqids2], names=["seqid_target", "seqid_query"]
        )
        self.result = SequenceDistanceMatrix(
            pd.DataFrame(distances_array, index=index, columns=[Metric.Uncorrected])
        )

    def _alignment_start(self) -> None:
        assert self.alignment is not None
        assert self.sequences is not None
        assert self.metrics
        if isinstance(self.sequences, SequenceData):
            sequences = self.sequences.get_dataframe()
            reference = sequences
            seqids1 = sequences.index.copy()
            seqids2 = sequences.index.copy()
        elif isinstance(self.sequences, SequencesPair):
            reference = self.sequences.target.get_dataframe()
            sequences = self.sequences.query.get_dataframe()
            seqids1 = reference.index.copy()
            seqids2 = sequences.index.copy()
        else:
            assert False
        if self.alignment is Alignment.AlreadyAligned:
            distances_array = calc.make_distance_array_aligned(
                reference["sequence"], sequences["sequence"]
            )
        elif self.alignment is Alignment.Pairwise:
            aligner = make_aligner()
            distances_array = calc.make_distance_array(
                aligner, reference["sequence"], sequences["sequence"]
            )
        else:
            assert self.alignment is not Alignment.AlignmentFree
        index = pd.MultiIndex.from_product(
            [seqids1, seqids2], names=["seqid_target", "seqid_query"]
        )
        self.result = SequenceDistanceMatrix(
            pd.DataFrame(distances_array, index=index, columns=list(Metric))[
                self.metrics
            ]
        )

    def start(self) -> None:
        if self.sequences is None:
            raise MissingArgument("sequences")
        if self.alignment is None:
            raise MissingArgument("alignment")
        if self.alignment is not Alignment.AlreadyAligned:
            self.sequences.normalize_sequences()
        if self.alignment is Alignment.AlignmentFree:
            self._alignment_free_start()
        elif not self.metrics:
            raise MissingArgument("metrics")
        else:
            self._alignment_start()


@dataclass
class VersusAllSummarizeArg:
    vouchers: Optional[VoucherPartition]
    species: Optional[GenusPartition]
    subspecies: Optional[SubsubspeciesPartition]

    @classmethod
    def from_path(
        cls, path: ValidFilePath, protocol: FileReader
    ) -> VersusAllSummarizeArg:
        self = cls(None, None, None)
        try:
            self.vouchers = VoucherPartition.from_path(path, protocol)
        except Exception:
            pass
        try:
            self.species = GenusPartition.from_path(path, protocol)
        except Exception:
            pass
        try:
            self.subspecies = SubsubspeciesPartition.from_path(path, protocol)
        except Exception:
            pass
        return self

    @classmethod
    def from_list(cls, data: List[DataType]) -> VersusAllSummarizeArg:
        self = cls(None, None, None)
        for table in data:
            if isinstance(table, VoucherPartition):
                self.vouchers = table
            if isinstance(table, GenusPartition):
                self.species = table
            if isinstance(table, SubsubspeciesPartition):
                self.subspecies = table
        return self

    def to_list(self) -> List[DataType]:
        return [
            table
            for table in (self.vouchers, self.species, self.subspecies)
            if table is not None
        ]


class VersusAllSummarize(Task[VersusAllSummary]):
    """
    Arguments:
        distances: SequenceDistanceMatrix
        data: VersusAllSummarizeArg
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.distances: Optional[SequenceDistanceMatrix] = None
        self.data: Optional[VersusAllSummarizeArg] = None

    def start(self) -> None:
        if self.distances is None:
            raise MissingArgument("distances")
        if self.data is None:
            raise MissingArgument("data")
        summary_table = (
            self.distances.get_dataframe()
            .reset_index()
            .rename(
                columns={
                    "seqid_target": SourcedColumn.query1("seqid"),
                    "seqid_query": SourcedColumn.query2("seqid"),
                }
            )
        )
        for table in self.data.to_list():
            for decorator in (SourcedColumn.query1, SourcedColumn.query2):
                query_dataframe = table.get_dataframe().rename(columns=decorator)
                query_dataframe.index.name = decorator(query_dataframe.index.name)
                summary_table = summary_table.join(
                    query_dataframe, on=decorator("seqid"), how="left"
                )

        summary_table["comparison_type"] = ""

        ones = pd.Series("1", index=summary_table.index)
        for taxon_level in ("genus", "species", "subspecies"):
            try:
                summary_table["comparison_type"] += ones.where(
                    summary_table[SourcedColumn.query1(taxon_level)]
                    == summary_table[SourcedColumn.query2(taxon_level)],
                    other="0",
                )
            except KeyError:
                summary_table["comparison_type"] += " "

        summary_table["comparison_type"] = summary_table["comparison_type"].map(
            defaultdict(
                lambda: "Undefined",
                {
                    "00 ": "inter-genus",
                    "000": "inter-genus",
                    "001": "inter-genus",
                    "10 ": "inter-species",
                    "100": "inter-species",
                    "101": "inter-species",
                    "11 ": "intra-species",
                    "111": "intra-species (same intraspecific lineage)",
                    "110": "intra-species (different intraspecific lineages)",
                },
            )
        )
        self.result = VersusAllSummary(summary_table)


class OutputSequenceDistances(Task[SequenceDistanceOutput]):
    """
    Arguments:
        sequence_matrix: SequenceDistanceMatrix
        metric: Metric
        ordering: RowOrdering
        species: Optional[SpeciesPartition]
        in_percent: bool
    If ordering is RowOrdering.Species, species is not None
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.sequence_matrix: Optional[SequenceDistanceMatrix] = None
        self.metric: Optional[Metric] = None
        self.ordering: Optional[RowOrdering] = None
        self.species: Optional[SpeciesPartition] = None
        self.in_percent = False

    def start(self) -> None:
        assert self.sequence_matrix is not None
        assert self.metric is not None
        assert self.ordering is not None
        dataframe = self.sequence_matrix.get_dataframe()[self.metric].copy()
        seqids = dataframe.index.to_frame()
        dataframe.loc[seqids["seqid_target"] == seqids["seqid_query"]] = float("nan")
        if self.in_percent:
            dataframe *= 100
        dataframe = dataframe.unstack(level="seqid_query")
        if self.ordering is RowOrdering.Species:
            assert self.species is not None
            species_index = (
                self.species.get_dataframe()
                .copy()
                .reindex(dataframe.index)
                .reset_index()
            )
            species_columns = (
                self.species.get_dataframe()
                .copy()
                .reindex(dataframe.columns)
                .reset_index()
            )
            species_index.rename(columns={"seqid_target": "seqid"}, inplace=True)
            species_columns.rename(columns={"seqid_query": "seqid"}, inplace=True)
            species_index = species_index[["species", "seqid"]]
            species_columns = species_columns[["species", "seqid"]]
            dataframe.set_axis(
                pd.MultiIndex.from_frame(species_index), axis="index", inplace=True
            )
            dataframe.set_axis(
                pd.MultiIndex.from_frame(species_columns), axis="columns", inplace=True
            )
            dataframe.sort_index(axis="index", level="species", inplace=True)
            dataframe.sort_index(axis="columns", level="species", inplace=True)
        else:
            dataframe.index.name = "seqid"
            dataframe.columns.name = "seqid"
        if self.ordering is RowOrdering.Alphabetical:
            dataframe.sort_index(axis="index", inplace=True)
            dataframe.sort_index(axis="columns", inplace=True)
        self.result = SequenceDistanceOutput(
            dataframe,
            metric=self.metric,
            ordering=self.ordering,
            in_percent=self.in_percent,
        )


@dataclass(frozen=True)
class SimpleStatisticResult:
    total: SimpleSequenceStatistic
    by_species: Optional[SimpleSpeciesStatistic]


class CalculateSimpleStatistic(Task):
    """
    Arguments:
        sequences: SequenceData
        species: Optional[SpeciesPartition]
    """

    def __init__(self, warn: WarningHandler):
        self.sequences: Optional[SequenceData] = None
        self.species: Optional[SpeciesPartition] = None

    def start(self) -> None:
        assert self.sequences is not None
        sequences_with_gaps = self.sequences.get_dataframe().copy()
        if self.species is not None:
            sequences_with_gaps.join(self.species.get_dataframe(), how="left")
        sequences = sequences_with_gaps.copy()
        sequences["sequence"] = sequences["sequence"].str.replace("-", "", regex=False)
        statistic_column = pd.concat(
            [
                sequences["sequence"].agg(sequence_statistics),
                sequences_with_gaps["sequence"].agg(sequence_statistics_with_gaps),
            ]
        )
        statistic_table = None
        if self.species is not None:
            statistic_table = pd.concat(
                [
                    sequences.groupby("species")["sequences"].agg(sequence_statistics),
                    sequences_with_gaps.groupby("species")["sequences"].agg(
                        sequence_statistics_with_gaps
                    ),
                ],
                axis="columns",
            )
        self.result = SimpleStatisticResult(
            total=statistic_column, by_species=statistic_table
        )


@dataclass(frozen=True)
class Dereplicated:
    """
    The type of the output of Dereplicate class
    """

    included: CompleteData
    excluded: CompleteData


class Dereplicate(Task[Iterator[Dereplicated]]):
    """
    Arguments:
        similarity: float
        length_threshold: Optional[int]
        keep_most_complete: Bool
        data: CompleteData
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.similarity = 0.07
        self.length_threshold: Optional[int] = None
        self.keep_most_complete = False
        self.data: Optional[CompleteData] = None
        self._calculate_distances = CalculateDistances(warn)
        self._calculate_distances.alignment = Alignment.AlignmentFree

    def start(self) -> None:
        if self.data is None:
            raise MissingArgument("data")
        self.result = self._dereplicate()

    def _dereplicate(self) -> Iterator[Dereplicated]:
        assert self.data is not None
        for chunk in self.data.get_chunks():
            assert chunk.dataframe is not None
            if self.length_threshold:
                chunk.dataframe = chunk.dataframe.loc[
                    chunk.dataframe["sequence"].str.len() >= self.length_threshold
                ]
            sequences = chunk.get_sequences()
            self._calculate_distances.sequences = sequences
            self._calculate_distances.start()
            distance_table = self._calculate_distances.result
            assert distance_table is not None
            distance_table.set_self_distance(float("inf"))

            nodes = chunk.dataframe.index

            # calculating components
            connected_table = distance_table.get_dataframe().loc[
                (distance_table.get_dataframe()[Metric.Uncorrected] <= self.similarity)
            ]
            graph = nx.from_pandas_edgelist(
                connected_table.index.to_frame(),
                source="seqid_target",
                target="seqid_query",
            )
            graph.add_nodes_from(nodes)
            components = nx.connected_components(graph)

            seqids_dereplicated: List[str] = []
            assert sequences.dataframe is not None

            for component in components:
                chosen_seqid = (
                    sequences.dataframe.loc[list(component), "sequence"]
                    .str.len()
                    .idxmax()
                )
                seqids_dereplicated.append(chosen_seqid)

            sequences.dataframe = sequences.dataframe.loc[seqids_dereplicated]

            excluded = chunk.split_sequences(sequences)

            yield Dereplicated(included=chunk, excluded=excluded)


@dataclass
class Decontaminated:
    """
    The type of output of Decontaminate task
    """

    decontaminated: CompleteData
    contaminates: CompleteData
    summary: DecontaminateSummary


class Decontaminate(Task[Iterator[Decontaminated]]):
    """
    Arguments:
        similarity: float
        alignment: Alignment
        reference: SequenceData
        data: CompleteData
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.similarity = 0.07
        self.alignment: Optional[Alignment] = None
        self.reference: Optional[SequenceData] = None
        self.data: Optional[CompleteData] = None
        self._calculate_distances = CalculateDistances(warn)
        self._calculate_distances.metrics = [Metric.Uncorrected]

    def start(self) -> None:
        if self.alignment is None:
            raise MissingArgument("alignment")
        if self.reference is None:
            raise MissingArgument("reference")
        if self.data is None:
            raise MissingArgument("data")
        self.result = self._decontaminate()

    def _decontaminate(self) -> Iterator[Decontaminated]:
        assert self.alignment is not None
        assert self.reference is not None
        assert self.data is not None
        self._calculate_distances.alignment = self.alignment

        for chunk in self.data.get_chunks():
            sequences = chunk.get_sequences()
            self._calculate_distances.sequences = SequencesPair(
                target=self.reference, query=sequences
            )
            self._calculate_distances.start()
            assert self._calculate_distances.result is not None
            distance_table = self._calculate_distances.result.get_dataframe()
            distance_table.reset_index(inplace=True)

            indices_closest = (
                distance_table[["seqid_query", Metric.Uncorrected]]
                .groupby("seqid_query")
                .idxmin()[Metric.Uncorrected]
                .squeeze()
                .dropna()
            )
            closest_table = distance_table.loc[indices_closest].rename(
                columns={"seqid_target": "closest possible contaminant"}
            )
            closest_table["is_contaminant"] = ""
            closest_table.loc[
                closest_table[Metric.Uncorrected] <= self.similarity, "is_contaminant"
            ] = "contaminant"
            decontaminates_seqids = closest_table.loc[
                closest_table[Metric.Uncorrected] > self.similarity,
                "seqid_query",
            ]
            closest_table.rename(columns={Metric.Uncorrected: "distance"}, inplace=True)
            closest_table.set_index("seqid_query", inplace=True)

            assert sequences.dataframe is not None
            sequences.dataframe = sequences.dataframe.loc[decontaminates_seqids]

            contaminates = chunk.split_sequences(sequences)
            yield Decontaminated(
                contaminates=contaminates,
                decontaminated=chunk,
                summary=DecontaminateSummary(closest_table),
            )