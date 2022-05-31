#!/usr/bin/env python3
from __future__ import annotations

from typing import (
    Callable,
    Generic,
    TypeVar,
    Optional,
    List,
    Iterator,
    Union,
    Type,
    Any,
)
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
    Decontaminate2Summary,
    VersusAllSummary,
    VersusReferenceSummary,
    SourcedColumn,
    VoucherPartition,
    GenusPartition,
    SpeciesPartition,
    SubspeciesPartition,
    SequenceDistanceOutput,
    RowOrdering,
    SimpleSequenceStatistic,
    SimpleSpeciesStatistic,
    MeanMinMaxDistances,
    TaxonRank,
    ValidFilePath,
)
from .rust_backend import calc, make_aligner
from .sequence_statistics import (
    sequence_statistics,
    sequence_statistics_with_gaps,
)

from .alfpy_distance import alfpy_distance_array, alfpy_distance_array2
from .config import Config

_Result = TypeVar("_Result")
_TaskClass = TypeVar("_TaskClass", bound="Task")


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


@dataclass(frozen=True)
class Progress:
    operation: str
    total_steps: Optional[int]  # None if total number of steps is unknown
    current_step: int


ProgressHandler = Callable[[Progress], None]


class Task(ABC, Generic[_Result]):
    """
    To run:
        * Set the argument attributes (depends on the subclass)
        * (optional) set self.config
        * call self.start()
        * self.result contains the output

    If self is instance of Task[T], self.result is instance of T
    """

    def __init__(self, warn: WarningHandler):
        self.warn = warn
        self.config: Optional[Config] = None
        self.result: Optional[_Result] = None
        self.progress_handler: Optional[ProgressHandler] = None
        self._total_steps: Optional[int] = None

    def subtask(self, task_class: Type[_TaskClass]) -> _TaskClass:
        """
        Create a subtask that inherits configuration, warning handler and progress handler
        """
        task = task_class(self.warn)
        task.config = self.config
        task.progress_handler = self.progress_handler
        return task

    def progress(self, *, operation: Any, step: int) -> None:
        if self.progress_handler:
            self.progress_handler(
                Progress(
                    operation=str(operation),
                    total_steps=self._total_steps,
                    current_step=step,
                )
            )

    @abstractmethod
    def start(self) -> None:
        pass


class CalculateDistances(Task[SequenceDistanceMatrix]):
    """
    Arguments:
        sequences: Union[SequenceData, SequencesPair]
        alignment: Alignment
        metrics: List[Metric]
    """

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
            if self.config is None:
                aligner = make_aligner()
            else:
                aligner = make_aligner(self.config.alignment_scores)
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
    subspecies: Optional[SubspeciesPartition]

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
            self.subspecies = SubspeciesPartition.from_path(path, protocol)
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
            if isinstance(table, SubspeciesPartition):
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
        in_percent: bool = False
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
        if self.sequence_matrix is None:
            raise MissingArgument("sequence_matrix")
        if self.metric is None:
            raise MissingArgument("metric")
        if self.ordering is None:
            raise MissingArgument("ordering")
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

    def make_file_name(self) -> str:
        return "Sequence summary statistic.txt"


class CalculateSimpleStatistic(Task[SimpleStatisticResult]):
    """
    Arguments:
        sequences: SequenceData
        species: Optional[SpeciesPartition]
    """

    def __init__(self, warn: WarningHandler):
        self.sequences: Optional[SequenceData] = None
        self.species: Optional[SpeciesPartition] = None

    def start(self) -> None:
        if self.sequences is None:
            raise MissingArgument("sequences")
        sequences_with_gaps = self.sequences.get_dataframe().copy()
        if self.species is not None:
            sequences_with_gaps = sequences_with_gaps.join(
                self.species.get_dataframe(), how="left"
            )
        sequences = sequences_with_gaps.copy()
        sequences["sequence"] = sequences["sequence"].str.replace("-", "", regex=False)
        statistic_column = pd.concat(
            [
                sequences["sequence"].agg(sequence_statistics),
                sequences_with_gaps["sequence"].agg(sequence_statistics_with_gaps),
            ]
        )
        statistic_column.name = None
        statistic_table = None
        if self.species is not None:
            statistic_table = pd.concat(
                [
                    sequences.groupby("species")["sequence"].agg(**sequence_statistics),
                    sequences_with_gaps.groupby("species")["sequence"].agg(
                        **sequence_statistics_with_gaps
                    ),
                ],
                axis="columns",
            )
        self.result = SimpleStatisticResult(
            total=SimpleSequenceStatistic(statistic_column),
            by_species=SimpleSpeciesStatistic(statistic_table),
        )


class Connect(Enum):
    """
    Describes whether calculations should be done between samples of the same group (Connect.Intra)
    or samples of different groups (Connect.Between)
    """

    Intra = auto()
    Between = auto()


class CalculateMeanMinMax(Task[MeanMinMaxDistances]):
    """
    Arguments:
        distances: SequenceDistanceMatrix
        partition: Union[SpeciesPartition, GenusPartition]
        connection: Connect
        metric: Metric
        in_percent: bool = False
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.distances: Optional[SequenceDistanceMatrix] = None
        self.partition: Union[None, SpeciesPartition, GenusPartition] = None
        self.connection: Optional[Connect] = None
        self.metric: Optional[Metric] = None
        self.in_percent: bool = False

    def start(self) -> None:
        if self.distances is None:
            raise MissingArgument("distances")
        if self.partition is None:
            raise MissingArgument("partition")
        if self.connection is None:
            raise MissingArgument("connection")
        if self.metric is None:
            raise MissingArgument("metric")
        if isinstance(self.partition, GenusPartition):
            taxon_rank = TaxonRank.Genus
            partition = self.partition.get_dataframe()[["genus"]].copy()
        elif isinstance(self.partition, SpeciesPartition):
            taxon_rank = TaxonRank.Species
            partition = self.partition.get_dataframe().copy()
        else:
            assert False
        dataframe = self.distances.get_dataframe()[[self.metric]].copy()
        if self.in_percent:
            dataframe = dataframe * 100
        seqids = dataframe.index.get_level_values(0).unique()
        same_sample = pd.MultiIndex.from_arrays([seqids, seqids])
        dataframe.loc[same_sample] = float("nan")
        if self.connection is Connect.Intra:
            partition = partition[str(taxon_rank)]
            partition_target = partition.reindex(dataframe.index, level="seqid_target")
            partition_query = partition.reindex(dataframe.index, level="seqid_query")
            same_taxon = partition_target.loc[partition_target == partition_query]
            mean_distances = dataframe.groupby(same_taxon).mean()
            min_distances = dataframe.groupby(same_taxon).min()
            max_distances = dataframe.groupby(same_taxon).max()
            mean_distances.index.name = str(taxon_rank)
            min_distances.index.name = str(taxon_rank)
            max_distances.index.name = str(taxon_rank)
        elif self.connection is Connect.Between:
            dataframe = dataframe.join(
                partition.rename(columns=lambda s: s + "_target"), on="seqid_target"
            ).join(partition.rename(columns=lambda s: s + "_query"), on="seqid_query")
            between_index_names = [
                str(taxon_rank) + "_target",
                str(taxon_rank) + "_query",
            ]
            mean_distances = dataframe.groupby(between_index_names).mean()
            min_distances = dataframe.groupby(between_index_names).min()
            max_distances = dataframe.groupby(between_index_names).max()
        else:
            assert False
        self.result = MeanMinMaxDistances(
            mean_distances,
            min_distances,
            max_distances,
            metric=self.metric,
            taxon_rank=taxon_rank,
            in_percent=self.in_percent,
        )


@dataclass
class VersusAllOutput:
    sequence_summary_statistic: SimpleStatisticResult
    distances: List[SequenceDistanceOutput]
    mean_min_max_distances: List[MeanMinMaxDistances]
    summary_statistics: VersusAllSummary


class VersusAll(Task[VersusAllOutput]):
    """
    Arguments:
        sequences: SequenceData
        alignment: Alignment
        metrics: List[Metric]
        species: Optional[SpeciesPartition]
        vouchers: Optional[VoucherPartition]
        subspecies: Optional[SubspeciesPartition]
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.sequences: Optional[SequenceData] = None
        self.alignment: Optional[Alignment] = None
        self.metrics: List[Metric] = []
        self.species: Optional[SpeciesPartition] = None
        self.vouchers: Optional[VoucherPartition] = None
        self.subspecies: Optional[SubspeciesPartition] = None

    def start(self) -> None:
        if self.sequences is None:
            raise MissingArgument("sequences")
        if self.alignment is None:
            raise MissingArgument("alignment")
        if self.alignment is not Alignment.AlignmentFree and self.metrics is None:
            raise MissingArgument("metrics")

        distances_task = self.subtask(CalculateDistances)
        distances_task.sequences = self.sequences
        distances_task.alignment = self.alignment
        distances_task.metrics = self.metrics
        distances_task.start()
        assert distances_task.result is not None

        simple_statistic_task = self.subtask(CalculateSimpleStatistic)
        simple_statistic_task.sequences = self.sequences
        simple_statistic_task.species = self.species
        simple_statistic_task.start()
        assert simple_statistic_task.result is not None

        distances: List[SequenceDistanceOutput] = []
        for ordering in RowOrdering:
            for in_percent in [False, True]:
                for metric in self.metrics:
                    if self.species is None and ordering is RowOrdering.Species:
                        continue
                    output_distances_task = self.subtask(OutputSequenceDistances)
                    output_distances_task.sequence_matrix = distances_task.result
                    output_distances_task.metric = metric
                    output_distances_task.ordering = ordering
                    output_distances_task.species = self.species
                    output_distances_task.in_percent = in_percent
                    output_distances_task.start()
                    assert output_distances_task.result is not None
                    distances.append(output_distances_task.result)

        mean_min_max_distances: List[MeanMinMaxDistances] = []

        try:
            if self.species:
                genera = GenusPartition.from_species(self.species)
            else:
                genera = None
        except Exception:
            genera = None

        for partition in [self.species, genera]:
            if partition is None:
                continue
            for in_percent in [False, True]:
                for connection in Connect:
                    for metric in self.metrics:
                        mean_min_max_task = self.subtask(CalculateMeanMinMax)
                        mean_min_max_task.distances = distances_task.result
                        mean_min_max_task.partition = partition
                        mean_min_max_task.connection = connection
                        mean_min_max_task.metric = metric
                        mean_min_max_task.in_percent = in_percent
                        mean_min_max_task.start()
                        assert mean_min_max_task.result is not None
                        mean_min_max_distances.append(mean_min_max_task.result)

        summarize_arg: VersusAllSummarizeArg = VersusAllSummarizeArg(
            self.vouchers, genera, self.subspecies
        )
        summarize_task = self.subtask(VersusAllSummarize)
        summarize_task.distances = distances_task.result
        summarize_task.data = summarize_arg
        summarize_task.start()
        assert summarize_task.result is not None
        self.result = VersusAllOutput(
            simple_statistic_task.result,
            distances,
            mean_min_max_distances,
            summarize_task.result,
        )


class VersusReference(Task[Iterator[VersusReferenceSummary]]):
    """
    Arguments:
        data: CompleteData
        reference: CompleteData
        alignment: Alignment
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.data: Optional[CompleteData] = None
        self.reference: Optional[CompleteData] = None
        self.alignment: Optional[Alignment] = None

    def start(self) -> None:
        if self.data is None:
            raise MissingArgument("data")
        if self.reference is None:
            raise MissingArgument("reference")
        if self.alignment is None:
            raise MissingArgument("alignment")
        self.result = self._process()

    def _process(self) -> Iterator[VersusReferenceSummary]:
        assert self.data is not None
        assert self.reference is not None
        assert self.alignment is not None

        for chunk in self.data.get_chunks():
            distance_task = self.subtask(CalculateDistances)
            distance_task.alignment = self.alignment
            distance_task.metrics = list(Metric)
            distance_task.sequences = SequencesPair(
                target=self.reference.get_sequences(), query=chunk.get_sequences()
            )
            distance_task.start()
            assert distance_task.result is not None

            distance_table = distance_task.result.get_dataframe()
            min_distance_table = distance_table.loc[
                distance_table.groupby(level="seqid_query")[Metric.Uncorrected].idxmin()
            ]
            min_distance_table.reset_index(inplace=True)
            min_distance_table.rename(
                columns={
                    "seqid_query": SourcedColumn.query("seqid"),
                    "seqid_target": SourcedColumn.reference("seqid"),
                },
                inplace=True,
            )
            summary = min_distance_table.join(
                chunk.get_dataframe().drop(columns="sequence").rename(columns=SourcedColumn.query),
                on=SourcedColumn.query("seqid"),
            ).join(
                self.reference.get_dataframe().drop(columns="sequence").rename(columns=SourcedColumn.reference),
                on=SourcedColumn.reference("seqid"),
            )
            yield VersusReferenceSummary(summary)


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

    `similarity` describes alignment-free distance between sequences.
    The correspondence to similarity percentage was empirically determined to be the following:
        0.07:   98%
        0.10:   95%
        0.25:   90%
        0.31:   80%
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.similarity = 0.07
        self.length_threshold: Optional[int] = None
        self.keep_most_complete = False
        self.data: Optional[CompleteData] = None
        self._calculate_distances = self.subtask(CalculateDistances)
        self._calculate_distances.alignment = Alignment.AlignmentFree

    def start(self) -> None:
        if self.data is None:
            raise MissingArgument("data")
        if self.progress_handler:
            self._total_steps = self.data.sequence_count()
        self.result = self._dereplicate()

    def _dereplicate(self) -> Iterator[Dereplicated]:
        assert self.data is not None
        step = 0
        for chunk in self.data.get_chunks():
            assert chunk.dataframe is not None
            step += len(chunk.dataframe)
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

            self.progress(operation="Dereplicate", step=step)
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

    `similarity` describes alignment-free distance between sequences.
    The correspondence to similarity percentage was empirically determined to be the following:
        0.07:   98%
        0.10:   95%
        0.25:   90%
        0.31:   80%
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.similarity = 0.07
        self.alignment: Optional[Alignment] = None
        self.reference: Optional[SequenceData] = None
        self.data: Optional[CompleteData] = None
        self._calculate_distances = self.subtask(CalculateDistances)
        self._calculate_distances.metrics = [Metric.Uncorrected]

    def start(self) -> None:
        if self.alignment is None:
            raise MissingArgument("alignment")
        if self.reference is None:
            raise MissingArgument("reference")
        if self.data is None:
            raise MissingArgument("data")
        if self.progress_handler:
            self._total_steps = self.data.sequence_count()
        self.result = self._decontaminate()

    def _decontaminate(self) -> Iterator[Decontaminated]:
        assert self.alignment is not None
        assert self.reference is not None
        assert self.data is not None
        self._calculate_distances.alignment = self.alignment

        step = 0
        for chunk in self.data.get_chunks():
            assert chunk.dataframe is not None
            step += len(chunk.dataframe)
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

            self.progress(operation="Decontaminate", step=step)
            yield Decontaminated(
                contaminates=contaminates,
                decontaminated=chunk,
                summary=DecontaminateSummary(closest_table),
            )


@dataclass
class Decontaminated2:
    """
    The type of output of Decontaminate task
    """

    decontaminated: CompleteData
    contaminates: CompleteData
    summary: Optional[Decontaminate2Summary]


class Decontaminate2(Task[Iterator[Decontaminated2]]):
    """
    Arguments:
        alignment: Alignment
        ingroup: SequenceData
        outgroup: SequenceData
        data: CompleteData
        outgroup_weight: float = 1

    `outgroup_weight` desribes how strong the distance to the outgroup is considered.
    A sequence is considered a contiminant if
    distance_to_ingroup > outgroup_weight * distances_to_outgroup
    """

    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.alignment: Optional[Alignment] = None
        self.ingroup: Optional[SequenceData] = None
        self.outgroup: Optional[SequenceData] = None
        self.data: Optional[CompleteData] = None
        self.outgroup_weight: float = 1
        self._calculate_distances = self.subtask(CalculateDistances)
        self._calculate_distances.metrics = [Metric.Uncorrected]

    def start(self) -> None:
        if self.alignment is None:
            raise MissingArgument("alignment")
        if self.ingroup is None:
            raise MissingArgument("ingroup")
        if self.outgroup is None:
            raise MissingArgument("outgroup")
        if self.data is None:
            raise MissingArgument("data")
        if self.progress_handler:
            self._total_steps = self.data.sequence_count()
        self.result = self._decontaminate2()

    def _decontaminate2(self) -> Iterator[Decontaminated2]:
        assert self.alignment is not None
        assert self.ingroup is not None
        assert self.outgroup is not None
        assert self.data is not None
        self._calculate_distances.alignment = self.alignment

        step = 0
        for chunk in self.data.get_chunks():
            assert chunk.dataframe is not None
            step += len(chunk.dataframe)
            sequences = chunk.get_sequences()
            self._calculate_distances.sequences = SequencesPair(
                target=self.outgroup, query=sequences
            )
            self._calculate_distances.start()
            assert self._calculate_distances.result is not None
            distance_table_outgroup = (
                self._calculate_distances.result.get_dataframe().reset_index()
            )
            self._calculate_distances.sequences = SequencesPair(
                target=self.ingroup, query=sequences
            )
            self._calculate_distances.start()
            distance_table_ingroup = (
                self._calculate_distances.result.get_dataframe().reset_index()
            )

            distances_to_ingroup = distance_table_ingroup.loc[
                distance_table_ingroup.groupby("seqid_query")[
                    Metric.Uncorrected
                ].idxmin()
            ].set_index("seqid_query")
            distances_to_outgroup = distance_table_outgroup.loc[
                distance_table_outgroup.groupby("seqid_query")[
                    Metric.Uncorrected
                ].idxmin()
            ].set_index("seqid_query")

            summary = distances_to_ingroup.rename(
                columns={
                    "seqid_target": "closest ingroup",
                    Metric.Uncorrected: "distance to ingroup",
                }
            ).join(
                distances_to_outgroup.rename(
                    columns={
                        "seqid_target": "closest outgroup",
                        Metric.Uncorrected: "distance to outgroup",
                    }
                ),
                how="inner",
            )

            summary["is contaminant"] = (
                summary["distance to ingroup"]
                > self.outgroup_weight * summary["distance to outgroup"]
            )

            assert sequences.dataframe is not None
            sequences.dataframe = sequences.dataframe.loc[~summary["is contaminant"]]

            contaminates = chunk.split_sequences(sequences)

            self.progress(operation="Decontaminate", step=step)
            yield Decontaminated2(
                contaminates=contaminates,
                decontaminated=chunk,
                summary=Decontaminate2Summary(summary),
            )
