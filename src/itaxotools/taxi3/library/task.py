#!/usr/bin/env python3

from typing import Callable, Generic, TypeVar, Optional, List, Iterator, Union
from abc import ABC, abstractmethod
from enum import Enum, auto
from dataclasses import dataclass

import pandas as pd
import networkx as nx

from .datatypes import (
    SequenceDistanceMatrix,
    SequenceData,
    Metric,
    CompleteData,
    DecontaminateSummary,
)
from .rust_backend import calc, make_aligner

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

            indices_closest = (
                distance_table[["seqid (query 1)", Metric.Uncorrected]]
                .groupby("seqid (query 1)")
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
            closest_table.rename(columns={Metric.Uncorrected, "distance"}, inplace=True)

            assert sequences.dataframe is not None
            sequences.dataframe = sequences.dataframe.loc[decontaminates_seqids]

            contaminates = chunk.split_sequences(sequences)
            yield Decontaminated(
                contaminates=contaminates,
                decontaminated=chunk,
                summary=DecontaminateSummary(closest_table),
            )
