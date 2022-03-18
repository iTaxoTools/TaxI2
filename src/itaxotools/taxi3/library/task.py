#!/usr/bin/env python3

from typing import Callable, Generic, TypeVar, Optional, List, Iterator, Union
from abc import ABC, abstractmethod
from enum import Enum, auto
from dataclasses import dataclass

import pandas as pd

from .datatypes import SequenceDistanceMatrix, SequenceData, Metric, CompleteData
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
    included: CompleteData
    excluded: CompleteData


class Dereplicate(Task[Iterator[CompleteData]]):
    def __init__(self, warn: WarningHandler):
        super().__init__(warn)
        self.similarity = 0.07
        self.length_threshold: Optional[int] = None
        self.keep_most_complete = False
