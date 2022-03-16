#!/usr/bin/env python3

from typing import Callable, Generic, TypeVar, Optional, List
from abc import ABC, abstractmethod
from enum import Enum, auto

import pandas as pd

from .datatypes import SequenceDistanceMatrix, SequenceData, Metric
from .rust_backend import calc, make_aligner

from .alfpy_distance import alfpy_distance_array, alfpy_distance_array2

_Result = TypeVar("_Result")


class Alignment(Enum):
    Pairwise = auto()
    AlignmentFree = auto()
    AlreadyAligned = auto()


class MissingArgument(Exception):
    pass


class Task(ABC, Generic[_Result]):
    def __init__(self, warn: Callable[[Warning], None]):
        self.warn = warn
        self.result: Optional[_Result] = None

    @abstractmethod
    def start(self) -> None:
        pass


class VersusAllComparison(Task[SequenceDistanceMatrix]):
    def __init__(self, warn: Callable[[Warning], None]):
        super().__init__(warn)
        self.sequences: Optional[SequenceData] = None
        self.alignment: Optional[Alignment] = None
        self.metrics: List[Metric] = []

    def _alignment_free_start(self) -> None:
        assert self.sequences is not None
        sequences = self.sequences.get_dataframe()
        distances_array = alfpy_distance_array(sequences["sequences"])
        seqids = sequences.index
        index = pd.MultiIndex.from_product(
            [seqids.copy(), seqids.copy()], names=["seqid1", "seqid2"]
        )
        self.result = SequenceDistanceMatrix(
            pd.DataFrame(distances_array, index=index, columns=[Metric.Uncorrected])
        )

    def _alignment_start(self) -> None:
        assert self.alignment is not None
        assert self.sequences is not None
        assert self.metrics
        sequences = self.sequences.get_dataframe()
        if self.alignment is Alignment.AlreadyAligned:
            distances_array = calc.make_distance_array_aligned(
                sequences["sequence"], sequences["sequence"]
            )
        elif self.alignment is Alignment.Pairwise:
            aligner = make_aligner()
            distances_array = calc.make_distance_array(
                aligner, sequences["sequence"], sequences["sequence"]
            )
        else:
            assert self.alignment is not Alignment.AlignmentFree
        seqids = sequences.index
        index = pd.MultiIndex.from_product(
            [seqids.copy(), seqids.copy()], names=["seqid1", "seqid2"]
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
