#!/usr/bin/env python3

from typing import Callable, Generic, TypeVar, Optional, List
from abc import ABC, abstractmethod
from enum import Enum, auto

from .datatypes import SequenceDistanceMatrix, SequenceData, Metric

_Result = TypeVar("_Result")


class Alignment(Enum):
    Pairwise = auto()
    AlignmentFree = auto()
    AlreadyAligned = auto()


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
