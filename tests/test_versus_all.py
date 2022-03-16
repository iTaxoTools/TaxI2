#!/usr/bin/env python3

from pathlib import Path

from itaxotools.taxi3.library.datatypes import (
    SequenceData,
    TabfileReader,
    ValidFilePath,
    Metric,
)
from itaxotools.taxi3.library.task import CalculateDistances, Alignment


def test_versus_all():
    input_path = ValidFilePath(Path(__file__).with_name("Scaphio_input_small.txt"))
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task = CalculateDistances(print)
    task.sequences = sequences
    task.alignment = Alignment.Pairwise
    task.metrics = list(Metric)
    task.start()
    print(task.result.get_dataframe())
