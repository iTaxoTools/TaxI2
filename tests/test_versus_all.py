#!/usr/bin/env python3

from pathlib import Path

from itaxotools.taxi3.library.datatypes import (
    SequenceData,
    TabfileReader,
    ValidFilePath,
    Metric,
)
from itaxotools.taxi3.library.task import (
    CalculateDistances,
    Alignment,
    VersusAllSummarize,
    VersusAllSummarizeArg,
)


def test_versus_all() -> None:
    input_path = ValidFilePath(Path(__file__).with_name("Scaphio_input_small.txt"))
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task = CalculateDistances(print)
    task.sequences = sequences
    task.alignment = Alignment.Pairwise
    task.metrics = list(Metric)
    task.start()
    print(task.result.get_dataframe())


def test_summary() -> None:
    input_path = ValidFilePath(Path(__file__).with_name("Scaphio_input_small.txt"))
    tested_path = Path(__file__).with_name("Summary_statistics.txt")
    output_path = (
        Path(__file__).parent / "temp_test_files" / "Summary_statistics_test.txt"
    )
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task_distances = CalculateDistances(print)
    task_distances.sequences = sequences
    task_distances.alignment = Alignment.Pairwise
    task_distances.metrics = list(Metric)
    task_distances.start()
    task = VersusAllSummarize(print)
    task.data = VersusAllSummarizeArg.from_path(input_path, TabfileReader())
    assert task.data.vouchers is not None
    assert task.data.species is not None
    assert task.data.subspecies is None
    task.distances = task_distances.result
    task.start()
    task.result.to_file(output_path)
    assert tested_path.read_text() == output_path.read_text()
    output_path.unlink()
