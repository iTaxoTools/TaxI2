#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import pytest

from itaxotools.taxi3.library.datatypes import (
    SequenceData,
    TabfileReader,
    ValidFilePath,
    Metric,
    RowOrdering,
    SpeciesPartition,
)
from itaxotools.taxi3.library.task import (
    CalculateDistances,
    Alignment,
    VersusAllSummarize,
    VersusAllSummarizeArg,
    OutputSequenceDistances,
)

TEST_DATA_DIR = Path(__file__).parent / "versus_all_data"
TMP_TEST_DIR = Path(__file__).parent / "temp_test_files"


def test_versus_all() -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task = CalculateDistances(print)
    task.sequences = sequences
    task.alignment = Alignment.Pairwise
    task.metrics = list(Metric)
    task.start()
    print(task.result.get_dataframe())


def test_summary() -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    tested_path = TEST_DATA_DIR / "Summary_statistics.txt"
    output_path = TMP_TEST_DIR / "Summary_statistics_test.txt"
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
    assert tested_path.read_text().split("\n") == output_path.read_text().split("\n")
    output_path.unlink()


@pytest.mark.parametrize("ordering", list(RowOrdering))
@pytest.mark.parametrize("in_percent", [False, True])
def test_sequence_distance_table_output(
    ordering: RowOrdering, in_percent: bool
) -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task_distances = CalculateDistances(print)
    task_distances.sequences = sequences
    task_distances.alignment = Alignment.Pairwise
    task_distances.metrics = [Metric.Uncorrected]
    task_distances.start()
    if ordering is RowOrdering.Species:
        species: Optional[SpeciesPartition] = SpeciesPartition.from_path(
            input_path, TabfileReader()
        )
    else:
        species = None
    task_table_output = OutputSequenceDistances(print)
    task_table_output.sequence_matrix = task_distances.result
    task_table_output.metric = Metric.Uncorrected
    task_table_output.ordering = ordering
    task_table_output.species = species
    task_table_output.in_percent = in_percent
    task_table_output.start()
    table = task_table_output.result
    assert table is not None
    tested_path = TEST_DATA_DIR / table.make_file_name()
    output_path = TMP_TEST_DIR / table.make_file_name()
    with open(output_path, mode="w") as output_file:
        print(table.description(), file=output_file)
    table.append_to_file(output_path)
    assert tested_path.read_text().split("\n") == output_path.read_text().split("\n")
    output_path.unlink()
