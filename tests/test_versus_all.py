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
    VoucherPartition,
    SubspeciesPartition,
    GenusPartition,
    TaxonRank,
    MeanMinMaxFileFormat,
    MeanMinMaxDistances,
)
from itaxotools.taxi3.library.task import (
    CalculateDistances,
    Alignment,
    VersusAllSummarize,
    VersusAllSummarizeArg,
    OutputSequenceDistances,
    CalculateSimpleStatistic,
    Connect,
    CalculateMeanMinMax,
    VersusAll,
)

TEST_DATA_DIR = Path(__file__).parent / "versus_all_data"
TMP_TEST_DIR = Path(__file__).parent / "temp_test_files"


@pytest.mark.legacy
def test_versus_all_distances() -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task = CalculateDistances(print)
    task.sequences = sequences
    task.alignment = Alignment.Pairwise
    task.metrics = list(Metric)
    task.start()
    print(task.result.get_dataframe())


@pytest.mark.legacy
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


@pytest.mark.legacy
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


@pytest.mark.legacy
def test_simple_statistics() -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    tested_path = TEST_DATA_DIR / "Sequence summary statistics.txt"
    output_path = TMP_TEST_DIR / "Sequence summary statistics.txt"
    sequences = SequenceData.from_path(input_path, TabfileReader())
    species = SpeciesPartition.from_path(input_path, TabfileReader())
    task = CalculateSimpleStatistic(print)
    task.sequences = sequences
    task.species = species
    task.start()
    with open(
        output_path, mode="w"
    ):  # Erase data that might have been left from a previous run
        pass
    task.result.total.append_to_file(output_path)
    with open(output_path, mode="a") as output_file:
        print("", file=output_file)
    task.result.by_species.append_to_file(output_path)
    assert tested_path.read_text().split("\n") == output_path.read_text().split("\n")
    output_path.unlink()


@pytest.mark.legacy
@pytest.mark.parametrize(
    "format, taxon_rank, connection, in_percent",
    [
        (format, TaxonRank.Species, Connect.Between, in_percent)
        for format in MeanMinMaxFileFormat
        for in_percent in [False, True]
    ]
    + [
        (MeanMinMaxFileFormat.MeanMinMax, TaxonRank.Genus, Connect.Between, False),
        (MeanMinMaxFileFormat.MeanMinMax, TaxonRank.Genus, Connect.Between, True),
        (MeanMinMaxFileFormat.MeanMinMax, TaxonRank.Species, Connect.Intra, False),
        (MeanMinMaxFileFormat.MeanMinMax, TaxonRank.Genus, Connect.Intra, False),
    ],
)


@pytest.mark.legacy
def test_mean_min_max(
    format: MeanMinMaxFileFormat,
    taxon_rank: TaxonRank,
    connection: Connect,
    in_percent: bool,
) -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    sequences = SequenceData.from_path(input_path, TabfileReader())
    task_distances = CalculateDistances(print)
    task_distances.sequences = sequences
    task_distances.alignment = Alignment.Pairwise
    task_distances.metrics = [Metric.Uncorrected]
    task_distances.start()
    if taxon_rank is TaxonRank.Species:
        taxons = SpeciesPartition.from_path(input_path, TabfileReader())
    else:
        taxons = GenusPartition.from_path(input_path, TabfileReader())
    task_mean_min_max = CalculateMeanMinMax(print)
    task_mean_min_max.distances = task_distances.result
    task_mean_min_max.partition = taxons
    task_mean_min_max.connection = connection
    task_mean_min_max.metric = Metric.Uncorrected
    task_mean_min_max.in_percent = in_percent
    task_mean_min_max.start()
    table = task_mean_min_max.result
    assert table is not None
    tested_path = TEST_DATA_DIR / table.make_file_name(format)
    output_path = TMP_TEST_DIR / table.make_file_name(format)
    with open(output_path, mode="w") as output_file:
        print(table.description(format), file=output_file)
    table.append_to_file(output_path, format)
    assert tested_path.read_text().split("\n") == output_path.read_text().split("\n")


@pytest.mark.legacy
def test_versus_all_mode() -> None:
    input_path = ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt")
    data = TabfileReader.read_data(input_path)
    versus_all_task = VersusAll(warn=print)
    for table in data:
        if isinstance(table, SequenceData):
            versus_all_task.sequences = table
        elif isinstance(table, VoucherPartition):
            versus_all_task.vouchers = table
        elif isinstance(table, SpeciesPartition):
            versus_all_task.species = table
        elif isinstance(table, SubspeciesPartition):
            versus_all_task.subspecies = table
    assert versus_all_task.sequences is not None
    versus_all_task.alignment = Alignment.Pairwise
    versus_all_task.metrics = [Metric.Uncorrected]
    versus_all_task.start()
    assert versus_all_task.result is not None
    tables = versus_all_task.result
    for table in (
        [
            tables.sequence_summary_statistic.total,
            tables.sequence_summary_statistic.by_species,
        ]
        + tables.distances
        + tables.mean_min_max_distances
        + [tables.summary_statistics]
    ):
        if table:
            print(table.get_dataframe().to_string())
