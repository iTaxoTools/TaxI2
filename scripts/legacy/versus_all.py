from pathlib import Path
from sys import argv
from time import perf_counter

from itaxotools.taxi3.library.config import AlignmentScores, Config
from itaxotools.taxi3.library.datatypes import (
    Metric, SequenceData, SpeciesPartition, SubspeciesPartition,
    TabfileReader, ValidFilePath, VoucherPartition)
from itaxotools.taxi3.library.task import (
    Alignment, CalculateDistances, VersusAll)

input = Path(argv[1])


metrics = [Metric.Uncorrected]

print('Calculating distances...')
task = CalculateDistances(warn=print)
task.sequences = SequenceData.from_path(ValidFilePath(input), TabfileReader)
task.alignment = Alignment.Pairwise
task.metrics = metrics
# task.config = config
task.start()

distances = task.result

print('Running Versus All analysis...')
task = VersusAll(warn=print)

data = TabfileReader.read_data(ValidFilePath(input))
for table in data:
    if isinstance(table, SequenceData):
        task.sequences = table
    elif isinstance(table, VoucherPartition):
        task.vouchers = table
    elif isinstance(table, SpeciesPartition):
        task.species = table
    elif isinstance(table, SubspeciesPartition):
        task.subspecies = table

task.distances = distances
task.alignment = Alignment.Pairwise
task.metrics = metrics
# task.config = config
task.start()

print('Printing results...')
tables = task.result
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
