from pathlib import Path
from sys import argv
from time import perf_counter

from itaxotools.taxi3.library.datatypes import (
    Metric, SequenceData, TabfileReader, ValidFilePath)
from itaxotools.taxi3.library.task import (
    Alignment, CalculateDistances, SequencesPair)

# Time taken for 50 vs 50 sample input: ~ 45s

path_data = Path(argv[1])
path_reference = Path(argv[2])
path_out = Path(argv[3])

def foo(obj):
    print(f"\r Loading... {obj.current_step}/{obj.total_steps} = {(round(float(obj.current_step) / float(obj.total_steps) * 100, 2))}", end="")

ts = perf_counter()

data = SequenceData.from_path(ValidFilePath(path_data), TabfileReader)
reference = SequenceData.from_path(ValidFilePath(path_reference), TabfileReader)

pairs = SequencesPair(target = data, query = reference)

task = CalculateDistances(warn=print)
task.sequences = pairs
task.progress_handler = foo
# task.alignment = Alignment.AlreadyAligned
task.alignment = Alignment.Pairwise
task.metrics = [
    Metric.Uncorrected,
    Metric.JukesCantor,
    Metric.Kimura2P,
    Metric.UncorrectedWithGaps,
]

task.start()

task.result.distance_matrix.to_csv(path_out, sep='\t')

tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
