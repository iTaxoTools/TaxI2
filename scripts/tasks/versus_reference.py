from pathlib import Path
from sys import argv
from time import perf_counter

from itaxotools.taxi3.library.datatypes import (
    CompleteData, TabfileReader, ValidFilePath)
from itaxotools.taxi3.library.task import Alignment, VersusReference

# Time taken for 50 vs 50 sample input: ~ 47s

path_data = Path(argv[1])
path_reference = Path(argv[2])
path_out = Path(argv[3])

ts = perf_counter()

data = CompleteData.from_path(ValidFilePath(path_data), TabfileReader)
reference = CompleteData.from_path(ValidFilePath(path_reference), TabfileReader)

task = VersusReference(warn=print)
task.alignment = Alignment.Pairwise
task.data = data
task.reference = reference

task.start()

for output in task.result:
    output.append_to_file(path_out)

tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
