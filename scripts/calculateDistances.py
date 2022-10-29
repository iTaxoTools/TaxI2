from itaxotools.taxi3.distances import *
from itaxotools.taxi3.sequences import *
from pathlib import Path
from sys import argv
from time import perf_counter


def calc(data, refrence, metric):
    for x in data.read(idHeader='seqid', seqHeader='sequence'):
        for y in refrence.read(idHeader='seqid', seqHeader='sequence'):
            yield metric.calculate(x, y)


path_data = Path(argv[1])
path_reference = Path(argv[2])
path_out = Path(argv[3])

ts = perf_counter()

data = SequenceFile.Tabfile(path_data)
refrence = SequenceFile.Tabfile(path_reference)

metric = DistanceMetric.Uncorrected()
distances = calc(data, refrence, metric)
outFile = DistanceFile.Linear(path_out)

outFile.write(distances)


tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
