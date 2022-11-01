from itaxotools.taxi3.distances import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.align import *
from itaxotools.taxi3.pairs import *
from pathlib import Path
from sys import argv
from time import perf_counter
from itertools import groupby


def calc(aligned_pairs):
    metrics = [
        DistanceMetric.Uncorrected(),
        DistanceMetric.UncorrectedWithGaps(),
        DistanceMetric.JukesCantor(),
        DistanceMetric.Kimura2P(),
    ]
    for x, y in aligned_pairs:
        for metric in metrics:
            yield metric.calculate(x, y)

def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.idx):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d


path_data = Path(argv[1])
path_reference = Path(argv[2])
path_out = Path(argv[3])

ts = perf_counter()

file_data = SequenceFile.Tabfile(path_data)
file_reference = SequenceFile.Tabfile(path_reference)

data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')
reference = Sequences.fromFile(file_reference, idHeader='seqid', seqHeader='sequence')

data = data.normalize()
reference = reference.normalize()

pairs = SequencePairs.fromProduct(data, reference)

aligner = PairwiseAligner.Biopython()
aligned_pairs = aligner.align_pairs(pairs)

distances = calc(aligned_pairs)

minimums = get_minimum(distances)

outFile = DistanceFile.Linear(path_out)
outFile.write(minimums)


tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
