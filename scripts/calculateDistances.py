from itaxotools.taxi3.distances import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.align import *
from itaxotools.taxi3.pairs import *
from pathlib import Path
from sys import argv
from time import perf_counter


def calc(aligned_pairs, metric):
    for x, y in aligned_pairs:
        yield metric.calculate(x, y)


def make_pairs(x, y):
    for a in x.read(idHeader='seqid', seqHeader='sequence'):
        for b in y.read(idHeader='seqid', seqHeader='sequence'):
            yield SequencePair(a, b)


path_data = Path(argv[1])
path_reference = Path(argv[2])
path_out = Path(argv[3])

ts = perf_counter()

data = SequenceFile.Tabfile(path_data)
refrence = SequenceFile.Tabfile(path_reference)

# pairs = SequencePairs.fromProduct(data, refrence)
pairs = make_pairs(data, refrence)

aligner = PairwiseAligner.Biopython()
aligned_pairs = aligner.align_pairs(pairs)

metric = DistanceMetric.Kimura2P()
distances = calc(aligned_pairs, metric)

outFile = DistanceFile.Linear(path_out)
outFile.write(distances)


tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
