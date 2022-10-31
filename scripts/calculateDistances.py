from itaxotools.taxi3.distances import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.align import *
from itaxotools.taxi3.pairs import *
from pathlib import Path
from sys import argv
from time import perf_counter


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


def make_pairs(a, b):
    yield from (SequencePair(x, y) for x in a for y in b)


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

# pairs = SequencePairs.fromProduct(data, reference)
pairs = make_pairs(data, reference)

aligner = PairwiseAligner.Biopython()
aligned_pairs = aligner.align_pairs(pairs)

distances = calc(aligned_pairs)

outFile = DistanceFile.Linear(path_out)
outFile.write(distances)


tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
