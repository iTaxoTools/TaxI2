from itaxotools.taxi3.distances import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.align import *
from itaxotools.taxi3.pairs import *
from pathlib import Path
from sys import argv
from time import perf_counter
import time


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

def progress(distances, total):

    for index, distance in enumerate(distances):
        print(f"\r Loading... {index}/{total} = {(round(float(index)/float(total)  * 100, 2))}", end="")
        yield distance


path_data = Path(argv[1])
path_reference = Path(argv[2])
path_out_linear = Path(argv[3])
path_out_matrixs = Path(argv[4])
path_out_pairs = Path(argv[5])

ts = perf_counter()

file_data = SequenceFile.Tabfile(path_data)
file_reference = SequenceFile.Tabfile(path_reference)

data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')
reference = Sequences.fromFile(file_reference, idHeader='seqid', seqHeader='sequence')

total = len(data) * len(reference) * 4

data = data.normalize()
reference = reference.normalize()

pairs = SequencePairs.fromProduct(data, reference)

aligner = PairwiseAligner.Biopython()
aligned_pairs = aligner.align_pairs_parallel(pairs)

outFile = SequencePairFile.Formatted(path_out_pairs)
aligned_pairs = outFile.iter_write(aligned_pairs)

distances = calc(aligned_pairs)

# outFile = DistanceFile.Matrix(path_out_matrixs)
# distances = outFile.iter_write(distances)

outFile = DistanceFile.LinearWithExtras(path_out_linear)
distances = outFile.iter_write(distances)

distances = progress(distances, total)

for _ in distances:
    pass

tf = perf_counter()

print(f'Time taken: {tf-ts:.4f}s')
