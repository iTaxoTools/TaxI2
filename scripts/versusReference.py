from itertools import groupby
from pathlib import Path
from sys import argv
from time import perf_counter

from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *


def calc(aligned_pairs, metric=DistanceMetric.Uncorrected()):
    for x, y in aligned_pairs:
        yield metric.calculate(x, y)

def calAll(disctances):
    metrics = [
        DistanceMetric.UncorrectedWithGaps(),
        DistanceMetric.JukesCantor(),
        DistanceMetric.Kimura2P(),
    ]
    for disctance in disctances:
        yield disctance
        for metric in metrics:
            yield metric.calculate(disctance.x, disctance.y)


def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d


def progress(distances, total):

    for index, distance in enumerate(distances):
        print(f"\r Loading... {index}/{total} = {(round(float(index)/float(total)  * 100, 2))}", end="")
        yield distance


def main():
    path_data = Path(argv[1])
    path_reference = Path(argv[2])
    path_out = Path(argv[3])

    ts = perf_counter()

    file_data = SequenceFile.Tabfile(path_data)
    file_reference = SequenceFile.Tabfile(path_reference)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')
    reference = Sequences.fromFile(file_reference, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = data.normalize()
    reference = reference.normalize()

    pairs = SequencePairs.fromProduct(data, reference)

    aligner = PairwiseAligner.Biopython()
    aligned_pairs = aligner.align_pairs(pairs)

    distances = calc(aligned_pairs)

    minimums = get_minimum(distances)

    alldistance = calAll(minimums)

    distances = progress(alldistance, total)

    outFile = DistanceFile.LinearWithExtras(path_out)
    outFile.write(distances)

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
