from itertools import groupby
from pathlib import Path
from sys import argv
from time import perf_counter
import os
from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *


def calc(aligned_pairs, metric=DistanceMetric.Uncorrected()):
    for x, y in aligned_pairs:
        if x.id == y.id:
            continue
        yield metric.calculate(x, y)


def progress(distances, total):

    for index, distance in enumerate(distances):
        index+=1
        print(f"\r Loading... {int(index)}/{total} = {(round(float(index)/float(total)  * 100, 2))} ", end="")
        yield distance


def decontaminateSeq(disctances, contaminatePath, decontaminatePath, summaryPath, similarity=0.09):

    decontFile = SequenceFile.Tabfile(decontaminatePath)
    contFile = SequenceFile.Tabfile(contaminatePath)

    with decontFile.open('w') as decontHandler, contFile.open('w') as contHandler, open(summaryPath, 'w') as summryFile:

        summryFile.write("seqid_query\tclosest possible contaminant\tdistance\tis_contaminant")

        for p, disctance in disctances:
            isContaminate = disctance.d <= similarity

            if isContaminate:
                decontHandler.write(p.x)
            else:
                contHandler.write(p.x)

            summryFile.write(f'\n{disctance.x.id}\t{disctance.y.id}\t{disctance.d}\t{isContaminate}')

            yield disctance


def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d


def multiply(g, n):
    return (x for x in g for i in range(n))


def normalizePairs(pairs):
    for pair in pairs:
        yield SequencePair(pair.x.normalize(), pair.y.normalize())


def main():

    path_data = Path(argv[1])
    path_reference = Path(argv[2])
    contaminatePath = Path(argv[3])
    decontaminatePath = Path(argv[4])
    summaryPath = Path(argv[5])

    ts = perf_counter()

    file_data = SequenceFile.Tabfile(path_data)
    file_reference = SequenceFile.Tabfile(path_reference)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')
    reference = Sequences.fromFile(file_reference, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    pairs = SequencePairs.fromProduct(data, reference)

    pairs = multiply(pairs, 2)

    normalizePair = normalizePairs(pairs)

    aligner = PairwiseAligner.Biopython()
    aligned_pairs = aligner.align_pairs(normalizePair)

    distances = calc(aligned_pairs)

    minimums = get_minimum(distances)

    allPairs = zip(pairs, minimums)

    d = decontaminateSeq(allPairs, contaminatePath, decontaminatePath, summaryPath)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
