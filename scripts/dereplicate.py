from itertools import groupby, chain
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
        yield metric.calculate(x, y)


def progress(distances, total):

    for index, distance in enumerate(distances):
        print(f"\r Processing... {int(index) + 1}/{total} = {(round(float(index)/float(total)  * 100 + 10, 2))} ", end="")
        yield distance


def dropShortest(data, lenTrashold):

    for d in data:
        if len(d.seq) >= lenTrashold:
            yield d


def getNotExcludedPair(pairs, excludedSet):
    for pair in pairs:
        if pair.x.id not in excludedSet and pair.y.id not in excludedSet:
            yield pair


def normalizePairs(pairs):
    for pair in pairs:
        yield SequencePair(pair.x.normalize(), pair.y.normalize())


def derepSeq(disctances, dereplicatedPath, excludedPath, headers, similarity=0.07):
    excludedSet = set()

    with open(excludedPath, 'w') as excludeFile, open(dereplicatedPath, 'w') as derepFile, open('summaryDrep.txt', 'w') as summryFile:

        excludeFile.write('\t'.join(headers))
        derepFile.write('\t'.join(headers))
        summryFile.write("seqid_query\tclosest possible calc\tdistance\tis_replicated")

        for disctance in disctances:
            isReplicated = disctance.d <= similarity
            extras = [v for v in disctance.x.extras.values()]
            extrasString = '\t'.join(extras)

            if isReplicated and disctance.y.id not in excludedSet:
                excludeFile.write(f"\n{disctance.x.id}\t{extrasString}\t{disctance.x.seq}")
                excludedSet.add(disctance.x.id)
            else:
                derepFile.write(f"\n{disctance.x.id}\t{extrasString}\t{disctance.x.seq}")

            summryFile.write(f'\n{disctance.x.id}\t{disctance.y.id}\t{disctance.d}\t{isReplicated}')

            yield disctance


def getSimilar(distances, similarityThreshold):

    for distance in distances:
        if distance.d <= similarityThreshold:
            yield True
        else:
            yield False


def groupSimilars(pairs):

    for k, g in groupby(pairs, lambda p: p[0].x.id):
        yield g


def multiply(g, n):
    return (x for x in g for i in range(n))


def removeIdentical(pairs):

    for pair in pairs:
        if pair.x.id != pair.y.id:
            yield pair


def excludeReplicate(groupSimilar, includeSet, excludedSet, summaryPath):
    with open(summaryPath, 'w') as summaryFile:
        summaryFile.write('seqid_query\tis_replicant\tseqid_replicant\tdistance\tlength\n')
        for g in groupSimilar:
            first = next(g)
            maxSeq = (first[0].x)
            maxlen = len(maxSeq.seq)
            querySequence = first[0].x
            for p, similar, distance in chain([first], g):
                if similar:
                    maxlen = max(maxlen, len(p.y.seq))
                    if maxlen == len(p.y.seq):
                        # we are excluding this
                        summaryFile.write(f'{maxSeq.id}\t{True}\t{p.y.id}\t{distance.d:.4f}\t{len(maxSeq.seq)}\n')
                        excludedSet.add(maxSeq.id)
                        maxSeq = p.y
            if querySequence.id not in excludedSet:
                # we are including this
                summaryFile.write(f'{querySequence.id}\t{False}\t{None}\t{None}\t{len(querySequence.seq)}\n')
            yield querySequence


def dereplicate(sequences, excludedSet, dereplicatedPath, excludedPath):
    derepFile = SequenceFile.Tabfile(dereplicatedPath)
    excFile = SequenceFile.Tabfile(excludedPath)
    with (
        derepFile.open('w') as derepHandler,
        excFile.open('w') as exclHandler,
    ):
        for sequence in sequences:

            if sequence.id in excludedSet:
                exclHandler.write(sequence)
            else:
                derepHandler.write(sequence)
            yield sequence


def main():
    dataPath = Path(argv[1])
    dereplicatedPath = Path(argv[2])
    excludedPath = Path(argv[3])
    summaryPath = Path(argv[4])

    lenTrashold = 10
    similarityThreshold = 0.07
    ts = perf_counter()

    excludedSet = set()
    includeSet = set()

    file_data = SequenceFile.Tabfile(dataPath)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = Sequences(lambda data, lenTrashold:dropShortest(data, lenTrashold), data, lenTrashold)

    pairs = SequencePairs.fromProduct(data, data)

    pairs = getNotExcludedPair(pairs, excludedSet)

    pairs = removeIdentical(pairs)

    pairs = multiply(pairs, 2)

    normalizePair = normalizePairs(pairs)

    aligner = PairwiseAligner.Biopython()

    aligned_pairs = aligner.align_pairs(normalizePair)

    distances = calc(aligned_pairs)
    distances = multiply(distances, 2)

    isSimilar = getSimilar(distances, similarityThreshold)

    allPairs = zip(pairs, isSimilar, distances)

    groupSimilar = groupSimilars(allPairs)

    sequences = excludeReplicate(groupSimilar, includeSet, excludedSet, summaryPath)

    sequences = dereplicate(sequences, excludedSet, dereplicatedPath, excludedPath)

    sequences = progress(sequences, total)

    for _ in sequences:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
