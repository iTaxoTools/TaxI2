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


def excludeReplicate(groupSimilar, includeSet, excludedSet):
    for g in groupSimilar:
        first = next(g)
        maxSeq = (first[0].x)
        maxlen = len(maxSeq.seq)
        for p, similar in chain([first], g):
            if similar:
                maxlen = max(maxlen, len(p.y.seq))
                if maxlen == len(p.y.seq):
                    if maxSeq.id in includeSet:
                        includeSet.remove(maxSeq.id)
                    excludedSet.add(maxSeq.id)
                    maxSeq = p.y
        yield first
        includeSet.add(maxSeq.id)


def dereplicate(generatorObject, excludedSet, dereplicatedPath, excludedPath, headers):
    with open(excludedPath, 'w') as excludeFile, open(dereplicatedPath, 'w') as derepFile:
        excludeFile.write('\t'.join(headers))
        derepFile.write('\t'.join(headers))
        for p, _ in generatorObject:
            extras = [v for v in p.x.extras.values()]
            extrasString = '\t'.join(extras)
            if p.x.id in excludedSet:
                excludeFile.write(f"\n{p.x.id}\t{extrasString}\t{p.x.seq}")
            else:
                derepFile.write(f"\n{p.x.id}\t{extrasString}\t{p.x.seq}")

            yield p


def main():
    global graph, vertices_no
    path_data = Path(argv[1])
    dereplicatedPath = 'dereplicated.txt'
    excludedPath = 'excluded.txt'
    lenTrashold = 10
    similarityThreshold = 0.07
    ts = perf_counter()

    excludedSet = set()
    includeSet = set()

    file_data = SequenceFile.Tabfile(path_data)

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

    isSimilar = getSimilar(distances, similarityThreshold)

    allPairs = zip(pairs, isSimilar)

    groupSimilar = groupSimilars(allPairs)

    pairData = excludeReplicate(groupSimilar, includeSet, excludedSet)

    headers = file_data.getHeader()

    d = dereplicate(pairData, excludedSet, dereplicatedPath, excludedPath, headers)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
