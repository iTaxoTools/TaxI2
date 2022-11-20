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
        yield metric.calculate(x, y)


def progress(distances, total):

    for index, distance in enumerate(distances):
        index+=1
        print(f"\r Loading... {int(index)}/{total} = {(round(float(index)/float(total)  * 100, 2))} ", end="")
        yield distance


def decontaminateSeq(allpairs, contaminatePath, decontaminatePath, summaryPath, outgroup_weight=1.0):

    decontFile = SequenceFile.Tabfile(decontaminatePath)
    contFile = SequenceFile.Tabfile(contaminatePath)

    with decontFile.open('w') as decontHandler, contFile.open('w') as contHandler, open(summaryPath, 'w') as summryFile:
        summryFile.write("seqid_query\tclosest possible inGroup id\tingroup distance\tclosest possible outGroup id\toutgroup distance\tis_contaminant")

        for inputData, ingroup_disctance, outgroup_disctance in allpairs:
            isContaminate = ingroup_disctance.d > (outgroup_weight * outgroup_disctance.d)

            if isContaminate:
                decontHandler.write(inputData)
            else:
                contHandler.write(inputData)

            summryFile.write(f'\n{ingroup_disctance.x.id}\t{ingroup_disctance.y.id}\t{ingroup_disctance.d}\t{outgroup_disctance.y.id}\t{outgroup_disctance.d}\t{isContaminate}')

            yield ingroup_disctance


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
    path_ingroup_reference = Path(argv[2])
    path_outgroup_reference = Path(argv[3])
    contaminatePath = Path(argv[4])
    decontaminatePath = Path(argv[5])
    summaryPath = Path(argv[6])

    ts = perf_counter()

    file_data = SequenceFile.Tabfile(path_data)
    file_ingroup_reference = SequenceFile.Tabfile(path_ingroup_reference)
    file_outgroup_reference = SequenceFile.Tabfile(path_outgroup_reference)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')
    ingroup_reference = Sequences.fromFile(file_ingroup_reference, idHeader='seqid', seqHeader='sequence')
    outgroup_reference = Sequences.fromFile(file_outgroup_reference, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    in_pairs = SequencePairs.fromProduct(data, ingroup_reference)
    out_pairs = SequencePairs.fromProduct(data, outgroup_reference)

    #Ingroup
    normalizeInPair = normalizePairs(in_pairs)
    aligner = PairwiseAligner.Biopython()
    aligned_in_pairs = aligner.align_pairs(normalizeInPair)

    in_distances = calc(aligned_in_pairs)

    in_minimums = get_minimum(in_distances)

    #outGroup
    normalizeOutPair = normalizePairs(out_pairs)
    aligned_out_pairs = aligner.align_pairs(normalizeOutPair)

    out_distances = calc(aligned_out_pairs)

    out_minimums = get_minimum(out_distances)

    allPairs = zip(data, in_minimums, out_minimums)

    d = decontaminateSeq(allPairs, contaminatePath, decontaminatePath, summaryPath)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
