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
        print(f"\r Loading... {index}/{total} = {(round(float(index)/float(total)  * 100, 2))} ", end="")
        yield distance


def decontaminateSeq(disctances, contaminatePath, decontaminatePath, summaryPath,similarity=0.09):
    excludedSet = set()
    includeSet = set()

    for disctance in disctances:
        extras = [v for k, v in disctance.y.extras.items()]
        extrasString = '\t'.join(extras)

        if disctance.d <= similarity:
            with open(decontaminatePath, 'a') as f:
                if disctance.x.id != disctance.y.id:
                    if disctance.x.id not in excludedSet:
                        f.write(f"\n{disctance.x.id}\t{extrasString}\t{disctance.x.seq}")
                        excludedSet.add(disctance.x.id)
        else:
            with open(contaminatePath, 'a') as f:
                if disctance.x.id not in excludedSet:
                    if disctance.x.id not in includeSet:
                        f.write(f"\n{disctance.x.id}\t{extrasString}\t{disctance.x.seq}")
                    includeSet.add(disctance.x.id)

        with open(summaryPath, 'a') as f:
            isContaminate = True if disctance.x.id in includeSet else False
            f.write(f'\n{disctance.x.id}\t{disctance.y.id}\t{disctance.d}\t{isContaminate}')

        yield disctance


def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d

def main():

    path_data = Path(argv[1])
    path_reference = Path(argv[2])
    contaminatePath = 'contaminates.txt'
    decontaminatePath = 'decontaminated.txt'
    summaryPath = 'summary.txt'

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

    headers = file_data.getHeader()

    if not os.path.exists(contaminatePath):
        with open(contaminatePath, 'w') as f:
            f.write('\t'.join(headers))

    if not os.path.exists(summaryPath):
        with open(summaryPath, 'w') as f:
            f.write("seqid_query\tclosest possible contaminant\tdistance\tis_contaminant")

    if not os.path.exists(decontaminatePath):
        with open(decontaminatePath, 'w') as f:
            f.write('\t'.join(headers))

    d = decontaminateSeq(minimums, contaminatePath, decontaminatePath, summaryPath)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
