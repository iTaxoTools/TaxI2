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

def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d

def main():

    path_data = Path(argv[1])
    dereplicatedPath = 'dereplicated.txt'
    excludedPath = 'excluded.txt'

    ts = perf_counter()

    file_data = SequenceFile.Tabfile(path_data)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = data.normalize()

    pairs = SequencePairs.fromProduct(data, data)

    aligner = PairwiseAligner.Biopython()
    aligned_pairs = aligner.align_pairs(pairs)

    distances = calc(aligned_pairs)

    minimums = get_minimum(distances)

    headers = file_data.getHeader()

    d = derepSeq(minimums, dereplicatedPath, excludedPath, headers)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
