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
        print(f"\r Loading... {index}/{total} = {(round(float(index)/float(total)  * 100, 2))} ", end="")
        yield distance


def derepSeq(disctances, dereplicatedPath, excludedPath, trashold=0.06):
    excludedSet = set()
    includeSet = set()

    for disctance in disctances:
        extras = [v for k, v in disctance.y.extras.items()]
        extrasString = '\t'.join(extras)

        if disctance.d <= trashold:
            with open(excludedPath, 'a') as f:
                if disctance.x.id != disctance.y.id and disctance.x.id not in excludedSet:
                    f.write(f"\n{disctance.y.id}\t{extrasString}\t{disctance.y.seq}")
                    excludedSet.add(disctance.y.id)
        else:
            with open(dereplicatedPath, 'a') as f:
                if disctance.x.id not in includeSet and disctance.x.id not in excludedSet:
                    f.write(f"\n{disctance.x.id}\t{extrasString}\t{disctance.x.seq}")
                    includeSet.add(disctance.x.id)

        yield disctance




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

    headers = file_data.getHeader()

    if not os.path.exists(excludedPath):
        with open(excludedPath, 'w') as f:
            f.write('\t'.join(headers))

    if not os.path.exists(dereplicatedPath):
        with open(dereplicatedPath, 'w') as f:
            f.write('\t'.join(headers))

    d = derepSeq(distances, dereplicatedPath, excludedPath)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
