from itertools import groupby
from pathlib import Path
from sys import argv
from time import perf_counter

from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.partitions import *
from typing import NamedTuple


class SubsetDistance(NamedTuple):
    subX: str
    subY: str
    d: float

class Stats(NamedTuple):
    minStat: float
    maxStat: float
    sumStat: float
    length: int
    @property
    def mean(self):
        return self.sumStat/self.length

def calc(aligned_pairs, metric=DistanceMetric.Uncorrected()):
    for x, y in aligned_pairs:
        if x.id == y.id:
            continue
        yield metric.calculate(x, y)


def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d


def get_maximum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = max(g, key = lambda x: x.d)
        yield d


def progress(distances, total):

    for index, distance in enumerate(distances):
        print(f"\r Loading... {index}/{total} = {(round(float(index)/float(total)  * 100, 2))}", end="")
        yield distance


def getSubsetDistances(distances, spartition):

    for distance in distances:
        yield SubsetDistance(spartition[distance.x.id], spartition[distance.y.id], distance.d)


def getSubsetPairs(subsetDistances, pairDict):

    for subsetDistance in subsetDistances:
        if not (subsetDistance.subX, subsetDistance.subY) in pairDict:
            pairDict[(subsetDistance.subX, subsetDistance.subY)] = []
        pairDict[(subsetDistance.subX, subsetDistance.subY)].append(subsetDistance.d)


        yield subsetDistance


def getStats(pairDict):
    minStat = float('inf')
    maxStat = float('-inf')
    sumPairs = 0
    for key, val in pairDict.items():
        length = len(val)
        for stat in val:
            sumPairs += float(stat)
            minStat = min(minStat, float(stat))
            maxStat = max(maxStat, float(stat))
        pairDict.update([(key, (minStat, maxStat, sumPairs/length))])
        minStat = float('inf')
        maxStat = float('-inf')
        sumPairs = 0


def multiply(g, n):
    return (x for x in g for i in range(n))

def main():
    path_data = Path(argv[1])

    ts = perf_counter()

    spartition_file = PartitionFile.Tabfile(path_data)
    gpartition_file = PartitionFile.Tabfile.Genus(path_data)
    spartition = Partition.fromFile(spartition_file, idHeader='seqid', subsetHeader='organism')
    gpartition = Partition.fromFile(gpartition_file, idHeader='seqid', subsetHeader='organism')




    file_data = SequenceFile.Tabfile(path_data)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = data.normalize()

    pairs = SequencePairs.fromProduct(data, data)

    aligner = PairwiseAligner.Biopython()
    aligned_pairs = aligner.align_pairs(pairs)

    distances = calc(aligned_pairs)
    distances = multiply(distances, 2)

    pairDict = {}
    subset_distances = getSubsetDistances(distances, spartition)

    subset_pairs = getSubsetPairs(subset_distances, pairDict)


    gpairDict = {}
    gdistances = getSubsetDistances(distances, gpartition)

    genus_pairs = getSubsetPairs(gdistances, gpairDict)

    for subset_pair in zip(subset_pairs, genus_pairs):
        pass


    getStats(pairDict)
    getStats(gpairDict)


    print(pairDict)
    print(gpairDict)


    return

    minimums = get_minimum(distances)


    distances = progress(alldistance, total)

    outFile = DistanceFile.LinearWithExtras(path_out)
    outFile.write(distances)

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
