from itertools import groupby
from pathlib import Path
from sys import argv
from time import perf_counter
from statistics import median, stdev
from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.partitions import *
from itaxotools.taxi3.handlers import *
from itaxotools.taxi3.handlers import *
from itaxotools.taxi3.statistics import *
from typing import NamedTuple
import numpy as np


class SubsetDistance(NamedTuple):
    subX: str
    subY: str
    d: float


class SummaryHandler(DistanceHandler.Linear.WithExtras):
    def _open(self, path, mode, spartitionDict, gpartitionDict):
        self.spartitionDict = spartitionDict
        self.gpartitionDict = gpartitionDict
        super()._open(path, mode, tagX=' (query 1)', tagY=' (query 2)')

    def _write_headers(self, file: FileHandler.Tabfile, line: list[Distance]):
        if self.wrote_headers:
            return
        idxHeader = self.idxHeader + self.tagX
        idyHeader = self.idyHeader + self.tagY
        extrasX = [key + self.tagX for key in line[0].x.extras.keys()]
        extrasY = [key + self.tagY for key in line[0].y.extras.keys()]
        metrics = [str(distance.metric) for distance in line]
        infoX = ('genus' + self.tagX, 'species' + self.tagX)
        infoY = ('genus' + self.tagY, 'species' + self.tagY)
        out = (idxHeader, idyHeader, *metrics, *extrasX, *extrasY, *infoX, *infoY, 'comparison_type')
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[Distance]):
        idx = line[0].x.id
        idy = line[0].y.id
        extrasX = line[0].x.extras.values()
        extrasY = line[0].y.extras.values()
        scores = [self.distanceToText(distance.d) for distance in line]
        genusX = self.gpartitionDict[idx]
        genusY = self.gpartitionDict[idy]
        speciesX = self.spartitionDict[idx]
        speciesY = self.spartitionDict[idy]
        comparison_type = self._get_comparison_type(genusX, genusY, speciesX, speciesY)
        out = (idx, idy, *scores, *extrasX, *extrasY, genusX, speciesX, genusY, speciesY, comparison_type)
        file.write(out)

    @staticmethod
    def _get_comparison_type(genusX, genusY, speciesX, speciesY) -> str:
        if genusX == genusY:
            if speciesX == speciesY:
                return 'intra-species'
            else:
                return 'inter-species'
        return 'inter-genus'


def calc(aligned_pairs, metric=DistanceMetric.Uncorrected()):
    for x, y in aligned_pairs:
        yield metric.calculate(x, y)


def calAll(aligned_pairs, metrics):
    for x, y in aligned_pairs:
        for metric in metrics:
            yield metric.calculate(x, y)


def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d


def progress(distances, total):

    for index, distance in enumerate(distances):
        print(f"\r Loading... {index}/{total} = {(round(float(index)/float(total)  * 100, 2))}", end="")
        yield distance


def getSubset(pairs, spartition):

    for pair in pairs:
        yield SubsetDistance(spartition[pair.x.id], spartition[pair.y.id], 0)


def getDictFromPartition(spartition, metrics):
    pairDict = {}
    allSubsets = set(spartition.values())
    for metric in metrics:
        pairDict[str(metric)] = {}
        for x in allSubsets:
            for y in allSubsets:
                pairDict[str(metric)][(x, y)] = []
    return pairDict


def calculate3Ms(pairDict):
    minStat = float('inf')
    maxStat = float('-inf')
    sumPairs = 0
    for key, val in pairDict.items():
        length = 0

        for stat in val:
            if stat == 0.0:
                continue
            length += 1
            sumPairs += float(stat)
            minStat = min(minStat, float(stat))
            maxStat = max(maxStat, float(stat))

        if sumPairs == 0:
            print(key, (float('nan'), float('nan'), float('nan')))
            pairDict.update([(key, (float('nan'), float('nan'), float('nan')))])
        else:
            print(key, (sumPairs/length, minStat, maxStat))
            pairDict.update([(key, (sumPairs/length, minStat, maxStat))])
        minStat = float('inf')
        maxStat = float('-inf')
        sumPairs = 0


def multiply(g, n):
    return (x for x in g for i in range(n))


def writeSubsetPairs(pairDict, path, metrics):

    headerList = []
    for metric in metrics:
        for m in ["mean", "minimum", "maximum"]:
            headerList.append(f'{m} {str(metric)}')

    with FileHandler.Tabfile(path, 'w', columns=('target', 'query', *headerList)) as f:
        species_target = '\r'
        for key, val in pairDict[str(metrics[0])].items():
            target, query = key[0], key[1]
            current_target = ''
            if species_target != target:
                current_target = target
                species_target = target
            buffer = []
            buffer.extend(str(x) for x in val)
            for metricIndx in range(1, len(metrics)):
                buffer.extend(str(x) for x in pairDict[str(metrics[metricIndx])][(target, query)])
            f.write((current_target, query, *buffer))


def addDistanceScore(distances, pairDict, genusPairDict, spartition, gpartition):

    for distance in distances:
        metric = str(distance.metric)
        pairDict[metric][spartition[distance.x.id], spartition[distance.y.id]].append(distance.d)
        genusPairDict[metric][gpartition[distance.x.id], gpartition[distance.y.id]].append(distance.d)

        yield distance


def writeSubsetAginstItself(pairDict, path, metrics):
    headerList = []
    for metric in metrics:
        for m in ["mean", "minimum", "maximum"]:
            headerList.append(f'{m} {str(metric)}')

    with FileHandler.Tabfile(path, 'w', columns=('target', *headerList)) as f:
        for key, val in pairDict[str(metrics[0])].items():
            sub1, sub2 = key[0], key[1]
            if sub1 == sub2:
                buffer = []
                buffer.extend(str(x) for x in val)
                for metricIndx in range(1, len(metrics)):
                    buffer.extend(str(x) for x in pairDict[str(metrics[metricIndx])][(sub1, sub2)])
                f.write((sub1, *buffer))

def iter_write_distances_linear(distances, path):
    with DistanceHandler.Linear.WithExtras(path, 'w') as file:
        for d in distances:
            file.write(d)
            yield d

def iter_write_distances_matrix(distances, path):
    with DistanceHandler.Matrix(path, 'w') as file:
        for d in distances:
            file.write(d)
            yield d


def calculateAll3Ms(pairDict):
    for metric in pairDict.keys():
        calculate3Ms(pairDict[metric])



def main():
    path_data = Path(argv[1])
    path_out = Path(argv[2])
    path_out_2 = Path(argv[3])
    pairDict = {}
    genusPairDict = {}
    metrics = [
        DistanceMetric.Uncorrected(),
        DistanceMetric.UncorrectedWithGaps(),
        DistanceMetric.JukesCantor(),
        DistanceMetric.Kimura2P(),
    ]

    ts = perf_counter()

    spartition_file = PartitionFile.Tabfile(path_data)
    gpartition_file = PartitionFile.Tabfile.Genus(path_data)
    spartitionDict = Partition.fromFile(spartition_file, idHeader='seqid', subsetHeader='organism')
    gpartitionDict = Partition.fromFile(gpartition_file, idHeader='seqid', subsetHeader='organism')

    pairDict = getDictFromPartition(spartitionDict, metrics)
    genusPairDict = getDictFromPartition(gpartitionDict, metrics)

    data = Sequences.fromPath(path_data, SequenceHandler.Tabfile, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = data.normalize()

    allStats = StatisticsCalculator()
    genusStats = dict()
    for genus in gpartitionDict.values():
        if genus not in genusStats:
            genusStats[genus] = StatisticsCalculator(group=genus)

    for seq in data:
        allStats.add(seq.seq)
        genus = gpartitionDict[seq.id]
        genusStats[genus].add(seq.seq)

    with StatisticsHandler.Single('stats.all', 'w') as file:
        stats = allStats.calculate()
        file.write(stats)

    with StatisticsHandler.Groups('stats.groups', 'w', group_name='genus') as file:
        for genus, calc in genusStats.items():
            stats = calc.calculate()
            file.write(stats)

    return
    # print(statsCalculator.calculateGenusStats())
    # print(statsCalculator.calculateAllStats())

    #write stats

    # calculateAllStatistics(data, stats)
    # writeStatistics(stats)

    # subsetStatistic = calculateStatistics(data, spartitionDict)
    # writeStatistics(subsetStatistic)

    #Create pairs

    pairs = SequencePairs.fromProduct(data, data)

    aligner = PairwiseAligner.Biopython()
    aligned_pairs = aligner.align_pairs_parallel(pairs)

    distances = calAll(aligned_pairs, metrics)

    distances = addDistanceScore(distances, pairDict, genusPairDict, spartitionDict, gpartitionDict)

    # write matrics file
    distances = iter_write_distances_matrix(distances, path_out)

    # write linear file
    distances = iter_write_distances_linear(distances, path_out_2)

    with SummaryHandler('summaryFile', 'w', spartitionDict, gpartitionDict) as file:
        for distance in distances:
            if 'organism' in distance.x.extras:
                del distance.x.extras['organism']
            if 'organism' in distance.y.extras:
                del distance.y.extras['organism']
            file.write(distance)

    #calculate mean, min, max

    calculateAll3Ms(pairDict)
    calculateAll3Ms(genusPairDict)

    print(genusPairDict)

    # #write subset pairs
    writeSubsetAginstItself(pairDict, 'subsetAgainstItself', metrics)
    writeSubsetPairs(pairDict, 'subsetPair', metrics)
    #
    #
    # #write genus pairs
    writeSubsetAginstItself(genusPairDict, 'genusAgainstItself', metrics)
    writeSubsetPairs(genusPairDict, 'genusPair', metrics)


    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
