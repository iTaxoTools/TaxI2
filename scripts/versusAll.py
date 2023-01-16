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
    def _open(self, path, mode, species_partition, genera_partition):
        self.species_partition = species_partition
        self.genera_partition = genera_partition
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
        genusX = self.genera_partition[idx]
        genusY = self.genera_partition[idy]
        speciesX = self.species_partition[idx]
        speciesY = self.species_partition[idy]
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


def create_metric_and_partition_combinations_lists(spartition, metrics):
    combinations = {}
    subsets = set(spartition.values())
    for metric in metrics:
        combinations[str(metric)] = {}
        for x in subsets:
            for y in subsets:
                combinations[str(metric)][(x, y)] = []
    return combinations


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


def calculateAll3Ms(pairDict):
    for metric in pairDict.keys():
        calculate3Ms(pairDict[metric])


def calculate_statistics_all(data: Sequences, path_stats: Path):
    allStats = StatisticsCalculator()

    for seq in data:
        allStats.add(seq.seq)
        yield seq

    with StatisticsHandler.Single(path_stats, 'w') as file:
        stats = allStats.calculate()
        file.write(stats)


def calculate_statistics_partition(data: Sequences, path: Path, partition: Partition, group_name: str):
    calculators = dict()
    for subset in partition.values():
        if subset not in calculators:
            calculators[subset] = StatisticsCalculator(group=subset)

    for seq in data:
        subset = partition[seq.id]
        calculators[subset].add(seq.seq)
        yield seq

    with StatisticsHandler.Groups(path, 'w', group_name=group_name) as file:
        for calc in calculators.values():
            stats = calc.calculate()
            file.write(stats)


def calculate_statistics_species(data: Sequences, path: Path, species_partition: Partition):
    yield from calculate_statistics_partition(data, path, species_partition, 'species')


def calculate_statistics_genera(data: Sequences, path: Path, genera_partition: Partition):
    yield from calculate_statistics_partition(data, path, genera_partition, 'genera')


def align_pairs(pairs: SequencePairs):
    aligner = PairwiseAligner.Biopython()
    return aligner.align_pairs_parallel(pairs)


def write_pairs(pairs: SequencePairs, path: Path):
    with SequencePairHandler.Formatted(path, 'w') as file:
        for pair in pairs:
            file.write(pair)
            yield pair


def calculate_distances(pairs: SequencePairs, metrics: list[DistanceMetric]):
    for x, y in pairs:
        for metric in metrics:
            yield metric.calculate(x, y)


def write_distances_linear(distances: Distances, path: Path):
    with DistanceHandler.Linear.WithExtras(path, 'w') as file:
        for d in distances:
            file.write(d)
            yield d


def write_distances_matrix(distances: Distances, metric: DistanceMetric, path: Path):
    with DistanceHandler.Matrix(path, 'w') as file:
        for d in distances:
            if d.metric.type == metric.type:
                file.write(d)
            yield d


def write_distances_multimatrix(distances: Distances, metrics: list[DistanceMetric], path: Path):
    for metric in metrics:
        distances = write_distances_matrix(distances, metric, path / f'matrix.{metric.label}.tsv')
    return distances


def main():
    path_data = Path(argv[1])
    path_out = Path(argv[2])
    path_out.mkdir(exist_ok=True)
    path_stats = path_out / 'stats'
    path_stats.mkdir(exist_ok=True)
    path_align = path_out / 'align'
    path_align.mkdir(exist_ok=True)
    path_distances = path_out / 'distances'
    path_distances.mkdir(exist_ok=True)

    metrics = [
        DistanceMetric.Uncorrected(),
        DistanceMetric.UncorrectedWithGaps(),
        DistanceMetric.JukesCantor(),
        DistanceMetric.Kimura2P(),
    ]

    ts = perf_counter()

    species_partition = Partition.fromPath(path_data, PartitionHandler.Tabfile, idHeader='seqid', subHeader='organism')
    genera_partition = Partition.fromPath(path_data, PartitionHandler.Tabfile, idHeader='seqid', subHeader='organism', filter=PartitionHandler.subset_first_word)

    data = Sequences.fromPath(path_data, SequenceHandler.Tabfile, idHeader='seqid', seqHeader='sequence')
    data = data.normalize()

    data_left = data
    data_left = calculate_statistics_all(data_left, path_stats / 'all.tsv')
    data_left = calculate_statistics_species(data_left, path_stats / 'species.tsv', species_partition)
    data_left = calculate_statistics_genera(data_left, path_stats / 'genera.tsv', genera_partition)
    total = len(data)

    species_combinations = create_metric_and_partition_combinations_lists(species_partition, metrics)
    genera_combinations = create_metric_and_partition_combinations_lists(genera_partition, metrics)

    pairs = SequencePairs.fromProduct(data_left, data)
    pairs = align_pairs(pairs)
    pairs = write_pairs(pairs, path_align / 'alignments.txt')

    distances = calculate_distances(pairs, metrics)
    distances = write_distances_linear(distances, path_distances / 'linear.tsv')
    distances = write_distances_multimatrix(distances, metrics, path_distances)

    for _ in distances:
        pass

    return


    distances = addDistanceScore(distances, species_combinations, genera_combinations, species_partition, genera_partition)

    # write matrics file
    distances = iter_write_distances_matrix(distances, path_out)

    # write linear file
    distances = iter_write_distances_linear(distances, path_out_2)

    with SummaryHandler('summaryFile', 'w', species_partition, genera_partition) as file:
        for distance in distances:
            if 'organism' in distance.x.extras:
                del distance.x.extras['organism']
            if 'organism' in distance.y.extras:
                del distance.y.extras['organism']
            file.write(distance)

    #calculate mean, min, max

    calculateAll3Ms(species_combinations)
    calculateAll3Ms(genera_combinations)

    print(genera_combinations)

    # #write subset pairs
    writeSubsetAginstItself(species_combinations, 'subsetAgainstItself', metrics)
    writeSubsetPairs(species_combinations, 'subsetPair', metrics)
    #
    #
    # #write genus pairs
    writeSubsetAginstItself(genera_combinations, 'genusAgainstItself', metrics)
    writeSubsetPairs(genera_combinations, 'genusPair', metrics)


    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
