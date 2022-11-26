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
from typing import NamedTuple
import numpy as np


class SubsetDistance(NamedTuple):
    subX: str
    subY: str
    d: float


class SequenceStatistics(NamedTuple):
    subset: str
    numberOfSequence: int
    lessThan100BP: int
    between100To300BP: int
    between301To1000BP: int
    greaterThan1000BP: int
    minimumLenght: int
    maximumLenght: int
    meanLenght: float
    medianLenght: float
    standardDeviation: float
    N50: float
    L50: float
    N90: float
    L90: float
    totalLength: int
    percentageOfA: float
    percentageOfC: float
    percentageOfG: float
    percentageOfT: float
    GCcontent: float
    percentageOfAmbiguity: float
    percentageOfMissingData: float
    percentageOfMissingDataWithGaps: float


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


def getDictFromPartition(spartition):
    pairDict = {}
    allSubsets = set(spartition.values())
    for x in allSubsets:
        for y in allSubsets:
            pairDict[(x, y)] = []
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


def writeSubsetPairs(pairDict, path):
    with open(path, 'w') as f:
        f.write(f'target\tquery\tmean pairwise uncorrected distance\tminimum pairwise uncorrected distance\tmaximum pairwise uncorrected distance\t\n')
        species_target = ''
        for key, val in pairDict.items():
            targetSubset, querySubset, mean, minn, maxx = key[0], key[1], val[0], val[1], val[2]

            if not species_target or species_target != targetSubset:
                f.write(f'{targetSubset}\t')
            while targetSubset == species_target:
                f.write(' ' * len(targetSubset) + '\t')
                break
            f.write(f'{querySubset}\t{mean}\t{minn}\t{maxx}\n')
            species_target = targetSubset


def calculateAllStatistics(data, stats):

    sequenceLengths = []
    for sequenceData in data:
        sequence = sequenceData.seq
        seqLength = len(sequence) - sequence.count('-')
        stats['totalSeq'] += 1
        stats['totalLengthOfSeq'] += seqLength
        sequenceLengths.append(seqLength)

        if seqLength <= 100:
            stats['lessThan100BP'] += 1
        elif seqLength <= 300:
            stats['between100-300BP'] += 1
        elif seqLength <= 1000:
            stats['between301-1000BP'] += 1
        else:
            stats['greaterThan1000BP'] += 1

        stats['minimumLenght'] = min(seqLength, stats['minimumLenght'])
        stats['maximumLenght'] = max(seqLength, stats['maximumLenght'])

        stats['percentageofA'] += sequence.count('A')
        stats['percentageofC'] += sequence.count('C')
        stats['percentageofG'] += sequence.count('G')
        stats['percentageofT'] += sequence.count('T')
        stats['GC content'] += sequence.count('G') + sequence.count('C')
        stats['percentageofAmbiguity'] += sequence.count('R') + sequence.count('Y')+ sequence.count('S')+ sequence.count('W')+ sequence.count('K')+ sequence.count('M')
        stats['percentageofMissingData'] += sequence.count('N') + sequence.count('?')
        stats['percentageofMissingDataWithGap'] += sequence.count('N') + sequence.count('?') + sequence.count('-')

    stats['meanLenght'] = stats['totalLengthOfSeq'] / stats['totalSeq']
    stats['medianLenght'] = median(sequenceLengths)
    stats['stdLenght'] = stdev(sequenceLengths)
    stats['N50'] = calculate_NL(sequenceLengths, 'N', 50)
    stats['L50'] = calculate_NL(sequenceLengths, 'L', 50)
    stats['N90'] = calculate_NL(sequenceLengths, 'N', 90)
    stats['L90'] = calculate_NL(sequenceLengths, 'L', 90)
    stats['percentageofA'] = stats['percentageofA'] / stats['totalLengthOfSeq']
    stats['percentageofC'] = stats['percentageofC'] / stats['totalLengthOfSeq']
    stats['percentageofG'] = stats['percentageofG'] / stats['totalLengthOfSeq']
    stats['percentageofT'] = stats['percentageofT'] / stats['totalLengthOfSeq']
    stats['GC content'] = stats['GC content'] / stats['totalLengthOfSeq']
    stats['percentageofAmbiguity'] = stats['percentageofAmbiguity'] / stats['totalLengthOfSeq']
    stats['percentageofMissingData'] = stats['percentageofMissingData'] / stats['totalLengthOfSeq']
    stats['percentageofMissingDataWithGap'] = stats['percentageofMissingDataWithGap'] / stats['totalLengthOfSeq']


def calculateStatistics(data, spartition):
    for sequenceData in data:
        seqID = spartition[sequenceData.id]
        sequence = sequenceData.seq
        seqLength = len(sequence) - sequence.count('-')
        total_seq = 1
        total_seq_length = seqLength
        lessThan100BP = 0
        between100_300BP = 0
        between301_1000BP = 0
        greaterThan1000BP = 0
        minimumLenght = seqLength
        meanLenght = total_seq_length / total_seq
        maximumLenght = seqLength
        medianLenght = median([seqLength])
        #stdLenght = stdev([seqLength])
        N50 = calculate_NL([seqLength], 'N', 50)
        L50 = calculate_NL([seqLength], 'L', 50)
        N90 = calculate_NL([seqLength], 'N', 90)
        L90 = calculate_NL([seqLength], 'L', 90)
        percentageofA = sequence.count('A') / total_seq_length
        percentageofT = sequence.count('T') / total_seq_length
        percentageofC = sequence.count('C') / total_seq_length
        percentageofG = sequence.count('G') / total_seq_length
        GC_content = (sequence.count('G') + sequence.count('C')) / total_seq_length
        percentageofAmbiguity = (sequence.count('R') + sequence.count('Y')+ sequence.count('S')+ sequence.count('W')+ sequence.count('K')+ sequence.count('M')) / total_seq_length
        percentageofMissingData = (sequence.count('N') + sequence.count('?')) / total_seq_length
        percentageofMissingDataWithGap = (sequence.count('N') + sequence.count('?') + sequence.count('-')) / total_seq_length

        if seqLength <= 100:
            lessThan100BP = 1
        elif seqLength <= 300:
            between100_300BP = 1
        elif seqLength <= 1000:
            between301_1000BP = 1
        else:
            greaterThan1000BP = 1

        yield SequenceStatistics(
            seqID,
            total_seq,
            lessThan100BP,
            between100_300BP,
            between301_1000BP,
            greaterThan1000BP,
            minimumLenght,
            maximumLenght,
            meanLenght,
            medianLenght,
            float('nan'),
            N50,
            L50,
            N90,
            L90,
            total_seq_length,
            percentageofA,
            percentageofC,
            percentageofG,
            percentageofT,
            GC_content,
            percentageofAmbiguity,
            percentageofMissingData,
            percentageofMissingDataWithGap
        )


def addDistanceScore(distances, pairDict, genusPairDict, spartition, gpartition):

    for distance in distances:
        pairDict[spartition[distance.x.id], spartition[distance.y.id]].append(distance.d)
        genusPairDict[gpartition[distance.x.id], gpartition[distance.y.id]].append(distance.d)

        yield distance


def calculate_NL(list_of_lengths, nOrL,arg):

    # tmp = []
    # for tmp_number in set(list_of_lengths):
    #     tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    # tmp.sort()
    #
    # if (len(tmp) % 2) == 0:
    #     median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    # else:
    #     median = tmp[int(len(tmp) / 2)]
    #
    # return median
    stats = {}
    seq_array = np.array(list_of_lengths)
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)

    nx = int(stats['total_bps'] * (arg / 100))
    csumn = min(csum[csum >= nx])
    l_level = int(np.where(csum == csumn)[0])
    n_level = int(sorted_lens[l_level])

    stats['L' + str(arg)] = l_level
    stats['N' + str(arg)] = n_level

    if nOrL.upper() == 'L':
        stats[nOrL.upper() + str(arg)] +=1
    return stats[nOrL.upper() + str(arg)]


def writeStatistics(stats):

    if type(stats) is dict:
        with open('AllStats.txt', 'w') as f:
            for k, v in stats.items():
                f.write(f"{k} : {v}\n")
    else:
        with open('individualStats.txt', 'w') as f:
            for stat in stats:
                f.write(''.join(str(stat)) + '\n')


def writeSubsetAginstItself(pairDict, path):
    with open(path, 'w') as f:
        f.write(f'species\tmean pairwise uncorrected distance\tminimum pairwise uncorrected distance\tmaximum pairwise uncorrected distance\n')
        for key, val in pairDict.items():
            sub1, sub2 = key[0], key[1]
            if sub1 == sub2:
                mean, mini, maxi = val[0], val[1], val[2]
                f.write(f'{sub1}\t{mean}\t{mini}\t{maxi}\n')

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

def main():
    path_data = Path(argv[1])
    path_out = Path(argv[2])
    path_out_2 = Path(argv[3])
    pairDict = {}
    genusPairDict = {}
    stats = {
        'totalSeq' : 0,
        'lessThan100BP' : 0,
        'between100-300BP' : 0,
        'between301-1000BP' : 0,
        'greaterThan1000BP' : 0,
        'minimumLenght' : float('inf'),
        'maximumLenght' : 0,
        'meanLenght' : 0,
        'medianLenght' : 0,
        'stdLenght' : 0,
        'N50': 0,
        'L50': 0,
        'N90': 0,
        'L90': 0,
        'totalLengthOfSeq': 0,
        'percentageofA' : 0,
        'percentageofC' : 0,
        'percentageofG' : 0,
        'percentageofT' : 0,
        'GC content' : 0,
        'percentageofAmbiguity': 0,
        'percentageofMissingData': 0,
        'percentageofMissingDataWithGap': 0}

    ts = perf_counter()

    spartition_file = PartitionFile.Tabfile(path_data)
    gpartition_file = PartitionFile.Tabfile.Genus(path_data)
    spartitionDict = Partition.fromFile(spartition_file, idHeader='seqid', subsetHeader='organism')
    gpartitionDict = Partition.fromFile(gpartition_file, idHeader='seqid', subsetHeader='organism')
    pairDict = getDictFromPartition(spartitionDict)
    genusPairDict = getDictFromPartition(gpartitionDict)

    data = Sequences.fromPath(path_data, SequenceHandler.Tabfile, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = data.normalize()

    #write stats

    calculateAllStatistics(data, stats)
    writeStatistics(stats)

    subsetStatistic = calculateStatistics(data, spartitionDict)
    writeStatistics(subsetStatistic)

    #Create pairs

    pairs = SequencePairs.fromProduct(data, data)

    aligner = PairwiseAligner.Biopython()
    aligned_pairs = aligner.align_pairs(pairs)

    distances = calc(aligned_pairs)

    distances = addDistanceScore(distances, pairDict, genusPairDict, spartitionDict, gpartitionDict)


    # write matrics file
    distances = iter_write_distances_matrix(distances, path_out)

    # write linear file
    distances = iter_write_distances_linear(distances, path_out_2)

    for distance in distances:
        pass

    #calculate mean, min, max
    calculate3Ms(pairDict)
    calculate3Ms(genusPairDict)

    #write subset pairs
    writeSubsetAginstItself(pairDict, 'subsetAgainstItself')
    writeSubsetPairs(pairDict, 'subsetPair')



    #write genus pairs
    writeSubsetAginstItself(genusPairDict, 'genusAgainstItself')
    writeSubsetPairs(genusPairDict, 'genusPair')



    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
