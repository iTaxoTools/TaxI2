from statistics import median, stdev
import numpy as np
from collections import Counter
import itertools


class StatisticsCalculator:
    def __init__(self):
        self.allStats = {
            'totalSeq': 0,
            'lessThan100BP': 0,
            'between100_300BP': 0,
            'between301_1000BP': 0,
            'greaterThan1000BP': 0,
            'minimumLength': float('inf'),
            'maximumLength': 0,
            'meanLength': 0,
            'medianLength': 0,
            'stdLength': 0,
            'N50': 0,
            'L50': 0,
            'N90': 0,
            'L90': 0,
            'total_seq_length': 0,
            'percentageA': 0,
            'percentageC': 0,
            'percentageG': 0,
            'percentageT': 0,
            'GC_content': 0,
            'percentageAmbiguity': 0,
            'percentageMissingData': 0,
            'percentageMissingDataWithGap': 0
        }
        self.genusStats = {}
        self.genusLengths = {}
        self.sequenceLengths = []

    def addSequences(self, sequenceData):
        self.seqStats = {}
        sequence = sequenceData.seq
        seqLength = len(sequence) - sequence.count('-')
        self.sequenceLengths.append(seqLength)

        self.seqStats['totalSeq'], self.allStats['totalSeq'] = 1, self.allStats['totalSeq'] + 1
        self.seqStats['total_seq_length'], self.allStats['total_seq_length'] = seqLength, self.allStats[
            'total_seq_length'] + seqLength

        self.seqStats['lessThan100BP'], self.allStats['lessThan100BP'] = 0, self.allStats['lessThan100BP']
        self.seqStats['between100_300BP'], self.allStats['between100_300BP'] = 0, self.allStats['between100_300BP']
        self.seqStats['between301_1000BP'], self.allStats['between301_1000BP'] = 0, self.allStats['between301_1000BP']
        self.seqStats['greaterThan1000BP'], self.allStats['greaterThan1000BP'] = 0, self.allStats['greaterThan1000BP']

        self.seqStats['minimumLength'], self.allStats['minimumLength'] = seqLength, min(self.allStats['minimumLength'],
                                                                                        seqLength)
        self.seqStats['meanLength'], self.allStats['meanLength'] = self.allStats['total_seq_length'] / self.allStats[
            'totalSeq'], 0
        self.seqStats['maximumLength'], self.allStats['maximumLength'] = seqLength, max(self.allStats['minimumLength'],
                                                                                        seqLength)
        self.seqStats['medianLength'], self.allStats['medianLength'] = median([seqLength]), 0
        self.seqStats['stdLength'], self.allStats['stdLength'] = float('nan'), 0

        self.seqStats['N50'], self.allStats['N50'] = self.calculate_NL([seqLength], 'N', 50), 0
        self.seqStats['L50'], self.allStats['L50'] = self.calculate_NL([seqLength], 'L', 50), 0
        self.seqStats['N90'], self.allStats['N90'] = self.calculate_NL([seqLength], 'N', 90), 0
        self.seqStats['L90'], self.allStats['L90'] = self.calculate_NL([seqLength], 'L', 90), 0

        self.seqStats['percentageA'], self.allStats['percentageA'] = sequence.count('A') / self.allStats[
            'total_seq_length'], \
                                                                     self.allStats['percentageA'] + sequence.count('A')
        self.seqStats['percentageT'], self.allStats['percentageT'] = sequence.count('T') / self.allStats[
            'total_seq_length'], \
                                                                     self.allStats['percentageT'] + sequence.count('T')
        self.seqStats['percentageC'], self.allStats['percentageC'] = sequence.count('C') / self.allStats[
            'total_seq_length'], \
                                                                     self.allStats['percentageC'] + sequence.count('C')
        self.seqStats['percentageG'], self.allStats['percentageG'] = sequence.count('G') / self.allStats[
            'total_seq_length'], \
                                                                     self.allStats['percentageG'] + sequence.count('G')

        self.seqStats['GC_content'], self.allStats['GC_content'] = (sequence.count('G') + sequence.count('C')) / \
                                                                   self.allStats['total_seq_length'], self.allStats[
                                                                       'GC_content'] + sequence.count(
            'G') + sequence.count('C')

        ambiguityCount = (sequence.count('R') + sequence.count('Y') + sequence.count('S') + sequence.count(
            'W') + sequence.count('K') + sequence.count('M'))
        self.seqStats['percentageAmbiguity'], self.allStats['percentageAmbiguity'] = ambiguityCount / self.allStats[
            'total_seq_length'], self.allStats['percentageAmbiguity'] + ambiguityCount

        missingDataCount = (sequence.count('N') + sequence.count('?'))
        self.seqStats['percentageMissingData'], self.allStats['percentageMissingData'] = missingDataCount / \
                                                                                         self.allStats[
                                                                                             'total_seq_length'], \
                                                                                         self.allStats[
                                                                                             'percentageMissingData'] + missingDataCount

        missingDataWithGapsCount = (sequence.count('N') + sequence.count('?') + sequence.count('-'))
        self.seqStats['percentageMissingDataWithGap'], self.allStats[
            'percentageMissingDataWithGap'] = missingDataWithGapsCount / self.allStats['total_seq_length'], \
                                              self.allStats[
                                                  'percentageMissingDataWithGap'] + missingDataWithGapsCount

        if seqLength <= 100:
            self.seqStats['lessThan100BP'], self.allStats['lessThan100BP'] = 1, self.allStats['lessThan100BP'] + 1
        elif seqLength <= 300:
            self.seqStats['between100_300BP'], self.allStats['between100_300BP'] = 1, self.allStats[
                'between100_300BP'] + 1
        elif seqLength <= 1000:
            self.seqStats['between301_1000BP'], self.allStats['between301_1000BP'] = 1, self.allStats[
                'between301_1000BP'] + 1
        else:
            self.seqStats['greaterThan1000BP'], self.allStats['greaterThan1000BP'] = 1, self.allStats[
                'greaterThan1000BP'] + 1

        return self.seqStats

    def addGenuses(self, seqStats, genusName):

        seqStats['meanLength'] = 0
        seqStats['medianLength'] = 0
        seqStats['stdLength'] = 0
        seqStats['N50'] = 0
        seqStats['L50'] = 0
        seqStats['N90'] = 0
        seqStats['L90'] = 0

        if genusName in self.genusStats and genusName in self.genusLengths:
            self.genusLengths[genusName].append(seqStats['total_seq_length'])
            for key, val in seqStats.items():

                if key == 'totalSeq':
                    self.genusStats[genusName]['totalSeq'] += 1
                elif key == 'total_seq_length':
                    self.genusStats[genusName]['total_seq_length'] += seqStats['total_seq_length']
                elif key == 'minimumLength':
                    self.genusStats[genusName]['minimumLength'] = min(self.genusStats[genusName]['minimumLength'],
                                                                      seqStats['minimumLength'])
                elif key == 'maximumLength':
                    self.genusStats[genusName]['maximumLength'] = min(self.genusStats[genusName]['maximumLength'],
                                                                      seqStats['maximumLength'])
                elif key == 'lessThan100BP':
                    self.genusStats[genusName]['lessThan100BP'] += seqStats['lessThan100BP']
                elif key == 'between100_300BP':
                    self.genusStats[genusName]['between100_300BP'] += seqStats['between100_300BP']
                elif key == 'between301_1000BP':
                    self.genusStats[genusName]['between301_1000BP'] += seqStats['between301_1000BP']
                elif key == 'greaterThan1000BP':
                    self.genusStats[genusName]['greaterThan1000BP'] += seqStats['greaterThan1000BP']
                elif key == 'percentageA':
                    self.genusStats[genusName]['percentageA'] += seqStats['percentageA'] * seqStats['total_seq_length']
                elif key == 'percentageC':
                    self.genusStats[genusName]['percentageC'] += seqStats['percentageC'] * seqStats['total_seq_length']
                elif key == 'percentageG':
                    self.genusStats[genusName]['percentageG'] += seqStats['percentageG'] * seqStats['total_seq_length']
                elif key == 'percentageT':
                    self.genusStats[genusName]['percentageT'] += seqStats['percentageT'] * seqStats['total_seq_length']
                elif key == 'GC_content':
                    self.genusStats[genusName]['GC_content'] += seqStats['GC_content'] * seqStats['total_seq_length']
                elif key == 'percentageAmbiguity':
                    self.genusStats[genusName]['percentageAmbiguity'] += seqStats['percentageAmbiguity'] * seqStats[
                        'total_seq_length']
                elif key == 'percentageMissingData':
                    self.genusStats[genusName]['percentageMissingData'] += seqStats['percentageMissingData'] * seqStats[
                        'total_seq_length']
                elif key == 'percentageMissingDataWithGap':
                    self.genusStats[genusName]['percentageMissingDataWithGap'] += seqStats[
                                                                                      'percentageMissingDataWithGap'] * \
                                                                                  seqStats['total_seq_length']
                else:
                    self.genusStats[genusName][key] += val
        else:
            self.genusStats[genusName] = seqStats
            self.genusLengths[genusName] = [seqStats['total_seq_length']]

    def calculateGenusStats(self):

        for key in self.genusStats.keys():
            self.genusStats[key]['meanLength'] = self.genusStats[key]['total_seq_length'] / self.genusStats[key][
                'totalSeq']
            self.genusStats[key]['medianLength'] = median(self.genusLengths[key])
            self.genusStats[key]['stdLength'] = stdev(self.genusLengths[key])

            self.genusStats[key]['N50'] = self.calculate_NL(self.genusLengths[key], 'N', 50)
            self.genusStats[key]['L50'] = self.calculate_NL(self.genusLengths[key], 'L', 50)
            self.genusStats[key]['N90'] = self.calculate_NL(self.genusLengths[key], 'N', 90)
            self.genusStats[key]['L90'] = self.calculate_NL(self.genusLengths[key], 'L', 90)

            self.genusStats[key]['percentageA'] = self.genusStats[key]['percentageA'] / self.genusStats[key][
                'total_seq_length']
            self.genusStats[key]['percentageC'] = self.genusStats[key]['percentageC'] / self.genusStats[key][
                'total_seq_length']
            self.genusStats[key]['percentageG'] = self.genusStats[key]['percentageG'] / self.genusStats[key][
                'total_seq_length']
            self.genusStats[key]['percentageT'] = self.genusStats[key]['percentageT'] / self.genusStats[key][
                'total_seq_length']

            self.genusStats[key]['GC_content'] = self.genusStats[key]['GC_content'] / self.genusStats[key][
                'total_seq_length']
            self.genusStats[key]['percentageAmbiguity'] = self.genusStats[key]['percentageAmbiguity'] / \
                                                          self.genusStats[key]['total_seq_length']
            self.genusStats[key]['percentageMissingData'] = self.genusStats[key]['percentageMissingData'] / \
                                                            self.genusStats[key]['total_seq_length']
            self.genusStats[key]['percentageMissingDataWithGap'] = self.genusStats[key][
                                                                       'percentageMissingDataWithGap'] / \
                                                                   self.genusStats[key][
                                                                       'total_seq_length']
        return self.genusStats

    def calculateAllStats(self):

        self.allStats['meanLength'] = self.allStats['total_seq_length'] / self.allStats['totalSeq']
        self.allStats['medianLength'] = median(self.sequenceLengths)
        self.allStats['stdLength'] = stdev(self.sequenceLengths)
        self.allStats['N50'] = self.calculate_NL(self.sequenceLengths, 'N', 50)
        self.allStats['L50'] = self.calculate_NL(self.sequenceLengths, 'L', 50)
        self.allStats['N90'] = self.calculate_NL(self.sequenceLengths, 'N', 90)
        self.allStats['L90'] = self.calculate_NL(self.sequenceLengths, 'L', 90)
        self.allStats['percentageA'] = self.allStats['percentageA'] / self.allStats['total_seq_length']
        self.allStats['percentageC'] = self.allStats['percentageC'] / self.allStats['total_seq_length']
        self.allStats['percentageG'] = self.allStats['percentageG'] / self.allStats['total_seq_length']
        self.allStats['percentageT'] = self.allStats['percentageT'] / self.allStats['total_seq_length']
        self.allStats['GC_content'] = self.allStats['GC_content'] / self.allStats['total_seq_length']
        self.allStats['percentageAmbiguity'] = self.allStats['percentageAmbiguity'] / self.allStats['total_seq_length']
        self.allStats['percentageMissingData'] = self.allStats['percentageMissingData'] / self.allStats[
            'total_seq_length']
        self.allStats['percentageMissingDataWithGap'] = self.allStats['percentageMissingDataWithGap'] / self.allStats[
            'total_seq_length']

        return self.allStats

    def calculate_NL(self, list_of_lengths, nOrL, arg):

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
            stats[nOrL.upper() + str(arg)] += 1
        return stats[nOrL.upper() + str(arg)]