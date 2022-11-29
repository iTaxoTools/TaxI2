from statistics import median, stdev
import numpy as np
from collections import Counter
import itertools


class StatisticsCalculator:
    def __init__(self):
        self.bufferStats = {
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
        self.bufferLengths = []

    def prepareStats(self, sequenceData):

        # todo: remove all mentions of allStats

        self.seqStats = {}
        sequence = sequenceData.seq
        seqLength = len(sequence) - sequence.count('-')

        self.seqStats['seqLength'] = seqLength

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

    def addSequence(self, sequence):
        seqStats = self.prepareStats(sequence)
        self.bufferLengths.append(seqStats['seqLength'])

        seqStats['meanLength'] = 0
        seqStats['medianLength'] = 0
        seqStats['stdLength'] = 0
        seqStats['N50'] = 0
        seqStats['L50'] = 0
        seqStats['N90'] = 0
        seqStats['L90'] = 0

        self.bufferLengths.append(seqStats['total_seq_length'])
        for key, val in seqStats.items():

            if key == 'totalSeq':
                self.bufferStats['totalSeq'] += 1
            elif key == 'total_seq_length':
                self.bufferStats['total_seq_length'] += seqStats['total_seq_length']
            elif key == 'minimumLength':
                self.bufferStats['minimumLength'] = min(self.bufferStats['minimumLength'],
                                                                  seqStats['minimumLength'])
            elif key == 'maximumLength':
                self.bufferStats['maximumLength'] = min(self.bufferStats['maximumLength'],
                                                                  seqStats['maximumLength'])
            elif key == 'lessThan100BP':
                self.bufferStats['lessThan100BP'] += seqStats['lessThan100BP']
            elif key == 'between100_300BP':
                self.bufferStats['between100_300BP'] += seqStats['between100_300BP']
            elif key == 'between301_1000BP':
                self.bufferStats['between301_1000BP'] += seqStats['between301_1000BP']
            elif key == 'greaterThan1000BP':
                self.bufferStats['greaterThan1000BP'] += seqStats['greaterThan1000BP']
            elif key == 'percentageA':
                self.bufferStats['percentageA'] += seqStats['percentageA'] * seqStats['total_seq_length']
            elif key == 'percentageC':
                self.bufferStats['percentageC'] += seqStats['percentageC'] * seqStats['total_seq_length']
            elif key == 'percentageG':
                self.bufferStats['percentageG'] += seqStats['percentageG'] * seqStats['total_seq_length']
            elif key == 'percentageT':
                self.bufferStats['percentageT'] += seqStats['percentageT'] * seqStats['total_seq_length']
            elif key == 'GC_content':
                self.bufferStats['GC_content'] += seqStats['GC_content'] * seqStats['total_seq_length']
            elif key == 'percentageAmbiguity':
                self.bufferStats['percentageAmbiguity'] += seqStats['percentageAmbiguity'] * seqStats[
                    'total_seq_length']
            elif key == 'percentageMissingData':
                self.bufferStats['percentageMissingData'] += seqStats['percentageMissingData'] * seqStats[
                    'total_seq_length']
            elif key == 'percentageMissingDataWithGap':
                self.bufferStats['percentageMissingDataWithGap'] += seqStats[
                                                                                  'percentageMissingDataWithGap'] * \
                                                                              seqStats['total_seq_length']
            # else:
            #     self.bufferStats[key] += val


    def calculate(self):

        self.bufferStats['meanLength'] = self.bufferStats['total_seq_length'] / self.bufferStats[
            'totalSeq']
        self.bufferStats['medianLength'] = median(self.genusLengths[key])
        self.bufferStats['stdLength'] = stdev(self.genusLengths[key])

        self.bufferStats['N50'] = self.calculate_NL(self.genusLengths[key], 'N', 50)
        self.bufferStats['L50'] = self.calculate_NL(self.genusLengths[key], 'L', 50)
        self.bufferStats['N90'] = self.calculate_NL(self.genusLengths[key], 'N', 90)
        self.bufferStats['L90'] = self.calculate_NL(self.genusLengths[key], 'L', 90)

        self.bufferStats['percentageA'] = self.bufferStats['percentageA'] / self.bufferStats[
            'total_seq_length']
        self.bufferStats['percentageC'] = self.bufferStats['percentageC'] / self.bufferStats[
            'total_seq_length']
        self.bufferStats['percentageG'] = self.bufferStats['percentageG'] / self.bufferStats[
            'total_seq_length']
        self.bufferStats['percentageT'] = self.bufferStats['percentageT'] / self.bufferStats[
            'total_seq_length']

        self.bufferStats['GC_content'] = self.bufferStats['GC_content'] / self.bufferStats[
            'total_seq_length']
        self.bufferStats['percentageAmbiguity'] = self.bufferStats['percentageAmbiguity'] / \
                                                      self.bufferStats['total_seq_length']
        self.bufferStats['percentageMissingData'] = self.bufferStats['percentageMissingData'] / \
                                                        self.bufferStats['total_seq_length']
        self.bufferStats['percentageMissingDataWithGap'] = self.bufferStats[
                                                                   'percentageMissingDataWithGap'] / \
                                                               self.bufferStats[
                                                                   'total_seq_length']
        return self.bufferStats

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
