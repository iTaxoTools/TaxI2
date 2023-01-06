from __future__ import annotations

import statistics
import numpy as np
from collections import Counter
from typing import NamedTuple
import itertools
from math import inf
from enum import Enum


class Counts(NamedTuple):

    total: int
    nucleotides: int
    missing: int
    gaps: int
    a: int
    t: int
    c: int
    g: int

    @classmethod
    def from_sequence(cls, seq: str) -> Counts:
        counter = Counter(seq)
        return cls(
            total = len(seq),
            nucleotides = len(seq) - counter['-'],
            missing = counter['N'],
            gaps = counter['-'],
            a = counter['A'],
            t = counter['T'],
            c = counter['C'],
            g = counter['G'],
        )


class Statistic(Enum):

    SequenceCount = 'totalSeq', int
    BP_0_100 = 'lessThan100BP', int
    BP_101_300 = 'between100_300BP', int
    BP_301_1000 = 'between301_1000BP', int
    BP_1001_plus = 'greaterThan1000BP', int
    Minimum = 'minimumLength', int
    Maximum = 'maximumLength', int
    Mean = 'meanLength', float
    Median = 'medianLength', float
    Stdev = 'stdLength', float
    N50 = 'N50', float
    L50 = 'L50', float
    N90 = 'N90', float
    L90 = 'L90', float
    Total = 'total_seq_length', float
    PercentA = 'percentageA', float
    PercentT = 'percentageT', float
    PercentC = 'percentageC', float
    PercentG = 'percentageG', float
    PercentGC = 'GC_content', float
    PercentAmbiguous = 'percentageAmbiguity', float
    PercentMissing = 'percentageMissingData', float
    PercentMissingGaps = 'percentageMissingDataWithGap', float

    def __init__(self, label, type):
        self.label = label
        self.type = type

    def __repr__(self):
        return f'<{type(self).__name__}.{self._name_}>'

    def __str__(self):
        return self.label


class Statistics(dict[Statistic, ...]):
    @classmethod
    def from_sequences(cls, sequences: iter[str]) -> Statistics:
        all_nucleotides = []
        length = 0

        bp_0_100 = 0
        bp_101_300 = 0
        bp_101_300 = 0
        bp_301_1000 = 0
        bp_1001_plus = 0
        minimum = inf
        maximum = -inf
        sum_total = 0
        sum_nucleotides = 0
        sum_missing = 0
        sum_gaps = 0
        sum_a = 0
        sum_t = 0
        sum_c = 0
        sum_g = 0

        counters = (Counts.from_sequence(seq) for seq in sequences)
        for counter in counters:
            length += 1
            all_nucleotides.append(counter.nucleotides)

            if counter.nucleotides <= 100:
                bp_0_100 += 1
            elif counter.nucleotides <= 300:
                bp_101_300 += 1
            elif counter.nucleotides <= 1000:
                bp_301_1000 += 1
            else:
                bp_1001_plus += 1

            minimum = min(minimum, counter.nucleotides)
            maximum = max(maximum, counter.nucleotides)
            sum_total += counter.total
            sum_nucleotides += counter.nucleotides
            sum_missing += counter.missing
            sum_gaps += counter.gaps
            sum_a += counter.a
            sum_t += counter.t
            sum_c += counter.c
            sum_g += counter.g

        mean = sum_nucleotides / length if length else 0
        median = statistics.median(all_nucleotides) if length else 0
        stdev = statistics.stdev(all_nucleotides) if len(all_nucleotides) > 1 else 0

        sum_cg = sum_c + sum_g
        sum_ambiguous = sum_nucleotides - sum_missing - sum_a - sum_t - sum_c - sum_g
        sum_missing_and_gaps = sum_missing + sum_gaps

        return cls({
            Statistic.SequenceCount: length,
            Statistic.BP_0_100: bp_0_100,
            Statistic.BP_101_300: bp_101_300,
            Statistic.BP_301_1000: bp_301_1000,
            Statistic.BP_1001_plus: bp_1001_plus,
            Statistic.Minimum: minimum,
            Statistic.Maximum: maximum,
            Statistic.Mean: mean,
            Statistic.Median: median,
            Statistic.Stdev: stdev,
            # Statistic.N50
            # Statistic.L50
            # Statistic.N90
            # Statistic.L90
            Statistic.Total: sum_total,
            Statistic.PercentA: 100 * sum_a / sum_total if sum_total else 0,
            Statistic.PercentT: 100 * sum_t / sum_total if sum_total else 0,
            Statistic.PercentC: 100 * sum_c / sum_total if sum_total else 0,
            Statistic.PercentG: 100 * sum_g / sum_total if sum_total else 0,
            Statistic.PercentGC: 100 * sum_cg / sum_total if sum_total else 0,
            Statistic.PercentAmbiguous: 100 * sum_ambiguous / sum_total if sum_total else 0,
            Statistic.PercentMissing: 100 * sum_missing / sum_total if sum_total else 0,
            Statistic.PercentMissingGaps: 100 * sum_missing_and_gaps / sum_total if sum_total else 0,
        })

    @staticmethod
    def calculate_NL(self, list_of_lengths, nOrL, arg):
        # needs adjusting

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
