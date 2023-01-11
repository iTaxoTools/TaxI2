from __future__ import annotations

import statistics
import numpy as np
from collections import Counter
from itertools import accumulate
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


class NL(NamedTuple):
    n: int
    l: int


class Statistic(Enum):

    SequenceCount = 'totalSeq', int
    BP_0 = 'zeroBP', int
    BP_1_100 = 'lessThan100BP', int
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
    PercentMissingGaps = 'percentageMissingDataWithGaps', float
    PercentGaps = 'percentageGaps', float

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
        nucleotide_counts = []
        length = 0

        bp_0 = 0
        bp_1_100 = 0
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

        counts = (Counts.from_sequence(seq) for seq in sequences)
        for count in counts:
            length += 1
            nucleotide_counts.append(count.nucleotides)

            if count.nucleotides == 0:
                bp_0 += 1
            elif count.nucleotides <= 100:
                bp_1_100 += 1
            elif count.nucleotides <= 300:
                bp_101_300 += 1
            elif count.nucleotides <= 1000:
                bp_301_1000 += 1
            else:
                bp_1001_plus += 1

            minimum = min(minimum, count.nucleotides)
            maximum = max(maximum, count.nucleotides)
            sum_total += count.total
            sum_nucleotides += count.nucleotides
            sum_missing += count.missing
            sum_gaps += count.gaps
            sum_a += count.a
            sum_t += count.t
            sum_c += count.c
            sum_g += count.g

        mean = sum_nucleotides / length if length else 0
        median = statistics.median(nucleotide_counts) if length else 0
        stdev = statistics.stdev(nucleotide_counts) if len(nucleotide_counts) > 1 else 0

        sum_cg = sum_c + sum_g
        sum_ambiguous = sum_nucleotides - sum_missing - sum_a - sum_t - sum_c - sum_g
        sum_missing_and_gaps = sum_missing + sum_gaps

        n_50, l_50 = cls._calculate_NL(nucleotide_counts, 50)
        n_90, l_90 = cls._calculate_NL(nucleotide_counts, 90)

        return cls({
            Statistic.SequenceCount: length,
            Statistic.BP_0: bp_0,
            Statistic.BP_1_100: bp_1_100,
            Statistic.BP_101_300: bp_101_300,
            Statistic.BP_301_1000: bp_301_1000,
            Statistic.BP_1001_plus: bp_1001_plus,
            Statistic.Minimum: minimum,
            Statistic.Maximum: maximum,
            Statistic.Mean: mean,
            Statistic.Median: median,
            Statistic.Stdev: stdev,
            Statistic.N50: n_50,
            Statistic.L50: l_50,
            Statistic.N90: n_90,
            Statistic.L90: l_90,
            Statistic.Total: sum_nucleotides,
            Statistic.PercentA: 100 * sum_a / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentT: 100 * sum_t / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentC: 100 * sum_c / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentG: 100 * sum_g / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentGC: 100 * sum_cg / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentAmbiguous: 100 * sum_ambiguous / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentMissing: 100 * sum_missing / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentMissingGaps: 100 * sum_missing_and_gaps / sum_total if sum_total else 0,
            Statistic.PercentGaps: 100 * sum_gaps / sum_total if sum_total else 0,
        })

    @staticmethod
    def _calculate_NL(counts: list[int], arg: int = 50) -> tuple[int, int]:
        if not any(counts):
            return NL(0, 0)
        counts = sorted(counts, reverse=True)
        target = sum(counts) * arg / 100
        sumsum = accumulate(counts)
        n = next((i for i, v in enumerate(sumsum) if v >= target), None)
        assert n is not None
        return NL(n + 1, counts[n])
