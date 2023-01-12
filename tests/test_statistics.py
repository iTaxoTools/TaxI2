from __future__ import annotations

from pathlib import Path
from sys import stderr
from typing import Callable, NamedTuple
from math import sqrt

import pytest

from itaxotools.taxi3.statistics import Counts, Statistic, Statistics


class CountTest(NamedTuple):
    count: str
    fixed: int
    sequence: str

    def validate(self):
        counts = Counts.from_sequence(self.sequence)
        tested = getattr(counts, self.count)
        assert tested == self.fixed


class StatisticTest(NamedTuple):
    stat: Statistic
    fixed: object
    sequences: list[str]
    precision: float = 0.00051

    def validate(self):
        stats = Statistics.from_sequences(self.sequences)
        tested = stats[self.stat]
        if isinstance(tested, float):
            assert abs(tested - self.fixed) <= self.precision
        else:
            assert tested == self.fixed


count_tests = [
    CountTest('total', 0, ''),
    CountTest('total', 4, 'ATCG'),
    CountTest('total', 1, 'N'),
    CountTest('total', 1, '-'),
    CountTest('nucleotides', 4, 'ATCG'),
    CountTest('nucleotides', 1, 'N'),
    CountTest('nucleotides', 0, '-'),
    CountTest('missing', 0, 'ATCG'),
    CountTest('missing', 1, 'N'),
    CountTest('missing', 0, '-'),
    CountTest('gaps', 0, 'ATCG'),
    CountTest('gaps', 0, 'N'),
    CountTest('gaps', 1, '-'),
    CountTest('a', 1, 'ATCG'),
    CountTest('t', 1, 'ATCG'),
    CountTest('c', 1, 'ATCG'),
    CountTest('g', 1, 'ATCG'),
]


statistic_tests = [
    StatisticTest(Statistic.SequenceCount, 0, []),
    StatisticTest(Statistic.SequenceCount, 1, ['ATCG']),
    StatisticTest(Statistic.SequenceCount, 7, ['ATCG'] * 7),

    StatisticTest(Statistic.NucleotideCount, 0, ['']),
    StatisticTest(Statistic.NucleotideCount, 0, ['-']),
    StatisticTest(Statistic.NucleotideCount, 4, ['ATCG']),
    StatisticTest(Statistic.NucleotideCount, 6, ['ATCG', 'AT']),
    StatisticTest(Statistic.NucleotideCount, 6, ['ATCG', 'AT--']),
    StatisticTest(Statistic.NucleotideCount, 6, ['ATCG', 'NN--']),

    StatisticTest(Statistic.BP_0, 1, ['']),
    StatisticTest(Statistic.BP_0, 7, [''] * 7),
    StatisticTest(Statistic.BP_0, 1, ['-']),
    StatisticTest(Statistic.BP_0, 0, ['N']),
    StatisticTest(Statistic.BP_0, 0, ['A']),

    StatisticTest(Statistic.BP_1_100, 0, ['A' * 0]),
    StatisticTest(Statistic.BP_1_100, 1, ['A' * 1]),
    StatisticTest(Statistic.BP_1_100, 1, ['A' * 100]),
    StatisticTest(Statistic.BP_1_100, 0, ['A' * 101]),

    StatisticTest(Statistic.BP_101_300, 0, ['A' * 100]),
    StatisticTest(Statistic.BP_101_300, 1, ['A' * 101]),
    StatisticTest(Statistic.BP_101_300, 1, ['A' * 300]),
    StatisticTest(Statistic.BP_101_300, 0, ['A' * 301]),

    StatisticTest(Statistic.BP_301_1000, 0, ['A' * 300]),
    StatisticTest(Statistic.BP_301_1000, 1, ['A' * 301]),
    StatisticTest(Statistic.BP_301_1000, 1, ['A' * 1000]),
    StatisticTest(Statistic.BP_301_1000, 0, ['A' * 1001]),

    StatisticTest(Statistic.BP_1001_plus, 0, ['A' * 1000]),
    StatisticTest(Statistic.BP_1001_plus, 1, ['A' * 1001]),
    StatisticTest(Statistic.BP_1001_plus, 1, ['A' * 9000]),

    StatisticTest(Statistic.BP_1_100, 0, ['-' * 100]),
    StatisticTest(Statistic.BP_101_300, 0, ['-' * 300]),
    StatisticTest(Statistic.BP_301_1000, 0, ['-' * 1000]),
    StatisticTest(Statistic.BP_1001_plus, 0, ['-' * 1001]),

    StatisticTest(Statistic.Minimum, 2, ['AT--', 'ATC', 'ATCG']),
    StatisticTest(Statistic.Maximum, 4, ['AT', 'ATC', 'ATCG', 'ATCG--']),

    StatisticTest(Statistic.Median, 3, ['AT', 'ATC', 'ATCG']),
    StatisticTest(Statistic.Median, 3, ['AT', 'ATCG']),

    StatisticTest(Statistic.Stdev, 0.0, ['ATCG', 'ATGC']),
    StatisticTest(Statistic.Stdev, 0.0, ['ATCG', 'ATGC--']),
    StatisticTest(Statistic.Stdev, sqrt(2), ['ATCG', 'ATATAT']),

    StatisticTest(Statistic.PercentA, 0.0, ['']),
    StatisticTest(Statistic.PercentA, 1.0, ['A']),
    StatisticTest(Statistic.PercentA, 1.0, ['A-']),
    StatisticTest(Statistic.PercentA, 0.5, ['AT']),
    StatisticTest(Statistic.PercentA, 0.5, ['AN']),

    StatisticTest(Statistic.PercentT, 0.0, ['']),
    StatisticTest(Statistic.PercentT, 1.0, ['T']),
    StatisticTest(Statistic.PercentT, 1.0, ['T-']),
    StatisticTest(Statistic.PercentT, 0.5, ['TA']),
    StatisticTest(Statistic.PercentT, 0.5, ['TN']),

    StatisticTest(Statistic.PercentC, 0.0, ['']),
    StatisticTest(Statistic.PercentC, 1.0, ['C']),
    StatisticTest(Statistic.PercentC, 1.0, ['C-']),
    StatisticTest(Statistic.PercentC, 0.5, ['CT']),
    StatisticTest(Statistic.PercentC, 0.5, ['CN']),

    StatisticTest(Statistic.PercentG, 0.0, ['']),
    StatisticTest(Statistic.PercentG, 1.0, ['G']),
    StatisticTest(Statistic.PercentG, 1.0, ['G-']),
    StatisticTest(Statistic.PercentG, 0.5, ['GT']),
    StatisticTest(Statistic.PercentG, 0.5, ['GN']),

    StatisticTest(Statistic.PercentGC, 0.0, ['']),
    StatisticTest(Statistic.PercentGC, 1.0, ['G']),
    StatisticTest(Statistic.PercentGC, 1.0, ['C']),
    StatisticTest(Statistic.PercentGC, 1.0, ['GC']),
    StatisticTest(Statistic.PercentGC, 1.0, ['GC-']),
    StatisticTest(Statistic.PercentGC, 0.5, ['ATCG']),
    StatisticTest(Statistic.PercentGC, 0.25, ['ATTG']),
    StatisticTest(Statistic.PercentGC, 0.25, ['ATTC']),

    StatisticTest(Statistic.PercentAmbiguous, 0.0, ['']),
    StatisticTest(Statistic.PercentAmbiguous, 1.0, ['RYSWKM']),
    StatisticTest(Statistic.PercentAmbiguous, 1.0, ['RYSWKM------']),
    StatisticTest(Statistic.PercentAmbiguous, 0.5, ['RYSWKMATATAT']),
    StatisticTest(Statistic.PercentAmbiguous, 0.5, ['RYSWKMNNNNNN']),
    StatisticTest(Statistic.PercentAmbiguous, 0.5, ['AY']),

    StatisticTest(Statistic.PercentMissing, 0.0, ['']),
    StatisticTest(Statistic.PercentMissing, 1.0, ['N']),
    StatisticTest(Statistic.PercentMissing, 1.0, ['N-']),
    StatisticTest(Statistic.PercentMissing, 0.5, ['NA']),
    StatisticTest(Statistic.PercentMissing, 0.0, ['ATCG']),
    StatisticTest(Statistic.PercentMissing, 0.0, ['ATCGY']),

    StatisticTest(Statistic.PercentMissingGaps, 0.0, ['']),
    StatisticTest(Statistic.PercentMissingGaps, 1.0, ['N']),
    StatisticTest(Statistic.PercentMissingGaps, 1.0, ['-']),
    StatisticTest(Statistic.PercentMissingGaps, 1.0, ['N-']),
    StatisticTest(Statistic.PercentMissingGaps, 0.5, ['NA']),
    StatisticTest(Statistic.PercentMissingGaps, 0.5, ['-A']),
    StatisticTest(Statistic.PercentMissingGaps, 0.5, ['ATCG--NN']),

    StatisticTest(Statistic.PercentGaps, 0.0, ['']),
    StatisticTest(Statistic.PercentGaps, 1.0, ['-']),
    StatisticTest(Statistic.PercentGaps, 0.5, ['A-']),
    StatisticTest(Statistic.PercentGaps, 0.5, ['N-']),
    StatisticTest(Statistic.PercentGaps, 0.5, ['Y-']),
    StatisticTest(Statistic.PercentGaps, 0.5, ['ATCG----']),

    StatisticTest(Statistic.N50, 0, ['']),
    StatisticTest(Statistic.L50, 0, ['']),
    StatisticTest(Statistic.N90, 0, ['']),
    StatisticTest(Statistic.L90, 0, ['']),

    StatisticTest(Statistic.N50, 4, ['AT', 'ATCG']),
    StatisticTest(Statistic.L50, 1, ['AT', 'ATCG']),

    StatisticTest(Statistic.N50, 4, ['ATCG', 'AT']),
    StatisticTest(Statistic.L50, 1, ['ATCG', 'AT']),

    StatisticTest(Statistic.N50, 4, ['ATCG', 'ATCG']),
    StatisticTest(Statistic.L50, 1, ['ATCG', 'ATCG']),

    StatisticTest(Statistic.N90, 4, ['ATCG', 'ATCG']),
    StatisticTest(Statistic.L90, 2, ['ATCG', 'ATCG']),

    StatisticTest(Statistic.N90, 4, ['ATCG', 'ATCG----']),
    StatisticTest(Statistic.L90, 2, ['ATCG', 'ATCG----']),

    # 9 contigs with lengths from 2 to 10
    StatisticTest(Statistic.N50, 8, list(map(lambda l: 'A' * l, range(2, 11)))),
    StatisticTest(Statistic.L50, 3, list(map(lambda l: 'A' * l, range(2, 11)))),
    StatisticTest(Statistic.N50, 8, list(map(lambda l: 'A' * l, range(2, 11))) + ['-']),
    StatisticTest(Statistic.L50, 3, list(map(lambda l: 'A' * l, range(2, 11))) + ['-']),
    StatisticTest(Statistic.N90, 4, list(map(lambda l: 'A' * l, range(2, 11)))),
    StatisticTest(Statistic.L90, 7, list(map(lambda l: 'A' * l, range(2, 11)))),

]


@pytest.mark.parametrize("test", count_tests)
def test_counts(test: CountTest) -> None:
    test.validate()


@pytest.mark.parametrize("test", statistic_tests)
def test_statistics(test: StatisticTest) -> None:
    test.validate()
