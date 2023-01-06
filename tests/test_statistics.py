from __future__ import annotations

from pathlib import Path
from sys import stderr
from typing import Callable, NamedTuple

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

    StatisticTest(Statistic.BP_0_100, 1, ['A' * 100]),
    StatisticTest(Statistic.BP_0_100, 0, ['A' * 300]),
    StatisticTest(Statistic.BP_0_100, 0, ['A' * 1000]),
    StatisticTest(Statistic.BP_0_100, 0, ['A' * 1001]),

    StatisticTest(Statistic.BP_101_300, 0, ['A' * 100]),
    StatisticTest(Statistic.BP_101_300, 1, ['A' * 300]),
    StatisticTest(Statistic.BP_101_300, 0, ['A' * 1000]),
    StatisticTest(Statistic.BP_101_300, 0, ['A' * 1001]),

    StatisticTest(Statistic.BP_301_1000, 0, ['A' * 100]),
    StatisticTest(Statistic.BP_301_1000, 0, ['A' * 300]),
    StatisticTest(Statistic.BP_301_1000, 1, ['A' * 1000]),
    StatisticTest(Statistic.BP_301_1000, 0, ['A' * 1001]),

    StatisticTest(Statistic.BP_1001_plus, 0, ['A' * 100]),
    StatisticTest(Statistic.BP_1001_plus, 0, ['A' * 300]),
    StatisticTest(Statistic.BP_1001_plus, 0, ['A' * 1000]),
    StatisticTest(Statistic.BP_1001_plus, 1, ['A' * 1001]),

    # StatisticTest(Statistic.BP_0_100, 0, ['-' * 100]),
    # StatisticTest(Statistic.BP_101_300, 0, ['-' * 300]),
    # StatisticTest(Statistic.BP_301_1000, 0, ['-' * 1000]),
    # StatisticTest(Statistic.BP_1001_plus, 0, ['-' * 1001]),
]


@pytest.mark.parametrize("test", count_tests)
def test_counts(test: CountTest) -> None:
    test.validate()


@pytest.mark.parametrize("test", statistic_tests)
def test_statistics(test: StatisticTest) -> None:
    test.validate()
