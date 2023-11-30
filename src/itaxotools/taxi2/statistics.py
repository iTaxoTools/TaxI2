from __future__ import annotations

import statistics
from collections import Counter
from enum import Enum
from itertools import accumulate
from math import inf, isinf
from pathlib import Path
from typing import Generator, Literal, NamedTuple

from .handlers import FileHandler, ReadHandle, WriteHandle
from .types import Percentage


class Counts(NamedTuple):
    total: int
    nucleotides: int
    missing: int
    gaps: int
    a: int
    c: int
    g: int
    t: int

    @classmethod
    def from_sequence(cls, seq: str) -> Counts:
        counter = Counter(seq)
        return cls(
            total=len(seq),
            nucleotides=len(seq) - counter["-"],
            missing=counter["N"],
            gaps=counter["-"],
            a=counter["A"],
            c=counter["C"],
            g=counter["G"],
            t=counter["T"],
        )


class NL(NamedTuple):
    N: int
    L: int


class Statistic(Enum):
    """Defines statistic labels & types. Order matters."""

    Group = "Group", str
    SequenceCount = "Total number of sequences", int
    NucleotideCount = "Total length of all sequences ", int
    BP_0 = "Number of sequences with 0 bp", int
    BP_1_100 = "Number of sequences with less than 100 bp", int
    BP_101_300 = "Number of sequences between 101-300 bp", int
    BP_301_1000 = "Number of sequences between 301-1000 bp", int
    BP_1001_plus = "Number of sequences with more than 1000 bp", int
    Minimum = "Minimum sequence length", int
    Maximum = "Maximum sequence length ", int
    Mean = "Mean sequence length  ", float
    Median = "Median sequence length  ", float
    Stdev = "Standard deviation of sequence length", float
    PercentA = "Percentage of base A", Percentage
    PercentC = "Percentage of base C", Percentage
    PercentG = "Percentage of base G", Percentage
    PercentT = "Percentage of base T", Percentage
    PercentGC = "GC content", Percentage
    PercentAmbiguous = "Percentage of ambiguity codes", Percentage
    PercentMissing = "Percentage of missing data ", Percentage
    PercentMissingGaps = "Percentage of missing data including gaps", Percentage
    PercentGaps = "Percentage of gaps", Percentage
    N50 = "N50 statistic", int
    L50 = "L50 statistic", int
    N90 = "N90 statistic", int
    L90 = "L90 statistic", int

    def __init__(self, label, type):
        self.label = label
        self.type = type

    def __repr__(self):
        return f"<{type(self).__name__}.{self._name_}>"

    def __str__(self):
        return self.label


class Statistics(dict[Statistic, object]):
    """Keep Enum order, convert values to the proper type"""

    def __init__(self, stats: dict[Statistic, object]):
        super().__init__({s: s.type(stats[s]) for s in Statistic if s in stats})

    @classmethod
    def from_sequences(self, sequences: iter[str], group: str = None) -> Statistics:
        calc = StatisticsCalculator(sequences, group)
        return calc.calculate()


class StatisticsCalculator:
    """Calculate statistics from sequences"""

    def __init__(self, sequences: iter[str] = [], group: str = None):
        self.group = group
        self.it = self.iter()
        next(self.it)

        for seq in sequences:
            self.add(seq)

    def add(self, seq: str) -> None:
        """Add a sequence to calculations"""
        self.it.send(seq)

    def calculate(self) -> Statistics:
        """Call only once to get final statistics"""
        result = self.it.send(None)
        if self.group:
            result[Statistic.Group] = self.group
        return Statistics(result)

    def iter(self) -> Generator[None | Statistics, str | None, None]:
        """Send sequences to update internal state,
        send None to get statistics"""
        nucleotide_counts = []

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

        sequence = yield
        while sequence is not None:
            count = Counts.from_sequence(sequence)
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

            sequence = yield

        length = len(nucleotide_counts)
        mean = sum_nucleotides / length if length else 0
        median = statistics.median(nucleotide_counts) if length else 0
        stdev = (
            statistics.pstdev(nucleotide_counts) if len(nucleotide_counts) > 1 else 0
        )

        sum_cg = sum_c + sum_g
        sum_ambiguous = sum_nucleotides - sum_missing - sum_a - sum_t - sum_c - sum_g
        sum_missing_and_gaps = sum_missing + sum_gaps

        n_50, l_50 = self._calculate_NL(nucleotide_counts, 50)
        n_90, l_90 = self._calculate_NL(nucleotide_counts, 90)

        yield {
            Statistic.SequenceCount: length,
            Statistic.NucleotideCount: sum_nucleotides,
            Statistic.BP_0: bp_0,
            Statistic.BP_1_100: bp_1_100,
            Statistic.BP_101_300: bp_101_300,
            Statistic.BP_301_1000: bp_301_1000,
            Statistic.BP_1001_plus: bp_1001_plus,
            Statistic.Minimum: minimum if not isinf(minimum) else 0,
            Statistic.Maximum: maximum if not isinf(maximum) else 0,
            Statistic.Mean: mean,
            Statistic.Median: median,
            Statistic.Stdev: stdev,
            Statistic.PercentA: sum_a / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentC: sum_c / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentG: sum_g / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentT: sum_t / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentGC: sum_cg / sum_nucleotides if sum_nucleotides else 0,
            Statistic.PercentAmbiguous: sum_ambiguous / sum_nucleotides
            if sum_nucleotides
            else 0,
            Statistic.PercentMissing: sum_missing / sum_nucleotides
            if sum_nucleotides
            else 0,
            Statistic.PercentMissingGaps: sum_missing_and_gaps / sum_total
            if sum_total
            else 0,
            Statistic.PercentGaps: sum_gaps / sum_total if sum_total else 0,
            Statistic.N50: n_50,
            Statistic.L50: l_50,
            Statistic.N90: n_90,
            Statistic.L90: l_90,
        }

    @staticmethod
    def _calculate_NL(counts: list[int], arg: int = 50) -> tuple[int, int]:
        if not any(counts):
            return NL(0, 0)
        counts = sorted(counts, reverse=True)
        target = sum(counts) * arg / 100
        sumsum = accumulate(counts)
        pos = next((i for i, v in enumerate(sumsum) if v >= target), None)
        assert pos is not None
        return NL(counts[pos], pos + 1)


class StatisticsHandler(FileHandler[Statistics]):
    def _open(
        self,
        path: Path,
        mode: Literal["r", "w"] = "w",
        float_formatter: str = "{:f}",
        percentage_formatter: str = "{:f}",
        percentage_multiply: bool = False,
        *args,
        **kwargs,
    ):
        self.formatters = {}
        self.formatters[float] = float_formatter
        self.formatters[Percentage] = percentage_formatter
        self.percentage_multiply = percentage_multiply
        super()._open(path, mode, *args, **kwargs)

    def _iter_read(self) -> ReadHandle[Statistics]:
        raise NotImplementedError()

    def statisticToText(self, value):
        if isinstance(value, Percentage) and self.percentage_multiply:
            value = Percentage(value * 100)
        formatter = self.formatters.get(type(value), "{}")
        return formatter.format(value)


class Single(StatisticsHandler):
    def _iter_write(self) -> WriteHandle[Statistics]:
        with FileHandler.Tabfile(self.path, "w") as file:
            try:
                stats = yield
                for stat, value in stats.items():
                    file.write((str(stat), self.statisticToText(value)))
                yield
                raise Exception("Can only write a single statistics instance")
            except GeneratorExit:
                return


class Groups(StatisticsHandler):
    def _open(
        self,
        path: Path,
        mode: Literal["r", "w"] = "w",
        group_name: str = "group",
        *args,
        **kwargs,
    ):
        self.group_name = group_name
        super()._open(path, mode, *args, **kwargs)

    def _iter_write(self) -> WriteHandle[Statistics]:
        self.wrote_headers = False

        with FileHandler.Tabfile(self.path, "w") as file:
            try:
                stats = yield
                self._check_stats(stats)
                self._write_headers(file, stats)
                self._write_stats(file, stats)
                while True:
                    stats = yield
                    self._check_stats(stats)
                    self._write_stats(file, stats)
            except GeneratorExit:
                return

    def _check_stats(self, stats: Statistics):
        if Statistic.Group not in stats:
            raise Exception("Statistics must contain a group name")

    def _write_headers(self, file: FileHandler, stats: Statistics):
        if self.wrote_headers:
            return
        stats = [str(stat) for stat in stats]
        out = (self.group_name, *stats[1:])
        file.write(out)
        self.wrote_headers = True

    def _write_stats(self, file: FileHandler, stats: Statistics):
        stats = [self.statisticToText(value) for value in stats.values()]
        file.write((*stats,))
