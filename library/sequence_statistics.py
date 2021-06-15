#!/usr/bin/env python3

from typing import Optional, Callable, Union, Dict

import pandas as pd


class LenBetween:
    """
    Provides count method

    LenBeetween(min, max).count(col) returns the number of of elements of col with length between min and max (inclusive)
    """

    def __init__(self, min: Optional[int] = None, max: Optional[int] = None):
        self.min = min
        self.max = max

    def count(self, col: pd.Series) -> int:
        if not self.min and not self.max:
            return col.count()
        lengths = col.str.len()
        if self.min:
            over_min = lengths >= self.min
        else:
            over_min = True
        if self.max:
            over_max = lengths <= self.max
        else:
            over_max = True
        return col.loc[over_min & over_max].count()


class OnLength:
    """
    OnLength('fun').apply(col) returns the result of aggregator 'fun' on the length of strings in col
    """

    def __init__(self, agg: str) -> None:
        self.agg = agg

    def apply(self, col: pd.Series) -> int:
        return col.str.len().agg(self.agg)


class ContigStat:
    """
    ContigStat(i).statN calculates N{i} statistic
    ContigStat(i).statL calculates L{i} statistic
    """

    def __init__(self, arg: float):
        self.arg = arg

    def statN(self, col: pd.Series) -> int:
        col = col.str.len()
        total_length = col.sum()
        col = col.sort_values(ascending=False)
        i_from_end = len(col.loc[col.cumsum() >= total_length * self.arg])
        return col.iat[-i_from_end]

    def statL(self, col: pd.Series) -> int:
        col = col.str.len()
        total_length = col.sum()
        col = col.sort_values(ascending=False)
        return len(col.loc[col.cumsum() < total_length * self.arg]) + 1


class Percentage:
    """Percentage(regexp).count(col) return the percentage of characters matching regexp in col"""

    def __init__(self, regexp: str):
        self.regexp = regexp

    def count(self, col: pd.Series) -> float:
        total_length = col.str.len().sum()
        matches = col.str.count(self.regexp).sum()
        return matches / total_length


Agg = Union[str, Callable[[pd.Series], int], Callable[[pd.Series], float]]

sequence_statistics: Dict[str, Agg] = {
    "Total number of sequences": "count",
    "Number of sequences of less than 100 bp": LenBetween(max=100).count,
    "Number of sequences of 101-300 bp": LenBetween(min=101, max=300).count,
    "Number of sequences of 301-1000 bp": LenBetween(min=301, max=1000).count,
    "Number of sequences of greater than 1000 bp": LenBetween(min=1001).count,
    "Minimum of sequence length": OnLength("min").apply,
    "Maximum of sequence length": OnLength("max").apply,
    "Mean of sequence length": OnLength("mean").apply,
    "Median of sequence length": OnLength("median").apply,
    "Standard deviation of sequence length": OnLength("std").apply,
    "N50 statistic": ContigStat(0.5).statN,
    "L50 statistic": ContigStat(0.5).statL,
    "N90 statistic": ContigStat(0.9).statN,
    "L90 statistic": ContigStat(0.9).statL,
    "Total length of all sequences": OnLength("sum").apply,
    "Percentage of base A": Percentage("[Aa]").count,
    "Percentage of base C": Percentage("[Cc]").count,
    "Percentage of base G": Percentage("[Gg]").count,
    "Percentage of base T": Percentage("[Tt]").count,
    "GC content": Percentage("[GCgc]").count,
    "Percentage of ambiguity codes": Percentage("[RYSWKMryswkm]").count,
    "Percentage of missing data": Percentage("[Nn?]").count,
}

sequence_statistics_with_gaps: Dict[str, Agg] = {
    "Percentage of missing data including gaps": Percentage("[Nn?-]").count
}
