#!/usr/bin/env python3

from typing import Optional, Callable, Union, List

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


Agg = Union[str, Callable[[pd.Series], int]]

sequence_statistics: List[Agg] = [
    "count",
    LenBetween(max=100).count,
    LenBetween(min=101, max=300).count,
    LenBetween(min=301, max=1000).count,
    LenBetween(min=1001).count,
    OnLength("min").apply,
    OnLength("max").apply,
    OnLength("mean").apply,
    OnLength("median").apply,
    OnLength("std").apply,
    ContigStat(0.5).statN,
    ContigStat(0.5).statL,
    ContigStat(0.9).statN,
    ContigStat(0.9).statL,
]