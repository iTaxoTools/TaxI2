#!/usr/bin/env python
from . import calculate_distances as calc
import math
import gc
import sys
import os
from typing import Any, TextIO

import pandas as pd
import numpy as np

from .seq import PDISTANCE, NDISTANCES
from .resources import get_resource


distances_short_names = [
    "p-distance",
    "JC distance",
    "K2P distance",
    "p-distance with gaps",
]


with open(get_resource("scores.tab")) as scores_file:
    scores_dict = {}
    for line in scores_file:
        score_name, _, val = line.partition("\t")
        try:
            scores_dict[score_name] = int(val)
        except ValueError as ex:
            raise ValueError(
                f"The value for '{score_name}' in data/scores.tab is not a number"
            ) from ex

try:
    GAP_PENALTY = scores_dict["gap penalty"]
    GAP_EXTEND_PENALTY = scores_dict["gap extend penalty"]
    END_GAP_PENALTY = scores_dict["end gap penalty"]
    END_GAP_EXTEND_PENALTY = scores_dict["end gap extend penalty"]
    MATCH_SCORE = scores_dict["match score"]
    MISMATCH_SCORE = scores_dict["mismatch score"]
except KeyError as ex:
    raise ValueError(f"'{ex.args[0]}' is missing in data/scores.tab") from ex


def make_aligner() -> Any:
    return calc.make_aligner(
        MATCH_SCORE,
        MISMATCH_SCORE,
        END_GAP_PENALTY,
        END_GAP_EXTEND_PENALTY,
        GAP_PENALTY,
        GAP_EXTEND_PENALTY,
    )


def show_alignment(aligner, target: str, query: str, file: TextIO) -> None:
    print(calc.show_alignment(aligner, target, query), file=file)


def make_distance_table(table: pd.DataFrame, already_aligned: bool) -> pd.DataFrame:
    """
    Makes the table of distances between sequences in table
    """
    if already_aligned:
        distance_array = calc.make_distance_array_aligned(
            table["sequence"], table["sequence"]
        )
    else:
        aligner = make_aligner()
        distance_array = calc.make_distance_array(
            aligner, table["sequence"], table["sequence"]
        )

    # prepare the index
    index = pd.MultiIndex.from_product(
        [table.index, table.index], names=["seqid (query 1)", "seqid (query 2)"]
    )
    # create distance table
    distance_table = pd.DataFrame(
        distance_array, index=index, columns=distances_short_names
    ).reset_index()
    distance_table = distance_table[
        distance_table["seqid (query 1)"] != distance_table["seqid (query 2)"]
    ]

    # add other columns
    table.drop(columns="sequence", inplace=True)
    distance_table = distance_table.join(table, on="seqid (query 1)")
    distance_table.rename(
        columns={col: (col + " (query 1)") for col in table.columns}, inplace=True
    )
    distance_table = distance_table.join(table, on="seqid (query 2)")
    distance_table.rename(
        columns={col: (col + " (query 2)") for col in table.columns}, inplace=True
    )

    # reorder columns
    col1 = [col for col in distance_table.columns if "query 1" in col]
    col2 = [col for col in distance_table.columns if "query 2" in col]
    distance_cols = [col for col in distance_table.columns if "query" not in col]
    distance_table = distance_table[col1 + col2 + distance_cols]
    gc.collect()
    return distance_table


def make_distance_table2(
    table: pd.DataFrame, reference_table: pd.DataFrame, already_aligned: bool
) -> pd.DataFrame:
    """
    Makes a table of distances between sequences in table and sequences in reference_table
    """
    if already_aligned:
        distance_array = calc.make_distance_array_aligned(
            table["sequence"], reference_table["sequence"]
        )
    else:
        aligner = make_aligner()
        distance_array = calc.make_distance_array(
            aligner, table["sequence"], reference_table["sequence"]
        )

    # prepare the index
    index = pd.MultiIndex.from_product(
        [table.index, reference_table.index],
        names=["seqid (query 1)", "seqid (query 2)"],
    )
    # create distance table
    distance_table = pd.DataFrame(
        distance_array, index=index, columns=distances_short_names
    ).reset_index()

    # add other columns
    table.drop(columns="sequence", inplace=True)
    reference_table = reference_table.drop(columns="sequence")
    distance_table = distance_table.join(table, on="seqid (query 1)")
    distance_table.rename(
        columns={col: (col + " (query 1)") for col in table.columns}, inplace=True
    )
    distance_table = distance_table.join(reference_table, on="seqid (query 2)")
    distance_table.rename(
        columns={col: (col + " (query 2)") for col in table.columns}, inplace=True
    )

    # reorder columns
    col1 = [col for col in distance_table.columns if "query 1" in col]
    col2 = [col for col in distance_table.columns if "query 2" in col]
    distance_cols = [col for col in distance_table.columns if "query" not in col]
    distance_table = distance_table[col1 + col2 + distance_cols]
    gc.collect()
    return distance_table
