#!/usr/bin/env python3

import numpy as np
import pandas as pd
import gc
import alfpy.bbc as bbc
from alfpy.utils import seqrecords
import itertools

distances_short_names = [
    "p-distance",
    "JC distance",
    "K2P distance",
    "p-distance with gaps",
]


def alfpy_distance_array(sequences: pd.Series) -> np.array:
    seqs = seqrecords.SeqRecords(
        id_list=sequences.index.to_list(), seq_list=sequences.to_list()
    )
    size = len(sequences)
    vector = bbc.create_vectors(seqs)
    dist = bbc.Distance(vector)
    arr = np.zeros(size * size)
    for i, j in itertools.combinations(range(size), 2):
        value = dist.pairwise_distance(i, j)
        arr[i * size + j] = value
        arr[j * size + i] = value
    return arr


def make_alfpy_distance_table(table: pd.DataFrame) -> pd.DataFrame:
    """
    Makes the table of distances between sequences in table
    """
    distance_array = alfpy_distance_array(table["sequence"])

    # prepare the index
    index = pd.MultiIndex.from_product(
        [table.index, table.index], names=["seqid (query 1)", "seqid (query 2)"]
    )
    # create distance table
    distance_data = {
        distance_name: distance_array for distance_name in distances_short_names
    }
    distance_table = pd.DataFrame(distance_data, index=index).reset_index()
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
