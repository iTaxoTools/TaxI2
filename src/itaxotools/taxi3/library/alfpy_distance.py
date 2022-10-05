#!/usr/bin/env python3

import numpy as np
import pandas as pd
import gc
import alfpy.bbc as bbc
import alfpy.ncd as ncd
from alfpy.utils import seqrecords
import itertools

distances_short_names = [
    "p-distance",
    "JC distance",
    "K2P distance",
    "p-distance with gaps",
]


def alfpy_distance_array_ncd(sequences: pd.Series) -> np.ndarray:
    """
    Calculate NCD distances between `sequences`
    """
    seqs = seqrecords.SeqRecords(
        id_list=sequences.index.to_list(), seq_list=sequences.to_list()
    )
    size = len(sequences)
    dist = ncd.Distance(seqs)
    arr = np.zeros(size * size)
    for i, j in itertools.combinations(range(size), 2):
        try:
            value = dist.pairwise_distance(i, j)
        except KeyError:
            value = 0
        arr[i * size + j] = value
        arr[j * size + i] = value
    return arr


def alfpy_distance_array2_ncd(sequences1: pd.Series, sequences2: pd.Series) -> np.ndarray:
    """
    Calculates NCD distances between `sequences1` and `sequences2`

    Returns `arr` where arr[i * len(sequences2) + j] is the distance between sequences1[i] and sequences2[j]
    """
    seqs = seqrecords.SeqRecords(
        id_list=sequences1.index.to_list() + sequences2.index.to_list(),
        seq_list=sequences1.to_list() + sequences2.to_list(),
    )
    size1 = len(sequences1)
    size2 = len(sequences2)
    dist = ncd.Distance(seqs)
    arr = np.zeros(size1 * size2)
    for i, j in itertools.product(range(size1), range(size2)):
        value = dist.pairwise_distance(i, size1 + j)
        arr[i * size2 + j] = value
    return arr


def alfpy_distance_array_bbc(sequences: pd.Series) -> np.ndarray:
    """
    Calculates BBC distances between `sequences`
    """
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


def alfpy_distance_array2_bbc(
    sequences1: pd.Series, sequences2: pd.Series
) -> np.ndarray:
    """
    Calculates BBC distances between `sequences1` and `sequences2`

    Returns `arr` where arr[i * len(sequences2) + j] is the distance between sequences1[i] and sequences2[j]
    """
    seqs = seqrecords.SeqRecords(
        id_list=sequences1.index.to_list() + sequences2.index.to_list(),
        seq_list=sequences1.to_list() + sequences2.to_list(),
    )
    size1 = len(sequences1)
    size2 = len(sequences2)
    vector = bbc.create_vectors(seqs)
    dist = bbc.Distance(vector)
    arr = np.zeros(size1 * size2)
    for i, j in itertools.product(range(size1), range(size2)):
        value = dist.pairwise_distance(i, size1 + j)
        arr[i * size2 + j] = value
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


def make_alfpy_distance_table2(
    table: pd.DataFrame, reference_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Makes a table of distances between sequences in table and sequences in reference_table
    """
    distance_array = alfpy_distance_array2(
        table["sequence"], reference_table["sequence"]
    )

    # prepare the index
    index = pd.MultiIndex.from_product(
        [table.index, reference_table.index],
        names=["seqid (query 1)", "seqid (query 2)"],
    )
    # create distance table
    distance_data = {
        distance_name: distance_array for distance_name in distances_short_names
    }
    distance_table = pd.DataFrame(distance_data, index=index).reset_index()

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
