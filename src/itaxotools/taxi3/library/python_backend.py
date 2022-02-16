#!/usr/bin/env python
import math
from multiprocessing import Pool, cpu_count
import gc

import pandas as pd
import numpy as np

from library.seq import seq_distances_ufunc, seq_distances_aligned_ufunc, PDISTANCE, NDISTANCES, make_aligner, show_alignment

distances_short_names = ['p-distance', 'JC distance',
                         'K2P distance', 'p-distance with gaps']


def _seq_par(x):
    return seq_distances_ufunc.outer(x[0], x[1])


def _seq_par_aligned(x):
    return seq_distances_aligned_ufunc.outer(x[0], x[1])


def _lines_for_data_len(tar, ref):
    _max_comparisons = 10000
    _cores = cpu_count()
    _lines = math.ceil(_max_comparisons / ref)
    _opt_split = math.ceil(tar / _cores)
    _lines = min(_lines, _opt_split)
    return _lines


def make_distance_table(table: pd.DataFrame, already_aligned: bool) -> pd.DataFrame:
    """
    Takes a series of sequences with a multi-index and returns a square dataframe

    Index and columns of the dataframe are the same as the index of the series

    The entries are arrays of pairwise distances
    """
    print('> Calculating distance table:')
    _tasks = 1
    _ufunc = _seq_par_aligned if already_aligned else _seq_par
    _len = len(table["sequence"])
    _lines = _lines_for_data_len(_len, _len)
    _sections = math.ceil(_len / _lines)
    print('> Number of Chunks:', _sections)
    chunks = np.array_split(np.asarray(table["sequence"]), _sections)
    ref_table = np.asarray(table["sequence"])
    x = [(i, ref_table) for i in chunks]
    res_list = []
    with Pool(maxtasksperchild=_tasks) as p:
        print('> {} / {}'.format(0, _sections))
        res = p.imap(_ufunc, x)
        for i, res in enumerate(res):
            print('> {} / {}'.format(i+1, _sections))
            res_list.append(res)
            # print(res) # Write this to a temp file
    print('> Merging Chunks:', _sections)
    distance_array = np.concatenate(res_list)

    # prepare indices
    seqid1 = table.index.copy()
    seqid1.name = "seqid (query 1)"
    seqid2 = table.index.copy()
    seqid2.name = "seqid (query 2)"
    # create distance table
    distance_table = pd.DataFrame(
        distance_array, index=seqid1, columns=seqid2).stack().reset_index(name="distances")
    distance_table = distance_table[distance_table["seqid (query 1)"]
                                    != distance_table["seqid (query 2)"]]

    # add other columns
    table.drop(columns="sequence", inplace=True)
    distance_table = distance_table.join(table, on="seqid (query 1)")
    distance_table.rename(columns={col: (col + " (query 1)")
                                   for col in table.columns}, inplace=True)
    distance_table = distance_table.join(table, on="seqid (query 2)")
    distance_table.rename(columns={col: (col + " (query 2)")
                                   for col in table.columns}, inplace=True)

    # reorder columns
    col1 = [col for col in distance_table.columns if "query 1" in col]
    col2 = [col for col in distance_table.columns if "query 2" in col]
    distance_table = distance_table[col1 + col2 + ["distances"]]
    distance_table[distances_short_names] = pd.DataFrame(
        distance_table["distances"].to_list(), index=distance_table.index)
    distance_table.drop(columns="distances", inplace=True)
    gc.collect()
    return distance_table


def make_distance_table2(table: pd.DataFrame, reference_table: pd.DataFrame, already_aligned: bool) -> pd.DataFrame:
    """
    Takes a series of sequences with a multi-index and returns a square dataframe

    Index and columns of the dataframe are the same as the index of the series

    The entries are arrays of pairwise distances
    """
    print('> Calculating distance table:')
    _tasks = 1
    _map_chunksize = 1
    _ufunc = _seq_par_aligned if already_aligned else _seq_par
    _len_tar = len(table["sequence"])
    _len_ref = len(reference_table["sequence"])
    _lines = _lines_for_data_len(_len_tar, _len_ref)
    _sections = math.ceil(_len_tar / _lines)
    print('> Number of Chunks:', _sections)
    chunks = np.array_split(np.asarray(table["sequence"]), _sections)
    ref_table = np.asarray(reference_table["sequence"])
    x = [(i, ref_table) for i in chunks]
    res_list = []
    with Pool(maxtasksperchild=_tasks) as p:
        print('> {} / {}'.format(0, _sections))
        res = p.imap(_ufunc, x, chunksize=_map_chunksize)
        for i, res in enumerate(res):
            print('> {} / {}'.format(i+1, _sections))
            res_list.append(res)
            # print(res) # Write this to a temp file
    print('> Merging Chunks:', _sections)
    distance_array = np.concatenate(res_list)

    # prepare indices
    seqid1 = table.index.copy()
    seqid1.name = "seqid (query 1)"
    seqid2 = reference_table.index.copy()
    seqid2.name = "seqid (query 2)"
    # create distance table
    distance_table = pd.DataFrame(
        distance_array, index=seqid1, columns=seqid2).stack().reset_index(name="distances")
    distance_table = distance_table[distance_table["seqid (query 1)"]
                                    != distance_table["seqid (query 2)"]]

    # add other columns
    table.drop(columns="sequence", inplace=True)
    reference_table = reference_table.drop(columns="sequence")
    distance_table = distance_table.join(table, on="seqid (query 1)")
    distance_table.rename(columns={col: (col + " (query 1)")
                                   for col in table.columns}, inplace=True)
    distance_table = distance_table.join(reference_table, on="seqid (query 2)")
    distance_table.rename(columns={col: (col + " (query 2)")
                                   for col in reference_table.columns}, inplace=True)

    # reorder columns
    col1 = [col for col in distance_table.columns if "query 1" in col]
    col2 = [col for col in distance_table.columns if "query 2" in col]
    distance_table = distance_table[col1 + col2 + ["distances"]]
    distance_table[distances_short_names] = pd.DataFrame(
        distance_table["distances"].to_list(), index=distance_table.index)
    distance_table.drop(columns="distances", inplace=True)
    gc.collect()
    return distance_table
