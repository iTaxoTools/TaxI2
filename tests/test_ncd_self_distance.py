#!/usr/bin/env python3

import pytest

import alfpy.bbc as bbc
import alfpy.ncd as ncd
from alfpy.utils import seqrecords

SEQUENCE = "ATAATCGGCGGATTTGGAAACTGATTAGTTCCCCTAATAATTGGAGCACCAGACATGGCCTTCCCACGAATAAATAACATAAGCTTCTGACTACTCCCACCATCCTTTCTCCTTCTATTAGCCTCTTCCGCAGTGGAAGCCGGGGCCGGAACAGGTTGAACTGTTTATCCCCCTCTAGCTGGCAACCTTGCACACGCTGGCCCATCCGTTGACCTAACAATTTTCTCCCTACACCTGGCTGGAGTCTCATCAATTTTAGGTGCAATTAATTTTATTACCACTATTATTAATATAAAACCCCCATCAGTCACCCAATACCAAACACCCCTTTTCGTCTGATCTGTGCTAATTACAGCAGTACTTCTACTCCTCTCCCTCCCTGTCC"


@pytest.mark.legacy
def test_bbc_self_distance() -> None:
    """
    Tests reflexivity of BBC distance
    """
    seqs = seqrecords.SeqRecords(id_list=[0], seq_list=[SEQUENCE])
    vectors = bbc.create_vectors(seqs)
    dist = bbc.Distance(vectors)
    assert dist.pairwise_distance(0, 0) == 0


@pytest.mark.legacy
def test_ncd_self_distance() -> None:
    """
    Asserts that the bug in NCD distance has not been fixed
    """
    # list of two copies of SEQUENCE
    seqs = seqrecords.SeqRecords(id_list=[0], seq_list=[SEQUENCE])
    dist = ncd.Distance(seqs)
    with pytest.raises(KeyError):
        # Distance between the sequence and itself raises an exception
        dist.pairwise_distance(0, 0) == 0
