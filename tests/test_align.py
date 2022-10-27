from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.align import PairwiseAligner
from itaxotools.taxi3.pairs import (
    SequencePair, SequencePairFile, SequencePairs)
from itaxotools.taxi3.sequences import Sequence, Sequences

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


@pytest.mark.xfail
def test_biopair() -> None:
    x = Sequence('idx', 'ATC')
    y = Sequence('idy', 'TCA')
    p = SequencePair(x, y)
    al = PairwiseAligner.Biopython()
    ap = al.align(p)
    print(ap)
    assert False
