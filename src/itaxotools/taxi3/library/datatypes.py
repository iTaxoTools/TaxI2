#!/usr/bin/env python3
from __future__ import annotations

from typing import List
from abc import ABC, abstractclassmethod
from pathlib import Path
from dataclasses import dataclass
import re

import pandas as pd


class _Header:

    RENAME_DICT = {
        "organism": "species",
        "specimen_identifier": "specimen_voucher",
        "specimen identifier": "specimen_voucher",
        "specimen voucher": "specimen_voucher",
    }

    def __init__(self, header: List[str]):
        self.header = list(map(_Header.rename, map(_Header.to_ascii, header)))

    @staticmethod
    def to_ascii(name: str) -> str:
        regex = r"([a-zA-Z0-9_-]+)"
        return "".join(re.findall(regex, name))

    @staticmethod
    def rename(name: str) -> str:
        if "sequence" in name:
            return "sequence"
        else:
            try:
                return _Header.RENAME_DICT[name]
            except KeyError:
                return name


class DataType(ABC):
    pass


@dataclass
class SequenceData(DataType):
    path: Path
    protocol: SequenceReader

    def get_dataframe(self) -> pd.DataFrame:
        return self.protocol.read(self.path, columns=["seqid", "sequence"])


class SequenceReader(ABC):
    @abstractclassmethod
    @staticmethod
    def read(path: Path, *, columns: List[str]) -> pd.DataFrame:
        pass
