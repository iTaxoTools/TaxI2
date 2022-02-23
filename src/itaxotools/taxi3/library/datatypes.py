#!/usr/bin/env python3
from __future__ import annotations

from typing import List, Set
from abc import ABC, abstractclassmethod
from pathlib import Path
from dataclasses import dataclass
import re

import pandas as pd
import openpyxl


class _Header:

    RENAME_DICT = {
        "organism": "species",
        "specimen_identifier": "specimen_voucher",
        "specimen identifier": "specimen_voucher",
        "specimen voucher": "specimen_voucher",
    }

    def __init__(self, header: List[str]):
        self.names = list(map(_Header.rename, map(_Header.to_ascii, header)))

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


class InvalidPath(Exception):
    pass


@dataclass
class SequenceData(DataType):
    def __init__(self, path: Path, protocol: SequenceReader):
        """
        raises `InvalidPath` if `path` does not exists
        """
        if not path.exists() or not path.is_file():
            raise InvalidPath
        self.path = path
        self.protocol = protocol

    def get_dataframe(self) -> pd.DataFrame:
        return self.protocol.read(self.path, columns=["seqid", "sequence"])


class SequenceReader(ABC):
    @abstractclassmethod
    @staticmethod
    def read(path: Path, *, columns: List[str]) -> pd.DataFrame:
        pass


class ColumnsNotFound(Exception):
    def __init__(self, columns: Set[str]):
        if not isinstance(columns, set):
            raise TypeError("ColumnsNotFound argument should be a set")
        self.args = (columns,)


class TabfileReader(SequenceReader):
    @staticmethod
    def read(path: Path, *, columns: List[str]) -> pd.DataFrame:
        with open(path, errors="replace") as file:
            header = _Header(file.readline().rstrip("\n").split("\t"))
        if not set(columns) in set(header.names):
            raise ColumnsNotFound(set(columns) - set(header.names))
        return pd.read_table(
            path,
            header=0,
            names=header.names,
            usecols=columns,
            na_filter=False,
            dtype=str,
            encoding_errors="replace",
        )


class XlsxReader(SequenceReader):
    @staticmethod
    def read(path: Path, *, columns: List[str]) -> pd.DataFrame:
        worksheet = openpyxl.load_workbook(path).worksheets[0]
        header = _Header(list(map(lambda cell: cell.value, next(worksheet.rows))))
        if not set(columns) in set(header.names):
            raise ColumnsNotFound(set(columns) - set(header.names))
        return pd.read_excel(
            path,
            sheet_name=0,
            header=0,
            names=header.names,
            usecols=columns,
            na_filter=False,
            dtype=str,
            engine="openpyxl",
            encoding_errors="replace",
        )
