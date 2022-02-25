#!/usr/bin/env python3
from __future__ import annotations

from typing import List, Set, Any
from abc import ABC, abstractclassmethod
from pathlib import Path
from dataclasses import dataclass
import re

import pandas as pd
import openpyxl

from itaxotools.DNAconvert.library.fasta import Fastafile
from itaxotools.DNAconvert.library.genbank import GenbankFile
from itaxotools.DNAconvert.library.record import Record


class InvalidPath(Exception):
    pass


class ValidFilePath(Path):
    """
    Subclass of Path that contains a path to an existing file
    """

    def __init__(self, *pathsegments: Any):
        super().__init__(*pathsegments)
        if not self.exists() or not self.is_file():
            raise InvalidPath


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


@dataclass
class SequenceData(DataType):
    def __init__(self, path: ValidFilePath, protocol: SequenceReader):
        """
        raises `InvalidPath` if `path` does not exists
        """
        self.path = path
        self.protocol = protocol

    def get_dataframe(self) -> pd.DataFrame:
        return self.protocol.read(self.path, columns=["seqid", "sequence"])


class SequenceReader(ABC):
    @abstractclassmethod
    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        pass


class ColumnsNotFound(Exception):
    def __init__(self, columns: Set[str]):
        if not isinstance(columns, set):
            raise TypeError("ColumnsNotFound argument should be a set")
        self.args = (columns,)


class TabfileReader(SequenceReader):
    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
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
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
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


class FastaReader(SequenceReader):
    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        if not set(columns) in {"seqid", "sequence"}:
            raise ColumnsNotFound(set(columns) - {"seqid", "sequence"})
        with open(path, errors="replace") as file:
            _, records = Fastafile.read(file)
            return pd.DataFrame(
                ([record["seqid"], record["sequence"]] for record in records()),
                columns=["seqid", "sequence"],
            )


class GenbankReader(SequenceReader):
    @staticmethod
    def select_columns(columns: List[str], record: Record) -> List[str]:
        return list(
            record[column] if column != "species" else record["organism"]
            for column in columns
        )

    @staticmethod
    def read(path: ValidFilePath, *, columns: List[str]) -> pd.DataFrame:
        with open(path, errors="replace") as file:
            names, records = GenbankFile.read(file)
        if not set(columns) in set(names):
            raise ColumnsNotFound(set(columns) - set(names))
        return pd.DataFrame(
            (GenbankReader.select_columns(columns, record) for record in records()),
            columns=columns,
        )
