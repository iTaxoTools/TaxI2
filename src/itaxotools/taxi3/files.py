from __future__ import annotations

from typing import Callable
from enum import Enum, auto
from pathlib import Path
from re import fullmatch


class Identifier:
    def isFasta(path: Path) -> bool:
        with path.open() as file:
            return bool(file.read(1) == '>')

    def isTabFile(path: Path) -> bool:
        with path.open() as file:
            line = file.readline()
            return bool(fullmatch(r'([^\t]+\t)+[^\t]+', line))


class FileFormat(Enum):
    Fasta = auto(), Identifier.isFasta
    Tabfile = auto(), Identifier.isTabFile
    Excel = auto(), lambda _: False
    Unknown = auto(), lambda _: True

    def __init__(self, _, identifier):
        self.identifier = identifier

    @classmethod
    def identify(cls, path: Path) -> FileFormat:
        for format in cls:
            if format.identifier(path):
                return format
        return cls.Unknown
