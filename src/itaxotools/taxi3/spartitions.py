from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

from .sequences import Sequence, Sequences
from .types import Container, Type


class Spartition(dict[str, str]):
    @classmethod
    def fromFile(cls, file: SpartitionFile, *args, **kwargs) -> Spartition:
        return file.read(*args, **kwargs)


class SpartitionFile(Type):
    """Handlers for spartition files"""

    def __init__(self, path: Path):
        self.path = path

    def read(self, *args, **kwargs) -> Spartition:
        raise NotImplementedError()


class Tabfile(SpartitionFile):
    pass


class Spart(SpartitionFile):
    pass
