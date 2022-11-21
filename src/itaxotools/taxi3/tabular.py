from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple


class Tabular:
    def __init__(
        self,
        path: Path,
        mode: 'r' | 'w' = 'r',
        columns: iter[int | str] = None,
        has_headers: bool = False,
    ):
        if columns is not None:
            has_headers = True
        self.path = path
        self.has_headers = has_headers
        self.columns = columns
        self.it = self.iter_read()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.it)

    @classmethod
    def open(cls, *args, **kwargs) -> ReadHandler:
        return cls(*args, **kwargs)

    def close(self) -> None:
        self.it.close()

    def iter_rows(self) -> iter[tuple[str, ...]]:
        with open(self.path) as file:
            for line in file:
                yield tuple(line.strip().split('\t'))

    def iter_read(self) -> iter[tuple[str, ...]]:
        rows = self.iter_rows()
        columns = self.columns
        if self.has_headers:
            header_row = next(rows)
        if columns is None:
            yield from rows
        else:
            columns = tuple(columns)
            if not len(columns):
                raise ValueError('Columns argument must contain at least one item')
            if isinstance(columns[0], str):
                try:
                    columns = tuple(header_row.index(x) for x in columns)
                except Exception as e:
                    missing = set(columns) - set(header_row)
                    raise ValueError(f'Column header(s) not found in file: {missing}') from e
            for row in rows:
                yield tuple(row[x] for x in columns)

    def read(self) -> tuple[str, ...] | None:
        try:
            return next(self.it)
        except StopIteration:
            return None
