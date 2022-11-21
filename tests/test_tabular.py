from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest

from itaxotools.taxi3.tabular import Tabular

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


Item = tuple[str]

class Items(NamedTuple):
    headers: Item
    items: list[Item]


class ReadTest(NamedTuple):
    fixture: Callable[[], Items]
    validator: Callable[[Item, Item], bool]
    input: str
    protocol: Tabular
    kwargs: dict = {}

    def validate(self):
        path = TEST_DATA_DIR / self.input
        for validate in (
            self.validate_context_iter,
            self.validate_context_read,
            self.validate_context_exhaust,
        ):
            fixture = self.fixture()
            validate(path, fixture)

    def validate_context_iter(self, path: Path, fixture: Items):
        with self.protocol.open(path, **self.kwargs) as file:
            for item, fixed in zip(file, fixture.items):
                assert self.validator(item, fixed)

    def validate_context_read(self, path: Path, fixture: Items):
        with self.protocol.open(path, **self.kwargs) as file:
            for fixed in fixture.items:
                item = file.read()
                assert self.validator(item, fixed)

    def validate_context_exhaust(self, path: Path, fixture: Items):
        with self.protocol.open(path, **self.kwargs) as file:
            fixed_iter = (item for item in fixture.items)
            while True:
                item = file.read()
                if item is None:
                    break
                fixed = next(fixed_iter)
                assert self.validator(item, fixed)


def val_all(item: Item, fixed: Item) -> bool:
    return item == fixed


def val_0_2(item: Item, fixed: Item) -> bool:
    return item == (fixed[0], fixed[2])


def items_simple() -> tuple:
    return Items(
        headers = ('header_1', 'header_2', 'header_3'),
        items = [
            ('item_1_1', 'item_1_2', 'item_1_3'),
            ('item_2_1', 'item_2_2', 'item_2_3'),
            ('item_3_1', 'item_3_2', 'item_3_3'),
        ]
    )


read_tests = [
    ReadTest(items_simple, val_all, 'simple.tsv', Tabular),
    ReadTest(items_simple, val_all, 'headers.tsv', Tabular, dict(has_headers=True)),
    ReadTest(items_simple, val_0_2, 'headers.tsv', Tabular, dict(columns=[0, 2])),
    ReadTest(items_simple, val_0_2, 'headers.tsv', Tabular, dict(columns=['header_1', 'header_3'])),
]


@pytest.mark.parametrize("test", read_tests)
def test_read_tabular(test: ReadTest) -> None:
    test.validate()


def test_read_tabular_missing_header() -> None:
    test = ReadTest(items_simple, val_all, 'headers.tsv', Tabular, dict(columns=['header_X']))
    with pytest.raises(ValueError):
        test.validate()


def test_read_tabular_zero_columns() -> None:
    test = ReadTest(items_simple, val_all, 'headers.tsv', Tabular, dict(columns=[]))
    with pytest.raises(ValueError):
        test.validate()
