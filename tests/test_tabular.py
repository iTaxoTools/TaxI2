from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.tabular import Tabular

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem


Item = tuple[str]

class Items(NamedTuple):
    headers: Item
    items: list[Item]


class ReadTest(NamedTuple):
    fixture: Callable[[], Items]
    asserter: Callable[[Item, Item], bool]
    input: str
    protocol: Tabular
    kwargs: dict = {}

    def validate_context_iter(self, path: Path, fixture: Items):
        with self.protocol.open(path, **self.kwargs) as file:
            for item, fixed in zip(file, fixture.items):
                self.asserter(item, fixed)

    def validate_context_read(self, path: Path, fixture: Items):
        with self.protocol.open(path, **self.kwargs) as file:
            for fixed in fixture.items:
                item = file.read()
                self.asserter(item, fixed)

    def validate_context_exhaust(self, path: Path, fixture: Items):
        with self.protocol.open(path, **self.kwargs) as file:
            fixed_iter = (item for item in fixture.items)
            while True:
                item = file.read()
                if item is None:
                    break
                fixed = next(fixed_iter)
                self.asserter(item, fixed)

    def validate_open_iter(self, path: Path, fixture: Items):
        file = self.protocol.open(path, **self.kwargs)
        for item, fixed in zip(file, fixture.items):
            self.asserter(item, fixed)
        file.close()

    @staticmethod
    def assert_all(item: Item, fixed: Item) -> bool:
        assert item == fixed

    @staticmethod
    def assert_0_2(item: Item, fixed: Item) -> bool:
        assert item == (fixed[0], fixed[2])


class WriteTest(NamedTuple):
    fixture: Callable[[], Items]
    output: str
    protocol: Tabular
    kwargs: dict = {}

    def validate_context(self, fixed_path: Path, output_path: Path, fixture: Items):
        with self.protocol.open(output_path, 'w', **self.kwargs) as file:
            for item in fixture.items:
                file.write(item)
        assert_eq_files(output_path, fixed_path)

    def validate_open(self, fixed_path: Path, output_path: Path, fixture: Items):
        file = self.protocol.open(output_path, 'w', **self.kwargs)
        for item in fixture.items:
            file.write(item)
        file.close()
        assert_eq_files(output_path, fixed_path)


def items_simple() -> tuple:
    return Items(
        headers = ('header_1', 'header_2', 'header_3'),
        items = [
            ('item_1_1', 'item_1_2', 'item_1_3'),
            ('item_2_1', 'item_2_2', 'item_2_3'),
            ('item_3_1', 'item_3_2', 'item_3_3'),
        ]
    )


@pytest.mark.parametrize(
    "test", [
    ReadTest(items_simple, ReadTest.assert_all, 'simple.tsv', Tabular),
    ReadTest(items_simple, ReadTest.assert_all, 'headers.tsv', Tabular, dict(has_headers=True)),
    ReadTest(items_simple, ReadTest.assert_0_2, 'headers.tsv', Tabular, dict(columns=[0, 2])),
    ReadTest(items_simple, ReadTest.assert_0_2, 'headers.tsv', Tabular, dict(columns=['header_1', 'header_3'])),
])
@pytest.mark.parametrize(
    "validator", [
    ReadTest.validate_context_iter,
    ReadTest.validate_context_read,
    ReadTest.validate_context_exhaust,
    ReadTest.validate_open_iter,
])
def test_read_tabular(test: ReadTest, validator: Callable) -> None:
    path = TEST_DATA_DIR / test.input
    fixture = test.fixture()
    validator(test, path, fixture)


def test_read_tabular_missing_header() -> None:
    test = ReadTest(items_simple, ReadTest.assert_all, 'headers.tsv', Tabular, dict(columns=['header_X']))
    path = TEST_DATA_DIR / test.input
    fixture = test.fixture()
    with pytest.raises(ValueError):
        test.validate_context_iter(path, fixture)


def test_read_tabular_zero_columns() -> None:
    test = ReadTest(items_simple, ReadTest.assert_all, 'headers.tsv', Tabular, dict(columns=[]))
    path = TEST_DATA_DIR / test.input
    fixture = test.fixture()
    with pytest.raises(ValueError):
        test.validate_context_iter(path, fixture)


def test_read_tabular_early_close() -> None:
    path = TEST_DATA_DIR / 'simple.tsv'
    file = Tabular.open(path)
    file.read()
    file.close()


@pytest.mark.parametrize(
    "test", [
        WriteTest(items_simple, 'simple.tsv', Tabular),
        WriteTest(items_simple, 'headers.tsv', Tabular, dict(columns=['header_1', 'header_2', 'header_3'])),
])
@pytest.mark.parametrize(
    "validator", [
    WriteTest.validate_context,
    WriteTest.validate_open,
])
def test_write_tabular(test: WriteTest, validator: Callable, tmp_path: Path) -> None:
    fixed_path = TEST_DATA_DIR / test.output
    output_path = tmp_path / test.output
    fixture = test.fixture()
    validator(test, fixed_path, output_path, fixture)
