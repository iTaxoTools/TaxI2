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

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> Path:
        return self.fixture()

    def validate_context_iter(self):
        with self.protocol.open(self.input_path, **self.kwargs) as file:
            for item, fixed in zip(file, self.fixed.items):
                self.asserter(item, fixed)

    def validate_context_read(self):
        with self.protocol.open(self.input_path, **self.kwargs) as file:
            for fixed in self.fixed.items:
                item = file.read()
                self.asserter(item, fixed)

    def validate_context_exhaust(self):
        with self.protocol.open(self.input_path, **self.kwargs) as file:
            fixed_iter = (item for item in self.fixed.items)
            while True:
                item = file.read()
                if item is None:
                    break
                fixed = next(fixed_iter)
                self.asserter(item, fixed)

    def validate_open_iter(self):
        file = self.protocol.open(self.input_path, **self.kwargs)
        for item, fixed in zip(file, self.fixed.items):
            self.asserter(item, fixed)
        file.close()

    @staticmethod
    def assert_all(item: Item, fixed: Item) -> bool:
        assert item == fixed

    @staticmethod
    def assert_0_2(item: Item, fixed: Item) -> bool:
        assert item == (fixed[0], fixed[2])


class HeaderTest(NamedTuple):
    fixture: Callable[[], Items]
    input: str
    protocol: Tabular

    @property
    def input_path(self) -> Path:
        return TEST_DATA_DIR / self.input

    @property
    def fixed(self) -> Path:
        return self.fixture()

    def validate(self):
        headers = self.protocol.headers(self.input_path)
        assert headers == self.fixed.headers


class WriteTest(NamedTuple):
    fixture: Callable[[], Items]
    output: str
    protocol: Tabular
    kwargs: dict = {}

    @property
    def fixed_path(self) -> Path:
        return TEST_DATA_DIR / self.output

    @property
    def fixed(self) -> Path:
        return self.fixture()

    def get_output_path(self, tmp_path) -> Path:
        return tmp_path / self.output

    def validate_context(self, output_path: Path):
        with self.protocol.open(output_path, 'w', **self.kwargs) as file:
            for item in self.fixed.items:
                file.write(item)
        assert_eq_files(output_path, self.fixed_path)

    def validate_open(self, output_path: Path):
        file = self.protocol.open(output_path, 'w', **self.kwargs)
        for item in self.fixed.items:
            file.write(item)
        file.close()
        assert_eq_files(output_path, self.fixed_path)


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
    validator(test)


def test_read_tabular_missing_header() -> None:
    test = ReadTest(items_simple, ReadTest.assert_all, 'headers.tsv', Tabular, dict(columns=['header_X']))
    with pytest.raises(ValueError):
        test.validate_context_iter()


def test_read_tabular_zero_columns() -> None:
    test = ReadTest(items_simple, ReadTest.assert_all, 'headers.tsv', Tabular, dict(columns=[]))
    with pytest.raises(ValueError):
        test.validate_context_iter()


def test_read_tabular_early_close() -> None:
    path = TEST_DATA_DIR / 'simple.tsv'
    file = Tabular.open(path)
    file.read()
    file.close()

@pytest.mark.parametrize(
    "test", [
    HeaderTest(items_simple, 'headers.tsv', Tabular),
])
def test_read_tabular_headers(test: HeaderTest) -> None:
    test.validate()


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
    output_path = test.get_output_path(tmp_path)
    validator(test, output_path)
