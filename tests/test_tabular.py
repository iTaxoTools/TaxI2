from __future__ import annotations

from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from utility import assert_eq_files

from itaxotools.taxi3.tabular import TabfileHandler, ExcelHandler, Row

TEST_DATA_DIR = Path(__file__).parent / Path(__file__).stem

Item = Row


class Items(NamedTuple):
    headers: Item | None
    items: list[Item]


class ReadTest(NamedTuple):
    fixture: Callable[[], Items]
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
        with self.protocol(self.input_path, **self.kwargs) as file:
            assert not file.closed
            for item, fixed in zip(file, self.fixed.items):
                assert item == fixed
        assert file.closed

    def validate_context_read(self):
        with self.protocol(self.input_path, **self.kwargs) as file:
            assert not file.closed
            for fixed in self.fixed.items:
                item = file.read()
                assert item == fixed
        assert file.closed

    def validate_context_exhaust(self):
        with self.protocol(self.input_path, **self.kwargs) as file:
            assert not file.closed
            fixed_iter = (item for item in self.fixed.items)
            while True:
                item = file.read()
                if item is None:
                    break
                fixed = next(fixed_iter)
                assert item == fixed
        assert file.closed

    def validate_open_iter(self):
        file = self.protocol(self.input_path, **self.kwargs)
        assert not file.closed
        for item, fixed in zip(file, self.fixed.items):
            assert item == fixed
        file.close()
        assert file.closed

    def validate_headers(self):
        with self.protocol(self.input_path, **self.kwargs) as file:
            assert not file.closed
            assert file.headers == self.fixed.headers
        assert file.closed


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
        headers = self.protocol.get_headers(self.input_path)
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
        with self.protocol(output_path, 'w', **self.kwargs) as file:
            assert not file.closed
            for item in self.fixed.items:
                file.write(item)
        assert file.closed
        assert_eq_files(output_path, self.fixed_path)

    def validate_open(self, output_path: Path):
        file = self.protocol(output_path, 'w', **self.kwargs)
        assert not file.closed
        for item in self.fixed.items:
            file.write(item)
        assert not file.closed
        file.close()
        assert file.closed
        assert_eq_files(output_path, self.fixed_path)


def items_simple_headers_all() -> tuple:
    return Items(
        headers = ('header_1', 'header_2', 'header_3'),
        items = [
            ('item_1_1', 'item_1_2', 'item_1_3'),
            ('item_2_1', 'item_2_2', 'item_2_3'),
            ('item_3_1', 'item_3_2', 'item_3_3'),
        ]
    )


def items_simple_plain_all() -> tuple:
    return Items(
        headers = None,
        items = [
            ('item_1_1', 'item_1_2', 'item_1_3'),
            ('item_2_1', 'item_2_2', 'item_2_3'),
            ('item_3_1', 'item_3_2', 'item_3_3'),
        ]
    )


def items_simple_headers_0_2() -> tuple:
    return Items(
        headers = ('header_1', 'header_3'),
        items = [
            ('item_1_1', 'item_1_3'),
            ('item_2_1', 'item_2_3'),
            ('item_3_1', 'item_3_3'),
        ]
    )


def items_simple_plain_0_2() -> tuple:
    return Items(
        headers = None,
        items = [
            ('item_1_1', 'item_1_3'),
            ('item_2_1', 'item_2_3'),
            ('item_3_1', 'item_3_3'),
        ]
    )


def items_simple_headers_0_2_1() -> tuple:
    return Items(
        headers = ('header_1', 'header_3', 'header_2'),
        items = [
            ('item_1_1', 'item_1_3', 'item_1_2'),
            ('item_2_1', 'item_2_3', 'item_2_2'),
            ('item_3_1', 'item_3_3', 'item_3_2'),
        ]
    )


def items_simple_plain_0_2_1() -> tuple:
    return Items(
        headers = None,
        items = [
            ('item_1_1', 'item_1_3', 'item_1_2'),
            ('item_2_1', 'item_2_3', 'item_2_2'),
            ('item_3_1', 'item_3_3', 'item_3_2'),
        ]
    )


@pytest.mark.parametrize(
    "test", [
    ReadTest(items_simple_plain_all, 'simple.tsv', TabfileHandler),
    ReadTest(items_simple_plain_0_2, 'simple.tsv', TabfileHandler, dict(columns=[0, 2])),
    ReadTest(items_simple_plain_0_2_1, 'simple.tsv', TabfileHandler, dict(columns=[0, 2], get_all_columns=True)),
    ReadTest(items_simple_headers_all, 'headers.tsv', TabfileHandler, dict(has_headers=True)),
    ReadTest(items_simple_headers_0_2, 'headers.tsv', TabfileHandler, dict(columns=[0, 2], has_headers=True)),
    ReadTest(items_simple_headers_0_2_1, 'headers.tsv', TabfileHandler, dict(columns=[0, 2], has_headers=True, get_all_columns=True)),
    ReadTest(items_simple_headers_0_2, 'headers.tsv', TabfileHandler, dict(columns=['header_1', 'header_3'])),
    ReadTest(items_simple_headers_0_2_1, 'headers.tsv', TabfileHandler, dict(columns=['header_1', 'header_3'], get_all_columns=True)),

    ReadTest(items_simple_plain_all, 'simple.xlsx', ExcelHandler),
    ReadTest(items_simple_plain_0_2, 'simple.xlsx', ExcelHandler, dict(columns=[0, 2])),
    ReadTest(items_simple_plain_0_2_1, 'simple.xlsx', ExcelHandler, dict(columns=[0, 2], get_all_columns=True)),
    ReadTest(items_simple_headers_all, 'headers.xlsx', ExcelHandler, dict(has_headers=True)),
    ReadTest(items_simple_headers_0_2, 'headers.xlsx', ExcelHandler, dict(columns=[0, 2], has_headers=True)),
    ReadTest(items_simple_headers_0_2_1, 'headers.xlsx', ExcelHandler, dict(columns=[0, 2], has_headers=True, get_all_columns=True)),
    ReadTest(items_simple_headers_0_2, 'headers.xlsx', ExcelHandler, dict(columns=['header_1', 'header_3'])),
    ReadTest(items_simple_headers_0_2_1, 'headers.xlsx', ExcelHandler, dict(columns=['header_1', 'header_3'], get_all_columns=True)),
])
@pytest.mark.parametrize(
    "validator", [
    ReadTest.validate_context_iter,
    ReadTest.validate_context_read,
    ReadTest.validate_context_exhaust,
    ReadTest.validate_open_iter,
    ReadTest.validate_headers,
])
def test_read_tabular(test: ReadTest, validator: Callable) -> None:
    validator(test)


def test_read_tabular_missing_header() -> None:
    test = ReadTest(items_simple_headers_all, 'headers.tsv', TabfileHandler, dict(columns=['header_X']))
    with pytest.raises(ValueError):
        test.validate_context_iter()


def test_read_tabular_zero_columns() -> None:
    test = ReadTest(items_simple_headers_all, 'headers.tsv', TabfileHandler, dict(columns=[]))
    with pytest.raises(ValueError):
        test.validate_context_iter()


def test_read_tabular_early_close() -> None:
    path = TEST_DATA_DIR / 'simple.tsv'
    file = TabfileHandler(path)
    file.read()
    assert not file.closed
    file.close()
    assert file.closed


@pytest.mark.parametrize(
    "test", [
    HeaderTest(items_simple_headers_all, 'headers.tsv', TabfileHandler),
    HeaderTest(items_simple_headers_all, 'headers.xlsx', ExcelHandler),
])
def test_read_tabular_headers(test: HeaderTest) -> None:
    test.validate()


@pytest.mark.parametrize(
    "test", [
    WriteTest(items_simple_plain_all, 'simple.tsv', TabfileHandler),
    WriteTest(items_simple_headers_all, 'headers.tsv', TabfileHandler, dict(columns=['header_1', 'header_2', 'header_3'])),
])
@pytest.mark.parametrize(
    "validator", [
    WriteTest.validate_context,
    WriteTest.validate_open,
])
def test_write_tabular(test: WriteTest, validator: Callable, tmp_path: Path) -> None:
    output_path = test.get_output_path(tmp_path)
    validator(test, output_path)
