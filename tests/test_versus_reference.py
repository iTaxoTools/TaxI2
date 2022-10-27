#!/usr/bin/env python3

from pathlib import Path
import pytest

from itaxotools.taxi3.library.datatypes import (
    TabfileReader,
    CompleteData,
    VersusReferenceSummary,
    ValidFilePath,
)
from itaxotools.taxi3.library.task import (
    Alignment,
    VersusReference,
)

TEST_DATA_DIR = Path(__file__).parent / "versus_reference_data"
TMP_TEST_DIR = Path(__file__).parent / "temp_test_files"


@pytest.mark.legacy
def test_versus_reference() -> None:
    reference = CompleteData.from_path(
        ValidFilePath(TEST_DATA_DIR / "Scaphio_input_small.txt"), TabfileReader
    )
    data = CompleteData.from_path(
        ValidFilePath(TEST_DATA_DIR / "Boophis_input.txt"), TabfileReader
    )
    task = VersusReference(print)
    task.data = data
    task.reference = reference
    task.alignment = Alignment.Pairwise
    task.start()

    output = TMP_TEST_DIR / "output"
    with open(output, mode="w"):  # clearing the file
        pass

    for chunk in task.result:
        chunk.append_to_file(output)

    print(output.read_text())
    output.unlink()
