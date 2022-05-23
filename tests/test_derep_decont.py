#!/usr/bin/env python3

from pathlib import Path

from itaxotools.taxi3.library.task import Decontaminate, Alignment
from itaxotools.taxi3.library.datatypes import (
    SequenceData,
    CompleteData,
    ValidFilePath,
    TabfileReader,
)

TEST_DATA_DIR = Path(__file__).parent / "tests/derep_decont_test_data"
TMP_TEST_DIR = Path(__file__).parent / "temp_test_files"


def test_decont2() -> None:
    data = CompleteData.from_path(
        ValidFilePath(TEST_DATA_DIR / "Scaphio_small_input.txt"), TabfileReader()
    )
    reference = SequenceData.from_path(
        ValidFilePath(TEST_DATA_DIR / "Scaphio_small_input_reference.txt"),
        TabfileReader(),
    )
    reference2 = SequenceData.from_path(
        ValidFilePath(TEST_DATA_DIR / "Scaphio_small_input_reference2.txt"),
        TabfileReader(),
    )
    decont2_task = Decontaminate(print)
    decont2_task.alignment = Alignment.AlignmentFree
    decont2_task.data = data
    decont2_task.reference = reference
    decont2_task.reference2 = reference2
    decont2_task.start()
    contaminates_output = TMP_TEST_DIR / "contaminates.txt"
    tested_contaminates = TEST_DATA_DIR / "Scaphio_small_input_contaminates.txt"
    decontaminated_output = TMP_TEST_DIR / "decontaminated.txt"
    tested_decontaminated = TEST_DATA_DIR / "Scaphio_small_input_decontaminated.txt"

    # clear output files
    with open(contaminates_output, mode="w"):
        pass
    with open(decontaminated_output, mode="w"):
        pass

    for chunk in decont2_task.result:
        chunk.contaminates.append_to_file(contaminates_output)
        chunk.decontaminated.append_to_file(decontaminated_output)

    assert (
        contaminates_output.read_text().splitlines()
        == tested_contaminates.read_text().splitlines()
    )
    assert (
        decontaminated_output.read_text().splitlines()
        == tested_decontaminated.read_text().splitlines()
    )

    contaminates_output.unlink()
    decontaminated_output.unlink()
