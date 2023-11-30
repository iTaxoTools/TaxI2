import re
from pathlib import Path


def assert_eq_files(
    file1: Path,
    file2: Path,
    case_sensitive=True,
    ignore=r"\s",
) -> None:
    regex = re.compile(r"[" + ignore + "]")
    text1 = file1.read_text()
    text2 = file2.read_text()
    text1 = regex.sub("", text1)
    text2 = regex.sub("", text2)
    if not case_sensitive:
        text1 = text1.lower()
        text2 = text2.lower()
    assert text1 == text2
