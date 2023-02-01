# PyInstaller entry points for setuptools
# https://pyinstaller.readthedocs.io/en/stable/hooks.html

# When a module is detected by PyInstaller, it will search
# for corresponding hooks and tests in this directory.

import itaxotools
from pathlib import Path


def get_namespace_dirs():
    return [str(Path(__file__).parent)]


def get_hook_dirs():
    return get_namespace_dirs()


def get_PyInstaller_tests():
    return get_namespace_dirs()
