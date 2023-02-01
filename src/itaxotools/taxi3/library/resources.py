#!/usr/bin/env python3

"""Get resource paths for Python3.8+"""


from typing import Any
import pathlib
import sys
import os


def _package_path_pyinstaller(package):
    if isinstance(package, str):
        parts = package.split('.')
    elif isinstance(package, type(sys)):
        parts = package.__name__.split('.')
    else:
        return None
    path = pathlib.Path(sys._MEIPASS)
    for part in parts:
        path /= part
    return path


def _package_path_file(package):
    if isinstance(package, str):
        file = sys.modules[package].__file__
    elif isinstance(package, type(sys)):
        file = package.__file__
    else:
        return None
    path = pathlib.Path(os.path.dirname(file))
    return path


def _package_path_importlib(package):
    return importlib.resources.files(package)


try:
    import importlib.resources
    from importlib.resources import files
    package_path = _package_path_importlib
except (ModuleNotFoundError, ImportError):
    if hasattr(sys, '_MEIPASS'):
        package_path = _package_path_pyinstaller
    else:
        package_path = _package_path_file

_resource_path = package_path("itaxotools.taxi3")


def get_resource(path: Any) -> str:
    return str(_resource_path / "resources" / path)
