
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-legacy", action="store_true", default=False, help="run legacy tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "legacy: mark test as legacy")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-legacy"):
        skip_legacy = pytest.mark.skip(reason="need --run-legacy option to run")
        for item in items:
            if "legacy" in item.keywords:
                item.add_marker(skip_legacy)
