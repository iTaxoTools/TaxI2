"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
from setuptools_rust import RustExtension, Binding
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="taxi2",
    version="2.1.0",
    description="TaxI2 - Calculation and analysis of pairwise sequence distances",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iTaxoTools/taxi2/",
    author="Stefanos Patmanidis, Vladimir Kharchev",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    package_dir={"": "src"},
    packages=find_namespace_packages(
        include=("itaxotools*",),
        where="src",
    ),
    rust_extensions=[
        RustExtension(
            "itaxotools.taxi2.library.calculate_distances",
            binding=Binding.PyO3,
            path=str(
                here
                / "src"
                / "itaxotools"
                / "taxi2"
                / "library"
                / "calculate_distances"
                / "Cargo.toml"
            ),
        )
    ],
    python_requires=">=3.9, <4",
    install_requires=[
        "itaxotools-common>=0.2.4",
        "DNAconvert==0.2.0",
        "spart_parser==0.1.1",
        "BioPython>=1.80",
        "alfpy",
        "appdirs",
        "numpy",
        "networkx",
        "openpyxl",
        "pandas",
        "seaborn",
        "scipy",
    ],
    extras_require={
        "dev": [
            "setuptools_rust",
            "pyinstaller",
            "pytest",
            "coverage",
            "flake8",
            "isort",
        ],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "taxi2=itaxotools.taxi2.taxi2:main",
        ],
        "pyinstaller40": [
            "hook-dirs = itaxotools.__pyinstaller:get_hook_dirs",
            "tests = itaxotools.__pyinstaller:get_pyinstaller_tests",
        ],
    },
)
