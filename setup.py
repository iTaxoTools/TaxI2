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
    name="taxi3",
    version="0.2.dev0",
    description="taxi3 description",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iTaxoTools/taxi3/",
    author="Vladimir Kharchev",
    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    package_dir={"": "src"},
    packages=find_namespace_packages(
        # exclude=('itaxotools.common*',),
        include=("itaxotools*",),
        where="src",
    ),
    rust_extensions=[
        RustExtension(
            "itaxotools.taxi3.library.calculate_distances",
            binding=Binding.PyO3,
            path=str(
                here
                / "src"
                / "itaxotools"
                / "taxi3"
                / "library"
                / "calculate_distances"
                / "Cargo.toml"
            ),
        )
    ],
    python_requires=">=3.9, <4",
    install_requires=[
        "alfpy",
        "appdirs",
        "BioPython",
        "DNAconvert",
        "numpy",
        "networkx",
        "openpyxl",
        "pandas",
        "seaborn",
    ],
    extras_require={
        "dev": ["pyinstaller"],
    },
    # Include all data from MANIFEST.in
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "taxi3=itaxotools.taxi3:main",
        ],
        "pyinstaller40": [
            "hook-dirs = itaxotools.__pyinstaller:get_hook_dirs",
            "tests = itaxotools.__pyinstaller:get_pyinstaller_tests",
        ],
    },
)
