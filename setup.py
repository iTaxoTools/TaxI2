"""A setuptools based setup module."""

from setuptools import find_namespace_packages, setup

setup(
    package_dir={"": "src"},
    include_package_data=True,
    packages=find_namespace_packages(
        include=("itaxotools*",),
        where="src",
    ),
)
