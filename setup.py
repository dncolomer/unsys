#!/usr/bin/env python3

import os
import io
from setuptools import find_packages, setup, Command

"""
git tag {VERSION}
git push --tags
python setup.py sdist
twine upload dist/*
"""

VERSION = "0.1.0"

here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = "\n" + f.read()

setup(
    name="ehsim",
    version=VERSION,
    author="D Colomer",
    author_email="dncolomer",
    description=("Entanglement Hypergraph Simulator"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords=["graph", "quantum"],
    url="https://github.com/dncolomer/ehsim_prototype/tarball/" + VERSION,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    classifiers=[],
    install_requires=[
        "numpy==1.20.1",
        "networkx==2.5",
        "hypernetx==0.3.7",
        "pandas==1.2.3",
    ],
    include_package_data=True,
)
