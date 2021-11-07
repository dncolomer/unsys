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

VERSION = "0.2.0"

here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = "\n" + f.read()

setup(
    name="unsys",
    version=VERSION,
    author="Daniel Colomer",
    author_email="uncertainsystems@gmail.com",
    description=("Uncertain Sytsems Toolkit"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords=["graph", "quantum"],
    url="https://github.com/dncolomer/_intent/tarball/" + VERSION,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    classifiers=[],
    install_requires=[
        "numpy==1.20.1",
        "sympy==1.7.1",
        "igraph==0.9.8",
        "celluloid==0.2.0",
        "networkx==1.2.0"
    ],
    include_package_data=True,
)
