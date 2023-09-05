#!/usr/bin/env python

# import mtpy

# Check for setuptools package:

from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()


requirements = [
    "numpy",
    "scipy",
    "matplotlib<=3.5.3",
    "pyyaml",
    "pyproj",
    "configparser",
    "mt_metadata",
    "mth5",
    "pandas",
    "geopandas",
    "contextily",
    "pyevtk",
    "loguru",
]

setup_requirements = [
    "pytest-runner",
]

test_requirements = [
    "pytest>=3",
]


setup(
    author="Jared Peacock, Alison Kirkby, Fei Zhang, Rakib Hassan, Lars Krieger, Stephan Thiel",
    author_email="jpeacock@usgs.gov",
    python_requires=">=3.5",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU GENERAL PUBLIC License v3",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    description="Python toolkit for standard magnetotelluric data processing.",
    install_requires=requirements,
    license="GNU GENERAL PUBLIC LICENSE v3",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="magnetotellurics",
    name="mtpy",
    packages=find_packages(include=["mtpy", "mtpy.*"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/MTgeophysics/mtpy/tree/v2",
    version="2.0.0",
    zip_safe=False,
    package_data={"": ["mtpy/utils/epsg.npy"]},
)
