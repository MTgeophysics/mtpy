#!/usr/bin/env python

# import mtpy

# Check for setuptools package:

from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()



LONG_DESC = """
MTPy is an open source Python package to assist with magnetotelluric (MT) data processing, analysis,
modelling, visualization and interpretation. 
"""

setup_kwargs = {}

# The advantage of setuptools is that EXE wrappers are created on Windows,
# which allows Tab-completion for the script names from the system Scripts
# folder.

# Add names of scripts here. You can also specify the function to call
# by adding :func_name after the module name, and the name of the script
# can be customized before the equals sign.

setup_kwargs["entry_points"] = {
    "console_scripts": [
        "ws2vtk = mtpy.utils.ws2vtk:main",
        "modem_pyqt = mtpy.gui.modem_pyqt:main",
        "modem_plot_response = mtpy.gui.modem_plot_response:main",
        "modem_plot_pt_maps = mtpy.gui.modem_plot_pt_maps:main",
        "modem_mesh_builder = mtpy.gui.modem_mesh_builder:main",
        "modem2vtk = mtpy.utils.modem2vtk:main",
        "occam1d_gui = mtpy.gui.occam1d_gui:main",
        "edi_editor = mtpy.gui.edi_editor:main",
    ]
}

# But many people will not have setuptools installed, so we need to handle
# the default Python installation, which only has Distutils:

if setuptools is False:
    # Different script specification style for ordinary Distutils:

    setup_kwargs["scripts"] = [
        s.split(" = ")[1].replace(".", "/").split(":")[0] + ".py"
        for s in setup_kwargs["entry_points"]["console_scripts"]
    ]
    del setup_kwargs["entry_points"]

    # "You must explicitly list all packages in packages: the Distutils will not
    # recursively scan your source tree looking for any directory with an
    # __init__.py file"

setup_kwargs["packages"] = [
    "mtpy",
    "mtpy.core",
    "mtpy.imaging",
    "mtpy.utils",
    "mtpy.modeling",
    "mtpy.modeling.modem",
    "mtpy.processing",
    "mtpy.analysis",
    "mtpy.gui",
]

setup_kwargs["install_requires"] = [
    "numpy>=1.8.1",
    "scipy>=0.14.0",
    "matplotlib",
    "pyyaml",
    "pyproj",
    "configparser",
    "mt_metadata",
    "mth5",
]

setup_kwargs["data_files"] = [("data", ["mtpy/utils/epsg.npy"])]

setup(
    name="mtpy",
    version=get_version("2.0.0"),
    author="Jared Peacock,Alison Kirkby,Fei Zhang,,Rakib Hassan, Jinming Duan",
    author_email="jpeacock@usgs.gov",
    description="Python toolkit for standard magnetotelluric data processing.",
    long_description=LONG_DESC,
    url="https://github.com/MTgeophysics/mtpy/tree/v2",
    include_package_data=True,
    license="GNU GENERAL PUBLIC LICENSE v3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    **setup_kwargs
)
