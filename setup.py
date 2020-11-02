#!/usr/bin/env python

#import mtpy

# Check for setuptools package:

try:
    from setuptools import setup
except ImportError:
    setuptools = False
    from distutils.core import setup
else:
    setuptools = True

import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

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

setup_kwargs['entry_points'] = {'console_scripts': 
                    ['ws2vtk = mtpy.utils.ws2vtk:main',
                     'modem_pyqt = mtpy.gui.modem_pyqt:main',
                     'modem_plot_response = mtpy.gui.modem_plot_response:main',
                     'modem_plot_pt_maps = mtpy.gui.modem_plot_pt_maps:main',
                     'modem_mesh_builder = mtpy.gui.modem_mesh_builder:main',
                     'modem2vtk = mtpy.utils.modem2vtk:main',
                     'occam1d_gui = mtpy.gui.occam1d_gui:main',
                     'edi_editor = mtpy.gui.edi_editor:main']}

# But many people will not have setuptools installed, so we need to handle
# the default Python installation, which only has Distutils:

if setuptools is False:
    # Different script specification style for ordinary Distutils:

    setup_kwargs['scripts'] = [
        s.split(' = ')[1].replace('.', '/').split(':')[0] + '.py' for s in 
        setup_kwargs['entry_points']['console_scripts']]
    del setup_kwargs['entry_points']

    # "You must explicitly list all packages in packages: the Distutils will not
    # recursively scan your source tree looking for any directory with an
    # __init__.py file"

setup_kwargs['packages'] = [ 
                            'mtpy',
                            'mtpy.core',
                            'mtpy.imaging',
                            'mtpy.utils',
                            'mtpy.modeling',
                            'mtpy.modeling.modem',
                            'mtpy.contrib',
                            'mtpy.contrib.netcdf',
                            'mtpy.processing',
                            'mtpy.analysis',
                            #'tests', 
                            #'mtpy.test',
                            'mtpy.uofa',
                            'mtpy.usgs',
                            'mtpy.gui']
     
setup_kwargs['install_requires'] = ['numpy>=1.8.1',
                                     'scipy>=0.14.0',
                                     'matplotlib',
                                     'pyyaml',
                                     'pyproj',
                                     'configparser']

setup_kwargs['data_files'] = [('data', ['mtpy/utils/epsg.npy'])]

setup(
	name="mtpy",
	version=get_version("mtpy/__init__.py") ,
	author="Alison Kirkby,Fei Zhang,Jared Peacock,Rakib Hassan, Jinming Duan",
    author_email="Fei.Zhang@ga.gov.au",
	description="Python toolkit for standard magnetotelluric data processing.",
	long_description=LONG_DESC,
    url="https://github.com/MTgeophysics/mtpy",
	#data_files=[('', ['mtpy/utils/epsg.npy',]),], #this will install datafiles in wearied palce such as ~/.local/
	include_package_data=True,
	license="GNU GENERAL PUBLIC LICENSE v3",
	classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        ],
	**setup_kwargs)
