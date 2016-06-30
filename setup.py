#!/usr/bin/env python

# Check for setuptools package:

try:
    from setuptools import setup
except ImportError:
    setuptools = False
    from distutils.core import setup
else:
    setuptools = True

setup_kwargs = {}

# The advantage of setuptools is that EXE wrappers are created on Windows,
# which allows Tab-completion for the script names from the system Scripts
# folder.

# Add names of scripts here. You can also specify the function to call
# by adding :func_name after the module name, and the name of the script
# can be customized before the equals sign.

setup_kwargs['entry_points'] = {'console_scripts': 
                    ['ws2vtk = mtpy.utils.ws2vtk:main',
                     'modem_plot_response = mtpy.gui.modem_plot_response:main',
                     'modem_plot_pt_maps = mtpy.gui.modem_plot_pt_maps:main',
					 'modem_mesh_builder = mtpy.gui.modem_mesh_builder:main',
                     'modem2vtk = mtpy.utils.modem2vtk:main',
                     'occam1d_gui = mtpy.gui.occam1d_gui:main']}

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

setup_kwargs['packages'] = ['mtpy',
                            'mtpy.core',
                            'mtpy.imaging',
                            'mtpy.utils',
                            'mtpy.modeling',
                            'mtpy.processing',
                            'mtpy.analysis',
                            'mtpy.test',
                            'mtpy.uofa',
                            'mtpy.usgs',
					'mtpy.gui']

	

setup(name = "mtpy", 
		version = '0.0.1',
		description = ("Collection of python tools for standard MT data processing."),
		license = "GNU GENERAL PUBLIC LICENSE v3",
		**setup_kwargs)
