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

setup_kwargs['entry_points'] = {'console_scripts': [
        'CombineEDIs = MTpy.utils.CombineEDIs:main',
        'runParalanaMT = MTpy.utils.runParalanaMT:main',
        'wsmt_pv = MTpy.utils.wsmt_pv:main',
        'occam2d_gui = MTpy.utils.gui.occam2d.v1.run1:main',
        'RunBIRRPSingleStation = MTpy.core.RunBIRRPSingleStation:main']}

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

    setup_kwargs['packages'] = ['MTpy',
                                'MTpy.core',
                                'MTpy.imaging',
                                'MTpy.utils']

setup(name="MTpy", **setup_kwargs)
