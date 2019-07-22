# -*- coding: utf-8 -*-
"""
Plot Penetration Depth from EDI files
Example how to use the module mtpy.imaging import penetration_depth1d as pd1d

Created on Wed Nov 08 12:04:38 2017

@author: u64125
"""
import os
import tempfile
from mtpy.imaging import penetration_depth1d as pd1d

#edipath = r'C:\Git\mtpy\examples\data\edi_files'  # avoid using \
# edipath = r'C:/mtpywin/mtpy/examples/data/edi_files_2/Synth00.edi'  # / is Unix and Win-Dos compatible
# or get this variable from the cmdline:  edipath = sys.argv[1]

try:
    # PACK_ROOT = os.environ['PACK_ROOT']
    # mtpy_path = os.path.join(PACK_ROOT, 'mtpy')
    mtpy_path = os.environ['MTPY_ROOT']
except:
    print("Warn: The environment variable MTPY_ROOT is not defined. We will guess")
    mtpy_path = os.path.abspath('../..')

if not os.path.isdir(mtpy_path):
    raise Exception("the guessed mtpy dir %s is not a folder!"% mtpy_path)

# change the variable below according to your edi files folder !!!
# edidir = r'C:/mtpywin/mtpy/data/edifiles'  # / is Unix and Win-Dos compatible
# or get this variable from the cmdline:  edidir = sys.argv[1]

edipath = os.path.join(mtpy_path,'examples/data/edi_files_2/Synth00.edi')

# savepath = r'C:\tmp'
temp_dir = tempfile.gettempdir()
print('Using temporary directory ' + temp_dir)
savepath = temp_dir

fig_dpi = 400 # change to your preferred file resolution

if os.path.isfile(edipath):
    pd1d.plot_edi_file(edipath, savefile=os.path.join(savepath, 'pen_depth.jpg'),
                       fig_dpi=fig_dpi)
    # rholist can be any of ['zxy','zyx','det'], default all of them
elif os.path.isdir(edipath):
    # plot_edi_dir(edi_path )
    pd1d.plot_edi_dir(edipath, rholist=['det'])
