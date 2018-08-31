# -*- coding: utf-8 -*-
"""
Plot Penetration Depth from EDI files
Example how to use the module mtpy.imaging import penetration_depth1d as pd1d

Created on Wed Nov 08 12:04:38 2017

@author: u64125
"""
import os, sys
from mtpy.imaging import penetration_depth1d as pd1d

#edipath = r'C:\Git\mtpy\examples\data\edi_files'  # avoid using \
edipath = r'C:/mtpywin/mtpy/examples/data/edi_files_2/Synth00.edi'  # / is Unix and Win-Dos compatible
# or get this variable from the cmdline:  edipath = sys.argv[1]

fig_dpi = 400 # change to your preferred file resolution

if os.path.isfile(edipath):
    pd1d.plot_edi_file(edipath, savefile='C:/tmp/pen_depth.jpg',
                       fig_dpi=fig_dpi)
    # rholist can be any of ['zxy','zyx','det'], default all of them
elif os.path.isdir(edipath):
    # plot_edi_dir(edi_path )
    pd1d.plot_edi_dir(edipath, rholist=['det'])
