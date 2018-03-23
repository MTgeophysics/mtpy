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
edipath = r'examples/data/edi_files'  # / is Unix and Win-Dos compatible

edi_path = edipath  # sys.argv[1]

if os.path.isfile(edi_path):
    pd1d.plot_edi_file(edi_path, savefile='C:/temp/pen_depth.jpg')
            # rholist can be any of ['zxy','zyx','det'], default all of them
elif os.path.isdir(edi_path):  # choose a suitable function below at run
    # plot_edi_dir(edi_path )
    pd1d.plot_edi_dir(edi_path, rholist=['det'])
