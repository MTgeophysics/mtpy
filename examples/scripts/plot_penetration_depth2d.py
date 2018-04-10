#! /usr/bin/env python
"""
Description:
    Example python script
    plot 3D penetration depth for a folder of EDI files

    input = path2edifolder

CreationDate:   23/03/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     23/03/2018   FZ

"""

import os, sys
from mtpy.imaging import penetration_depth2d as pen2d

# change the variable below according to your edi files folder !!!
edidir = r'C:/mtpywin/mtpy/data/edifiles'  # / is Unix and Win-Dos compatible
# or get this variable from the cmdline:  edidir = sys.argv[1]

if not os.path.isdir(edidir):
    print ("please provide the path to edi folder")
    sys.exit(1)

period_index_list = [0, 1, 10, 20, 30, 40, 50, 59]  # user to customise

# show three different kind of calculated pen-depth
pen2d.plot2Dprofile(edidir, period_index_list, 'det')

pen2d.plot2Dprofile(edidir, period_index_list, 'zxy')

pen2d.plot2Dprofile(edidir, period_index_list, 'zyx')




