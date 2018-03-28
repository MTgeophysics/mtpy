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
from mtpy.imaging import penetration_depth3d as pen3d

# change the variable below according to your edi files folder !!!
edidir = r'E:/mtpywin/mtpy/examples/data/edi2'  # / is Unix and Win-Dos compatible
# or get this variable from the cmdline:  edidir = sys.argv[1]

if not os.path.isdir(edidir):
    print ("please provide the path to edi folder")
    sys.exit(1)

# provide the index of period, eg 10
pen3d.plot_latlon_depth_profile(edidir, 10, 'det', showfig=True, savefig=False)

# OR provide a period value, which must be identified by user according to the EDI files.
pen3d.plot_latlon_depth_profile(edidir, 0.009846, savefig=False)


