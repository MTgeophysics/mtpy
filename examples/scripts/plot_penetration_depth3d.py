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
import glob
from mtpy.core.edi_collection import EdiCollection
from mtpy.imaging import penetration_depth3d as pen3d

# change the variable below according to your edi files folder !!!
# USE / for Unix and Win-Dos compatible
edidir = r'C:\mtpywin\mtpy\Alison_penetrationDepth3D\EDI_files'
#edidir = r'C:/mtpywin/mtpy/examples/data/edi2'
# or get this variable from the cmdline:  edidir = sys.argv[1]
savepath = r'C:\tmp'


if not os.path.isdir(edidir):
    print ("please provide the path to edi folder")
    sys.exit(1)

edifiles = glob.glob(os.path.join(edidir, "*.edi"))
edis_obj = EdiCollection(edilist= edifiles)

# edis_obj.select_periods(percentage=5.0)

#per_freq = edis_obj.get_periods_by_stats(percentage=80.0)
#print(per_freq)

# for aper in edis_obj.all_unique_periods[:2]:
#     print(aper, edis_obj.get_period_occurance(aper))
#     pen3d.plot_latlon_depth_profile(edidir, aper, savefig=True, savepath=savepath,
#                                     fig_dpi=400)  # change to your preferred file resolution

# provide the index of period - must be an integer (e.g. 10)
# change to your preferred file resolution
pen3d.plot_latlon_depth_profile(edidir, 4, 'det', showfig=True, savefig=True, savepath=savepath, fig_dpi=400)
# OR
# provide a period value, which must be identified by user according to the EDI files.
# must be a float e.g. 10.0 or 1.5

#pen3d.plot_latlon_depth_profile(edidir, 25.2, savefig=True, savepath=savepath, fig_dpi=400) # change to your preferred file resolution
#pen3d.plot_latlon_depth_profile(edidir, 100.0, savefig=True, savepath=savepath, fig_dpi=400) # change to your preferred file resolution


