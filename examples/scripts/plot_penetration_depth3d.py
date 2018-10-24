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
    LastUpdate:     23/10/2018   FZ

"""
import glob
import os
import sys

from mtpy.core.edi_collection import EdiCollection
from mtpy.imaging import penetration_depth3d as pen3d

# change the variable below according to your edi files folder
# USE / for Unix and Win-Dos compatible
edidir = r'C:/mtpywin/mtpy/Alison_penetrationDepth3D/EDI_files'  # unequal periods  across EDI files
# edidir = r"E:/Data/MT_Datasets/Cloncurry"  # unequal periods across EDI files
# edidir = r'C:/mtpywin/mtpy/examples/data/edi2' # equal periods  across EDI files

savepath = r'C:/tmp'

if not os.path.isdir(edidir):
    print ("please provide the path to edi folder")
    sys.exit(1)

edifiles = glob.glob(os.path.join(edidir, "*.edi"))

# Create plot for a period index number 1,2,3 ..., 10 for determinant. This may not make sense sometimes.
# change to your preferred file resolution
pen3d.plot_latlon_depth_profile(edidir, 4, 'det', showfig=True, savefig=True, savepath=savepath, fig_dpi=400)

# The recommended way is to use a float value for the period.
# provide a float value (e.g. 100.0) of the period, which should be in the EDI files.

pen3d.plot_latlon_depth_profile(edidir, 25.2, savefig=True, savepath=savepath, fig_dpi=400)
# pen3d.plot_latlon_depth_profile(edidir, 100.0, savefig=True, savepath=savepath, fig_dpi=400)

#############################################################
# More Stats analysis on the EDI files periods

edis_obj = EdiCollection(edilist=edifiles)
edis_obj.select_periods(percentage=5.0)

per_freq = edis_obj.get_periods_by_stats(percentage=80.0)
print(per_freq)

for aper in edis_obj.all_unique_periods[:3]:
    print(aper, edis_obj.get_period_occurance(aper))
    pen3d.plot_latlon_depth_profile(edidir, aper, showfig=True, savefig=True, savepath=savepath,
                                    fig_dpi=400)  # change to your preferred file resolution
