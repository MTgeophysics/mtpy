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
import tempfile

from mtpy.core.edi_collection import EdiCollection
from mtpy.imaging import penetration_depth3d as pen3d


# https://www.geeksforgeeks.org/python-handling-recursion-limit/
# sys.setrecursionlimit(1000000)

try:
    # PACK_ROOT = os.environ['PACK_ROOT']
    # mtpy_path = os.path.join(PACK_ROOT, 'mtpy')
    mtpy_path = os.environ['MTPY_ROOT']
except:
    print("Warn: The environment variable MTPY_ROOT is not defined. We will guess")
    mtpy_path = os.path.abspath('../..')

if not os.path.isdir(mtpy_path):
    raise Exception("the guessed mtpy dir %s is not a folder!"% mtpy_path)

edidir = os.path.join(mtpy_path,'examples/data/edi2')
edidir = r"C:\Users\u25656\Desktop\Wenping_EDI132\MT086_Edited_EDIs"
edidir = "/Datasets/MTWorkflow/Wenping_EDI132/MT086_Edited_EDIs"
# edidir = r"C:\Githubz\mtpy\examples\data\edi2"

# If you have own edi change the variable below according to your edi files folder
# USE / for Unix and Win-Dos compatible
# edidir = r'C:/mtpywin/mtpy/Alison_penetrationDepth3D/EDI_files'  # unequal periods  across EDI files
# edidir = '/g/data/ha3/fxz547/Data/3D_MT_data_edited_fromDuanJM'


temp_dir = tempfile.gettempdir()
print('Using temporary directory ' + temp_dir)
savepath = temp_dir
# savepath = r'C:/tmp'

if not os.path.isdir(edidir):
    print ("please provide the path to edi folder")
    sys.exit(1)


edifiles = glob.glob(os.path.join(edidir, "*.edi"))


print("***************** Stats analysis on the EDI files [periods list]  *******************************")

edis_obj = EdiCollection(edilist=edifiles)
edis_obj.select_periods(percentage=5.0)

print("*****************  the periods with high occurance in EDI files **************************")
per_freq = edis_obj.get_periods_by_stats(percentage=99.0)  # occur in 99% of EDI files
print(per_freq)
for aper in per_freq:
    print (aper, edis_obj.get_period_occurance(aper))

print("***************** For each uniqu period, show it's occurance in the EDI files ************")
for aper in edis_obj.all_unique_periods:
    print(aper, edis_obj.get_period_occurance(aper))

print("****** To do the pen3d plotting, pick a value from the above periods, which have high occurance (~100%) ******************")
# pen3d plot had issues introduced in penetration.py  https://github.com/MTgeophysics/mtpy/commit/817c9f4a6384460d57974fb7a6f80e04a8a4ce97

# Create plot for a period index number 1,2,3 ..., 10 for determinant. This may not make sense sometimes.
# change to your preferred file resolution
#pen3d.plot_latlon_depth_profile(edidir, 4, 'det', showfig=True, savefig=True, savepath=savepath, fig_dpi=400)

# The recommended way is to use a float value for the second argument, which is frequency (NOT period).
# provide a float value (e.g. 100.0), which should be in the 1/FREQ  of the EDI files.

#pen3d.plot_latlon_depth_profile(edidir, 25.2, savefig=True, savepath=savepath, fig_dpi=400)

#pen3d.plot_latlon_depth_profile(edidir,  9.372, savefig=True, savepath=savepath, fig_dpi=400)  # period =0.1067s ec
# limit hit pen3d.plot_latlon_depth_profile(edidir,  0.01, savefig=True, savepath=savepath, fig_dpi=400)


# User MUST provide period value in seconds ("/Datasets/MTWorkflow/Wenping_EDI132/MT086_Edited_EDIs")
#pen3d.plot_latlon_depth_profile(edidir,  0.1067, savefig=True, savepath=savepath, fig_dpi=400)  # period =0.1067s ec
pen3d.plot_latlon_depth_profile(edidir,  95.33, savefig=True, savepath=savepath, fig_dpi=400)
pen3d.plot_latlon_depth_profile(edidir,  48.9956, savefig=True, savepath=savepath, fig_dpi=400)
pen3d.plot_latlon_depth_profile(edidir,  10.839, savefig=True, savepath=savepath, fig_dpi=400)


