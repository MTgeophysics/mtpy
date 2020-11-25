# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:51:02 2015

@author: Alison Kirkby

get dimensionality/strike angle from an edi file

"""
import os

os.chdir(r"C:\mtpywin\mtpy")  # change to path where mtpy is installed

from mtpy.core.mt import MT
from mtpy.core.z import Z, Tipper
from mtpy.analysis.geometry import dimensionality


# directory containing edis
edi_path = r"C:\mtpywin\mtpy\examples\data\edi_files_2"
savepath = r"C:/tmp"

# edi file name
edi_file = os.path.join(edi_path, "Synth00.edi")

# read edi file into an MT object
mtObj = MT(edi_file)

# use the phase tensor to determine which frequencies are 1D/2D/3D
dim = dimensionality(
    z_object=mtObj.Z,
    skew_threshold=5,  # threshold in skew angle (degrees) to determine if data are 3d
    eccentricity_threshold=0.1,  # threshold in phase ellipse eccentricity to determine if data are 2d (vs 1d)
)

# create a True/False array to mask with
mask = dim < 3

new_Z_object = Z(
    z_array=mtObj.Z.z[mask], z_err_array=mtObj.Z.z_err[mask], freq=mtObj.Z.freq[mask]
)

new_Tipper_object = Tipper(
    tipper_array=mtObj.Tipper.tipper[mask],
    tipper_err_array=mtObj.Tipper.tipper_err[mask],
    freq=mtObj.Tipper.freq[mask],
)

mtObj.write_mt_file(
    save_dir=savepath,
    fn_basename="Synth00_new",
    file_type="edi",  # edi or xml format
    new_Z_obj=new_Z_object,  # provide a z object to update the data
    new_Tipper_obj=new_Tipper_object,  # provide a tipper object to update the data
    longitude_format="LONG",  # write longitudes as 'LON' or 'LONG'
    latlon_format="dd"  # write as decimal degrees (any other input
    # will write as degrees minutes seconds
)
