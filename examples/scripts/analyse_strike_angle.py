# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:51:02 2015

@author: Alison Kirkby

get dimensionality/strike angle from an edi file

"""
import os.path as op
import os
os.chdir(r'C:\mtpywin\mtpy') # change to path where mtpy is installed
import numpy as np

from mtpy.core.mt import MT
from mtpy.analysis.geometry import dimensionality,strike_angle


# directory containing edis
edi_path = r'C:\mtpywin\mtpy\examples\data\edi_files_2'

# edi file name
edi_file = os.path.join(edi_path,'Synth00.edi')

# read edi file into an MT object
mtObj = MT(edi_file)

# determine strike angle for the 2d parts of impedance tensor
strike = strike_angle(z_object=mtObj.Z,
                     skew_threshold=5, # threshold in skew angle (degrees) to determine if data are 3d
                     eccentricity_threshold=0.1 # threshold in phase ellipse eccentricity to determine if data are 2d (vs 1d)
                      )

print(strike)