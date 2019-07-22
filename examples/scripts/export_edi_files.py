#!/bin/env python
"""
Description:
    Export interpolated edi files
References:
 
CreationDate:   4/12/18
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     4/12/18   RH
    LastUpdate:     4/13/18   Alison Kirkby
"""


import os

import glob
import numpy as np

from mtpy.core.edi_collection import EdiCollection


# directory format for linux users
#edi_path = "../data/edi_files" 
#savepath '/tmp'

# directory format for windows users
edi_path = r'C:\mtpywin\mtpy\examples\data\edi_files' 
savepath = r'C:\tmp'


edi_files = glob.glob(os.path.normpath(os.path.abspath(os.path.join(edi_path, "*.edi"))))
edi_collection = EdiCollection(edi_files)

# period_list to interpolate onto
period_list = np.logspace(-4, 3, 43)

edi_collection.export_edi_files(savepath, 
                                period_list=period_list, # if not provided, will search edi files and find
                                                         # periods present in at least 10% of edis
                                period_buffer=2 # factor to stretch interpolation by. For example: if period_buffer=2
                                                 # then interpolated data points will only be included if they are
                                                 # within a factor of 2 of a true data point.
                                )
