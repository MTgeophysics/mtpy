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

edi_collection.create_phase_tensor_csv(savepath,
                                       period_list=period_list # list of periods, if provided it 
                                                        # interpolates onto those frequencies, 
                                                        # if blank it exports all periods in file.
                                       )
