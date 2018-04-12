#!/bin/env python
"""
Description:
    Export interpolated edi files
References:
 
CreationDate:   4/12/18
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     4/12/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mtpy.core.edi_collection import is_num_in_seq, EdiCollection
from mtpy.core.mt import MT
import os
import glob
import numpy as np

edi_path = "../data/edi_files"

edi_files = glob.glob(os.path.normpath(os.path.abspath(os.path.join(edi_path, "*.edi"))))
edi_collection = EdiCollection(edi_files)

period_list = np.logspace(-2, 3, 15)

edi_collection.export_edi_files('/tmp', period_list=period_list)