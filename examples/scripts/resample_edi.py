#!/bin/env python
"""
Description:
    Read and write and edi file
References:
 
CreationDate:   4/12/18
Developer:      Alison.Kirkby@ga.gov.au
 
Revision History:
    LastUpdate:     5/15/18   Alison Kirkby
"""


import os
os.chdir(r'C:/mtpywin/mtpy') # change to path where mtpy is installed


from mtpy.core.mt import MT
from mtpy.utils.calculator import get_period_list


# directory format for windows users
edi_path = r'C:\mtpywin\mtpy\examples\data\edi_files_2' 
savepath = r'C:\tmp'


edi_file = os.path.join(edi_path,'Synth00.edi')

mtObj = MT(edi_file)

#new_freq_list = mtObj.Z.freq # takes every second frequency
new_freq_list = 1./get_period_list(1e-4,1e3,5) # 5 periods per decade from 0.0001 to 100000 s

# create new Z and Tipper objects containing interpolated data
new_Z_obj, new_Tipper_obj = mtObj.interpolate(new_freq_list)

# write a new edi file using the new data
mtObj.write_mt_file(save_dir=savepath, 
                    fn_basename='Synth00_new', 
                    file_type='edi', # edi or xml format
                    new_Z_obj=new_Z_obj, # provide a z object to update the data
                    new_Tipper_obj=new_Tipper_obj, # provide a tipper object to update the data
                    longitude_format='LONG', # write longitudes as 'LON' or 'LONG'
                    latlon_format='dd' # write as decimal degrees (any other input
                                       # will write as degrees minutes seconds
                    )