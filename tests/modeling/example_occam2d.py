# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 19:13:24 2023

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

from mtpy import MT
from mtpy.core.mt_data import MTData

# =============================================================================

## load in data
edi_path = Path(r"c:\Users\jpeacock\OneDrive - DOI\EDI_FILES")
edi_list = edi_path.glob("LV*.edi")

md = MTData()

for edi_fn in edi_list:
    mt_obj = MT()
    mt_obj.read_tf_file(edi_fn)

    md.add_station(mt_obj)


## generate a profile based on station locations
x1, y1, x2, y2, profile_from_stations = md.generate_profile()

## get stations within a 2000 meter radius of profile line
profile_md = md.get_profile(x1, y1, x2, y2, 2000)


## generate a profile based on geoelectrical strike
x1, y1, x2, y2, profile_from_strike = profile_md.generate_profile_from_strike(
    70
)

## get stations
strike_md = profile_md.get_profile(x1, y1, x2, y2, None)
