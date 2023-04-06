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

edi_path = Path(r"c:\Users\jpeacock\OneDrive - DOI\EDI_FILES")

edi_list = edi_path.glob("LV*.edi")

md = MTData()

for edi_fn in edi_list:
    mt_obj = MT()
    mt_obj.read_tf_file(edi_fn)

    md.add_station(mt_obj)
