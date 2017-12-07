#!/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Compute Phase Tensors from ModEM Dat File and output to CSV files

Usage Examples:
    python scripts/modem_data_to_phase_tensor.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat [OutDir]
    python scripts/modem_data_to_phase_tensor.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat [OutDir]


Developer:      fei.zhang@ga.gov.au
LastUpdate:     08/09/2017
LastUpdate:     05/12/2017 FZ moved the function into the module mtpy.modeling.modem.Data
"""

import sys
from mtpy.modeling.modem import Data
from mtpy.mtpy_globals import NEW_TEMP_DIR

if __name__ == "__main__":

    file_dat = sys.argv[1]
    if len(sys.argv)>2:
        outdir = sys.argv[2]
    else:
        outdir=NEW_TEMP_DIR

    obj = Data()
    
    obj.compute_phase_tensor(file_dat, outdir)
