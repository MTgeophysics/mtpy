#! /usr/bin/env python
"""
Description:
    Example python script to create grid formats (geotiff and ascii) from ModEM model file

References:


CreationDate:   6/11/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     6/11/2019   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""
import sys
from mtpy.utils.convert_modem_data_to_geogrid import create_geogrid

# =============================================
# quick test of this script
#  python examples/scripts/create_geogrid_from_resistivity_model.py  /c/Data/JinMing_GridData_sample/JM_model_002/EFTF_NLCG_002.dat /c/Data/JinMing_GridData_sample/JM_model_002/EFTF_NLCG_002.rho /c/temp
#  python examples/scripts/create_geogrid_from_resistivity_model.py  /c/Data/JinMing_GridData_sample/EFTF_MT_model/EF_NLCG_001.dat  /c/Data/JinMing_GridData_sample/EFTF_MT_model/EF_NLCG_001.rho /c/temp
#  python examples/scripts/create_geogrid_from_resistivity_model.py  /c/Data/Alison_201910/Alison_ModEM_Grid/MT075_ModEM_files/ModEM_Data.dat  /c/Data/Alison_201910/Alison_ModEM_Grid/MT075_ModEM_files/Modular_MPI_NLCG_004.rho /c/temp/
# ---------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("USAGE: python %s dat_file rho_file out_dir" % sys.argv[0])
        sys.exit(1)

    dat_file = sys.argv[1]
    rho_file = sys.argv[2]
    out_dir = sys.argv[3]
    # Before calling the function create_geogrid(), a user should
    # provide the right optional parameters.
    # Otherwise, default parameters will be used, which may not make sense
    kwargs = {
        "xpad": 6,
        "ypad": 6,
        "zpad": 10,
        "grid_size": 7500,
        "center_lat": None,  # to override what is found in the file.dat, which may be incorrect.
        "center_lon": None,  # to override what is found in the file.dat, which may be incorrect.
        "source_proj": None,
        "depth_index": None,  # slices [0,1,2,10] to be output
    }

    print("User Options:", kwargs)

    create_geogrid(dat_file, rho_file, out_dir, **kwargs)  # ,depth_index=[0,1,2,10])
