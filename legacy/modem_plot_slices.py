#!/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Visualize Horizontal and Vertical Slices of the ModEM's output Model: *.dat and *.rho (same as *.ws) files
References:
    https://gajira.atlassian.net/browse/ALAMP-31
Usage Examples:
    python mtpy/imaging/modem_plot_slices.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho ns 300
    python mtpy/imaging/modem_plot_slices.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat  /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Model.ws z -1000 1000
CreationDate:   20/09/2017
Developer:      fei.zhang@ga.gov.au

LastUpdate:     13/10/2017   FZ
"""

import glob
import os
import sys

from mtpy.modeling.modem.data_model_analysis import DataModelAnalysis

#######################
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("********************************************")

        print(
            "USAGE: python %s %s %s %s %s"
            % (sys.argv[0], "file.dat", "file.rho", "[z|ns|ew]", "[list_of_location]")
        )

        print("********************************************")

        sys, exit(1)

    # Take commandline input
    if len(sys.argv) == 2:  # A model dir provided
        modeldir = sys.argv[1]
        datf = os.path.join(modeldir, "ModEM_Data.dat")
        rhofiles = glob.glob(os.path.join(modeldir, "*.rho"))

        # print(rhofiles)

        if len(rhofiles) < 1:
            print("No rho files found in the dir %s", modeldir)
            sys.exit(1)
        else:
            # the file with highest numbers in the last 3 numbers before *.rho
            rhof = sorted(rhofiles)[-1]

        print("Effective Files Used in the Plot...: ", datf, rhof)

    # dat and rho file both provided
    if len(sys.argv) >= 3:
        datf = sys.argv[1]
        rhof = sys.argv[2]

    slice_locs = []
    if len(sys.argv) >= 4:
        slice_orient = sys.argv[3]
        slice_locs = sys.argv[4:]  # a list of depth where h-slice to be visualized
        # slice_locs=[-2000, -1900, -1700, -1500, -1200, -1000, -800, -600, -400, -200,
        #             0, 20, 50, 80, 100,150, 200, 400, 600,800,1000,
        #             2000,3000,4000,5000,6000,7000,8000,9000,10000]

    if len(slice_locs) < 1:
        slice_locs = None

    # construct plot object
    # self = DataModelAnalysis(datf, rhof)  # default  map_scale='m')
    myObj = DataModelAnalysis(datf, rhof, map_scale="km")

    myObj.set_plot_orientation(slice_orient)  # horizontal at a given depth z
    # myObj.set_plot_orientation('ew')
    # myObj.set_plot_orientation('ns')

    # --------------------- visualize slices
    myObj.plot_multi_slices(slice_list=slice_locs)
