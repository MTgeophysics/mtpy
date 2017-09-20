# -*- coding: utf-8 -*-
# !/bin/env python
"""
Description:
    Visualize Horizontal and Vertical Slices of the ModEM's output Model: *.dat and *.rho (same as *.ws) files

Usage:
    python mtpy/imaging/modem_plot_slices.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho 300
    python mtpy/imaging/modem_plot_slices.py /e/tmp/GA_UA_edited_10s-10000s_16/ModEM_Data.dat  /e/tmp/GA_UA_edited_10s-10000s_16/ModEM_Model.ws -1000 1000

CreationDate:   20/09/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     8/09/2017   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import sys
import glob
from mtpy.modeling.modem_output_to_views import ModemSlices

#######################
if __name__ == "__main__":

    # Take commandline input
    if (len(sys.argv) == 2):  # A model dir provided
        modeldir = sys.argv[1]
        datf = os.path.join(modeldir, 'ModEM_Data.dat')
        rhofiles = glob.glob(os.path.join(modeldir, '*.rho'))

        print(rhofiles)

        if len(rhofiles) < 1:
            print ("No rho files found in the dir %s", modeldir)
            sys.exit(1)
        else:
            # the file with highest numbers in the last 3 numbers before *.rho
            rhof = sorted(rhofiles)[-1]

        print("Effective Files Used in Plot: ", datf, rhof)

    # dat and rho file both provided
    if len(sys.argv) >= 3:
        datf = sys.argv[1]
        rhof = sys.argv[2]

    # construct plot object
    #self = ModemSlices(datf, rhof)  # default  map_scale='m')
    myObj = ModemSlices(datf, rhof, map_scale='km')

#--------------------- visualize slices:
    # myObj.set_plot_orientation('ew')
    # myObj.set_plot_orientation('ns')
    # horizontal at a given depth z
    myObj.set_plot_orientation('z')


    if len(sys.argv) >= 4:
        slice_locs = sys.argv[3:] # a list of depth where h-slice to be visualized
        # slice_locs=[-2000, -1900, -1700, -1500, -1200, -1000, -800, -600, -400, -200,
        #             0, 20, 50, 80, 100,150, 200, 400, 600,800,1000,
        #             2000,3000,4000,5000,6000,7000,8000,9000,10000]
    else:
        slice_locs= None

    myObj.plot_multi_slices(slice_list=slice_locs)
