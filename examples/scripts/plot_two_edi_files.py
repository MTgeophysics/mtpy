#! /usr/bin/env python
"""
Description:
    Example script to plot the data of two edi files to compare them

CreationDate:   15/11/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     15/11/2019   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""
import sys
from mtpy.core.mt import MT
from mtpy.imaging.plot_mt_response import PlotMTResponse

# ======================================================================================================================
# How to Run this script:
# python examples/scripts/plot_two_edi_files.py /c/Githubz/mtpy/data/edifiles/15125A.edi
# OR provide edi files
# python examples/scripts/plot_two_edi_files.py /c/Githubz/mtpy/data/edifiles/15125A.edi /c/Githubz/mtpy/data/edifiles/15129A.edi
# ======================================================================================================================
if __name__ == "__main__":

    if len(sys.argv)<2:
        print("USAGE: python %s edifile1 [edifile2]"%sys.argv[0])
        sys.exit(1)
    elif len(sys.argv) == 2: # one edi file provided
        edifile = sys.argv[1]  # /c/Githubz/mtpy/data/edifiles/15125A.edi
        mt_obj = MT(edifile)

        rp1 = PlotMTResponse(z_object=mt_obj.Z,  # this is mandatory
                             # t_object=mt_obj.Tipper,
                             # pt_obj=mt_obj.pt,
                             station=mt_obj.station,
                             #plot_tipper='yr',  # plots the real part of the tipper
                             plot_num=3)  # plot_num =1 xy + yx; 2 all; 3 xy yx det

        # rp1.xy_color = (.5,.5,.9)
        # rp1.xy_marker = '*'
        # rp1.redraw_plot()
    elif(len(sys.argv)==3):    # overlay 2 edi files provided
        edifile = sys.argv[1]  #
        mt_obj = MT(edifile)
        edifile2 = sys.argv[2]  # /c/Githubz/mtpy/data/edifiles/15126A.edi
        mt_obj2 = MT(edifile2)

        rp1 = PlotMTResponse(z_object=mt_obj.Z,  # this is mandatory
                             station=mt_obj.station,
                             plot_yn='n',plot_num=2, # plot_num =1 xy + yx; =2 all 4 components; =3 xy yx det
                             # pt_obj=mt_obj.pt,  # ellispses
                             # t_object=mt_obj.Tipper,
                             #plot_tipper='yr'  # plots the real part of the tipper
                             )


        rp1.station = rp1.station + " and " + mt_obj2.station
        rp1.plot(overlay_mt_obj=mt_obj2)

