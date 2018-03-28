#! /usr/bin/env python
"""
Description:
    A command line python script to show how to plot phase tensor map

How2Run quick test this script in command line
# examples:
# python mtpy/imaging/phase_tensor_maps.py  /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 /e/MTPY2_Outputs
# OR Windows Prompt, where path begin with a drive letter such as C:
# python mtpy/imaging/phase_tensor_maps.py  C:/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 C:/MTPY2_Outputs

CreationDate:   29/03/2018
Developer:      fei.zhang@ga.gov.au

Revision History:

"""

# import section
import sys
import glob
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps

#  =======================================================================
if __name__ == "__main__":

# 1 get a folder of edi files
    edidir = sys.argv[1]  #'E:/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/'  # Windows
    edi_file_list = glob.glob(edidir + '/*.edi')
    print ("Edi files: ", edi_file_list)

# 2 get a frequency value (from the edi files)
    plot_freq = float(sys.argv[2])  # 159.0

# 3 get the output dir
    savedir = sys.argv[3]


    ptm_obj = PlotPhaseTensorMaps(fn_list=edi_file_list,
                        plot_freq=plot_freq,
                        ftol=.2,
                        xpad=0.02,
                        plot_tipper='yr',
                        edgecolor='k',
                        lw=0.1,
                        alpha=1,
                        minorticks_on=False,
                        ellipse_size=.02,
                        ellipse_range=[-10, 10, 2],
                        ellipse_colorby='skew',
                        arrow_size = 0.05,
                        ellipse_cmap='mt_seg_bl2wh2rd')

    ptm_obj.export_params_to_file(save_path=savedir)
