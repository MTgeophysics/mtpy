#! /usr/bin/env python
"""
Description:
    plot phase tensor ellipses and tipers into one figure

CreationDate:   1/11/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     1/11/2018   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os

from mtpy.utils.shapefiles_creator import plot_phase_tensor_ellipses_and_tippers

# =============================================
# Begin to run this script

# Define a folder of edi files.
# edi_dir = r'C:\mtpywin\mtpy\examples\data\edi_files_2'
edi_dir = 'E:/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM'
# edi_dir = r"E:\Data\MT_Datasets\GA_UA_edited_10s-10000s"
# edi_dir = r"E:\Data\MT_Datasets\Isa_EDI_edited_10Hz_1000s"
# edi_dir = r"E:\Data\MT_Datasets\728889\EDI_files" # narrow area, not shown
# edi_dir = r"E:\Data\MT_Datasets\75098\EDI_files"  # cross UTM zones 51 and 52, No tippers?
# edi_dir =r"E:\Data\MT_Datasets\75099_Youanmi\EDI_Files_edited"
# edi_dir = sys.argv[1]

# set output dir
out_dir='C:/temp2'


# call the function over the periods index list
for pindex in range(0, 3):  # how many periods to do?
    pngfile = os.path.join(out_dir, "phase_tensor_tipper_%s.png" % pindex)
    plot_phase_tensor_ellipses_and_tippers(edi_dir, outfile=pngfile, iperiod=pindex)
