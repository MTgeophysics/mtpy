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
from mtpy.utils.shapefiles_creator import plot_phase_tensor_ellipses_and_tippers

edi_dir = '/home/bren/mtpy/examples/data/edi_files_2'

out_dir = '/tmp'

# call the function over the periods index list
for pindex in range(0, 3):  # how many periods to do?
    plot_phase_tensor_ellipses_and_tippers(edi_dir, out_dir=out_dir, iperiod=pindex)
