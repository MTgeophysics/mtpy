#! /usr/bin/env python
"""
Description:
    Generate shapefiles for phase tensors and tippers from EDI data for
    display in a GIS viewer.

CreationDate:   1/11/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     1/11/2018   FZ
    LastUpdate:     27/03/2020  BM - create interface for drawing
                                     shapefiles from EDI data.
"""
# Change working directory if running on windows and mtpy is not
# installed as a package.
# import os
# os.chdir(r'C:\mtpywin')

from mtpy.utils.shapefiles_creator import create_tensor_tipper_shapefiles

edi_dir = '/path/to/edi/data'
out_dir = '/path/to/output/directory'
# EPSG code of the EDI data.
src_epsg = 4326
# EPSG code of the output (i.e. same CRS as the tiff you will be
# displaying on.
dst_epsg = 4326
# List of indicies for periods to plot
periods = [0, 1, 2]

# Create and save shapefiles.
create_tensor_tipper_shapefiles(edi_dir, out_dir, src_epsg, dst_epsg, period_indicies=periods)

# Plots the shapefiles as .png. Currently not working due to absence
# of 'descartes' library.
# from mtpy.utils.shapefiles_creator import plot_phase_tensor_ellipses_and_tippers
# for p in periods:
#    plot_phase_tensor_ellipses_and_tippers(edi_dir, out_dir=out_dir, iperiod=p)
