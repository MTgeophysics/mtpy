#! /usr/bin/env python
"""
Description:
    Generate shapefiles for phase tensors and tippers from EDI data for
    display in a GIS viewer.

CreationDate:   1/11/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     1/11/2018   FZ
    brenainn.moushall@ga.gov.au 27-03-2020 17:34:34 AEDT:  
        Rename script, change description and add interface for creating 
        shapefiles from EDI data. Plotting as image not working but left
        as commented code.
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
# displaying on).
dst_epsg = 4326
# List of periods in seconds to plot. Can also provide a single value.
# The nearest available period will be selected.
periods = [0., 100.]

# Create and save shapefiles.
create_tensor_tipper_shapefiles(edi_dir, out_dir, periods, src_epsg, dst_epsg)

# Plots the shapefiles as .png. Currently not working due to absence
# of 'descartes' library.
# from mtpy.utils.shapefiles_creator import plot_phase_tensor_ellipses_and_tippers
# for p in periods:
#    plot_phase_tensor_ellipses_and_tippers(edi_dir, out_dir=out_dir, iperiod=p)
