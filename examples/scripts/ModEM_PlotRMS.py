# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plots data/model misfit for a given period, at all sites, separately for each
of the impedance tensor modes + tipper

Revision History:
    brenainn.moushall@ga.gov.au 31-03-2020 13:42:41 AEDT:
        - Add option for plotting impedance or tippers
        - Add shapefile creation
        - Allow specifying period in seconds
        - Add plotting on geotiff background
"""
import os

os.chdir(r'C:\mtpywin\mtpy')

import os.path as op

import numpy as np

from mtpy.modeling.modem import PlotRMSMaps

wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'
savepath = r'C:\tmp'

filestem = op.join(wd,'Modular_MPI_NLCG_004')
resid_fn=op.join(wd,filestem + '.res')

# Parameter explanations (TODO: add to user guide):

# plot_elements: can plot only impedance or tippers by setting to
# 'impedance' or 'tippers'. Set as 'both' or leave out to plot impedance
# and tippers.

# bimg: path to a geotiff to use as map background image. If the CRS of
# the geotiff and the model differ, then 'PlotRMSMaps.model_epsg' must
# also be provided so model coordinates can be correctly projected onto
# the geotiff.

# period: can choose period by providing 'period_index' or by providing
# period in seconds to 'period'. If 'period' is provided it will take
# priority over 'period_index'. The closest available period will be
# selected, so the plotted period may be different from what was chosen.

probj = PlotRMSMaps(resid_fn,
                    # period=100.,  # can specify a period in seconds
                    period_index='all',
                    rms_cmap='jet', # choose matplotlib colormap or set to None
                    rms_max=5,
                    plot_elements='both',
                    # bimg=r'C:\path\to\a\background_image.tif'
                    )

# Can write RMS map as shapefiles by calling 'create_shapefiles'. This
# will use the period, plot_elements etc. attributes provided to the
# PlotRMSMaps class above.

# dst_epsg: the CRS of the shapefile output. This should be set to the
# EPSG code of the geotiff the shapefiles are intended to be displayed
# on.

# save_path: an optional save_path. If not provided, will save to the
# PlotRMSMaps.save_path atribute. Shapefiles will be saved within
# subdirectories labelled with component and period within this
# directory.

probj.create_shapefiles(dst_epsg=4326, save_path=savepath)

probj.save_figure(save_path=savepath,
                  save_fig_dpi = 400 # change to your preferred figure resolution
                  )
