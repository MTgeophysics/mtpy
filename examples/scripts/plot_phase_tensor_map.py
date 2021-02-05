# -*- coding: utf-8 -*-
"""
Plots phase tensor ellipses as a map for a given frequency.

Created on Fri Nov 22 07:29:58 2013
@author: Alison Kirkby

Revision History:
    brenainn.moushall@ga.gov.au 26-03-2020 15:05:46 AEDT:
        Add plotting of geotiff as basemap background.
"""
import os
import os.path as op

import mtpy.imaging.phase_tensor_maps as pptmaps

# directory containing edis

#edipath = r'C:\mtpywin\mtpy\examples\data\edi2'
edipath = '/mtpy/examples/data/ET_edi'

# whether or not to save the figure to file
save = False

# full path to file to save to
savepath = '/tmp'

# frequency to plot
plot_freq = 1e-2

# value to color ellipses by, options are phimin,phimax,skew
colorby = 'strike'
ellipse_range = [-9, 9]
rotation_angle = 0

image_fn = 'phase_tensor_map%1is_' % (int(1. / plot_freq)) + colorby + '.png'

# gets edi file names as a list
elst = [op.join(edipath, f) for f in os.listdir(edipath) if f.endswith('.edi')]

# Notes about background image:
#  bimg_band: setting as 'None' will read all bands. If the raster is
#  RGB/A, then leave this as 'None' (cmap will also be ignored). If the
#  data is singleband, then 'None' will also retrieve the only available
#  band. If you need to display a specific band, provide the band number
#  (with '1' being the first band.)
m = pptmaps.PlotPhaseTensorMaps(fn_list=elst,
                                plot_freq=plot_freq,
                                fig_size=(6, 4),
                                ftol=.1,
                                xpad=0.02,
                                plot_tipper='yr',
                                edgecolor='k',
                                lw=0.5,
                                minorticks_on=False,
                                ellipse_colorby=colorby,
                                ellipse_range=ellipse_range,
                                ellipse_size=0.01,
                                arrow_head_width=0.002,
                                arrow_head_length=0.002,
                                #ellipse_cmap='mt_seg_bl2wh2rd'
                                station_dict={'id': (5, 7)},
                                background_image='/mtpy/examples/data/gravity/ET_gravity.tif'
 #                               bimg_band=1, #  Optional, set to None by default
                                #bimg_cmap='viridis'  # Optional, set to 'viridis' by default
                                rotation_angle=rotation_angle)


if save:
    m.save_figure(op.join(savepath, image_fn),
                  fig_dpi=400)  # change to your preferred file resolution
    print("Complete! Figure saved to {}".format(op.join(savepath, image_fn)))
else:
    print("Complete!")
