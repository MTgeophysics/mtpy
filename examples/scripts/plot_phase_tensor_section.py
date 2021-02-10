# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period) 
"""

from mtpy.imaging.phase_tensor_pseudosection import PlotPhaseTensorPseudoSection
import os.path as op
import os

# path to edis
edi_path = r"c:\Users\jpeacock\Documents\GitHub\mtpy\examples\data\edi_files_2"

# save path
savepath = r"C:\tmp"

# edi list
elst = [
    op.join(edi_path, edi) for edi in os.listdir(edi_path) if ((edi.endswith(".edi")))
]  # and edi.startswith('GB')

# create a plot object
plotObj = PlotPhaseTensorPseudoSection(fn_list = elst,
                                 linedir='ns', # 'ns' if the line is closer to north-south, 'ew' if line is closer to east-west
                                 stretch=(17,8), # determines (x,y) aspect ratio of plot
                                 station_id=(0,10), # indices for showing station names
                                 plot_tipper = 'yri', # plot tipper ('y') + 'ri' means real+imag
                                 font_size=5,
                                 lw=0.5,
                                 rotation_angle=20,
                                 ellipse_dict = {'ellipse_colorby':'skew_seg',# option to colour by phimin, phimax, skew, skew_seg
                                                 'ellipse_range':[-12, 12, 3]} # set color limits - default 0,90 for phimin or max,
                                                                         # [-12,12] for skew. If plotting skew_seg need to provide
                                                                         # 3 numbers, the 3rd indicates interval, e.g. [-12,12,3]
                                 )

# update ellipse size (tweak for your dataset)
plotObj.ellipse_size = 2.5

plotObj.plot()

# plotObj.save_figure(save_fn = op.join(savepath,'PhaseTensorSection.png'),
#                    fig_dpi=400) # change to your preferred file resolution
