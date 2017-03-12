# coding: utf-8
"""
Demo: Plot Phase Tensor Pseudo Section

run this script in spyder (ipython) will show a series of consistent plots
But run this script in commpandline will have varying sized figures
	python phase_tensor_pseudo_section_plot.py
"""
# In[ ]:

import os, sys
import glob

#1) set env variables PYTHONPATH and MTPYPATH before starting jupyther notebook,
# OR alternatively
#2) do the following two statements. Then it will all works fine.
#sys.path.insert(0,'/Softlab/Githubz/mtpy2')


import mtpy.imaging.mtplot as mtplot
from mtpy.imaging.mtplottools import MTArrows, MTEllipse
# from mtpy.tests.common import MTPY_DEVELOPMENT_TEST_DATA
from mtpy.imaging.phase_tensor_pseudo_section_plot import PlotPhaseTensorPseudoSection


# We have a list of edi files for this demo in the repo
edi_dir=sys.argv[1]
edi_file_list = glob.glob(edi_dir+'/*.edi')



print(edi_file_list)


# Define the ellipse and arrow properties
ellipse_dict = {'range': (20, 70), 'cmap': 'mt_bl2gr2rd',
                'colorby': 'phimin', 'size': 10}
ellipse = MTEllipse(ellipse_dict=ellipse_dict)
arrow = MTArrows({'size': 60, 'head_length': 4})



# Plot the phase tensor pseudo section
pt1 = PlotPhaseTensorPseudoSection(
    edi_file_list,
    data_type='z',
    ellipse=ellipse,
    arrow=arrow,
    fn_list=edi_file_list,
    tscale='frequency',
    ellipse_freq=1,  # plot an ellipse at every frequency value
    plot_tipper='yri',
    stretch=(1500, 35),
    scale_arrow=False,
    fig_size=[10, 12]
    )


# Change some properties and replot
pt1.ellipse_freq = 2  # plot ellipse at every second frequency
pt1.arrow.arrow_size = 50  # change arrow size to 50
pt1.fig_size=[10, 12]
pt1.font_size=14
pt1.redraw_plot()


# In[ ]:

# change arrow properties
arrow.arrow_size = 100
arrow.arrow_head_width = 3
arrow.arrow_head_length = 4

# plot every 3rd ellipse
pt1 = PlotPhaseTensorPseudoSection(
    edi_file_list,
    data_type='z',
    ellipse=ellipse,
    arrow=arrow,
    fn_list=edi_file_list,
    tscale='frequency',
    ellipse_freq=5,  #=3 plot every 3rd ellipse
    plot_tipper='yri', # plot real and imaginary tipper arrows
    stretch=(1500, 35),
    scale_arrow=False,
    fig_size=[10, 12],
    font_size=14,
    )



# Colorby 'skew'
ellipse.ellipse_colorby = 'skew'
arrow.arrow_size = 50
# plot every 3rd ellipse
pt1 = PlotPhaseTensorPseudoSection(
    edi_file_list,
    data_type='z',
    ellipse=ellipse,
    arrow=arrow,
    fn_list=edi_file_list,
    tscale='frequency',
    ellipse_freq=3,  # plot every 3rd ellipse
    plot_tipper='yri', # plot real and imaginary tipper arrows
    stretch=(1500, 35),
    scale_arrow=False,
    fig_size=[10, 12],
    font_size=14
   )


# Colorby 'normalized_skew'
ellipse.ellipse_colorby = 'normalized_skew'
# change arrow size
arrow.arrow_size = 40
# plot every 4th ellipse
pt1 = PlotPhaseTensorPseudoSection(
    edi_file_list,
    data_type='z',
    ellipse=ellipse,
    arrow=arrow,
    fn_list=edi_file_list,
    tscale='frequency',
    ellipse_freq=4,  # plot every 4th ellipse
    plot_tipper='yri', # plot real and imaginary tipper arrows
    stretch=(1500, 35),
    scale_arrow=False,
    fig_size=[10, 12],
    font_size=14
    )


# Colorby 'ellipticity'
ellipse.ellipse_colorby = 'ellipticity'

# plot every 4th ellipse
pt1 = PlotPhaseTensorPseudoSection(
    edi_file_list,
    data_type='z',
    ellipse=ellipse,
    arrow=arrow,
    fn_list=edi_file_list,
    tscale='frequency',
    ellipse_freq=4,  # plot every 4th ellipse
    plot_tipper='yri', # plot real and imaginary tipper arrows
    stretch=(1500, 35),
    scale_arrow=False,
    fig_size=[10, 12],
    font_size=14
    )

