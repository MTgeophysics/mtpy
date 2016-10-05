# coding: utf-8
"""
Example
=========

	* plotting phase tensor maps using mtpy.imaging.mtplot.plot_pt_maps

JP 2016
"""
# Import necessary modules
#----------------------------------------
import os
import mtpy.imaging.mtplot as mtplot
#----------------------------------------

# create a list of .edi files
edi_path = os.path.join(os.getcwd(), 'data', 'HalfSpaceSQC')
edi_list = [os.path.join(edi_path, edi_fn) for edi_fn in os.listdir(edi_path) if edi_fn.endswith('.edi')]

# create a pt map object
ptm = mtplot.plot_pt_map(fn_list=edi_list)

# to learn more about plot_pt_map object
# type in an python console >>> help(ptm)

# Adjust parameters to make it look nice
#----------------------------------------
ptm.fig_size = [6, 6]      # figure size so it fits in the window
ptm.ellipse_size = .005    # ellipse size
ptm.xpad = -.035           # padding between last ellipse and plot border in x-direction
ptm.ypad = -.035           # padding between last ellipse and plot border in y-direction
ptm.plot_freq = 10         # change the plot frequency to a value that shows something interesting
ptm.fig_num = 1
ptm.redraw_plot()

# Add in tipper data
#----------------------------------------
ptm.plot_tipper = 'yri'              # plot real 'r' and imaginary 'i' parts
ptm.arrow_head_length = .001         # arrow head length
ptm.arrow_head_width = .0005         # arrow head width
ptm.arrow_head_height = .0005        # arrow head height
ptm.arrow_lw = .5                    # arrow line width 
ptm.arrow_size = ptm.ellipse_size    # arrow size
ptm.arrow_legend_xborderpad = .0075  # distance from border to arrow legend in x-direction
ptm.arrow_legend_yborderpad = .0075  # distance from border to arrow legend in y-direction 
ptm.arrow_legend_fontpad = .005      # distance between arrow and text in legend
ptm.fig_num = 2       
ptm.redraw_plot()

# Plot phase tensor skew
#----------------------------------------
ptm.ellipse_colorby = 'skew_seg'       # skew angle that is segmented into color bins
ptm.ellipse_cmap = 'mt_seg_bl2wh2rd'   # blue to white to red segmented colors
ptm.ellipse_range = (-6, 6, 2)         # range of skew values
ptm.fig_num = 3        
ptm.redraw_plot()

