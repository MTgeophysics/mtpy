# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:51:02 2015

@author: Alison Kirkby

strike analysis using roseplot function

"""
import os.path as op
import os
# os.chdir(r'C:\mtpywin\mtpy') # change to path where mtpy is installed
os.chdir(r'c:\Users\jpeacock\Documents\GitHub\mtpy\data\edifiles2') # change to path where mtpy is installed

from mtpy.imaging.plotstrike import PlotStrike


# directory containing edis
edipath = r'c:\Users\jpeacock\Documents\GitHub\mtpy\data\edifiles2'


# # full path to file to save to
# savepath = r'C:\mtpywin\mtpy\examples\plots\edi_plots'


# gets edi file names as a list
elst = [op.join(edipath,f) for f in os.listdir(edipath) if (f.endswith('.edi'))]# and f.startswith('GL')


### this will plot the estimated strike duplicated across quadrants and 
# plot the orthogonal component.  Plot type 2 will plot all estimates of 
# strike into one ploe 
strike_plot = PlotStrike(fn_list=elst,
                         fold=False,
                         show_ptphimin=False,
                         plot_type=2)

# If you want to remove the orthogonal component
strike_plot.plot_orthogonal = False
strike_plot.fig_num = 2
strike_plot.plot()

# if you want to plot Tipper strike estimates as well
strike_plot.plot_tipper = True
strike_plot.fig_num = 3
strike_plot.plot()

# if you want to rotate the data
# note this will rotate the data N30W
strike_plot.rotation_angle = -30
strike_plot.fig_num = 4
strike_plot.plot()

# to rotate back
# strike_plot.rotation_angle = 30

# # if you want to plot per decade plots
# strike_plot.plot_type = 1
# strike_plot.fig_num = 5
# strike_plot.text_pad = 1.4
# strike_plot.plot()

# # if you want to only look at a few period ranges
# # not the range is given in log10 of the period
# strike_plot.plot_range = [-2, 0] 
# strike_plot.fig_num = 6
# strike_plot.plot()

# # if you want a vertical orientation instead of horizontal
# strike_plot.plot_orientation = 'v'
# strike_plot.plot_range = 'data'
# strike_plot.fig_num = 7
# strike_plot.plot()




