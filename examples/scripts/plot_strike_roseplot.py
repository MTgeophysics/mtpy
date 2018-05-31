# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:51:02 2015

@author: Alison Kirkby

strike analysis using roseplot function

"""
import os.path as op
import os
os.chdir(r'C:\mtpywin\mtpy') # change to path where mtpy is installed


from mtpy.imaging.plotstrike import PlotStrike


# directory containing edis
edipath = r'C:\mtpywin\mtpy\examples\data\edi_files_2'


# full path to file to save to
savepath = r'C:\mtpywin\mtpy\examples\plots\edi_plots'


# gets edi file names as a list
elst = [op.join(edipath,f) for f in os.listdir(edipath) if (f.endswith('.edi'))]# and f.startswith('GL')


strikeplot = PlotStrike(fn_list=elst,
                        plot_type=1 # 1 means divide into separate plots for different decades
                                    # 2 means combine all data into one rose plot
                                    )
strikeplot.save_plot(savepath,
                     file_format='png',
                     fig_dpi=400 
                     )