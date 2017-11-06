
"""
Created on Fri Sep 20 14:58:51 2013

@author: Alison Kirkby

plots the resistivity model and/or the responses from 2D occam inversions

"""
import os
import os.path as op

import mtpy.modeling.occam2d as o2d
#import mtpy.modeling.occamtools as ot
from tests.beta import *


# import matplotlib.pyplot as plt
# plt.ion() # make figure disappear automatically
#plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc


def test_fun():
    """
    test function
    :return: T/F
    """

    # path to directory containing inversion files
    idir = os.path.join(SAMPLE_DIR,'Occam2d')

    # save path, to save plots to
    savepath = TEMP_OUT_DIR
    offset = 0


    #go to model results directory and find the latest iteration file
    iterfile = max([f for f in os.listdir(idir) if f[-5:]=='.iter'])
    respfile = max([f for f in os.listdir(idir) if f[-5:]=='.resp'])

    datafn = 'OccamDataFile.dat'
    # get the iteration number
    iterno = iterfile[-7:-5]
    outfilename = iterfile[:-5]

    plotmodel = True # set to True to plot the resistivity model
    plotresponses = True # set to True to plot the responses
    save = True

    #plot the model
    figure = plt.Figure()

    # horizontal padding on the edges for plotting, in km
    xpad = 1

    # plot the model
    if plotmodel:
        plotm = o2d.PlotModel(iter_fn=op.join(idir,iterfile),
                          data_fn=op.join(idir,datafn),
                          station_font_pad = 0.5,
                          station_font_size = 6,
                          station_font_rotation = 75,
                          climits=(0.,2.5),# colour scale limits
                          xpad=xpad,
                          dpi=300,# resolution of figure
                          fig_aspect = 0.5,   # aspect ratio between horizontal and vertical scale
                          ylimits = (0,10), # depth limits
                          stationid=(-1,3)) # index of station name to plot
        if save:
            plotm.save_figure(op.join(savepath, outfilename+'_resmodel.png') ) # this will produce 1 figure .png


    # plot the responses
    if plotresponses:
        plotresponse = o2d.PlotResponse(op.join(idir,datafn),
                                        resp_fn = op.join(idir,respfile),
                                        plot_type = ['pb35','pb40']
                                        )
        if save:
            plotresponse.save_figures(savepath) # this will produce 2 .pdf file


############################################
if __name__ == "__main__":
    test_fun()


"""
/usr/bin/python2.7 /Softlab/Githubz/mtpy/tests/beta/test_occam2d_plotmodel_and_responses.py
Note: need scipy version 0.14.0 or higher or interpolation might not work.
Reading from /Softlab/Githubz/mtpy/examples/model_files/Occam2d/OccamDataFile.dat
    profile_angle = 97.0
    geoelectric_strike = 7.0
    number of sites = 15
    number of frequencies = 19
Saved figure to: /Softlab/Githubz/mtpy/temp/ITER05_resmodel.png
Reading from /Softlab/Githubz/mtpy/examples/model_files/Occam2d/OccamDataFile.dat
    profile_angle = 97.0
    geoelectric_strike = 7.0
    number of sites = 15
    number of frequencies = 19
Traceback (most recent call last):
  File "/Softlab/Githubz/mtpy/tests/beta/test_occam2d_plotmodel_and_responses.py", line 85, in <module>
    test_fun()
  File "/Softlab/Githubz/mtpy/tests/beta/test_occam2d_plotmodel_and_responses.py", line 77, in test_fun
    plot_type = ['pb35','pb40']
  File "/Softlab/Githubz/mtpy/mtpy/modeling/occam2d.py", line 3389, in __init__
    self.plot()
    """
