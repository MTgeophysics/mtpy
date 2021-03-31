"""
Created on Fri Sep 20 14:58:51 2013

@author: Alison Kirkby

plots the resistivity model and/or the responses from 2D occam inversions

"""
import os

os.chdir(r"C:/Git/mtpy")
import matplotlib.pyplot as plt
import mtpy.modeling.occam2d as o2d

# import mtpy.modeling.occamtools as ot
import os.path as op

# path to directory containing inversion files
idir = r"C:\Git\mtpy\examples\model_files\Occam2d"

# save path, to save plots to
savepath = r"C:\Git\mtpy\examples\plots\Occam2d"
offset = 0


# go to model results directory and find the latest iteration file
iterfile = max([f for f in os.listdir(idir) if f[-5:] == ".iter"])
respfile = max([f for f in os.listdir(idir) if f[-5:] == ".resp"])

datafn = "OccamDataFile.dat"
# get the iteration number
iterno = iterfile[-7:-5]
outfilename = iterfile[:-5]

plotmodel = True  # set to True to plot the resistivity model
plotresponses = True  # set to True to plot the responses
save = True

# plot the model
figure = plt.Figure()

# horizontal padding on the edges for plotting, in km
xpad = 1

# plot the model
if plotmodel:
    o2d.PlotModel(
        iter_fn=op.join(idir, iterfile),
        data_fn=op.join(idir, datafn),
        station_font_pad=0.5,
        station_font_size=6,
        station_font_rotation=75,
        climits=(0.0, 2.5),  # colour scale limits
        xpad=xpad,
        dpi=300,  # resolution of figure
        fig_aspect=0.5,  # aspect ratio between horizontal and vertical scale
        ylimits=(0, 10),  # depth limits
        stationid=(-1, 3),
    )  # index of station name to plot
    if save:
        plt.savefig(op.join(savepath, outfilename + "_resmodel.png"), dpi=300)

# plot the responses
if plotresponses:
    plotresponse = o2d.PlotResponse(
        op.join(idir, datafn),
        resp_fn=op.join(idir, respfile),
        plot_type=["pb35", "pb40"],
    )
    if save:
        plotresponse.save_figures(savepath)
