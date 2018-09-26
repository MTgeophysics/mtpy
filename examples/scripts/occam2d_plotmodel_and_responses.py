"""

plots the resistivity model and/or the responses from 2D occam inversions

"""
import os
import matplotlib.pyplot as plt
import mtpy.modeling.occam2d as o2d
import os.path as op

# path to directory containing inversion files
idir = r'C:/mtpywin/mtpy/examples/model_files/Occam2d'

# save path, to save plots to
savepath = r'C:/tmp'

# specify iteration file to plot (will search for the latest iteration if not provided)
iterfile = None       # User can specify the iteration file to plot, e.g.'smooth25.iter'

#go to model results directory and find the latest iteration file
if iterfile is None:
	iterfile = max([f for f in os.listdir(idir) if f[-5:]=='.iter'])
	respfile = max([f for f in os.listdir(idir) if f[-5:]=='.resp'])
else:
	respfile = iterfile.replace('.iter','.resp')


datafn = 'OccamDataFile.dat'
# get the iteration number
iterno = iterfile[-7:-5]

plotmodel = True # set to True to plot the resistivity model
plotresponses = False # set to True to plot the responses
save = False

#plot the model
figure = plt.Figure()


if plotmodel:
    pltm = o2d.PlotModel(iter_fn=op.join(idir,iterfile),
                      data_fn=op.join(idir,datafn),
                      station_font_pad = 2.5,
                      station_font_size = 4,
                      station_font_rotation = 90,
		              font_size=4,
                      climits=(0,4),# colour scale limits
                      xpad=0.1,
                      yscale = 'm',
                      yminorticks = 0.250,
                      dpi=300,# resolution of figure
                      fig_aspect = 0.25,   # aspect ratio between horizontal and vertical scale
                      ylimits = (-2,20), # depth limits
                      stationid=(-1,3)) # index of station name to plot
    if save:
        pltm.save_figure(op.join(savepath,'model.png'))
        
if plotresponses:
    plts = o2d.PlotResponse(op.join(idir,datafn),
                                resp_fn = op.join(idir,respfile)
                                )
    plts.save_figures(op.join(idir,'Copna_response'))   

