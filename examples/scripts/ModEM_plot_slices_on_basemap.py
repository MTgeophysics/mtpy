"""
Plot ModEM depth slices on a basemap
"""


import os.path as op
import numpy as np

from mtpy.utils import gis_tools
from mtpy.modeling.modem import Data,Model
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import colors


pad=0.2

wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'

savepath = r'C:/tmp'


# get epsg and centre position of model
epsg = 28353 # epsg code for projection the model was projected to when creating the grid

# define model and data files
model_fn = op.join(wd,'Modular_MPI_NLCG_004.rho')
data_fn = op.join(wd,'ModEM_Data.dat')


mObj = Model()
mObj.read_model_file(model_fn = model_fn)

dObj = Data()
dObj.read_data_file(data_fn = data_fn)


# get easting and northing of model grid
east = mObj.grid_east + dObj.center_point['east']
north = mObj.grid_north + dObj.center_point['north']

# grid centres
gcx,gcy = [[np.mean(arr[i:i+2]) for i in range(len(arr)-1)] for arr in [east,north]]

# make a meshgrid, save the shape
east_grid,north_grid = np.meshgrid(east,north)
shape = east_grid.shape

# project grid to lat, lon
lonr,latr = gis_tools.epsg_project(east_grid,north_grid,epsg,4326)

# get boundaries for plot
minLon, maxLon, minLat, maxLat = lonr.min()-pad, lonr.max()+pad, latr.min()-pad, latr.max()+pad

# get corresponding resistivity values & station locations
resvals = mObj.res_model.copy()
sloc = dObj.station_locations




plt.figure(figsize=(20,15))


subplot = 1
for i in range(42,90,4): # list of depth slices to plot. This example plots 
                         # every 4th depth starting from 46 and finishing at 90
    plt.subplot(3,4,subplot) # change to nrows, ncols desired

    # make basemap
    bm = Basemap(resolution='l',projection='tmerc',llcrnrlon=minLon-pad,lat_0=(minLat+maxLat)/2.,lon_0=(minLon+maxLon)/2., 
                 llcrnrlat=minLat-pad, urcrnrlon=maxLon+pad, urcrnrlat=maxLat+pad)
    x,y = bm(lonr,latr)
    bm.drawmeridians(np.arange(np.floor(bm.lonmin),np.ceil(bm.lonmax),2),labels=[False,False,False,True])
    bm.drawparallels(np.arange(np.floor(bm.latmin),np.ceil(bm.latmax),2),labels=[True,False,False,False])
    bm.drawcoastlines()
    bm.drawstates()

    # set some parameters
    mpldict={}
    mpldict['cmap'] = 'jet_r'
    mpldict['norm'] = colors.LogNorm()
    mpldict['vmin'] = 1
    mpldict['vmax'] = 1e4

    plt.pcolormesh(x,y,resvals[:,:,i],**mpldict)
    plt.title('Depth slice at %.1f km depth'%(mObj.grid_z[i]/1e3))

    xp,yp=bm(sloc.lon,sloc.lat)
    plt.plot(xp,yp,'k+')
    subplot += 1

    plt.show()
#    plt.savefig(op.join(savepath,'DepthSlice_%.1fkmDepth.png'%(mObj.grid_z[i]/1e3)),dpi=400)