# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: u64125
"""

import os.path as op

import matplotlib.pyplot as plt
import numpy as np

import mtpy.modeling.modem as mtmn

# INPUTS #
# define a workdir for your environ
workdir = r'V:\Geology\conductivity_modelling'
workdir = r'E:\Githubz\mtpy2\examples\data\ModEM_files'
# workdir = r'/Softlab/Githubz/mtpy2/examples/data/ModEM_files'
# workdir = r'/g/data/ha3/fxz547/Githubz/mtpy2/examples/data/ModEM_files'

modeldir = op.join(workdir, 'VicSynthetic07')

# plot orientation ('ns' (north-south),'ew' (east-west) or 'z' (horizontal
# slice))
plotdir = 'z'
# slice location, in local grid coordinates (if it is a z slice, this is
# slice depth)
slice_location = 10000
# maximum distance in metres from vertical slice location and station
stationdist = 50000
# z limits (positive down so order is reversed)
zlim = (1e5, -5e3)
# colour limits
clim = [0.3, 3.7]

iterfn = 'Modular_MPI_NLCG_019.rho'
datafn = 'ModEM_Data_noise10inv.dat'

# END INPUTS #


read_data = True
if read_data:
    doo = mtmn.Data()
    doo.read_data_file(op.join(modeldir, datafn))
    moo = mtmn.Model(model_fn=op.join(modeldir, iterfn))
    moo.read_model_file()

# get grid centres
gcz = np.mean([moo.grid_z[:-1], moo.grid_z[1:]], axis=0)
gceast, gcnorth = [np.mean([arr[:-1], arr[1:]], axis=0)
                   for arr in [moo.grid_east, moo.grid_north]]

# distance from slice to grid centre locations
if plotdir == 'ew':
    sdist = np.abs(gcnorth - slice_location)
elif plotdir == 'ns':
    sdist = np.abs(gceast - slice_location)
elif plotdir == 'z':
    sdist = np.abs(gcz - slice_location)

# find closest slice index to specified location
sno = np.where(sdist == np.amin(sdist))[0][0]

# get data for plotting
if plotdir == 'ew':
    X, Y, res = moo.grid_east, moo.grid_z, np.log10(moo.res_model[sno, :, :].T)
    ss = np.where(
        np.abs(
            doo.station_locations['rel_north'] -
            np.median(gcnorth)) < stationdist)[0]

    sX, sY = doo.station_locations['rel_east'][
        ss], doo.station_locations['elev'][ss]
    xlim = (moo.grid_east[moo.pad_east], moo.grid_east[-moo.pad_east - 1])
    ylim = zlim
    title = 'East-west slice at {}km north'.format(gcnorth[sno])
elif plotdir == 'ns':
    X, Y, res = moo.grid_north, moo.grid_z, np.log10(
        moo.res_model[:, sno, :].T)
    # indices for selecting stations close to profile
    ss = np.where(
        np.abs(
            doo.station_locations['rel_east'] -
            np.median(gceast)) < stationdist)[0]

    sX, sY = doo.station_locations['rel_north'][
        ss], doo.station_locations['elev'][ss]
    xlim = (moo.grid_north[moo.pad_north], moo.grid_north[-moo.pad_north - 1])
    ylim = zlim
    title = 'North-south slice at {}km east'.format(gceast[sno])
elif plotdir == 'z':
    X, Y, res = moo.grid_east, moo.grid_north, np.log10(
        moo.res_model[:, :, sno])
    sX, sY = doo.station_locations[
        'rel_east'], doo.station_locations['rel_north']
    xlim = (moo.grid_east[moo.pad_east], moo.grid_east[-moo.pad_east - 1])
    ylim = (moo.grid_north[moo.pad_north], moo.grid_north[-moo.pad_north - 1])
    title = 'Depth slice at {}km'.format(gcz[sno])

# make the plot
plt.figure()
plt.pcolormesh(X, Y, res, cmap='bwr_r')
plt.xlim(*xlim)
plt.ylim(*ylim)

# plot station locations
plt.plot(sX, sY, 'kv')

# set title
plt.title(title)

if plotdir == 'z':
    plt.gca().set_aspect('equal')
plt.clim(*clim)
plt.colorbar()

plt.show()
