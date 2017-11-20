# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: u64125
"""

import os
import os.path as op

import matplotlib.pyplot as plt
import numpy as np

import legacy.modeling.modem_data as md
import mtpy.modeling.modem_model as mm

workdir = r'C:\Git\mtpy\examples\data'
modeldir = op.join(workdir, 'ModEM_files')

read_data = True
iterfn = max([ff for ff in os.listdir(modeldir) if ff.endswith('.rho')])

if read_data:
    doo = md.Data()
    doo.read_data_file(op.join(modeldir, 'ModEM_Data.dat'))
    moo = mm.Model(model_fn=op.join(modeldir, iterfn))
    moo.read_model_file()

snoew = 10
snons = 10
snoz = np.where(moo.grid_z > 80000)[0][0]
gcz = np.mean([moo.grid_z[:-1], moo.grid_z[1:]], axis=0)
plotdir = 'ew'

if plotdir == 'ew':
    X, Y, res = moo.grid_east, moo.grid_z, np.log10(
        moo.res_model[snoew, :, :].T)
    xlim = (-25000, 25000)
    ylim = (1e4, 0)
    sliceinfo = ''
elif plotdir == 'ns':
    X, Y, resf, resi = moo.grid_north, moo.grid_z, np.log10(moo.res_model[
                                                            :, snons, :].T)
    xlim = (-162500, 162500)
    ylim = (1e4, 0)
    sliceinfo = ''
elif plotdir == 'z':
    X, Y, resf, resi = moo.grid_east, moo.grid_north, np.log10(moo.res_model[
                                                               :, :, snoz])
    xlim = (-25000, 25000)
    ylim = (-162500, 162500)
    sliceinfo = ' depth {}km'.format(gcz[snoz])

titles = ['Forward model', 'Recovered model']


plt.figure()
plt.pcolormesh(X, Y, res, cmap='bwr_r')
plt.xlim(*xlim)
plt.ylim(*ylim)
# plt.title(titles[i]+sliceinfo)

if plotdir == 'z':
    plt.gca().set_aspect('equal')
    for ss in range(len(doo.station_locations)):
        plt.plot(
            doo.station_locations['rel_east'][ss],
            doo.station_locations['rel_north'][ss],
            'k.')
        plt.text(doo.station_locations['rel_east'][ss], doo.station_locations['rel_north'][ss],
                 doo.station_locations['station'][ss], fontsize=8)
# plt.clim(0.3,3.7)
plt.colorbar()
