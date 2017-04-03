# -*- coding: utf-8 -*-
"""
Created on Fri May 20 15:24:26 2016

@author: a1193899
"""
import mtpy.modeling.modem_new as mtmn
import os
import os.path as op
import numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt

import mayavi.mlab as mlab

workdir = r'M:\RegionalSurvey\MT056_Cloncurry\3Dmodels\AMTband\1000hz_topo'
workdir = r'E:\Data\Modeling\Isa\100hs_flat_BB'

# read the model and data files
moo = mtmn.Model()
moo.read_model_file(op.join(workdir,'Isa_run3_NLCG_049.rho'))

doo=mtmn.Data()
doo.read_data_file(op.join(workdir,'IIsa_run3_NLCG_049.dat'))

# eastings and northings of end points of slice
point1xy,point2xy = (0,0),(1000,1000)

x0,y0 = np.median(doo.station_locations['east']-doo.station_locations['rel_east']),\
        np.median(doo.station_locations['north']-doo.station_locations['rel_north'])

slicex = np.arange(point1xy[0],point2xy[0],1000) + x0
slicey = np.arange(point1xy[1],point2xy[1],1000) + y0


# turn x y station coordinates into 2d arrays with depth
xi,zi = np.meshgrid(slicex,moo.grid_z)
yi,zi = np.meshgrid(slicey,moo.grid_z)

# get x, y and z coordinates of each resistivity point in the model
mox,moy,moz = np.meshgrid(moo.grid_east,moo.grid_north,moo.grid_z)

# x y z points of res model in correct format for griddata
moxyz = np.vstack([mox.flatten(),moy.flatten(),moz.flatten()]).T

# x y z interpolation points in correct format
xyzi = np.vstack([xi.flatten(),yi.flatten(),zi.flatten()]).T

# log resistivity in correct format
logres = np.log10(moo.res_model.flatten())

# resistivity along a snake connecting the points
snakeres = si.griddata(moxyz,logres,xyzi,method='nearest').reshape(*list(xi.shape))


# plot the model
#plt.pcolormesh(stxi,stz,snakeres,cmap='jet_r')
#plt.ylim(100000,0)

# plot the model (3d)
k = np.where(stz>100000)[0][0]
image = mlab.mesh(xi[:k],yi[:k],-zi[:k],scalars=-snakeres[:k])
mlab.points3d(stx,sty,np.zeros_like(stx),color=(0,0,0),mode='sphere',scale_factor=2500,vmin=-1000,vmax=0)
#mlab.savefig(op.join(workdir,'model_iter052.png'))
mlab.savefig(op.join(workdir,'model_iter046.pdf'))
