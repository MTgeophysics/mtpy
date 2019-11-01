# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:03:22 2019

@author: u64125
"""

from mtpy.modeling.modem import Model,Data
from mtpy.utils import gis_tools,convert_modem_data_to_geogrid
from mtpy.utils.calculator import nearest_index
from pyproj import Proj
import numpy as np
import os
from scipy.interpolate import RegularGridInterpolator


# wd = r'M:\AusLAMP\AusLAMP_NSW\Release\Model_release\MT075_DepthSlice_ArcGIS_ascii_grids'
# wdmod = r'C:\Users\u64125\OneDrive - Geoscience Australia\AusLAMP_NSW\Modelling\ModEM\NSWinv141'

wd = r'C:\Data\Alison_201910\MyOutput'
wdmod = r'C:\Data\Alison_201910\Alison_ModEM_Grid\MT075_ModEM_files'
filestem = 'Modular_MPI_NLCG_004'

mObj = Model()
mObj.read_model_file(os.path.join(wdmod,filestem+'.rho'))

dObj = Data()
dObj.read_data_file(os.path.join(wdmod,'ModEM_Data.dat'))

gce,gcn,gcz = [np.mean([arr[:-1],arr[1:]],axis=0) for arr in [mObj.grid_east,mObj.grid_north,mObj.grid_z]]
gce,gcn = gce[6:-6],gcn[6:-6]  # padding big-sized edge cells
# ge,gn = mObj.grid_east[6:-6],mObj.grid_north[6:-6]

print(gce)
print(gcn)
print(gcz)

print("Shapes E, N Z =", gce.shape, gcn.shape, gcz.shape)

fileext = '.asc'
ascfilelist = [ff for ff in os.listdir(wd) if ff.endswith(fileext)]

resgrid_nopad = mObj.res_model[::-1][6:-6,6:-6]

cs=7500  # grid size takes as a middle/medium value

newgridx,newgridy = np.meshgrid(np.arange(gce[0],gce[-1]+cs,cs),
                                np.arange(gcn[0],gcn[-1]+cs,cs))

header = """ncols        146
nrows        147
xllcorner    -164848.1035642
yllcorner    5611364.73539792
cellsize     7500"""
# May need to shift by half cellsize -cs/2
# [1]: -164848.1035642 -3750
# Out[1]: -168598.1035642
#
# In [2]: 5611364.73539792 - 3750
# Out[2]: 5607614.73539792


for ascfn in ascfilelist:
    depth = float(ascfn[10:-5])  #
    di = nearest_index(depth,gcz)
    # define interpolation function (interpolate in log10 measure-space)
    # See https://docs.scipy.org/doc/scipy-0.16.0/reference/interpolate.html
    interpfunc = RegularGridInterpolator((gce,gcn),np.log10(resgrid_nopad[:,:,di].T))
    # evaluate on the regular grid points, which to be output into geogrid formatted files
    newgridres = 10**interpfunc(np.vstack([newgridx.flatten(),newgridy.flatten()]).T).reshape(newgridx.shape)

    print("resistivity grid shape: ", newgridres.shape)
#    for i in range(len(gce)):
#        for j in range(len(gcn)):
#            newgridres[np.where(np.all([newgridx>ge[i],
#                                        newgridx <ge[i+1],
#                                        newgridy>gn[j],
#                                        newgridy <gn[j+1]],axis=0))] = resgrid_nopad[j,i,di]
#            print("assigned resistivity",i,j)
#    with open(os.path.join(wd,ascfn)) as ascfile:
#        header = ''
#        for line in ascfile:
#            if str.isalnum(line[:1]):
#                header += line
#            else:
#                break
    np.savetxt(os.path.join(wd,ascfn),newgridres,header=header,comments='',fmt='%.3e')
    print("Saved depth",depth)
