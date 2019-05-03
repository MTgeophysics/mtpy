# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Plot root-mean-square misfit (RMS across all periods) at each site


"""
import os.path as op
import os
os.chdir(r'C:/mtpywin/mtpy')
#from mtpy.imaging.plot_response import PlotResponse
from mtpy.modeling.modem import Residual
import matplotlib.pyplot as plt

#### Inputs ####
wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM'
savepath = r'C:\tmp'
filestem = 'Modular_MPI_NLCG_004'

datafn = op.join(wd,'ModEM_Data.dat')
respfn = op.join(wd,filestem+'.dat')


# read residual file into a residual object
residObj = Residual(residual_fn=op.join(wd,filestem+'.res'))
residObj.read_residual_file()
residObj.get_rms()


# get some parameters as attributes
lat,lon,east, north,rel_east, rel_north, rms,station = [residObj.rms_array[key] for key in ['lat','lon','east','north','rel_east','rel_north','rms','station']]

# create the figure
plt.figure()
plt.scatter(east,north,c=rms,cmap='bwr')
for i in range(len(station)):
    plt.text(east[i],north[i],station[i],fontsize=8)

plt.colorbar()
plt.clim(1,4)