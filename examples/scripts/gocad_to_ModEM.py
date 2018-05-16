# -*- coding: utf-8 -*-
"""
Created on Mon Feb 02 11:38:58 2015

@author: Alison Kirkby

create modem input files

"""

import os.path as op
import os

os.chdir(r'C:/mtpywin/mtpy')

import numpy as np

from mtpy.modeling.modem import Model,Data,Covariance

# gocad file
gocad_sgrid_file = r'C:\mtpywin\mtpy\examples\model_files\gocad\sgrid_from_gocad'


## Define some other parameters ##
epsg = 28355 # epsg code your model was projected to. See http://spatialreference.org/
             # to find the epsg code for your projection

# workdir (path to save to)
workdir = 'C:\tmp'


## Create or define data file ##
# two options here. If we are modifying an old model, can provide the old data
# file to get the locations. Alternatively, provide an array of x,y locations 
# (longitude and latitude) in the case of a new forward model
#### option 1 ####
data_fn = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2\ModEM_Data.dat'
dObj = Data()
dObj.read_data_file(data_fn=data_fn)

### option 2 ####
stationxyfile = r'C:\mtpywin\mtpy\examples\model_files\gocad\stationxy.txt'
xy = np.loadtxt(stationxyfile, usecols=(1,2))
station_names = np.loadtxt(stationxyfile,usecols=(0,),dtype='S10')
period_list = np.logspace(1,4,19) # change to your period list

# create data file
dObj = Data(epsg=epsg)
dObj._initialise_empty_data_array(xy,period_list,location_type='LL',stationnames=station_names)
dObj.write_data_file(fill=False,save_path=workdir)

# Create model file
mObj = Model(Data = dObj,save_path=workdir)
mObj.read_gocad_sgrid_file(gocad_sgrid_file)
mObj.model_fn_basename = op.join(workdir,'ModEM_Model_forward.ws')
mObj.write_model_file()

# create covariance file
cov = Covariance(save_path=workdir,
                      smoothing_east=0.3,
                      smoothing_north=0.3,
                      smoothing_z=0.3)
cov.write_covariance_file(model_fn=mObj.model_fn)

