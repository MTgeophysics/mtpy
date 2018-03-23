# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Save depth slice data as x, y, res, for input to GMT.



"""

import os
import os.path as op


os.chdir(r'C:\mtpywin\mtpy') # change to path to your mtpy installation to ensure correct version is used
mtpy2 = False
from mtpy.modeling.modem import Model
from mtpy.modeling.modem import Data

######## inputs ##############
wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'
savepath = r'C:\test'


model_fn = op.join(wd,'Modular_MPI_NLCG_004.rho')
data_fn = op.join(wd,'ModEM_Data.dat')
location_type ='LL'# 'EN' to save eastings/northings, 'LL' to save longitude/latitude, need to provide model_epsg if using 'LL'
model_epsg = 28355 # epsg number model was projected in. common epsg numbers:
                  # 28351 (GDA94, mga zone 51), 28352 (mga zone 52), 28353 (mga zone 53),
                  # 28354 (mga zone 54), 28355 (mga zone 55), 28356 (mga zone 56)
                  # 3112 (Geoscience Australia Lambert)
                  # go to http://spatialreference.org/ref/epsg/?search=&srtext=Search for more info
model_utm_zone = '55S'# alternative to epsg, can provide utm zone
depth_indices = [44,45,46] # indices that define depth

##############################
# get the real-world origin from the data file
dataObj = Data(data_fn = data_fn, 
               model_epsg=model_epsg)
dataObj.read_data_file()
origin = [dataObj.center_point['east'][0],dataObj.center_point['north'][0]]

modObj = Model()
modObj.read_model_file(model_fn = op.join(wd,model_fn))
modObj.write_xyres(origin=origin,
                   savepath=savepath,
                   location_type=location_type,
                   model_epsg=model_epsg,
                   model_utm_zone=None,
                   log_res=True,
                   outfile_basename='GMTtest',
                   depth_index = depth_indices
                   )