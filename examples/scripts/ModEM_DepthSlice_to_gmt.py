# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

Save depth slice data as x, y, res, for input to GMT.



"""

from mtpy.modeling.modem import Model, Data
from tests import MODEM_DIR, TEST_TEMP_DIR

######## inputs ##############
savepath = TEST_TEMP_DIR.joinpath("ModEM")
if not savepath.exists():
    savepath.mkdir()


model_fn = MODEM_DIR.joinpath("Modular_MPI_NLCG_004.rho")
data_fn = MODEM_DIR.joinpath("ModEM_Data.dat")
location_type = "LL"  # 'EN' to save eastings/northings, 'LL' to save longitude/latitude, need to provide model_epsg if using 'LL'
model_epsg = 28355  # epsg number model was projected in. common epsg numbers:
# 28351 (GDA94, mga zone 51), 28352 (mga zone 52), 28353 (mga zone 53),
# 28354 (mga zone 54), 28355 (mga zone 55), 28356 (mga zone 56)
# 3112 (Geoscience Australia Lambert)
# go to http://spatialreference.org/ref/epsg/?search=&srtext=Search for more info
model_utm_zone = "55S"  # alternative to epsg, can provide utm zone
depth_indices = [44, 45, 46]  # indices that define depth

##############################
# get the real-world origin from the data file
dataObj = Data(model_epsg=model_epsg)
dataObj.read_data_file(data_fn)
origin = [dataObj.center_point["east"][0], dataObj.center_point["north"][0]]

modObj = Model()
modObj.read_model_file(model_fn=model_fn)
modObj.write_xyres(
    origin=origin,
    savepath=savepath,
    location_type=location_type,
    model_epsg=model_epsg,
    model_utm_zone=None,
    log_res=True,
    outfile_basename="GMTtest",
    depth_index=depth_indices,
)
