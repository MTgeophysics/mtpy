# -*- coding: utf-8 -*-
"""
Created on Mon Feb 02 11:38:58 2015

@author: Alison Kirkby

create modem input files
This script includes topography in the model. To not include topography, 
set number of air layers to zero (recommended) or comment out add_topography 
line. Note: setting number of air layers to zero will add bathymetry but not 
topography.

"""

import mtpy.modeling.modem_new as mtmn
import mtpy.core.edi as mtedi
import os.path as op
import os

# path to save to
workdir = r'V:\Software\mtpy\template_scripts\ModEM_example_build'
workdir = r'E:/ModEM_example_build'

# path where topography is located, if using
#sdir= r'V:\Software\mtpy\template_scripts\example_topo'
sdir= r'E:\Githubz\mtpy2\examples\etop1.xyz'

if not op.exists(workdir):
    os.mkdir(workdir)

# epsg to project to. Google epsg 'your projection'
epsg = 28354
edipath = r'E:/Datasets/MT_Datasets/GA_UA_edited_10s-10000s/'

# list of edi files, this line searches for all files in edipath ending with '.edi'
edi_list = [op.join(edipath,ff) for ff in os.listdir(edipath) if ff.endswith('.edi')]

# period list (can take periods from one of the edi files, or just specify 
# periods directly using the logspace function (commented out))
eo = mtedi.Edi(edi_list[0])
period_list = 1./eo.Z.freq
#period_list = np.logspace(-3,3)


do = mtmn.Data(edi_list=edi_list,
               inv_mode = '2',
               period_list=period_list,
               epsg=epsg
               )

do.write_data_file(save_path=workdir)


# create model file
mo = mtmn.Model(Data=do,
                cell_size_east=2000,
                cell_size_north=2000,
                pad_north=7, # number of padding cells in each of the north and south directions
                pad_east=7,# number of east and west padding cells
                pad_z=6, # number of vertical padding cells
                pad_stretch_v=3, # factor to increase by in padding cells (vertical)
                pad_stretch_h=3, # factor to increase by in padding cells (vertical)
                n_airlayers = 10, #number of air layers
                res_model=100, # halfspace resistivity value for reference model
                n_layers=80, # total number of z layers, including air
                z1_layer=100, # first layer thickness
                epsg=epsg, # epsg
                z_target_depth=120000)

mo.make_mesh()
# write a model file to initialise a resistivity model
mo.write_model_file(save_path=workdir)

# add topography to res model
mo.add_topography(op.join(sdir,'etopo1.asc'),interp_method='nearest')



# make covariance file
cov = mtmn.Covariance(mask_arr=mo.covariance_mask,
                      save_path=workdir,
                      smoothing_east=0.3,
                      smoothing_north=0.3,
                      smoothing_z=0.3
                      )
cov.write_covariance_file(model_fn=mo.model_fn)
