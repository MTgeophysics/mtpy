# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 10:21:05 2014

@author: Alison Kirkby

"""

import os
import mtpy.utils.filehandling as fh
import pek1dclasses as pek1dc


def generate_inputfiles(**input_parameters):
    
    data_kwds = ['wd', 'datafile', 'errorfloor_z', 
                 'errorfloor_type', 'epath', 'mode']
    control_kwds = ['wd', 'type_struct', 'type_aniso',
                    'value_struct', 'value_aniso', 'imax']
    inmodel_kwds = ['wd', 'modeldir', 'inmodel_vals']

    data_inputs = {}
    control_inputs = {}
    inmodel_inputs = {}
    
    build_inmodel = False
    for key in input_parameters.keys():
        if key in data_kwds:
            data_inputs[key] = input_parameters[key]
        if key in control_kwds:
            control_inputs[key] = input_parameters[key]
        if key in inmodel_kwds:
            inmodel_inputs = input_parameters[key]
        if key == 'build_inmodel':
            build_inmodel = input_parameters[key]
            
        
    Data = pek1dc.Data(**data_inputs)
    savepath = fh.make_unique_folder(Data.wd,os.path.basename(Data.epath)[:5]+Data.mode)
    Data.write_datafile(wd = savepath)
    
    Ctl = pek1dc.Control(**input_parameters)
    Ctl.write_ctlfile()
    
    if build_inmodel:
        Inmodel = pek1dc.Inmodel(inmodel_inputs)
        Inmodel.build_inmodel()
    



    

    