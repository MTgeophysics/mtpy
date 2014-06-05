# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 10:21:05 2014

@author: Alison Kirkby

"""

import os
import mtpy.utils.filehandling as fh
import pek1dclasses as pek1dc


def generate_inputfiles(epath, **input_parameters):
    
    """
    generate input files for a model. 
    
    -----------------------Compulsory parameter--------------------------------
    **epath** the full path to the edi file.
    
    -----------------------Recommended parameter-------------------------------
    **wd** working directory, default is the edi directory. A new directory
           is created under this directory to put all the input files into

    ------------------------Optional Parameters--------------------------------
    **datafile** name for the input file, if not specified, name is taken from
                 the edi file
    **errorfloor_z** error floor for the input z values, can be an absolute
                     value or relative (e.g. 0.1 means 10%)
                     default is 0.1
    **errorfloor_type** type of error floor, either 'relative' or 'absolute'
                        default is relative.
    **type_struct** type of structure penalty, default is 6
    **type_aniso** type of anisotropy penalty, default is 2
    **value_struct** structural penalty weights to apply, default is [1,10,100]
    **value_aniso** anisotropy penalty weights to apply, default is [1,10,100]
    **imax** maximum number of iterations to run, default is 100


    to generate an a priori (inmodel) file, need to put keyword
    **build_inmodel** = True, default is False
    
    also need to specify the following parameters:
    **inmodel_vals**
    

    
    inmodel_modeldir = string, folder containing previous model run with same 
    resolution, necessary for constructing the layer depths in the inmodel file.
    inmodel_vals = dictionary structured as follows:
    {layer top depth:[minimum_resistivity, maximum_resistivity, strike]}
    
    """
    
    data_kwds = ['wd','datafile', 'errorfloor_z', 
                 'errorfloor_type', 'epath', 'mode']
    control_kwds = ['type_struct', 'type_aniso',
                    'value_struct', 'value_aniso', 'imax']
    inmodel_kwds = ['inmodel_vals']

    data_inputs = {'epath':epath}
    control_inputs = {}
    inmodel_inputs = {}
    
    build_inmodel = False
    for key in input_parameters.keys():
        print key
        if key in data_kwds:
            data_inputs[key] = input_parameters[key]
        if key in control_kwds:
            control_inputs[key] = input_parameters[key]
        if key in inmodel_kwds:
            inmodel_inputs = input_parameters[key]
        if key == 'build_inmodel':
            build_inmodel = input_parameters[key]


    Data = pek1dc.Data(**data_inputs)
    # make a save path to match the edi file
    savepath = fh.make_unique_folder(Data.wd,os.path.basename(Data.epath)[:5]+Data.mode)
    Data.write_datafile(wd = savepath)
    
    # update the working directory to the new savepath
    control_inputs['wd'] = savepath
    inmodel_inputs['wd'] = savepath
    
    Ctl = pek1dc.Control(**control_inputs)
    Ctl.write_ctlfile()
    
    if build_inmodel:
        if 'inmodel_modeldir' in input_parameters.keys():
            Inmodel = pek1dc.Inmodel(input_parameters['inmodel_modeldir'],
                                     **inmodel_inputs)
            Inmodel.write_inmodel()

def parse_arguments(arguments):
    """
    parse command line arguments
    
    """
    return


def create_filelist(wd,subfolder_list = None):
    """
    create a list of full paths to edi files    
    
    """
    return
    
    
def build_run(**input_parameters):
    """
        
    
    """
    return
    

    