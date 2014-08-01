# -*- coding: utf-8 -*-
"""
Created on Fri Aug 01 15:52:47 2014

@author: Alison Kirkby
"""
import mtpy.modeling.occam2d as o2d
import mtpy.modeling.pek1dclasses as p1dc
import numpy as np
import os
import pek2dforward as p2d

class Model():
    """
    class for creating and reading model files
    
    """
   
    
    def __init__(self,**input_parameters):
        self.working_directory = '.'
        self.edi_directory = None
        self.occam_configfile = None
        self.parameters_model = {}
        self.parameters_model['model_depth'] = 100
        self.parameters_model['no_layers'] = 25
        self.parameters_model['max_blockwidth'] = 1000
        self.mesh = None
        self.meshlocations_x = None
        self.meshlocations_z = None
        self.meshblockwidths_x = None
        self.meshblockthicknesses_z = None
        self.profile_easts = None
        self.profile_norths = None
        self.inversion1d_dirdict = {}
        self.inversion1d_masterdir = '.'
        self.inversion1d_modelno
        self.inversion1d_bins = np.logspace(0,3,4)

        self.edifiles = []

        self.Data = None

        self.modelfile = 'model'


        update_dict = {}

        #correcting dictionary for upper case keys
        input_parameters_nocase = {}
        for key in input_parameters.keys():
            input_parameters_nocase[key.lower()] = input_parameters[key]

        update_dict.update(input_parameters_nocase)

        for dictionary in [self.parameters_startup, self.parameters_inmodel, 
                                    self.parameters_mesh, self.parameters_data]:
            for key in dictionary.keys():
                if key in update_dict:
                    #check if entry exists:
                    try:
                        value = float(update_dict[key])
                        dictionary[key] = value
                    except:
                        value = update_dict[key]
                        dictionary[key] = value
                        if type(value) in [str]:
                            if value.strip().lower()=='none':
                                 dictionary[key] = None
                                 
        for key in update_dict:
            try:
                value = getattr(self,key)
                if len(update_dict[key]) > 0:
                    try:
                        value = float(update_dict[key])
                        setattr(self,key,value)
                    except:
                        value = update_dict[key]
                        setattr(self,key,value)
                        if type(value) in [str]:
                            if value.strip().lower()=='none':
                                setattr(self,key,None)
            except:
                continue

        self.input_parameters = update_dict
        
    
    def build_model(self):
        """
        build model file
        """
        pass
    
    
    def get_station_meshblock_numbers(self):
        """
        
        """

        try:
            ivals = []
            ii = 0
            for sl in self.stationlocations:
                for j in range(len(self.meshlocations_x[:-1])):
                    if (sl>self.meshlocations_x[j])&(sl<=self.meshlocations_x[j+1]):
                        ivals.append(ii)
                    ii += 1
        except AttributeError:
            print "no stationlocations, please build mesh first"
         
                     
    
    def build_mesh(self):
        """
        create a mesh using occam2d
        """
        
        # create an occam2d setup object
        so = o2d.Setup(wd=self.working_directory,
                       edi_directory=self.edi_directory,
                       configfile=self.occam_configfile,)
        so.read_edifiles(edi_dir=self.edi_directory)
        # create an occam2d data object
        so.Data = o2d.Data(edilist=so.edifiles,wd=so.wd,**so.parameters_data)
        # set up meshlocations
        so.setup_mesh_and_model()
        self.stationlocations = np.array(so.Data.stationlocations)/1000.
        
        # set occam mesh attributes to pek2d object
        for attribute in ['meshlocations_x','meshlocations_z',
                          'meshblockwidths_x','meshblockthicknesses_z',
                          'profile_easts','profile_norths','Data']:
                              if 'mesh' in attribute:
                                  attvalue = np.array(getattr(so,attribute))/1000.
                              else:
                                  attvalue = getattr(so,attribute)
                              setattr(self,attribute,attvalue)
        
    def get_1d_results(self):
        """
        get 1d inversion results to apply to inputs of 2d model
        
        """
        if not hasattr(self,'Data'):
            self.build_mesh()
            
        self.inversion1d_dirdict = \
        p2d.find_directory(self.Data.stations,
                           self.inversion_masterdir,
                           start_list=self.inversion1d_dirdict)
                           
        models1d = {}
        
        for key in self.inversion1d_dirdict:
            idir = self.inversion1d_dirdict[key]
            mod = p1dc.Model(idir)
            mod.read_model()
            models1d[key] = mod.models[self.inversion1d_modelno-1]

        self.models1d = models1d

def bin_results(in_array)


def find_directory(search_string,masterdir,start_dict = {}):

    for s in search_string:
        for d in os.listdir(masterdir):
            append = False
            if s in os.path.basename(d):
                append = True
            else:
                slst = s.strip().split('_')
                for ss in slst:
                    if ss in d:
                        append = True
            if append:
                start_dict[s] = d
        
        
                            
    
        
    
    

#class Response():
#    """
#    class to contain outputs of forward modelling
#    
#    """