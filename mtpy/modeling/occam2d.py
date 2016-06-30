# -*- coding: utf-8 -*-
"""
Spin-off from 'occamtools'
(Created August 2011, re-written August 2013)

Tools for Occam2D

authors: JP/LK


Classes:
    - Data
    - Model
    - Setup
    - Run
    - Plot
    - Mask


Functions:
    - getdatetime
    - makestartfiles
    - writemeshfile
    - writemodelfile
    - writestartupfile
    - read_datafile
    - get_model_setup
    - blocks_elements_setup


"""
#==============================================================================
import numpy as np
import scipy as sp
from scipy.stats import mode
import sys
import os
import os.path as op
import subprocess
import shutil
import fnmatch
import datetime
from operator import itemgetter
import time
import matplotlib.colorbar as mcb
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

import mtpy.core.edi as MTedi
import mtpy.modeling.winglinktools as MTwl
import mtpy.utils.conversions as MTcv
import mtpy.utils.filehandling as MTfh
import mtpy.utils.configfile as MTcf
import mtpy.analysis.geometry as MTgy
import mtpy.utils.exceptions as MTex
import scipy.interpolate as si

reload(MTcv)
reload(MTcf)
reload(MTedi)
reload(MTex)

#==============================================================================

occamdict = {'1':'resxy','2':'phasexy','3':'realtip','4':'imagtip','5':'resyx',
             '6':'phaseyx'}

#------------------------------------------------------------------------------

class Setup():
    """
    Dealing with the setup  for an Occam2D run. Generate 'startup', 'inmodel', 
    'mesh' files. Calling Data() for generating a suitable input data file.

    Setting up those files within one (pre-determined) folder, so Occam can be 
    run there straight away.

    """


    def __init__(self, configfile = None, **input_parameters):


        self.parameters_startup = {}
        self.parameters_inmodel = {}
        self.parameters_data = {}
        self.parameters_mesh = {}
        self.parameters_prejudice = {}

        self.parameters_startup['description'] = 'generic MTpy setup'

        self.parameters_startup['iter_format'] = 'OCCAM_ITER'
        self.parameters_startup['datetime_string'] = datetime.datetime.now().strftime(
                                                             '%Y/%m/%d %H:%M:%S')

        self.parameters_startup['no_iteration'] = 0
        self.parameters_startup['roughness_start'] = 1.0E+07
        self.parameters_startup['reached_misfit'] = 0        
        self.parameters_startup['roughness_type'] = 1
        self.parameters_startup['debug_level'] = 1
        self.parameters_startup['mu_start'] = 5.0
        self.parameters_startup['max_no_iterations'] = 30
        self.parameters_startup['target_rms'] = 1.5
        self.parameters_startup['rms_start'] = 1000
        self.parameters_startup['halfspace_resistivity'] = 100.


        self.parameters_inmodel['no_sideblockelements'] = 5
        self.parameters_inmodel['no_bottomlayerelements'] = 4
        self.parameters_inmodel['firstlayer_thickness'] = 100
        self.parameters_inmodel['no_mergedsideblockelements'] = 5
        self.parameters_inmodel['sidepadding_increasefactor'] = 2.
        #model depth is in km!
        self.parameters_inmodel['model_depth'] = 100
        self.parameters_inmodel['no_layers'] = 25
        self.parameters_inmodel['max_blockwidth'] = 1000

        self.parameters_inmodel['model_name'] = 'Modelfile generated with MTpy'
        self.parameters_inmodel['block_merge_threshold'] = 0.75
        # factor determining amount of increase with each bottom padding
        self.parameters_inmodel['bottompadding_scalefactor'] = 1.
        
        self.parameters_inmodel['roughness_exceptions'] = None
        self.parameters_inmodel['build_roughness_exceptions'] = False
        self.parameters_inmodel['interface_dir'] = None
        self.parameters_inmodel['interface_filelist'] = None
        self.parameters_inmodel['elevation_filename'] = None

       
        self.parameters_data['strike'] = None

        self.parameters_data['rho_errorfloor'] = 0.
        self.parameters_data['phase_errorfloor'] = 0.
        self.parameters_data['tipper_errorfloor'] = 0.
        self.parameters_data['tipper_errorfloor_abs'] = 0.
        self.parameters_data['azimuth'] = 0

        self.parameters_data['mode'] = 'tetm'
        self.parameters_data['edi_type'] = 'z'

        self.parameters_data['min_frequency'] = None
        self.parameters_data['max_frequency'] = None
        self.parameters_data['max_no_frequencies'] = None

        self.parameters_mesh['mesh_title'] = 'Mesh file generated with MTpy'
        
        self.parameters_prejudice = {}
        self.parameters_prejudice['constraints_type'] = 'interfaces'
        self.parameters_prejudice['build_prejudice_file'] = False
        self.parameters_prejudice['interfaces'] = [] # list of arrays containing 2 columns, distance x along profile and depth z of interface
        self.parameters_prejudice['interface_resistivity_values'] = None
        self.parameters_prejudice['interface_prejudice_weights'] = None# 
        self.parameters_prejudice['interface_filelist'] = None
        self.parameters_prejudice['welldata'] = None        
        self.parameters_prejudice['prejudice_weight'] = 1.
 
        self.mesh = None
        self.meshlocations_x = None
        self.meshlocations_z = None
        self.meshblockwidths_x = None
        self.meshblockthicknesses_z = None
        self.profile_easts = None
        self.profile_norths = None
        self.modelblocklocations_x = None
        self.modelblocklocations_z = None

        self.inmodel = None
 
        self.no_parameters = None

        self.welldata = None
        self.well_xy = None

        self.edifiles = []

        self.Data = None
        self.Prejudice = None
        self.datafile = 'occaminputdata.dat'
        self.meshfile = 'mesh'
        self.inmodelfile = 'inmodel'
        self.startupfile = 'startup'
        self.staticsfile = None
        self.prejudicefile = None

        self.edi_directory = None
        #working directory
        self.wd = '.'

        update_dict = {}
        if configfile is not None:
            try:
                configfile = op.join(os.curdir,configfile)
                if not op.isfile(configfile):
                    raise
            except:
                raise MTex.MTpyError_inputarguments('Error - Configuration '\
                    'file {0} does not exist!'.format(configfile))
            try:
                update_dict = {}
                raw_configfile_content = MTcf.read_configfile(configfile)
                for k in raw_configfile_content.keys():

                    temp_dict = raw_configfile_content[k]
                    update_dict.update(temp_dict)
            except:
                raise MTex.MTpyError_inputarguments('Error - Configuration '\
                    'file {0} cannot be read!'.format(configfile))

        #correcting dictionary for upper case keys
        input_parameters_nocase = {}
        for key in input_parameters.keys():
            input_parameters_nocase[key.lower()] = input_parameters[key]

        update_dict.update(input_parameters_nocase)

        for dictionary in [self.parameters_startup, self.parameters_inmodel, 
                           self.parameters_mesh, self.parameters_data,
                           self.parameters_prejudice]:
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
        
        if self.parameters_prejudice['build_prejudice_file']:
            self.prejudicefile = 'prejudice'


    def read_configfile(self, configfile):

        cf = op.abspath(configfile)
        if not op.isfile(cf):
            print 'Warning - config file not found {0}'.format(cf)
            return

        config_dictionary = MTcf.read_configfile(cf)

        no_p = self.update_parameters(config_dictionary)
        no_a = self.update_attributes(config_dictionary)
        print '{0} parameters and attributes updated'.format(no_a + no_p)


    def update_parameters(self, **parameters_dictionary):

        input_parameters_nocase = {}
        for key in parameters_dictionary.keys():
            input_parameters_nocase[key.lower()] = parameters_dictionary[key]

        if self.validate_parameters(input_parameters_nocase) is False:
            print 'Error - parameters invalid \n'
            return

        counter = 0
        for dictionary in [self.parameters_startup, self.parameters_inmodel, 
                                    self.parameters_mesh, self.parameters_data]:
            for key in dictionary.keys():
                if key in input_parameters_nocase:
                    dictionary[key] = input_parameters_nocase[key]
                    counter += 1

        return counter



    def validate_parameters(self, **parameters_dictionary):
        
        valid = True

        return valid


    def update_attributes(self, **attributes_dictionary):
        input_attributes_nocase = {}
        for key in attributes_dictionary.keys():
            input_attributes_nocase[key.lower()] = attributes_dictionary[key]

        if self.validate_attributes(input_attributes_nocase) is False:
            print 'Error - attributes invalid \n'
            return
        counter = 0
        for attr in dir(self):
            if attr in input_attributes_nocase:
                self.setattr(attr,input_attributes_nocase[attr])
                counter += 1

        return counter 


    def validate_attributes(self, **attributes_dictionary):
        
        valid = True

        return valid


    def add_edifiles_directory(self, directory = None):

        if directory is None:
            print 'Error - provide directory name'
            return

        if op.isdir(directory) is False:
            print 'Warning - not a valid directory - cannot browse for EDI files: {0}'.format(directory)
            return

        edilist_raw = fnmatch.filter(os.listdir(directory),'*.[Ee][Dd][Ii]')
        edilist_full = [op.abspath(op.join(edi_dir,i)) for i in edilist_raw]
        
        counter = 0
        for edi in edilist_full:
            try:
                if op.isfile(edi):
                    self.edifiles.append(edi)
                    counter += 1
            except:
                continue
        
        print 'Added {0} Edi files'.format(counter)


    def add_edifiles(self, edilist):
        
        if not np.iterable(edilist):
            print 'Error - provide valid file list'
            return

        counter = 0
        for edi in edilist:
            try:
                fn = op.abspath(op.join(os.curdir,edi))
                if op.isfile(fn):
                    self.edifiles.append(fn)
                    counter += 1                
            except:
                continue

        print 'Added {0} Edi files'.format(counter)


    def remove_edifiles(self, edilist):
        if not np.iterable(edilist):
            print 'Error - provide valid file list'
            return

        counter = 0
        for edi in edilist:
            try:
                fn = op.abspath(op.join(os.curdir,edi))
                if fn in self.edifiles:
                    self.edifiles.remove(fn)
                counter += 1                
            except:
                continue

        print 'Removed {0} Edi files'.format(counter)


    def read_edifiles(self, edi_dir = None):
        
        if self.edi_directory is None:
            self.edi_directory = op.abspath(os.curdir)

        if (edi_dir is not None):
            edi_dir = op.abspath(op.join(os.curdir,edi_dir))
            if (edi_dir):
                self.edi_directory = edi_dir
            else:
                print 'Warning - given directory not found: {0} \n\t-'\
                    ' using current directory instead:{1}'.format(edi_dir,os.curdir)


        edilist_raw = fnmatch.filter(os.listdir(self.edi_directory),'*.[Ee][Dd][Ii]')
        edilist_full = [op.abspath(op.join(self.edi_directory,i)) for i in edilist_raw]
        edilist = []
        for edi in edilist_full:
            if op.isfile(edi):
                edilist.append(edi)
            else: 
                continue

        self.edifiles = edilist
       

    def write_datafile(self):

        try:
            data_object = Data(edilist = self.edifiles, wd = self.wd, **self.parameters_data)
        except:
            print 'cannot write data file'
            raise
        self.stationlocations = data_object.stationlocations
        data_object.writefile(self.datafile)


        #self.datafile = data_object.datafile
        self.Data = data_object
        

    def setup_mesh_and_model(self):
        """
        Build the mesh and inmodel blocks from given data and parameters.

        Attributes required: 

        - self.stationlocations
        - self.parameters_inmodel

        """
        #given as offset on the profile line
        lo_sites = self.Data.stationlocations
        n_sites  = len(lo_sites)

        #maximum width of MODEL block - implicitely defines finiteness of the mesh
        maxblockwidth = float(self.parameters_inmodel['max_blockwidth'])

        #define vertical setup
        #old version:
        #number of layers per depth meters decade
        #layers_per_decade     = float(self.parameters_inmodel['no_layersperdecade'])
        #new:
        model_depth_m         = float(self.parameters_inmodel['model_depth'])*1000.
        #depth of first layer
        first_layer_thickness = float(self.parameters_inmodel['firstlayer_thickness'])
        #so "layers_per_decade" layers between "first_layer_thickness" and 10* "first_layer_thickness"
        #then the same number between 10* "first_layer_thickness" and 100* "first_layer_thickness"

        #padding layer does NOT count in here
        n_layers = int(float(self.parameters_inmodel['no_layers']))

        #number of padding mesh layers at the bottom 
        n_bottompadding = int(float(self.parameters_inmodel['no_bottomlayerelements']))
        bottompadding_scalefactor = self.parameters_inmodel['bottompadding_scalefactor']

        #number of padding mesh columns to the sides
        n_sidepadding    = int(float(self.parameters_inmodel['no_sideblockelements']))
        nmerge_sidepadding  = int(float(self.parameters_inmodel['no_mergedsideblockelements']))
        
        if (n_sidepadding-nmerge_sidepadding)%2 == 1:
            nmerge_sidepadding += 1
        sidepadding_increasefactor = float(self.parameters_inmodel['sidepadding_increasefactor'])

        #1. check, if inter-station spacing is smaller than the allowed max block size
        #if not, add dummy station locations

        lo_allsites = []
        lo_distances = []
        lo_real_station_distances = []
        no_dummys = 0


        print '\nlength of station profile: {0:.1f} km '.format((lo_sites[-1]-lo_sites[0])/1000.)
        print 'Azimuth of profile: {0:.1f} degrees'.format(self.Data.azimuth)
        if self.Data.strike is not None:
            print 'Assumed strike: {0:.1f} degrees'.format(self.Data.strike)
        else:
            print 'Strike orientation unknown'

        for idx_site,location in enumerate(lo_sites):
            lo_allsites.append(location)
            if idx_site == len(lo_sites)-1:
                break

            distance = np.abs(lo_sites[idx_site+1] - location)
            lo_real_station_distances.append(distance)
            if distance >= maxblockwidth:
                dummys = int(distance/maxblockwidth) 
                smallblockwidth = distance/float(dummys+1)
                no_dummys += dummys
                for d in range(dummys):
                    lo_allsites.append(location + (d+1)* smallblockwidth)
                
                    lo_distances.append(smallblockwidth)
            
            else:
                lo_distances.append(distance)

        print '\nadded {0} dummy stations'.format(no_dummys)
        totalstations = no_dummys+len(lo_sites)
        totalmeshblocknumber = 2*n_sidepadding+4+2*(totalstations)
        totalmodelblocknumber = 4+totalstations+(n_sidepadding - nmerge_sidepadding)/2
        print '{0} stations in total => \n\t{1} meshblocks and {2} modelblocks expected in top layer'.format(
                                    totalstations,totalmeshblocknumber, totalmodelblocknumber)

        #2. determine padding column widths:
        paddingwidth = 0.5 * np.max(lo_real_station_distances)
        meshnodelocations = []
        #add left half of block under first station 
        leftedge = lo_allsites[0] - lo_distances[0]/2.
        meshnodelocations.insert(0,leftedge)

        #add extra block column on the left of first station
        #consists of two mesh cells
        leftedge -= paddingwidth
        meshnodelocations.insert(0,leftedge)
        leftedge -= paddingwidth
        meshnodelocations.insert(0,leftedge)

            
        #3. split the inner station gaps into 2 mesh blocks each 
        print lo_allsites
        for idx,station in enumerate(lo_allsites):
            meshnodelocations.append(station)
            if idx == len(lo_allsites)-1:
                break
            #include point in the middle between here and next station
            meshnodelocations.append(station+(lo_allsites[idx+1]-station)/2.)

        #add right half of block under last station 
        rightedge = lo_allsites[-1] + lo_distances[-1]/2.
        meshnodelocations.append(rightedge)

        #add extra block column on the right of last station
        #consists of two mesh cells
        rightedge += paddingwidth
        meshnodelocations.append(rightedge)
        rightedge += paddingwidth
        meshnodelocations.append(rightedge)

        
        #add 2 side padding blocks with expon. increasing width of N mesh cells
        padding_absolute = 0
        padding_nmerge = 0
        
        for p in range(n_sidepadding):
            current_padding = sidepadding_increasefactor**(p+1)*paddingwidth
            if current_padding > 100000:
                current_padding = 100000

            padding_absolute+=current_padding
            if p < nmerge_sidepadding:
                padding_nmerge += current_padding
            
            rightedge += current_padding
            meshnodelocations.append(rightedge)

            leftedge -= current_padding
            meshnodelocations.insert(0,leftedge) 
        
        print '\nlength of model profile: {0:.1f} km (from {1:.1f} to {2:.1f})'.format(
                                        (meshnodelocations[-1]-meshnodelocations[0])/1000.,
                                        meshnodelocations[0]/1000., meshnodelocations[-1]/1000.)
        if len(meshnodelocations) > 1000:
            raise MTex.MTpyError_occam('Error - Occam cannot handle more than 1000'\
                    ' lateral mesh blocks - increase size of "max_blockwidth" !')        

        #4.determine the overall width of mesh blocks
        lo_meshblockwidths = []
        for loc in range(len(meshnodelocations)-1):            
            lo_meshblockwidths.append( meshnodelocations[loc+1] - meshnodelocations[loc] )

        #5. build top layer modelblocks by merging paddings and then 2 blocks each:
        lo_columns_to_merge = []
        lo_modelblockwidths = [] 
        current_meshblock_index = 0

        lo_columns_to_merge.append(nmerge_sidepadding)
        lo_modelblockwidths.append(padding_nmerge)
        current_meshblock_index += nmerge_sidepadding

        # merge the extra padding cells that haven't been merged into a big block at the edge
        for ii in range((n_sidepadding - nmerge_sidepadding)/2):
            lo_columns_to_merge.append(2)
            lo_modelblockwidths.append(lo_meshblockwidths[current_meshblock_index] 
                                + lo_meshblockwidths[current_meshblock_index+1] )
            current_meshblock_index += 2
            
        #merge the extra column at the left:
        lo_columns_to_merge.append(2)
        
        lo_modelblockwidths.append( lo_meshblockwidths[current_meshblock_index] 
                                + lo_meshblockwidths[current_meshblock_index+1] )
        current_meshblock_index += 2

        for idx,location in enumerate(lo_allsites):
            #each site is on top of a block, consisting of 2 mesh cells each
            lo_columns_to_merge.append(2)
            lo_modelblockwidths.append(lo_meshblockwidths[current_meshblock_index] 
                                + lo_meshblockwidths[current_meshblock_index+1] )
            current_meshblock_index += 2
        
        #merge right extra column
        lo_columns_to_merge.append(2)
        lo_modelblockwidths.append( lo_meshblockwidths[current_meshblock_index] 
                                + lo_meshblockwidths[current_meshblock_index+1] )
        current_meshblock_index += 2

        # merge the extra padding cells that haven't been merged into a big block at the edge
        for ii in range((n_sidepadding - nmerge_sidepadding)/2):
            lo_columns_to_merge.append(2)
            lo_modelblockwidths.append(lo_meshblockwidths[current_meshblock_index] 
                                + lo_meshblockwidths[current_meshblock_index+1] )
            current_meshblock_index += 2


        #merge the side padding columns on the right
        lo_columns_to_merge.append(nmerge_sidepadding)
        lo_modelblockwidths.append(np.sum(padding_absolute))
        current_meshblock_index += nmerge_sidepadding


        #6.right side of left most model block - effectively the edge of the padding
        #given with resspect to location of the first station
        #idx of station 1 is n_sidepadding + 2(extra column) + 1 (half the block under the station)
        #binding_offset =  meshnodelocations[n_sidepadding+3] - meshnodelocations[n_sidepadding]
        
        binding_offset =  meshnodelocations[nmerge_sidepadding]

        #should be identical!
        no_x_nodes = current_meshblock_index + 1
        nodey = len(lo_meshblockwidths) + 1 #vertical nodes

        ncol0 = len(lo_columns_to_merge) # number of blocks in the first layer
        
        #------
        #7. now turn to depths - set up the z axis for the mesh:
        # the model layers increas in thickness with depth
        # the increase is in decadic logarithm, starting with the given depth 
        # of the first layer

        #part to be scaled logarithmically:
        log_part_thickness = model_depth_m - (n_layers-1) * first_layer_thickness
        depths = np.logspace(np.log10(first_layer_thickness),
                                np.log10(log_part_thickness), 
                                n_layers) + np.arange(n_layers) * first_layer_thickness
       
        #no_decades = int(n_layers/layers_per_decade)+1
        #no_depthpoints_max = layers_per_decade * no_decades
        #depthscale = 10**np.linspace(0,no_decades,no_depthpoints_max + 1) 

        #lo_model_depths = list((depthscale[:n_layers-1] * first_layer_thickness))
        lo_model_depths = list(depths)
        
        lo_mesh_depths = []
        lo_rows_to_merge = []

        #2 mesh-blocks for a model layer only in the topmost layer
        for idx, depth in enumerate(lo_model_depths):
            if idx == 0:
                newdepth = depth/2.
            #else:
            #    newdepth = depth - (depth -  lo_model_depths[idx-1])/2.
                lo_mesh_depths.append(newdepth)
                lo_rows_to_merge.append(2)
                lo_mesh_depths.append(depth)
                continue
            lo_mesh_depths.append(depth)
            lo_rows_to_merge.append(1)
      

        lo_mesh_thicknesses = []

        for idx,depth in enumerate(lo_mesh_depths):
            if idx == 0:
                thickness = depth
            else:
                thickness = depth - lo_mesh_depths[idx-1]
            lo_mesh_thicknesses.append(thickness)


        max_thickness = np.max(lo_mesh_thicknesses)
        maxdepth = lo_mesh_depths[-1]

        for i in range(n_bottompadding):
            lo_mesh_thicknesses.append(max_thickness)
            lo_mesh_depths.append(maxdepth+max_thickness)
            maxdepth += max_thickness
            max_thickness *= bottompadding_scalefactor
        
        lo_model_depths.append(lo_model_depths[-1]+n_bottompadding*max_thickness)

        #just to be safe!
#        self.parameters_inmodel['no_layers'] = len(lo_model_depths)

        lo_model_thicknesses = []
        for idx,depth in enumerate(lo_model_depths):
            if idx == 0:
                thickness = depth
            else:
                thickness = depth - lo_model_depths[idx-1] 
            lo_model_thicknesses.append(thickness)


        lo_rows_to_merge.append(n_bottompadding)


        no_z_nodes = len(lo_mesh_depths) +1

        
        self.parameters_inmodel['bindingoffset']      = binding_offset
        self.parameters_inmodel['max_number_columns'] = ncol0
        self.parameters_inmodel['lo_merged_lines']    = lo_rows_to_merge
        self.meshblockwidths_x                        = lo_meshblockwidths
        self.meshblockthicknesses_z                   = lo_mesh_thicknesses
         
        self.meshlocations_z                          = lo_mesh_depths
        self.meshlocations_x                          = meshnodelocations
        self.parameters_mesh['no_nodes_x']            = no_x_nodes
        self.parameters_mesh['no_nodes_z']            = no_z_nodes

        #mesh DONE
        #-----------------------------------------------------
        #defining the actual blocks:

        trigger    = self.parameters_inmodel['block_merge_threshold']

        modelblockstrings = []
        lo_column_numbers = []
        num_params = 0
        

        ncol = ncol0
        #loop over all model layers
        #lo_modelblockwidths include the generaic merging of side padding mesh cells
        #as well as the merging of each 2 neighbouring cells
        #
        for layer_idx, thickness in enumerate(lo_model_thicknesses):
            #start with second column, because side paddings are never merged
            block_idx = 1
            #sweep columns
            while block_idx+1 < ncol-1 :
                #if two neighbouring blocks are wider than the thickness of the 
                #layer (times a factor), nothing happens
                if thickness < (trigger*(lo_modelblockwidths[block_idx]+
                                        lo_modelblockwidths[block_idx+1])):
                    block_idx += 1
                    continue
                #otherwise: avoid the vertical over-exaggeration by merging 2 
                #neighboring blocks:
                else:
                    #concatenate/merge blocks

                    lo_modelblockwidths[block_idx] += lo_modelblockwidths[block_idx+1]
                    lo_columns_to_merge[block_idx]   += lo_columns_to_merge[block_idx+1]
                    lo_modelblockwidths.pop(block_idx+1)
                    lo_columns_to_merge.pop(block_idx+1)

                    ncol -=1

            lo_column_numbers.append(ncol)

            tempstring = ""
            for j in range(ncol):
                tempstring += "%i "%(lo_columns_to_merge[j])
            tempstring += "\n"
            modelblockstrings.append(tempstring)

            num_params += ncol
        print 'depth of model (incl.padding): {0:.1f} km'.format(lo_model_depths[-1]/1000.)
        print '\nnumber of mesh layers: {0} ({1} model layers + 1 split top layer'\
                        ' + {2} bottom-padding)'.format(len(lo_mesh_thicknesses),
                                                    n_layers,n_bottompadding)
        print 'number of model blocks: {0}\n'.format(num_params)
        self.no_parameters = num_params
        self.parameters_inmodel['lo_modelblockstrings'] = modelblockstrings
        self.parameters_inmodel['lo_column_numbers']    = lo_column_numbers
        self.parameters_inmodel['lo_merged_columns']    = lo_columns_to_merge
        self.parameters_inmodel['lo_model_thicknesses'] = lo_model_thicknesses



    def write_meshfile(self):

        """
        Create the mesh file

        Attributes required:

        - self.meshlocations_x
        - self.meshlocations_z
        - self.meshnodes_x
        - self.meshnodes_z
        - self.meshfile
        - self.wd

        """
        mesh_positions_vert = self.meshlocations_z
        mesh_positions_hor  = self.meshlocations_x
        mesh_widths         = self.meshblockwidths_x
        mesh_depths         = self.meshblockthicknesses_z

        n_nodes_x           = self.parameters_mesh['no_nodes_x'] 
        n_nodes_z           = self.parameters_mesh['no_nodes_z']

        mesh_outstring =''

        temptext = '{0}\n'.format(self.parameters_mesh['mesh_title'])
        mesh_outstring += temptext

        temptext = "{0} {1} {2} {0} {0} {3}\n".format(0,n_nodes_x,n_nodes_z,2)
        mesh_outstring += temptext

        temptext = ""
        counter = 0 
        for i in range(n_nodes_x-1):
            temptext += "%.1f "%(mesh_widths[i])
            counter +=1 
            if counter == 10:
                temptext += '\n'
                counter = 0
        temptext +="\n"
        mesh_outstring += temptext

        temptext = ""
        counter = 0 
        for i in range(n_nodes_z-1):
            temptext += "%.1f "%(mesh_depths[i])
            counter +=1 
            if counter == 10:
                temptext += '\n'
                counter = 0
        temptext +="\n"
        mesh_outstring += temptext

        mesh_outstring +="%i\n"%(0)

        #occam source code has hard coded limit for reading not more than 1000
        #characters per line:
        if n_nodes_x < 1000:
            # -1 since the surface is counted as uppermost node
            for j in range(4*(n_nodes_z-1)):
                tempstring=''
                tempstring += (n_nodes_x-1)*"?"
                tempstring += '\n'
                mesh_outstring += tempstring
        else:
            #hope that occam fortran code can handle line breaks:
            #most likely not....!!!
            for j in range(4*(n_nodes_z-1)):
                tempstring=''
                counter = 0
                for k in range(n_nodes_x-1):
                    tempstring += "?"
                    counter += 1
                    if counter == 1000:
                        #tempstring += '\n'
                        counter = 0
                if counter != 0 :
                    tempstring += '\n'
                mesh_outstring += tempstring


        self.mesh = mesh_outstring

        fn = op.join(self.wd,self.meshfile)       
        F_mesh = open(fn,'w')
        F_mesh.write(mesh_outstring)
        F_mesh.close()


    def write_inmodelfile(self):
        """
        Generate inmodel file.

        Require attributes:
        - self.parameters_inmodel['lo_modelblockstrings']
        - self.parameters_inmodel['lo_columnnumbers']
        - self.parameters_inmodel['lo_merged_columns']
        - self.parameters_inmodel['bindingoffset']
        - self.no_layers
        """

        modelblockstrings = self.parameters_inmodel['lo_modelblockstrings']
        lo_merged_columns = self.parameters_inmodel['lo_merged_columns']
        lo_merged_lines   = self.parameters_inmodel['lo_merged_lines']
        lo_column_numbers = self.parameters_inmodel['lo_column_numbers']
        boffset           = self.parameters_inmodel['bindingoffset']
        n_layers          = int(self.parameters_inmodel['no_layers'])

        model_outstring =''

        temptext = "Format:           {0}\n".format("OCCAM2MTMOD_1.0")
        model_outstring += temptext
        temptext = "Model Name:       {0}\n".format(self.parameters_inmodel['model_name'])
        model_outstring += temptext
        temptext = "Description:      {0}\n".format("Random Text")
        model_outstring += temptext
        temptext = "Mesh File:        {0}\n".format(self.meshfile)
        model_outstring += temptext
        temptext = "Mesh Type:        {0}\n".format("PW2D")
        model_outstring += temptext
        if self.staticsfile is not None:
            temptext = "Statics File:     {0}\n".format(self.staticsfile)
        else:
            temptext = "Statics File:     none\n"
        model_outstring += temptext
        if self.prejudicefile is not None:
            temptext = "Prejudice File:   {0}\n".format(self.prejudicefile)
        else:
            temptext = "Prejudice File:   none\n"
        model_outstring += temptext
        temptext = "Binding Offset:   {0:.1f}\n".format(boffset)
        model_outstring += temptext
        temptext = "Num Layers:       {0}\n".format(n_layers + 1 )
        model_outstring += temptext

        for k in range(n_layers+1):
            n_meshlayers  = lo_merged_lines[k]
            n_meshcolumns = lo_column_numbers[k]
            temptext="{0} {1}\n".format(n_meshlayers, n_meshcolumns)
            model_outstring += temptext
            temptext = modelblockstrings[k]
            model_outstring += temptext
            #model_outstring += "\n"
            

        if self.parameters_inmodel['build_roughness_exceptions']:
            for interface in self.parameters_inmodel['interface_filelist']:
                self.project_interface(interface)
            if self.parameters_inmodel['elevation_filename'] is not None:
                self.subtract_elevation(self.parameters_inmodel['elevation_filename'])
            self.build_roughness_exceptions()
            re = self.parameters_inmodel['roughness_exceptions']
            temptext = "Number Exceptions:{0}\n".format(-len(re))
            temptext += '\n'.join([' '.join([str(rr) for rr in r]+['0']) for r in re])
        else:
            temptext = "Number Exceptions:{0}\n".format(0)

        model_outstring += temptext
        
        self.inmodel = model_outstring
        fn = op.join(self.wd,self.inmodelfile)        
        F_model = open(fn,'w')
        F_model.write(model_outstring)
        F_model.close()



    def write_startupfile(self):
        """
        Generate startup file

        Require attributes:

        -


        """

        startup_outstring =''

        temptext = "Format:           {0}\n".format(self.parameters_startup['iter_format'])
        startup_outstring += temptext

        temptext = "Description:      {0}\n".format(self.parameters_startup['description'])
        startup_outstring += temptext

        temptext = "Model File:       {0}\n".format(self.inmodelfile)
        startup_outstring += temptext

        temptext = "Data File:        {0}\n".format(self.datafile)
        startup_outstring += temptext

        temptext = "Date/Time:        {0}\n".format(self.parameters_startup['datetime_string'])
        startup_outstring += temptext

        temptext = "Max Iter:         {0}\n".format(int(float(self.parameters_startup['max_no_iterations'])))
        startup_outstring += temptext

        temptext = "Target Misfit:    {0:.2f}\n".format(float(self.parameters_startup['target_rms']))
        startup_outstring += temptext

        temptext = "Roughness Type:   {0}\n".format(int(self.parameters_startup['roughness_type']))
        startup_outstring += temptext
    
        temptext = "Debug Level:      {0}\n".format(int(self.parameters_startup['debug_level']))
        startup_outstring += temptext

        temptext = "Iteration:        {0}\n".format(int(self.parameters_startup['no_iteration']))
        startup_outstring += temptext
    
        temptext = "Lagrange Value:   {0}\n".format(self.parameters_startup['mu_start'])
        startup_outstring += temptext
        
        temptext = "Roughness Value:  {0}\n".format(self.parameters_startup['roughness_start'])
        startup_outstring += temptext
        
        temptext = "Misfit Value:     {0}\n".format(float(self.parameters_startup['rms_start']))
        startup_outstring += temptext
        
        temptext = "Misfit Reached:   {0}\n".format(int(self.parameters_startup['reached_misfit']))
        startup_outstring += temptext
        
        temptext = "Param Count:      {0}\n".format(int(self.no_parameters))
        startup_outstring += temptext
        
        temptext = ""
        counter = 0 

        if type(self.parameters_startup['halfspace_resistivity']) not in [list,np.ndarray]:
            for l in range(self.no_parameters):
                temptext += "{0:.1g}  ".format(np.log10(float(self.parameters_startup['halfspace_resistivity'])))
                counter += 1
                if counter == 20:
                    temptext += '\n'
                    counter = 0
        else:
            temptext +='  '.join([str(res) for res in self.parameters_startup['halfspace_resistivity']])

        temptext += "\n"
        startup_outstring += temptext
        
        
     

        fn =  op.join(self.wd,self.startupfile)
        F_startup = open(fn,'w')
        F_startup.write(startup_outstring)
        F_startup.close()



    def generate_inputfiles(self, edi_dir=None):

        edi_directory = self.edi_directory
        if edi_dir is not None:
            if op.isdir(edi_dir):
                edi_directory = edi_dir 

        if not op.isdir(self.wd):
            os.makedirs(self.wd)

        self.read_edifiles(edi_directory)
        try:
            self.write_datafile()
        except:
            raise


        self.setup_mesh_and_model()
        self.write_meshfile()
        self.write_inmodelfile()
        self.write_startupfile()
        self.write_configfile()
        if self.prejudicefile is not None:
            self.write_prejudicefile()

        print '\nInput files in working directory {0}: \n'.format(op.abspath(self.wd))
        print '{0}'.format(op.basename(self.datafile))
        print '{0}'.format(op.basename(self.meshfile))
        print '{0}'.format(op.basename(self.inmodelfile))
        print '{0}'.format((self.startupfile))
        print '{0}'.format(op.basename(self.configfile))


        print '\n\t\t DONE !\n\n'



    def write_configfile(self):
        wd = op.abspath(self.wd)
        fn = 'occam2d_configuration.cfg'

        self.configfile = op.join(wd,fn)

        all_configs_dict = {}
        for dictionary in [self.parameters_startup, self.parameters_inmodel, 
                                    self.parameters_mesh, self.parameters_data]:

            for key,value in dictionary.items():
                if type(value) in [str,float,int]:
                    all_configs_dict[key] = value 

        for key,value in vars(self).items():
            if value is None:
                all_configs_dict[key] = value
                continue
            if type(value) in [float,int]:
                all_configs_dict[key] = value
                continue
            if type(value) in [str]:
                if len(value.split())>3:
                    continue
                all_configs_dict[key] = value

        for key, value in  vars(vars(self)['Data']).items():
            if value is None:
                all_configs_dict[key] = value
                continue
            if type(value) in [str,float,int]:
                all_configs_dict[key] = value
                continue

        occam_run_dict = {}
        occam_run_dict['OCCAM_run'] = all_configs_dict
        # print all_configs_dict
        # sys.exit()
        MTcf.write_dict_to_configfile(occam_run_dict,self.configfile)

        return 

    def make_savepath(self,basename='inv',start=0,suffix=''):
        """
        make a directory name in which to save results
        ensures that each new inversion is saved in a separate directory
        
        Author: Alison Kirkby
        """
        
        
        while os.path.exists(os.path.join(self.wd,basename + '%03i'%start + suffix)):
            start += 1
        
        svpath = basename + '%03i'%start +suffix
        
        return svpath
        

    def get_modelblock_info(self):
        """
        get model block locations and model block numbers and store them in the
        setup object.        
        
        Author: Alison Kirkby (2013)
        """
        blockmerges_x = [[int(val) for val in line.strip('\n').split()] for line in self.parameters_inmodel['lo_modelblockstrings']]
        print blockmerges_x[0]
        blockmerges_z = self.parameters_inmodel['lo_merged_lines']
        meshlocations_z = np.array([0]+self.meshlocations_z)
        
#        meshlocations_z.insert(0,0)
        num=1
        block_nums,block_x,block_z = [],[],[self.meshlocations_z[0]]
        block_nums_array = np.zeros((len(self.meshlocations_z)-1,len(self.meshlocations_x)-1))
        checkerbox = np.zeros_like(block_nums_array)
        iz = 0
        for j,z in enumerate(blockmerges_z):
            ix = 0
            block_z.append(meshlocations_z[iz+z-1])
            tmp_lst = []
            # list to contain x block locations for this row
            tmp_lst_xpos = [self.meshlocations_x[0]]
            for i,x in enumerate(blockmerges_x[j]):
#                if ((x > 2) and (j==0)):
#                    print x,ix,self.meshlocations_x[ix:ix+x]
#                xpos = np.mean(np.array(self.meshlocations_x[ix:ix+x]))
                tmp_lst.append(num)
                # append next x location
                tmp_lst_xpos.append(self.meshlocations_x[ix+x-1])
                block_nums_array[iz:iz+z,ix:ix+x] = num
                checkerbox[iz:iz+z,ix:ix+x] = np.random.random()
                ix += x
                num += 1
            iz += z
            block_nums.append(np.array(tmp_lst))
            block_x.append(np.array(tmp_lst_xpos))
        self.modelblocklocations_x = np.array(block_x)
        self.modelblocklocations_z = np.array(block_z)
        self.block_nums = np.array(block_nums)
        self.block_nums_array = np.array(block_nums_array)
        self.checkerbox = checkerbox

    
    def read_2dinterfacefile(self,fn,skiprows=1):
        

        # get origin of profile
        if self.Data.profile_origin is None:
            self.Data.get_profile_origin()
        xp0,yp0 = self.Data.profile_origin
        
        # get origin for interface from the first line
        ff = open(fn)
        xd0,yd0 = [float(val) for val in ff.readline().strip().split()[-2:]]
        # get length along profile and depth
        ld,zd = np.loadtxt(fn,skiprows=skiprows,unpack = True)

        # shift l so it aligns with the profile
        shift = ((xd0-xp0)**2. + (yd0-yp0)**2.)**0.5

        ld += shift
        
        if np.median(zd) < 0: zd = -zd
        
        return ld,zd

    def project_interface(self,fn,skiprows=1,fmt='xy'):
        """
        project an interface (from an xyz text file) onto the profile.
        adds a numpy array of z locations (corresponding to x model block locations)
        to the variable interfaces

        intdir = directory where interface text files are held
        fn = filename of interface text file
        skiprows = number of header rows to skip when reading file  
        fmt = xy (eastings, northings, depth) or l (length along profile, depth) with origin on first line       
        
        Author: Alison Kirkby (2013)
        """
        
        self.Data.get_profile_origin()
        self.get_modelblock_info()
            
        az = np.deg2rad(90-self.Data.azimuth)
        xp0,yp0 = self.Data.profile_origin
        
        xp = self.Data.profile_origin[0]+self.modelblocklocations_x[0]*np.cos(az)
        yp = self.Data.profile_origin[1]+self.modelblocklocations_x[0]*np.sin(az)
        
        if fmt == 'xy':
            data = np.loadtxt(fn,skiprows=skiprows)
            f = si.CloughTocher2DInterpolator(data[:,0:2],data[:,2])
            zp = np.array([f(xp[i],yp[i]) for i in range(len(xp))])
        else:
            ld,zd = self.read_2dinterfacefile(fn)

            f = si.interp1d(ld,zd,bounds_error=False)
            zp = f(self.modelblocklocations_x[0])

        
        
        # check for any null values
        for i,z in enumerate(zp):
            if np.isnan(z):
                ii = 0
                while np.isnan(zp[i+ii]):
                    if i < len(zp)/2:
                        ii += 1
                    else:
                        ii -= 1
                zp[i] = zp[i+ii]

        if not hasattr(self,'interfaces'):
            self.interfaces = []
        self.interfaces.append(zp)
            
        
    def subtract_elevation(self,intdir,elevfn):
        """
        subtract elevation from all projected interfaces in interface variable
        
        Author: Alison Kirkby (2013)
        """
        if hasattr(self,'interfaces'):    
            self.project_interface(elevfn) 
            
            for i in self.interfaces[:-1]:
                i -= self.interfaces[-1]
                
            self.interfaces = self.interfaces[:-1]
        else:
            print "nothing to subtract; please use project_interface to define interfaces first"
            return
    
    def build_roughness_exceptions(self):
        """
        generate a list of model block pairs between which to build a roughness
        exception and add these to the inmodel file
        
        Author: Alison Kirkby (2013)
        """
        
        if not hasattr(self,'interfaces'):
            print "nothing to build roughness exceptions from; please use project_interface to define interfaces first"
            return
        else:
            # get block centres
            block_z = (self.modelblocklocations_z[1:] + self.modelblocklocations_z[:-1])/2.
            block_x = [(self.modelblocklocations_x[i][1:] + self.modelblocklocations_z[i][:-1])/2. \
            for i in range(len(self.modelblocklocations_x))]
            x_vals = block_x[0]
            block_nums = self.block_nums

            blockmerges_x = [[int(val) for val in line.strip('\n').split()] for line in self.parameters_inmodel['lo_modelblockstrings']]
            
            values = []          
            cell_x = []
            
            for interface in self.interfaces:
                z_msl = -interface
                values_tmp = []
                cell_x_tmp = []
                values_tmp_v = []
                
                for i in range(len(interface)):
                    z = block_z[block_z-z_msl[i]<0][-1] # get out first block for which z_msl value is greater than z location of block
                    j = list(block_z).index(z) # get the index j for the z value
                    
                    crit = abs(block_x[j]-x_vals[i])==min(abs(block_x[j]-x_vals[i])) #  define criteria to find closest x block
                    bn1 = block_nums[j][crit][0]
                    ii = list(block_nums[j]).index(bn1)
                    iii = 0
                    while sum(blockmerges_x[j+1][:iii]) <= sum(blockmerges_x[j][:ii]):
                        iii += 1
                    iii -= 1
                    bn2 = block_nums[j+1][iii]
                    if not [block_x[j][ii],block_x[j+1][iii]] in cell_x_tmp:
                        values_tmp.append([bn1,bn2])
                        cell_x_tmp.append([block_x[j][ii],block_x[j+1][iii]])
                    if i != 0:
                        bmin1 = np.min(values_tmp[-1][0])
                        bmin2 = np.min(values_tmp[-2][0])
                        bmax1 = np.max(values_tmp[-1][0])
                        bmax2 = np.max(values_tmp[-2][0])
                        if bmin1 - bmin2 != 1:
                            if bmin1 > bmin2:                                                                 
                                bn3 = [bmin1-1,bmin1]
                                bn4 = [bmax2,bmax2+1]
                            else: 
                                bn3 = [bmin2,bmin2+1] 
                                bn4 = [bmax1-1,bmax1]
                            for bval in [bn3,bn4]:
                                if not bn3 in values_tmp_v:
                                    values_tmp_v.append(bval)
                                
                
                values += values_tmp+values_tmp_v
                cell_x.append(cell_x_tmp)
                
            if self.parameters_inmodel['roughness_exceptions'] is None:
                self.parameters_inmodel['roughness_exceptions'] = []  

            self.parameters_inmodel['roughness_exceptions'] += values

        
    def project_well_loc(self):
        """
        project a well location onto the profile.
        stores the result in the setup variable well_locations.
        x,y = x,y locations of well
        
        Author: Alison Kirkby
        """
        distances = []

        for well_xy in self.well_xy:
        
            [xp,yp] = well_xy       
            
            # get profile origin
            self.Data.get_profile_origin()
            
            # project drillhole location onto profile and get distance along profile
            [m,c1] = self.Data.profile
            [x0,y0] = self.Data.profile_origin
            x = (yp+(1.0/m)*xp-c1)/(m+(1.0/m))
            y = m*x+c1
            x -= x0
            y -= y0
            distances.append((x**2.+y**2.)**(0.5))
        self.well_locations = distances
        
        
    def setup_prejudice(self):
        
        for key in ['modelblocklocations_x','modelblocklocations_z',
                    'block_nums','block_nums_array']:
            if not hasattr(self,key):
                self.get_modelblock_info()
                break
        if self.parameters_prejudice['interface_filelist'] is not None:
            for intfile in self.parameters_prejudice['interface_filelist']:
                ld,zd = self.read_2dinterfacefile(intfile)
                self.parameters_prejudice['interfaces'].append(np.vstack([ld,zd]).T)
        self.Prejudice = Prejudice(self,**self.parameters_prejudice)
             
             
    def write_prejudicefile(self):
        """
        build a prejudice file
        """
        
        self.setup_prejudice()
        
        header  = 'FORMAT:           OCCAM2MTPREJ_2.0\n'
        
        prejvals = np.hstack(self.Prejudice.prejudice)
        weights = np.hstack(self.Prejudice.weights)
        data = np.vstack([np.hstack(self.block_nums),np.log10(prejvals),weights]).T
        header += 'NO. PARMS:        {}'.format(len(data))
        
        # assign prejudice values also as starting resistivities
        self.parameters_startup['halfspace_resistivity'] = np.log10(prejvals)
        self.write_startupfile()

        np.savetxt(op.join(self.wd,self.prejudicefile),data,fmt=['%1i','%.3f','%.3f'],header=header,comments='')
        

class Prejudice():
    def __init__(self,Setup,**input_parameters):
        

        self.constraints_type = 'interfaces'
        self.interfaces = None # list (in order of increasing depth) of arrays containing 2 columns, 
                               # distance x along profile and depth z of interface
        self.interface_resistivity_values = None # list of prejudice resistivities, 1 longer than interfaces. in same order as interfaces
        self.interface_prejudice_weights = None
        self.welldata = None

        # transfer necessary parameters from setup object
        for key in ['modelblocklocations_x','modelblocklocations_z',
                    'block_nums','block_nums_array',
#                    'interface_prejudice_weights','interface_resistivity_values',
                    'meshlocations_x','meshlocations_z']:
            if hasattr(Setup,key):
                setattr(self,key,getattr(Setup,key)) 
            else:
                print "Can't build prejudice, need model block information from setup object"
    
        for key in input_parameters.keys():
            if hasattr(self,key):
                setattr(self,key,input_parameters[key])
        
        if self.constraints_type == 'interfaces':
            self.assign_interface_prejudice()
    

    def assign_interface_prejudice(self):
        """
        Need to have interfaces and interface resistivities defined to use this option
        
        """
        if ((self.interfaces is None) or (self.interface_resistivity_values is None)):
            print "Can't build prejudices, need to define interfaces"
            return
        else:
            # remove zero and negative resistivities
            self.interface_resistivity_values = [max(1e-6,rv) for rv in self.interface_resistivity_values]
            # initialise a prejudice array. Starting resistivity is the value for the 1st layer then successively move down the model
            prejudice_array = np.ones_like(self.block_nums_array)*self.interface_resistivity_values[0]
            weights_array = np.zeros_like(self.block_nums_array)
#            blockx,blockz = self.modelblocklocations_x,self.modelblocklocations_z
            meshx,meshz = np.meshgrid(self.meshlocations_x,self.meshlocations_z)
            # get centres of mesh blocks
            meshx = (meshx[1:,1:] + meshx[1:,:-1])/2.
            meshz = (meshz[1:,1:] + meshz[:-1,1:])/2.
#            print meshx[0]
            for i in range(len(self.interfaces)):
                ld,zd = self.interfaces[i].T
                # create interpolation function to use if needed
                f = si.interp1d(ld,zd,bounds_error=False)
                # empty 1d array to contain depths
                zline = np.zeros(len(meshx[0]))
                for ii in range(len(meshx[0])):

                    if ii == 0:
                        zc0 = meshx[0][ii]-10.
                        zc1 = (meshx[0][ii]+meshx[0][ii+1])/2.
                    elif ii == len(meshx[0])-1:
                        zc1 = meshx[0][ii]+10.
                        zc0 = (meshx[0][ii]+meshx[0][ii-1])/2.
                    else:
                        zc0 = (meshx[0][ii]+meshx[0][ii-1])/2.
                        zc1 = (meshx[0][ii]+meshx[0][ii+1])/2.                        
#                    print zc0,zc1
                    zvalues = zd[(ld>zc0)&(ld<zc1)]
                    
                    # if there are >2 z values, take the median
                    if len(zvalues) > 2:
                        zval = np.median(zvalues)
                    # otherwise use the interpolation function to get value at the mesh node
                    else:
                        zval = f(meshx[0][ii])
#                        print blockx[jj][ii],blockx[jj][ii+1]

                    zline[ii] = zval
                    
                    # deal with edges
                    zline[(meshx[0]<=ld[0])] = zd[0]
                    zline[(meshx[0]>=ld[-1])] = zd[-1]

                # assign all values below zval to the value for the next layer
                prejudice_array[meshz > zline] = self.interface_resistivity_values[i+1]
                weights_array[meshz > zline] = self.interface_prejudice_weights[i+1]
                
            # now assign mesh values to block numbers
            prejudice = self.block_nums.copy()*0.
            weights = self.block_nums.copy()*0.
            for jj in range(len(prejudice)):
                for ii in range(len(prejudice[jj])):
#                    print prejudice[jj]
#                    print self.block_nums[jj][ii],prejudice_array[self.block_nums_array==self.block_nums[jj][ii]]
                    prejudice[jj][ii] = np.median(prejudice_array[self.block_nums_array==self.block_nums[jj][ii]])
                    weights[jj][ii] = np.median(weights_array[self.block_nums_array==self.block_nums[jj][ii]])
            self.prejudice = prejudice
            self.prejudice_array = prejudice_array
            self.weights = weights
            self.weights_array = weights_array                                     
    
    def subsample_welldata(self):
        """
        subsample a resistivity array by averaging
        input: 1-D array or list of 1-D arrays containing resistivity and depth data
        result: resitivity values and block indices stored in variables 
        prejudice_resitivity and prejudice_indices_z
        
        Author: Alison Kirkby
        """
        
        blockmerges_z = self.parameters_inmodel['lo_merged_lines']
        blockedges_z = [0]+[self.meshlocations_z[i-1] for i in [int(sum(blockmerges_z[:ii])) for ii in range(1,len(blockmerges_z)+1)]]

        
        if not hasattr(self,'prejudice_resisitivity'):
            self.prejudice_resistivity = []
        
        if not hasattr(self,'prejudice_indices_z'):
            self.prejudice_indices_z = []
                
        for w in range(len(self.welldata)):
            self.prejudice_resistivity.append([])
            self.prejudice_indices_z.append([])
                        
            for i in range(1,len(blockedges_z)):
                temp_res = []
                for j in range(len(self.welldata[w])):
                    if blockedges_z[i-1]<self.welldata[w][j,1]<blockedges_z[i]:
                        temp_res.append(self.welldata[w][j,0])
                if len(temp_res)>0:
                    self.prejudice_resistivity[-1].append(np.around(np.median(temp_res),1))
                    self.prejudice_indices_z[-1].append(i)                



        
 
       
    def get_blocknums_from_wells(self):
        """
        get block numbers for prejudice file
        requires prejudice_resistivity and prejudice_indices_z to be defined
        
        prejudice_weights = weighting (range 0.0 to 1.0) for each prejudice
        value (see occam documentation for definition)
        
        Author: Alison Kirkby        
        """

        
        res_values = self.prejudice_resistivity
        num_values = sum([len(l) for l in res_values])        
        
        if not hasattr(self,'prejudice_blocknums'):
            self.prejudice_blocknums = []
        if not hasattr(self,'prejudice_fstrings'):
            self.prejudice_fstrings = ["FORMAT:           OCCAM2MTPREJ_2.0\n",
                                       "NO. PARAMS:       "+str(num_values)+'\n']
     
        i = 0
        
        for w in range(len(res_values)):
            self.prejudice_blocknums.append([])
            distance = self.well_locations[w]

            for i in range(len(res_values[w])):
                row = self.prejudice_indices_z[w][i]
                # define model block locations for given row
                modc = self.modelblocklocations_x[row]
                # find index of closest model block to x location of drillhole:
                col = list(abs(modc - distance)).index(np.min(abs(modc - distance)))
                # append block nums
                blocknum = self.block_nums[row][col]
                self.prejudice_blocknums[w].append(blocknum)
                self.prejudice_fstrings.append(' '.join(['%01i'%blocknum,
                                                         '%.03f'%np.log10(self.prejudice_resistivity[w][i]),
                                                         '%.01f'%self.prejudice_weight])+'\n')

    def build_prejudicefile_wells(self):
        """
        build model prejudices from well data and write to a file.
        
        first need to define:
        self.well_xy = list of lists containing xy coordinates
        self.welldata = list of 2-D arrays containing res, depth 
        (shape = (num_wells,num_datapoints,2)) for each well_xy

        
        Author: Alison Kirkby        
        """    
        if self.well_xy is None:
            print "please first define well xy location"
            return
            
        if self.welldata is None:
            print "please first define well data"
            return

        if not hasattr(self,'block_nums'):
            self.get_modelblock_info() 
        
        self.project_well_loc()
        self.prejudice_subsample_welldata()
        self.prejudice_get_blocknums()
        
    
        pf = open(os.path.join(self.wd,self.prejudicefile),'wb')
        pf.writelines(self.prejudice_fstrings)
        pf.close()




        
#------------------------------------------------------------------------------


class Data():
    """
    Handling input data.

    Generation of suitable Occam data file(s) from Edi files/directories.
    Reading and writing data files.
    Allow merging of data files.
    """
    def __init__(self, edilist = None, wd = None, **data_parameters):

        self.wd = op.abspath(os.curdir)
        self.datafile = 'OccamDataFile.dat'
        self.edilist = []

        if edilist is not None:
            if np.iterable(edilist):
                self.edilist = edilist
        
        if wd is not None:
            if op.isdir(wd):
                self.wd = op.abspath(wd)

        self._strike_set = False
        self.strike = None
        self.azimuth = 0.
        self.profile = None
        self.frequencies = None
        self.stations = []
        self.stationlocations = []
        self.data = []
        self.mode = 'tetm'
        self.profile_offset = 0.
        self.format = 'OCCAM2MTDATA_1.0'
        self.title = 'MTpy-OccamDatafile'
        self.edi_type = 'z'
        self.profile_origin = None

        self.phase_errorfloor = 5
        self.rho_errorfloor = 10
        self.tipper_errorfloor = 5
        self.tipper_errorfloor_abs = None

        self.min_frequency = None
        self.max_frequency = None
        self.max_no_frequencies = None


        for key in data_parameters:
            setattr(self,key,data_parameters[key])
        
        if ('strike' in data_parameters) and (data_parameters['strike'] is not None):
            self._strike_set = True

        try:
            self.generate_profile()
        except:
            print 'cannot generate profile'
            raise
        
        try:
            self.build_data()
        except:
            print 'cannot build data file'
            raise

       

    def readfile(self,fn):
        if not op.isfile(fn):
            print 'Error - not a valid file: {0}'.fn

        self.filename = op.basename(fn)
        self.wd = op.split(fn)[0]

        F_in = file(fn,'r')
        datafile_raw = F_in.read()
        F_in.close()

        #string is reduced each step, i.e. cut off the sections, 
        #which are already read in
        reduced_string = self._read_format(datafile_raw)
        reduced_string = self._read_title(datafile_raw)
        reduced_string = self._read_sites(datafile_raw)
        reduced_string = self._read_offsets(datafile_raw)
        reduced_string = self._read_frequencies(datafile_raw)
        
        self._read_data(reduced_string)

    def _find_string(key,datastring):
        
        index = datastring.lower().find('key')
        return index

    def _read_format(self,datastring):
        idx = _find_string('format',datastring)
        reduced_string = datastring[idx:]
        data_list = datastring.split('\n')
        line = data_list[0]
        line = line.strip().split(':')
        self.format = line[1].strip().lower()

        return reduced_string 

    def _read_title(self,datastring):
        idx = _find_string('title',datastring)
        reduced_string = datastring[idx:]
        data_list = datastring.split('\n')
        line = data_list[0]
        line = line.strip().split(':')
        self.title = line[1].strip().lower()

        return reduced_string 
        

    def _read_sites(self,datastring):

        idx = _find_string('sites',datastring)
        reduced_string = datastring[idx:]
        data_list = datastring.split('\n')
        line = data_list[0]
        line = line.strip().split(':')
        no_sites = int(float(line[1].strip().lower()))
        lo_stations = []
        for idx in range(no_sites):
            sta = data_list[idx+1].strip()
            lo_stations.append(sta)

        self.stations = lo_stations

        return reduced_string 
        

    def _read_offsets(self,datastring):
        idx = _find_string('offsets',datastring)
        reduced_string = datastring[idx:]
        data_list = datastring.split('\n')
        line = data_list[0]
        line = line.strip().split(':')
        no_sites = len(self.stations)
        lo_offsets = []
        for idx in range(no_sites):
            offset = float(data_list[idx+1].strip())
            lo_offsets.append(offset)

        self.stationlocations = lo_offsets

        return reduced_string 
        

    def _read_frequencies(self,datastring):
        idx = _find_string('frequencies',datastring)
        reduced_string = datastring[idx:]
        data_list = datastring.split('\n')
        line = data_list[0]
        line = line.strip().split(':')
        no_freqs = int(float(line[1]))

        lo_freqs = []
        for idx in range(no_freqs):
            freq = float(data_list[idx+1].strip())
            lo_freqs.append(freq)

        self.frequencies = lo_freqs

        return reduced_string 

    def _read_data(self,datastring):
        idx = _find_string('data',datastring)
        reduced_string = datastring[idx:]
        data_list = datastring.split('\n')
        line = data_list[0]
        line = line.strip().split(':')
        no_data = int(float(line[1]))

        lo_data = []
        idx = 0
        row_idx = 2
        while idx < no_data:
            row = data_list[row_idx].strip().split()
            if row[0][0] == '#':
                row_idx += 1
                continue
            rowlist = [float(i) for i in row]
            lo_data.append(rowlist)
            row_idx += 1
            idx += 1


    def build_data(self):
        """Data file Generation

        Read all Edi files. 
        Extract frequencies. 
        Read in strike. If strike = None: find average strike over all stations and frequencies. 
        90 degree strike ambiguity leads to choice of strike: Larger angle with profile line is chosen.
        Rotate Z and Tipper: X components are along strike, Y orthogonal.  
        Extract off-diagonal data from Z. Extract Tipper x-component (along profile).

        Collect all information sorted according to occam specifications.

        Data of Z given in muV/m/nT = km/s
        Error is assumed to be 1 stddev.
        """ 
        

        #set data modes
        lo_modes = []
        modes = self.mode.lower().strip()
        
        # adding functionality to do log resistivity
        if 'log' in modes:
            rescode_te, rescode_tm = 1,5
        else:
            rescode_te, rescode_tm = 9,10
        
        if 'both' in modes :
            lo_modes.extend([rescode_te, rescode_tm,2,6])  
        if 'te' in modes:
            lo_modes.extend([rescode_te,2])
        if 'tm' in modes:
            lo_modes.extend([rescode_tm,6])
        if ('tipper' in modes): 
            lo_modes.extend([3,4])
        if 'all' in modes :
            lo_modes.extend([rescode_te, rescode_tm,2,6,3,4])


        lo_modes = sorted(list(set(lo_modes))) 

        #set data frequencies
        min_freq = self.min_frequency
        max_freq = self.max_frequency
        no_freqs_max = self.max_no_frequencies
 

        lo_all_freqs = []
        for lo_f in self.station_frequencies:
            lo_all_freqs.extend(list(lo_f))
        lo_all_freqs = sorted(list(set(lo_all_freqs)),reverse=True)

        if (min_freq is None) or (min_freq < min(lo_all_freqs) ) or (min_freq > max(lo_all_freqs) ) :
            min_freq = min(lo_all_freqs)
        if (max_freq is None) or (max_freq > max(lo_all_freqs) ) or (max_freq < min(lo_all_freqs) ) :
            max_freq = max(lo_all_freqs)
        
        lo_all_freqs_tmp = []
        for f in  lo_all_freqs:
            if min_freq <= f <= max_freq :
                lo_all_freqs_tmp.append(f)
            else:
                continue

        if len(lo_all_freqs_tmp) == 0 :
            print 'No frequencies in user-defined interval [{0},{1}]'.format(min_freq, max_freq)
            sys.exit()


        #check, if frequency list is longer than given max value
        if no_freqs_max is not None:
            no_freqs_max = int(float(no_freqs_max))
            if no_freqs_max < len(lo_all_freqs_tmp):
                lo_all_freqs_tmp2 = []
                excess = len(lo_all_freqs_tmp)/float(no_freqs_max)
                if excess < 2:
                    offset = 0

                else:
                    stepsize = (len(lo_all_freqs_tmp)-1)/no_freqs_max
                    offset = stepsize/2.
                indices = np.array(np.around(np.linspace(offset,len(lo_all_freqs_tmp)-1-offset,no_freqs_max),0))
                if indices[0]>(len(lo_all_freqs_tmp)-1-indices[-1]):
                    indices -= 1
                for idx in indices:
                    index = int(np.round(idx,0))+1
                    lo_all_freqs_tmp2.insert(0,lo_all_freqs_tmp[-index])
                

                lo_all_freqs_tmp = lo_all_freqs_tmp2

        self.frequencies = np.array(lo_all_freqs_tmp)


        #collect data 
        self.data = []

        for idx_s, station in enumerate(self.stations):
            station_number = idx_s + 1
            Z = self.Z[idx_s]
            T = self.Tipper[idx_s]

            rho = Z.resistivity
            phi = Z.phase
            rho_err = Z.resistivity_err
            phi_err = Z.phase_err
            z_array = Z.z
            zerr_array = Z.zerr

            for freq_num,freq in enumerate(self.frequencies):

                frequency_number = freq_num + 1 #OCCAM indices start with 1 

                #extract the freqs available for the respective station
                station_freqs = self.station_frequencies[idx_s]
                #skip, if the listed frequency is not available for the station
                if not (freq in station_freqs):
                    continue

                #find the respective frequency index for the station     
                idx_f = np.abs(station_freqs-freq).argmin()

                for mode in lo_modes:
                    append = True
                    if mode in [9,1,2] :
                        raw_rho_value = rho[idx_f][0,1]
                        value = raw_rho_value
                        #value = np.log10(raw_rho_value)
                        absolute_rho_error = rho_err[idx_f][0,1]
                        try:
                            if raw_rho_value == 0:
                                raise
                            relative_rho_error = np.abs(absolute_rho_error/raw_rho_value)
                        except:
                            relative_rho_error = 0.

                        # can't have zero errors, occam crashes
                        if relative_rho_error == 0:
                            if self.rho_errorfloor is not None:
                                relative_rho_error = self.rho_errorfloor/100.

                        if mode == 9 :
                            if self.rho_errorfloor is not None:
                                if self.rho_errorfloor/100. > relative_rho_error:
                                    relative_rho_error = self.rho_errorfloor/100.
                            error = np.abs(relative_rho_error * raw_rho_value)   #relative_error/np.log(10.)
                            #error = np.abs(relative_rho_error/np.log(10.))

                        elif mode == 1 :
                            if self.rho_errorfloor is not None:
                                if self.rho_errorfloor/100. > relative_rho_error:
                                    relative_rho_error = self.rho_errorfloor/100.
                            error = np.abs(relative_rho_error)/np.log(10)
                            value = np.log10(value)

                        elif mode == 2 :
                            raw_phi_value = phi[idx_f][0,1]
                            if raw_phi_value >=180:
                                raw_phi_value -= 180
                            value = raw_phi_value %180
                            # phase needs to be between 0 and 90 for occam
                            if value > 90:
                                append = False
                            if self.phase_errorfloor is not None:
                                if self.phase_errorfloor/100. > relative_rho_error:
                                    relative_rho_error = self.phase_errorfloor/100.
                            if relative_rho_error >= 2.:
                                error = 180.
                            else:
                                error = np.degrees(np.arcsin(0.5*relative_rho_error))#relative_error*100.*0.285
                            # maximum possible error +/- 90 degrees 
                            error = min(error,90.)
                        
                            
                    if mode in [10,5,6] :
                        raw_rho_value = rho[idx_f][1,0]
                        value = raw_rho_value
                        #value = np.log10(raw_rho_value)
                        absolute_rho_error = rho_err[idx_f][1,0]
                        
                        try:
                            if raw_rho_value == 0:
                                raise
                            relative_rho_error = np.abs(absolute_rho_error/raw_rho_value)
                        except:
                            relative_rho_error = 0.

                        # can't have zero errors, occam crashes
                        if relative_rho_error == 0:
                            if self.rho_errorfloor is not None:
                                relative_rho_error = self.rho_errorfloor/100.                       
                            
                        if mode == 10 :
                            if self.rho_errorfloor is not None:
                                if self.rho_errorfloor/100. > relative_rho_error:
                                    relative_rho_error = self.rho_errorfloor/100.
                            error = np.abs(relative_rho_error * raw_rho_value)   #relative_error/np.log(10.)
                            #error = np.abs(relative_rho_error /np.log(10.))

                        elif mode == 5 :
                            if self.rho_errorfloor is not None:
                                if self.rho_errorfloor/100. > relative_rho_error:
                                    relative_rho_error = self.rho_errorfloor/100.
                            error = np.abs(relative_rho_error)/np.log(10)
                            value = np.log10(value)

                        elif mode == 6 :
                            raw_phi_value = phi[idx_f][1,0]
                            if raw_phi_value >=180:
                                raw_phi_value -= 180
                            value = raw_phi_value %180
                            # phase needs to be between 0 and 90 for occam
                            if value > 90:
                                append = False
                            if self.phase_errorfloor is not None:
                                if self.phase_errorfloor/100. > relative_rho_error:
                                    relative_rho_error = self.phase_errorfloor/100.
                            if relative_rho_error >= 2.:
                                error = 180.
                            else:
                                error = np.degrees(np.arcsin(0.5*relative_rho_error))#relative_error*100.*0.285
                            # maximum possible error +/- 90 degrees
                            error = min(error,90.)
                            

                    elif mode in [3,4] :
                        if T.tipper is None:
                            #print 'no Tipper data for {0} Hz at station {1}'.format(freq, station_number) 
                            continue

                        tipper = T.tipper[idx_f]
                        try: 
                            tippererr = T.tippererr[idx_f]
                        except:
                            #print 'no Tipper error for station {0}/frequency {1}'.format(station_number,frequency_number)
                            tippererr = None


                        if mode == 3 :
                            value = np.real(tipper[0,1])

                        if mode == 4 :
                            value = np.imag(tipper[0,1])

                        # get tipper error if it exists
                        if tippererr is None:
                            raw_error = 0
                            if self.tipper_errorfloor is not None:
                                raw_error = np.abs((self.tipper_errorfloor/100.)*value)
                        else:
                            raw_error = tippererr[0,1] 

                        error = raw_error
                        
                        # set error floor
                        if self.tipper_errorfloor is not None:
                            error = max(error,np.abs((self.tipper_errorfloor/100.)*value))
                        
                        # set an absolute minimum tipper error to apply when value is close to zero
                        if self.tipper_errorfloor_abs is not None:
                            error = max(error, self.tipper_errorfloor_abs)
                    
                    if append:
                        self.data.append([station_number,frequency_number,mode,value,np.abs(error)])


    def generate_profile(self):
        """
            Generate linear profile by regression of station locations.

            Stations are projected orthogonally onto the profile. Calculate 
            orientation of profile (azimuth) and position of stations on the 
            profile.

            Sorting along the profile is always West->East.
            (In unlikely/synthetic case of azimuth=0, it's North->South)


            (self.stationlocations, self.azimuth, self.stations)

        """


        self.station_coords = []
        self.stations = []
        self.station_frequencies = []
        
        self.Z = []
        self.Tipper = []

        lo_strike_angles = []

        lo_easts = []
        lo_norths = []
        utmzones = []

        lo_wrong_edifiles = []

        for edifile in self.edilist:
            edi = MTedi.Edi()
            try:
                edi.readfile(edifile,datatype=self.edi_type)
            except:
                lo_wrong_edifiles.append(edifile)
                continue

            if self.strike is None:
                try:
                    lo_strike_angles.extend(list(MTgy.strike_angle(edi.Z.z[np.where(MTgy.dimensionality(edi.Z.z)!=1)])[:,0]%90))
                except:
                    pass
            self.station_coords.append([edi.lat,edi.lon,edi.elev])
            self.stations.append(edi.station)
            self.station_frequencies.append(np.around(edi.freq,5))
            try:
                self.Tipper.append(edi.Tipper)
            except:
                self.Tipper.append(None)
                
            self.Z.append(edi.Z)
            utm = MTcv.LLtoUTM(23,edi.lat,edi.lon)
            lo_easts.append(utm[1])
            lo_norths.append(utm[2])
            utmzones.append(int(utm[0][:-1]))
       

        for i in lo_wrong_edifiles:
            self.edilist.remove(i)

        if len(self.edilist) == 0:
            raise

        if self.strike is None:
            try:
                self.strike = np.mean(lo_strike_angles)
            except:
                #empty list or so....
                #can happen, if everyhing is just 1D
                self.strike = 0.

        main_utmzone = mode(utmzones)[0][0]


        for idx, zone in enumerate(utmzones):
            if zone == main_utmzone:
                continue
            utm = MTcv.LLtoUTM(23,self.station_coords[idx][0],self.station_coords[idx][1],main_utmzone)

            lo_easts[idx] = utm[1]
            lo_norths[idx] = utm[2]
        
        lo_easts = np.array(lo_easts)
        lo_norths = np.array(lo_norths)

        # check regression for 2 profile orientations:
        # horizontal (N=N(E)) or vertical(E=E(N))
        # use the one with the lower standard deviation
        profile1 = sp.stats.linregress(lo_easts, lo_norths)
        profile2 = sp.stats.linregress(lo_norths, lo_easts)
        profile_line = profile1[:2]
        #if the profile is rather E=E(N), the parameters have to converted 
        # into N=N(E) form:
        if profile2[4]<profile1[4]:
            profile_line = (1./profile2[0], -profile2[1]/profile2[0])

        #profile_line = sp.polyfit(lo_easts, lo_norths, 1) 
        self.azimuth = (90-(np.arctan(profile_line[0])*180/np.pi))%180


        
        #rotate Z according to strike angle, 

        #if strike was explicitely given, use that value!

        #otherwise:
        #have 90 degree ambiguity in strike determination
        #choose strike which offers larger angle with profile
        #if profile azimuth is in [0,90].

        if self._strike_set is False:
            if 0 <= self.azimuth < 90:
                if np.abs(self.azimuth - self.strike) < 45:
                    self.strike += 90
            elif 90 <= self.azimuth < 135:
                if self.azimuth - self.strike < 45:
                    self.strike -= 90
            else:
                if self.azimuth - self.strike >= 135:
                    self.strike += 90
         

        self.strike = self.strike%180


        rotation_angle = self.strike
        
        for old_z in self.Z:
            original_rotation_angle = np.array(old_z.rotation_angle)
            effective_rot_angle = rotation_angle - original_rotation_angle
            old_z.rotate(effective_rot_angle)
        
        # rotate tipper to profile azimuth, not strike. Need angle to be between

        rotation_angle = (self.azimuth - 90) % 180
            
        
        for old_tipper in self.Tipper:
            try:
                original_rotation_angle = np.array(old_tipper.rotation_angle)
                effective_rot_angle = rotation_angle - original_rotation_angle
                old_tipper.rotate(effective_rot_angle)
            except:
                pass


        projected_stations = []
        lo_offsets = []
        profile_vector = np.array([1,profile_line[0]])
        profile_vector /= np.linalg.norm(profile_vector)

        for idx,sta in enumerate(self.stations):
            station_vector = np.array([lo_easts[idx],lo_norths[idx]-profile_line[1]])
            position = np.dot(profile_vector,station_vector) * profile_vector 
            lo_offsets.append(np.linalg.norm(position))
            projected_stations.append([position[0],position[1]+profile_line[1]])

        lo_offsets -= min(lo_offsets)


        #Sort from West to East:
        profile_idxs = np.argsort(lo_offsets)
        if self.azimuth == 0:
            #Exception: sort from North to South
            profile_idxs = np.argsort(lo_norths)


        #sorting along the profile
        projected_stations = [projected_stations[i] for i in profile_idxs]
        projected_stations =  np.array(projected_stations)
        lo_offsets = np.array([lo_offsets[i] for i in profile_idxs])
        lo_offsets -= min(lo_offsets)

        self.station_coords = [self.station_coords[i] for i in profile_idxs]
        self.stations = [self.stations[i] for i in profile_idxs]
        self.station_frequencies = [self.station_frequencies[i] for i in profile_idxs]
        self.Z = [self.Z[i] for i in profile_idxs]
        self.Tipper = [self.Tipper[i] for i in profile_idxs]
        lo_easts = np.array([lo_easts[i] for i in profile_idxs])
        lo_norths = np.array([lo_norths[i] for i in profile_idxs])
              

        self.profile = profile_line
        self.stationlocations = lo_offsets
        self.easts = lo_easts
        self.norths = lo_norths
        #print self.stationlocations

        #plot profile and stations:
        if 0:
            lo_all_easts = list(lo_easts)
            lo_all_easts.extend(list(projected_stations[:,0]))
            lo_all_norths = list(lo_norths)
            lo_all_norths.extend(list(projected_stations[:,1]))
            x_extent = max(lo_all_easts) - min(lo_all_easts)
            y_extent = max(lo_all_norths) - min(lo_all_norths)
            plt.close('all')
            lfig = plt.figure(4, dpi=200)#, figsize=(2,2))
            plt.clf()
            ploty = sp.polyval(profile_line, sorted(lo_all_easts))
            lax = lfig.add_subplot(1, 1, 1,aspect='equal')
            lax.plot(sorted(lo_all_easts), ploty, '-k', lw=1)
            lax.scatter(lo_easts,lo_norths,color='b',marker='+')
            lax.scatter(projected_stations[:,0], projected_stations[:,1],color='r',marker='x')
            lax.set_title('Original/Projected Stations')
            lax.set_ylim(np.min([lo_norths.min(),projected_stations[:,1].min()])-0.2*y_extent, 
                                            np.max([lo_norths.max(),projected_stations[:,1].max()])+0.2*y_extent)
            lax.set_xlim(np.min([lo_easts.min(),projected_stations[:,0].min()])-0.2*x_extent, 
                                            np.max([lo_easts.max(),projected_stations[:,0].max()])+0.2*x_extent)
            lax.set_xlabel('Easting (m)', 
                           fontdict={'size':4, 'weight':'bold'})
            lax.set_ylabel('Northing (m)',
                           fontdict={'size':4, 'weight':'bold'})
            plt.show()
            #raw_input()




    def writefile(self, filename = None):

        if filename is not None:
            try:
                fn = op.abspath(op.join(self.wd,filename))
                self.datafile = op.split(fn)[1]
                #self.wd = op.abspath(op.split(fn)[0])
            except:
                self.datafile = op.abspath(op.join(self.wd,'OccamDataFile.dat')) 

        outstring = ''

        outstring += 'FORMAT:'+11*' '+self.format+'\n'
        outstring += 'TITLE:'+12*' '+'{0} - profile azimuth {1:.1f} deg -'\
                    ' strike {2:.1f} deg\n'.format(self.title,self.azimuth,
                                                                 self.strike)
        outstring += 'SITES:'+12*' '+'{0}\n'.format(len(self.stations))
        for s in self.stations:
            outstring += '    {0}\n'.format(s)
        outstring += 'OFFSETS (M):\n'
        for l in self.stationlocations:
            outstring += '    {0}\n'.format(l + self.profile_offset)
        outstring += 'FREQUENCIES:      {0}\n'.format(len(self.frequencies))
        for f in self.frequencies:
            outstring += '    {0}\n'.format(f)
        outstring += 'DATA BLOCKS:      {0}\n'.format(len(self.data))

        outstring += 'SITE    FREQ    TYPE    DATUM    ERROR\n'
        for d in self.data:
            outstring += '{0}    {1}    {2}    {3}    {4}\n'.format(*d)

        outfn = op.abspath(op.join(self.wd,self.datafile))

        F = open(outfn,'w')
        F.write(outstring)
        F.close()

    def get_profile_origin(self):
        """
        get the origin of the profile in real world coordinates
        
        Author: Alison Kirkby (2013)
        """

        x,y = self.easts,self.norths
        x1,y1 = x[0],y[0]
        [m,c1] = self.profile
        x0 = (y1+(1.0/m)*x1-c1)/(m+(1.0/m))
        y0 = m*x0+c1 
        self.profile_origin = [x0,y0]
    


class Model():
    """
    Handling of Occam output files.

    Reading, writing, renaming,... of 'ITER' and 'RESP' files. 
    """
    

class Plot():
    """
    Graphical representations of in- and output data.

    Provide gui for masking points.
    Represent output models.
    """

class Run():
    """
    Run Occam2D by system call.

    Future plan: implement Occam in Python and call it from here directly.
    """


class Mask(Data):
    """
    Allow masking of points from data file (effectively commenting them out, 
    so the process is reversable). Inheriting from Data class.
    """



