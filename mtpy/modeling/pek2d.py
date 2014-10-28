# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 14:06:45 2014

@author: a1655681
"""
import mtpy.modeling.occam2d as o2d
import mtpy.modeling.pek1dclasses as p1dc
import numpy as np
import os
import os.path as op
import pek2dforward as p2d
import string
import scipy.interpolate as si
import mtpy.utils.filehandling as fh
import mtpy.core.edi as mtedi
import mtpy.modeling.pek2dforward as p2d



class Model():
    """
    class for creating and reading model files
    
    """
   
    
    def __init__(self, working_directory, **input_parameters):
        self.working_directory = working_directory
        self.edi_directory = None
        self.occam_configfile = None
        self.parameters_model = {}
        self.parameters_model['no_sideblockelements'] = 5
        self.parameters_model['no_bottomlayerelements'] = 4
        self.parameters_model['firstlayer_thickness'] = 100
        
        #model depth is in km!
        self.parameters_model['model_depth'] = 100
        self.parameters_model['no_layers'] = 25
        self.parameters_model['max_blockwidth'] = 1000
        self.n_airlayers = 5

        self.mesh = None
        self.meshlocations_x = None
        self.meshlocations_z = None
        self.meshblockwidths_x = None
        self.meshblockthicknesses_z = None
        self.profile_easts = None
        self.profile_norths = None
        self.inversion1d_dirdict = {}
        self.inversion1d_masterdir = '.'
        self.inversion1d_modelno = 0
        self.inversion1d_imethod = 'nearest'
        self.binsize_resistivitylog10 = 1.
        self.binsize_strike = 20.
        self.build_from_1d = True
        self.rotation = 0.
        self.modelfile = 'model.dat'
        self.anisotropy_min_depth = 0.
        self.strike = 0.

        self.edifiles = []

        self.Data = None

        self.modelfile = 'model'
        

        update_dict = {}

        #correcting dictionary for upper case keys
        input_parameters_nocase = {}
        for key in input_parameters.keys():
            input_parameters_nocase[key.lower()] = input_parameters[key]

        update_dict.update(input_parameters_nocase)

        for dictionary in [self.parameters_model]:
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
                if update_dict[key] is not None:
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

        if self.edifiles == []:
            if self.edi_directory is not None:
                try:
                    self.edifiles = os.listdir(self.edi_directory)
                except IOError:
                    print("failed to find edi directory")
                    pass
                
                
    def build_model(self):
        """
        build model file string
        """
        # build a forward model object
        ro = p2d.Model(self.working_directory,**self.input_parameters)
        ro.build_mesh()
        ro.build_aircells()


        # assign relavent parameters to pek 2d inverse object
        for at in ['stationlocations','parameters_model',
                   'meshlocations_x','meshlocations_z',
                   'meshblockwidths_x','meshblockthicknesses_z',
                   'profile_easts','profile_norths','Data',
                   'meshblockthicknesses_zair','meshlocations_zair']:
                       attvalue = getattr(ro,at)
                       setattr(self,at,attvalue)


        ro.get_station_meshblock_numbers()
        if ro.build_from_1d:
            ro.get_1d_results()
            ro.interpolate_1d_results()
            for at in ['inversion1d_dirdict','inversion1d_modelno',
                       'models1d','resistivity','stationlocations',
                       'blockcentres_x','blockcentres_z']:
                       attvalue = getattr(ro,at)
                       setattr(self,at,attvalue)                           
        else:
            if not hasattr(ro,'resistivity'):
                print "Cannot build model, please provide resistivity or a 1d"+\
                " model to build resistivity from first"

        ro.get_station_meshblock_numbers()
        self.stationblocknums=ro.stationblocknums
        
        self.build_modelfilestring()
        self.write_modelfile()
        
        
    def write_modelfile(self):
        outfile = open(op.join(self.working_directory,self.modelfile),'w')
        outfile.write(self.modelfilestring)
        outfile.close()

        
       
    def build_modelfilestring(self):    
        # initialise a list containing info for model file
        modelfilestring = []
        # add header info
        modelfilestring.append('NEW')
        modelfilestring.append('    1')
        modelfilestring.append('     1.000')     

        # add string giving number of cells:
        modelfilestring.append(''.join(['%5i'%i for i in [len(self.meshlocations_x),
                                                         len(self.meshlocations_z)+self.n_airlayers,
                                                             self.n_airlayers]]))

        # add strings giving horizontal and vertical mesh steps
        meshz = list(self.meshblockthicknesses_zair)+list(self.meshblockthicknesses_z)
        for meshstep in [self.meshblockwidths_x,meshz]:
            modelfilestring.append\
            (p2d.create_multiple_line_string(meshstep,
                                             10,'%10.3f'))

        
        # add resistivity map
        rmap = ('%5i'%0*len(self.resistivity[0])+'\n')*self.n_airlayers
        rmap += '\n'.join([''.join('%5i'%ii for ii in i) for i in \
        np.arange(np.size(self.resistivity[:,:,0])).reshape(np.shape(self.resistivity)[:2])])
        modelfilestring.append(rmap)
        
        # add number of resistivity domains (+1 to include air)
        modelfilestring.append('%5i'%(np.size(self.resistivity[:,:,0])))
        
        # add dictionary contents, assuming rvertical = rmax, slant and dip zero
        # first, air layer, properties always the same
        modelfilestring.append('    0   0     -1.00      0.00      0.00      0.00      0.00      0.00 0 0 0 0 0 0')
        # second, dictionary contents
        no = 1
        for j in range(len(self.resistivity)):
            for i in range(len(self.resistivity[j])):
                # initialise a list containing resx,resy,strike
                rlist = list(self.resistivity[j,i])
                # insert resz (assumed to be same as resy)
                rlist.insert(2,rlist[1])
                # insert dip and slant (assumed zero)
                rlist += [0.,0.]
                if rlist[1]/rlist[0] == 1.:
                    aniso = '   0'
                    invert_key = ' 1 1 1 0 0 0'
                else:
                    aniso = '   1'
                    invert_key = ' 1 1 1 1 1 0'
                modelfilestring.append(''.join(['%5i'%no,aniso]+['%10.2f'%i for i in rlist]+[invert_key]))
                no += 1
        # append bathymetry index, at this stage only 0 allowed:
        modelfilestring.append('%5i'%0)
        
        # append number of calculation points (stations):
        modelfilestring.append('%5i'%len(self.stationblocknums))
        
        # append rotation
        modelfilestring.append('%10.2f'%self.rotation)
    
        # append station blocknums
        modelfilestring.append(p2d.create_multiple_line_string(self.stationblocknums,
                                                               5,'  %03i'))
                                                               
        modelfilestring.append('%5i'%0)
        self.modelfilestring = '\n'.join(modelfilestring)


    def build_datafiles(self):

        imethod = 'nearest'        
        
        self.read_edifiles()
        if self.eos is None:
            self.get_edis()
        
        eos = self.eos
        eos_sorted = []
        num_freq = int(self.parameters_data['max_no_frequencies'])
        
        strike = self.strike
     
        nx = len(self.meshblockwidths_x)

        for i in range(nx):
             for jj, j in enumerate(self.stationlocations):
                 if self.meshlocations_x[i]<=j<self.meshlocations_x[i+1]:
                     self.stn_blocknums[i+1] = self.stations[jj]        
        

            
        min_val = max([min([np.log10(i) for i in 1./eo.freq]) for eo in self.eos])
        max_val = min([max([np.log10(i) for i in 1./eo.freq]) for eo in self.eos])
        
        
        
        periodlst = []
        
        for eo in self.eos:
            for period in 1./eo.freq:
                if len(periodlst) > 0:
                    closest_period_diff = np.amin(np.abs(np.array(periodlst)-period))
                else:
                    closest_period_diff = 99999
#                print closest_period_diff,period
                if closest_period_diff/period > ftol:
                    if 10**min_val <= period <= 10**max_val:
                        periodlst.append(period)
        periodlst.sort()


        if len(periodlst) > num_freq:
            sub_factor = int(np.ceil(float(len(periodlst))/num_freq))
            print "sub_factor",sub_factor
            periodlst = [periodlst[i] for i in range(0,len(periodlst),sub_factor)]


        mode = self.parameters_data['mode']
        if type(mode) in [str]:
            mode = mode.split(',')
            self.parameters_data['mode'] = mode
        datafiles = {}
        
        for ee,eo in enumerate(eos_sorted):
            eo.rotate(round(strike))
            datfn = str(stn_blocknums_keys[ee])+'_'+eo.station+'.dat'

            
            datfstr = ''
            

            zerr = eo.Z.zerr
            z = eo.Z.z
            ze_rel = zerr/np.abs(z)
            
            # set error floors
            ef_arr = np.array([np.float(i) for i in self.parameters_ctl['errorfloors']])
            if max(ef_arr) > 0.:
                ef = np.amin(ef_arr[ef_arr>0.])
            else:
                ef = 0.
            ze_rel[ze_rel<ef] = ef
            zerr = ze_rel * np.abs(z)
            zvar = zerr**2
            
            # set diagonals errors to minimum of off current value and off diagonal value.
            # because off diags are often very small values            
            for ze_sub in zerr:
                min_offdiags = min(ze_sub[0,1],ze_sub[1,0])
                ze_sub[ze_sub<min_offdiags] = min_offdiags
                
                
            f_pr = si.interp1d(np.log10(1./eo.freq),np.real(z),axis = 0,kind = imethod)
            f_pi = si.interp1d(np.log10(1./eo.freq),np.imag(z),axis = 0,kind = imethod)
            f_err = si.interp1d(np.log10(1./eo.freq),zvar,axis = 0,kind = imethod)
            
            self.freq = 1./(np.array(periodlst))


            for pv in periodlst:
                datlst = '{0:>12}'.format('%.06f'%(pv))
                for ii in range(4):
                    datlst += '{0:>12}'.format('%.06f'%f_pr(np.log10(pv)).flatten()[ii])
                    datlst += '{0:>12}'.format('%.06f'%f_pi(np.log10(pv)).flatten()[ii])
                    datlst += '{0:>12}'.format('%.06f'%f_err(np.log10(pv)).flatten()[ii])
                datlst += 4*'{0:>12}'.format('%.06f'%0.0)
                datfstr += ''.join(datlst)+'\n'
            
            datafiles[datfn] = datfstr
            
        self.datafile_strings = datafiles

    def get_edis(self):
        """
        get edi files as edi objects
        """

        self.read_edifiles()
        
        self.eos = [mtedi.Edi(filename = os.path.join(self.edi_directory,edi)) for edi in self.edifiles]
        self.stations = [eo.station for eo in self.eos]            