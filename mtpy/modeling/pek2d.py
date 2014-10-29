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
        
        self.parameters_ctl = {}
        self.parameters_ctl['ctl_string'] = 'TAB'
        self.parameters_ctl['units_string'] = 'PR'
        self.parameters_ctl['quadrants'] = '++--'
        self.parameters_ctl['orientation_string'] = '0  0.d0  0.d0'
        self.parameters_ctl['convergence_string'] = '1  6  1.d-4'
        self.parameters_ctl['roughness_string'] = '2  1000.0d0  1000.0d0  0.d0'
        self.parameters_ctl['anisotropy_penalty_string'] = '2  1000.d0  0.d0'
        self.parameters_ctl['anisotropy_ctl_string'] = '1.d0  1.d0  1.d0'


        self.parameters_model = {}
        self.parameters_model['no_sideblockelements'] = 5
        self.parameters_model['no_bottomlayerelements'] = 4
        self.parameters_model['firstlayer_thickness'] = 100
        
        #model depth is in km!
        self.parameters_model['model_depth'] = 100
        self.parameters_model['no_layers'] = 25
        self.parameters_model['max_blockwidth'] = 1000
        self.parameters_data = {}
        self.parameters_data['strike'] = 0.
        self.parameters_data['errorfloor'] = dict(z=np.array([[0.05,0.05],
                                                              [0.05,0.05]]),
                                                  tipper=np.array([0.01,0.01]))
        self.parameters_data['max_no_frequencies'] = 80
        self.parameters_data['mode'] = [1,1,1,1,1,1]
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
        self.resfile = 'pb.res'
        self.cvgfile = 'pb.cvg'
        self.outfile = 'pb.out'
        self.pexfile = 'pb.pex'
        self.andfile = 'pb.and'
        self.exlfile = 'pb.exl'        

        update_dict = {}

        #correcting dictionary for upper case keys
        input_parameters_nocase = {}
        for key in input_parameters.keys():
            input_parameters_nocase[key.lower()] = input_parameters[key]

        update_dict.update(input_parameters_nocase)

        for dictionary in [self.parameters_model,self.parameters_data]:
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
           
    def build_inputfiles(self):
        self.inversiondir = fh.make_unique_folder(self.working_directory,basename='aniso')
        os.mkdir(op.join(self.working_directory,self.inversiondir))
        self.build_model()
        self.write_modelfile()
        self.write_datafiles()
        self.write_ctlfile()
                
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
        
        
    def write_modelfile(self):
        
        if not hasattr(self,'modelfilestring'):
            self.build_model()
        outfile = open(op.join(self.working_directory,
                               self.inversiondir,
                               self.modelfile),'w')
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
                                                             self.n_airlayers+1]]))

        # add strings giving horizontal and vertical mesh steps
        meshz = list(self.meshblockthicknesses_zair)+list(self.meshblockthicknesses_z)
        for meshstep in [self.meshblockwidths_x,meshz]:
            modelfilestring.append\
            (p2d.create_multiple_line_string(meshstep,
                                             10,'%10.3f'))

        
        # add resistivity map
        rmap = ('%5i'%0*len(self.resistivity[0])+'\n')*self.n_airlayers
        rmap += '\n'.join([''.join('%5i'%ii for ii in i) for i in \
        np.arange(1,np.size(self.resistivity[:,:,0])+1).reshape(np.shape(self.resistivity)[:2])])
        modelfilestring.append(rmap)
        
        # add number of resistivity domains (+1 to include air)
        modelfilestring.append('%5i'%(np.size(self.resistivity[:,:,0])+1))
        
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
        self.modelfilestring = '\n'.join(modelfilestring)+'\n'


    def build_data(self):

        imethod = 'nearest'
        ftol = 0.000001
        
        num_freq = int(self.parameters_data['max_no_frequencies'])

        # get minimum and maximum periods
        min_val = max([min(1./zo.freq) for zo in self.Data.Z])
        max_val = min([max(1./zo.freq) for zo in self.Data.Z])
        
        
        periodlst = []
        
        for period in 1./self.Data.frequencies:
            if len(periodlst) > 0:
                # find the difference between the period and the closest period already in the list
                closest_period_diff = np.amin(np.abs(np.array(periodlst)-period))
            else:
                # otherwise set period difference to a large number
                closest_period_diff = 99999

            # check whether the fractional difference is bigger than the tolerance set
#            print closest_period_diff,closest_period_diff/period,
            if closest_period_diff/period > ftol:
                if min_val <= period <= max_val:
                    periodlst.append(period)
        periodlst.sort()
#        print periodlst
        # if number of periods still too long based on the number of frequencies set
        # then take out some frequencies
        n = 1
        while len(periodlst) > num_freq:
            to_remove = int(float(len(periodlst))/num_freq*n)-n
            periodlst = periodlst[:to_remove]+periodlst[to_remove+1:]
            n += 1
#        print periodlst

        mode = self.parameters_data['mode']
        if type(mode) in [str]:
            mode = mode.split(',')
            self.parameters_data['mode'] = mode
        datafile_data = {}
        
        for ee,zo in enumerate(self.Data.Z):
            to=self.Data.Tipper[ee]
            datfn = str(self.stationblocknums[ee])+'_'+self.Data.stations[ee]+'.dat'


            zerr = zo.zerr
            z = zo.z
            ze_rel = zerr/np.abs(z)
            
            terr = to.tippererr
            t = to.tipper
            te_rel = terr/np.abs(t)
            
            # set error floors
            efz = self.parameters_data['errorfloor']['z']
            eft = self.parameters_data['errorfloor']['tipper']
            for i in range(2):
                for j in range(2):
                    ze_rel[ze_rel<efz[i,j]] = efz[i,j]
                te_rel[te_rel<eft[i]] = eft[i]
            zerr = ze_rel * np.abs(z)
            terr = te_rel * np.abs(t)
            terr[terr<0.02]=0.02
            zvar = zerr**2
      
            # create interpolation functions to interpolate z and tipper values
            properties = dict(z_real=np.real(z),z_imag=np.imag(z),
                              z_var=zvar,tipper_real=np.real(t),
                              tipper_imag=np.imag(t),tipper_err=terr)
            properties_interp = {}
            for key in properties.keys():
                f = si.interp1d(np.log10(1./zo.freq),properties[key],
                                axis = 0,kind = imethod)
                properties_interp[key] = f(np.log10(periodlst))

            datafile_data[datfn] = properties_interp
            
            
        self.datafile_data = datafile_data
        self.freq = 1./(np.array(periodlst))
            
    def build_datafiles(self):
        
        if not hasattr(self,'datafile_data'):
            self.build_data()
        dfstrings = {}
        
        for dfile in self.datafile_data.keys():
            datfstr = '{:<3} '.format(len(self.freq))+\
            ' '.join([str(i) for i in self.parameters_data['mode']])+'\n'
            for pv in range(len(self.freq)):
                datlst = '{0:>12}'.format('%.06f'%(1./(self.freq[pv])))
                for ii in range(2):
                    for jj in range(2):
                        for pval in [k for k in self.datafile_data[dfile].keys() if 'z' in k]:
#                            print self.datafile_data[dfile][pval][pv][ii,jj]
                            datlst += '{0:>12}'.format('%.06f'%self.datafile_data[dfile][pval][pv][ii,jj])
                for ii in range(2):
                    for pval in [k for k in self.datafile_data[dfile].keys() if 't' in k]:
                        datlst += '{0:>12}'.format('%.06f'%self.datafile_data[dfile][pval][pv][0,ii])

                datfstr += ''.join(datlst)+'\n'
            dfstrings[dfile] = datfstr
        
        self.datafile_strings = dfstrings

    def write_datafiles(self):
        
        if not hasattr(self,'datafile_strings'):
            self.build_datafiles()
            
        exlf = open(os.path.join(self.working_directory,self.inversiondir,self.exlfile),'w')
        dfkeys=self.datafile_strings.keys()
        dfkeys.sort()
        for dfile in dfkeys:
            f = open(op.join(self.working_directory,self.inversiondir,dfile),'w')
            f.write(self.datafile_strings[dfile])
            f.close()
            exlf.write(dfile+'\n')
        exlf.close()

    def write_ctlfile(self):
        ctrf = open(op.join(self.working_directory,self.inversiondir,'pb.ctr'),'w')
        
        if type(self.parameters_data['mode']) == str:
            self.parameters_data['mode'] = self.parameters_data['mode'].split(',')
        ef = np.hstack([self.parameters_data['errorfloor'][l].flatten() for l in ['z','tipper']])

        
        clist = []
        clist.append(self.exlfile)
        clist.append(self.parameters_ctl['ctl_string']+self.parameters_ctl['units_string']+self.parameters_ctl['quadrants'])
        clist.append(' '.join([str(i) for i in self.parameters_data['mode']]))
        clist.append(' '.join([str(i) for i in ef]))
        clist.append(self.parameters_ctl['orientation_string'])
        clist.append(self.modelfile)
        clist.append(self.resfile)
        clist.append(self.parameters_ctl['convergence_string'])
        clist.append(self.parameters_ctl['roughness_string'])
        clist.append(self.parameters_ctl['anisotropy_penalty_string'])
        clist.append(self.parameters_ctl['anisotropy_ctl_string'])
        clist.append(self.cvgfile)
        clist.append(self.outfile)
        clist.append(self.pexfile)
        clist.append(self.andfile)       
        
        self.controlfile_string = '\n'.join(clist)
        ctrf.write(self.controlfile_string)

        ctrf.close()