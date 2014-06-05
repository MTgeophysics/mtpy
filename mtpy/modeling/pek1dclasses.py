# -*- coding: utf-8 -*-
"""
Created on Thu Jun 05 13:55:13 2014

@author: Alison Kirkby

"""

import mtpy.core.edi as mtedi
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import mtpy.utils.exceptions as MTex


    

class Control():    
    
    def __init__(self, **input_parameters):
        
        self.run_input = [1,0,0.1,40,1.05,1,0]
        # define control file parameters
        self.imax = 100 # max number of iterations
        self.type_struct = 6
        self.type_aniso = 2 # type of structure and anisotropy penalties
        self.value_struct = [0.1,1.0,10.0] # values for the structure penalty weights
        self.value_aniso = [0.1,1.0,10.0] # values for the anisotropy penalty weights   
        self.wd = '.'

        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])
        
        if not os.path.exists(self.wd):
            os.mkdir(self.wd) 
            
            
    def write_ctlfile(self):
        """
        write control file
        """

        
        # create control file
        ctlfile = open(os.path.join(self.wd,'inregulm.dat'),'wb') # control file name is hardcoded into software!
        
        # define number of weights
        nw_struct = len(self.value_struct)
        nw_aniso = len(self.value_aniso)  
        
        for thing in [(2,self.imax),(nw_struct,nw_aniso),(self.type_struct,self.type_aniso)]:
            ctlfile.write('%1i%6i\n'%thing)
        for thing in [self.value_struct,self.value_aniso]:
            ctlfile.write('  '.join([str(i) for i in thing])+'\n')
        ctlfile.close()
        print "written control file to {}".format(self.wd)



class Inmodel():
    

    def __init__(self, inmodel_modeldir, **input_parameters):
        self.wd = '.'
        self.inmodel_modeldir = inmodel_modeldir
        self.inmodel_vals = {0:[100,100,0]} # dictionary containing values for 
                                            # inmodel file, in format topdepth: [minres,maxres,strike]
        self.elevation_file = None
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])


    def build_inmodel(self):
        """
        build an inmodel file to be used as a constraint
        need to give it a dictionary containing values (list of rmin,rmax and strike) and bottom depths
        depths are the keys, resistivities are the values
        and a modeldir - needs to have the same steps as
        the model planned to run.
        """

        modelf = open(os.path.join(self.inmodel_modeldir,'ai1mod.dat'))
        modelf.readline()
        
        flag = True
        model = []        
        
        ii = 1
        while flag:
            try:
                m = [float(i) for i in modelf.readline().strip().split()]
                if len(m) > 0:
                    if ii % 2 == 0:
                        model.append(m)
                    ii += 1
            except:
                flag = False
        
        model = np.array(model)
        model[:,2:] = 0.
        mvals = model[:,2:]
        mi = model[:,0]
        mdepths = [0.]+list(model[:,1])
        mthick = np.array([mdepths[i+1]-mdepths[i] for i in range(len(mi))])
        
        keys = self.inmodel_vals.keys()
        keys.sort()
        for key in keys:
            cond = model[:,1] >= key
            mvals[cond] = np.array(self.inmodel_vals[key])
        
        self.inmodel = np.vstack([mi,mthick,mvals.T]).T
        
    
    def write_inmodel(self, wd=None):
        """
        """
        
        if wd is not None:
            self.wd = wd
        
        if not hasattr(self,'inmodel'):
            self.build_inmodel()
        
        np.savetxt(os.path.join(self.wd,'inmodel.dat'),
                   self.inmodel,
                   fmt=['%5i','%11.4e','%11.4e','%11.4e','%11.4e'])
        print "written inmodel file to {}".format(self.wd)


class Data():
    """
    deals with input data from 1d inversions, including creating a data file
    and reading a data file afterwards to compare with inversion responses
    
    """    
    def __init__(self, **input_parameters):
        self.wd = None
        self.respfile = 'ai1dat.dat'
        self.datafile = None
        self.errorfloor_z = 0.1
        self.errorfloor_type = 'relative' # choose whether to set an absolute or relative error floor
        self.epath = None
        self.mode = 'I'

        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])   

        # default working directory is epath if it is specified, otherwise
        # current directory
        if self.wd is None:
            if self.epath is not None:
                self.wd = os.path.dirname(self.epath)
            else:
                self.wd = '.'

    def build_data(self):
        """
        create data to write to an input file
        
        """
        
        # read edi file to edi object
        eo = mtedi.Edi()
        eo.readfile(self.epath)
        self.edi_object = eo
        
        # define z
        zr = np.real(eo.Z.z)
        # sign of imaginary component needs to be reversed for this particular code
        zi = -np.imag(eo.Z.z)
        ze = eo.Z.zerr
        z = zr + 1j*zi
        if self.errorfloor_type == 'relative':
            zer = ze/np.abs(z)
            zer[(zer<self.errorfloor_z)] = self.errorfloor_z
            ze = np.abs(z)*zer
            
            #  set errors in off-diagonals to the minimum error of on-diagonal components
            for ze_sub in ze:
                min_offdiags = min(ze_sub[0,1],ze_sub[1,0])
                ze_sub[ze_sub<min_offdiags] = min_offdiags   
        
        elif self.errorfloor_type == 'absolute':
            ze[ze<self.errorfloor_z] = self.errorfloor_z
            
        # define header info for data file
        header = '{:>5}\n{:>5}'.format(self.mode,len(eo.Z.resistivity))
        
        # create data array
        data_list = [1./eo.freq]
        for i in range(2):
            for j in range(2):
                if self.mode == 'I':
                    dd = [zr,ze,zi,ze]

                for d in dd:
                    data_list.append(d[:,i,j])
        
        self.header = header
        self.data = np.vstack(data_list).T


    def write_datafile(self, wd = None):
        """
        write data to file
        
        """
        
        if wd is not None:
            self.wd = wd
        
        self.build_data()
        
        # define format list for writing data file
        fmt = ['%14.5f']+['%12.5e']*16
        
        # define file name and save data file
        fname_bas = self.edi_object.station[:5]
        self.datafile = fname_bas+'.dat'
        fname = os.path.join(self.wd,self.datafile)

        np.savetxt(fname,self.data,fmt=fmt,header=self.header,comments='')    

        
    def read_datafile(self):
        """
        read data file into the data object.
        calculate resistivity and phase
        
        """
        if self.datafile is None:
            dlst = [i for i in os.listdir(self.wd) if i[-4:] == '.dat']
            default_files = ['ai1dat.dat','ai1mod.dat','ai1fit.dat',
                             'inmodel.dat','inregulm.dat']
            for dd in default_files:
                try:
                    dlst.remove(dd)
                except:
                    pass
            if len(dlst) == 1:
                self.datafile = dlst[0]
            else:
                print "please define datafile"
                return           

        # define path to file
        datafpath = os.path.join(self.wd,self.datafile)
        self.mode = open(datafpath).readline().strip().split()[0]
        data = np.loadtxt(datafpath,skiprows = 2)
        self.freq = 1./data[:,0]
   
        if self.mode == 'I':   
            zr = np.vstack([data[:,i] for i in range(len(data[0])) if (i-1)%4 == 0])
            ze = np.vstack([data[:,i] for i in range(len(data[0])) if (i-2)%4 == 0])
            zi = np.vstack([data[:,i] for i in range(len(data[0])) if (i-3)%4 == 0])
            z = zr + 1j*zi
            self.z = z.T.reshape(len(z[0]),2,2)
            self.zerr = ze.T.reshape(len(z[0]),2,2)
            zi = np.imag(self.z)
            zr = np.real(self.z)
            
            freq2 = np.zeros(np.shape(self.z))
            for i in range(len(freq2)):
                freq2[i,:,:] = 1./data[:,0][i]             
            
#            self.resistivity = 
            self.resistivity = 0.2*np.abs(self.z)**2/freq2
            self.resistivity_err = (self.zerr/np.abs(self.z))*self.resistivity
            
            q = np.zeros(np.shape(self.resistivity))
            q[(zr<0)&(zi<0)] = np.pi
            q[(zr<0)&(zi>0)] = -np.pi
            self.q = q

            phase = np.rad2deg(np.arctan(np.imag(self.z)/np.real(self.z))+q)
            phase[phase<-180] += 360
            self.phase = phase
            self.phase_err = np.rad2deg(np.arcsin(self.zerr/np.abs(self.z)))
            
        elif self.mode == 'R':
            res = np.vstack([data[:,i] for i in range(len(data[0])) if (i-1)%4 == 0])
            self.resistivity = res.T.reshape(len(res[0]),2,2)
            res_err = np.vstack([data[:,i] for i in range(len(data[0])) if (i-2)%4 == 0])
            self.resistivity_err = res_err.T.reshape(len(res_err[0]),2,2)
            
            phs = np.vstack([data[:,i] for i in range(len(data[0])) if (i-3)%4 == 0])
            self.phase = phs.T.reshape(len(phs[0]),2,2)
            phs_err = np.vstack([data[:,i] for i in range(len(data[0])) if (i-4)%4 == 0])
            self.phase_err = phs_err.T.reshape(len(phs_err[0]),2,2)



class Model():
    """
    deals with outputs from 1d inversions
    """    
    
    def __init__(self,wkdir,**input_parameters):
        
        self.wd = wkdir
        self.modelfile = 'ai1mod.dat'
        self.respfile = 'ai1dat.dat'
        self.fitfile = 'ai1fit.dat'
        self.inmodelfile = 'inmodel.dat'
        self.misfit_threshold = 1.1
        self.station = None
        
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])       
        
        
    def find_bestmodel(self):
        """
        find the smoothest model that fits the data within self.misfit_threshold
        """
        # read fitfile
        self.read_fit()
        fit = self.fit
        
        # define parameters
        mis = self.misfit
        s = self.penalty_structure
        a = self.penalty_anisotropy
        
        # define function to minimise
        f = a*s*np.abs(a-s)/(a+s)
        
        # define the parameters relating to the best model
        self.params_bestmodel = fit[f==min(f[mis<min(mis)*self.misfit_threshold])][0]
        self.params_fittingmodels = fit[mis<min(mis)*self.misfit_threshold]

        
    def read_model(self):
        """
        read all models into an array
        """
        
        fpath = os.path.join(self.wd,self.modelfile)
        
        nlayers = 0
        flag = True
        modelf = open(fpath)
        modelf.readline()
        
        while flag:
            try:
                nlayers = int(modelf.readline().strip().split()[0])
            except:
                flag = False        
        
        models = np.genfromtxt(fpath,skiprows=1,invalid_raise=False)
        self.models = models.reshape(0.5*len(models)/nlayers,2*nlayers,5)
        
    def read_inmodel(self):
        """
        read the inmodel file
        """
        
        # read in file
        inmodel = np.loadtxt(os.path.join(self.wd,self.inmodelfile))
        
        # convert layer thicknesses to depths
        depths = np.array([[sum(inmodel[:i,1]),sum(inmodel[:i+1,1])] for i in range(len(inmodel))]).flatten()
        values = np.zeros((len(inmodel)*2,5))
        ii = 0
        for val in inmodel:
            for i in range(2):
                values[ii] = val
                ii += 1
        
        self.inmodel = np.vstack([values[:,0],depths,values[:,2],values[:,3],values[:,4]]).T
        
        

class Response():
    """
    deals with responses from 1d inversions
    """    
    
    def __init__(self,wkdir,**input_parameters):
        
        self.wd = wkdir
        self.respfile = 'ai1dat.dat'
        self.misfit_threshold = 1.1
        self.station = None
        
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])    


    def read_respfile(self):
        """
        read respfile into a data object
        """
        
        # define path to file
        respfpath = os.path.join(self.wd,self.respfile)
        respf = open(respfpath)

        # find out number of models
        n = 0
        for line in respf.readlines():
            if 'REG' in line:
                n += 1

        # load model responses into an array
        resp = np.genfromtxt(respfpath,skiprows=1,invalid_raise=False)
        resmod = np.vstack([resp[:,i] for i in range(len(resp[0])) if (i-1)%2 == 0])
        phsmod = np.vstack([resp[:,i] for i in range(len(resp[0])) if i!= 0 and (i-2)%2 == 0])
        period = resp[:len(resp)/n,0]

        self.resistivity = resmod.T.reshape(n,len(resp)/n,2,2)
        self.phase = phsmod.T.reshape(n,len(resp)/n,2,2)
        self.freq = 1./period
        zabs = np.zeros((n,len(resp)/n,2,2))
        for m in range(n):
            for f in range(len(self.freq)):
                zabs[m,f] = (self.resistivity[m,f]/(0.2*period[f]))**0.5
        zr = zabs*np.cos(np.deg2rad(self.phase))
        zi = zabs*np.sin(np.deg2rad(self.phase))
        self.z = zr + 1j*zi
               


class Fit():
    """
    deals with outputs from 1d inversions
    """    
    
    def __init__(self,wkdir,**input_parameters):
        
        self.wd = wkdir
        self.fitfile = 'ai1fit.dat'
        self.misfit_threshold = 1.1
        self.station = None
        
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])    
        
        
    def find_nperiods(self):
        """
        find number of periods used in inversion
        """
        
        # find out number of periods
        respfpath = os.path.join(self.wd,self.respfile)
        respf = open(respfpath)
        respf.readline()
        
        n = 0
        line = respf.readline()
        while 'REG' not in line:
            line = respf.readline()
            n += 1
        self.n_periods = n-1
        

    def read_fit(self):
        """
        read fit file to give structure and anisotropy penalties and penalty weights
        """
        # load the file with fit values in it
        fit = np.loadtxt(os.path.join(self.wd,self.fitfile))
        
        # find number of periods
        self.find_nperiods()        
        
        # total misfit
        self.misfit = (fit[:,5]/(self.n_periods*8.))**0.5
        
        # structure and anisotropy penalty
        self.penalty_structure = fit[:,6]
        self.penalty_anisotropy = fit[:,7]
        self.weight_structure = fit[:,2]
        self.weight_anisotropy = fit[:,4]
        self.modelno = fit[:,0]
        self.fit = fit
