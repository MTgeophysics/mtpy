# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 10:21:05 2014

@author: a1655681
"""

import mtpy.core.edi as e
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si


class Setup():
    """
    Dealing with the setup  for a 1d anisotropic inversion run. Generate model,
    and control files.

    Setting up those files within one (pre-determined) folder, so code can be 
    run there straight away.

    """
    def __init__(self,wkdir,epath,efile, **input_parameters):
        
        self.run_input = [1,0,0.1,40,1.05,1,0]
        # define control file parameters
        self.imax = 100 # max number of iterations
        self.type_struct = 6
        self.type_aniso = 2 # type of structure and anisotropy penalties
        self.value_struct = [0.1,1.0,10.0] # values for the structure penalty weights
        self.value_aniso = [0.1,1.0,10.0] # values for the anisotropy penalty weights   
        self.errorfloor_z = 0.1
        self.errorfloor_type = 'relative' # choose whether to set an absolute or relative error floor
        self.wd = wkdir
        self.epath = epath
        self.efile = efile
        self.mode = 'I'
        self.strike_fmax = 0.01 # maximum frequency to calculate strike
        self.lat = 0.0
        self.lon = 0.0
        self.inmodel_vals = {0:[100,100,0]} # dictionary containing values for 
                                            # inmodel file, in format topdepth: [minres,maxres,strike]
        
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])
        
        if not os.path.exists(self.wd):
            os.mkdir(self.wd)
    
    def generate_inputfiles(self,modeldir = None):

        os.chdir(self.wd)    
        self.write_datafile()
        self.write_ctlfile()
        if self.run_input[-1] == 1:
            if modeldir is not None:
                self.write_inmodel(self.inmodel_vals,modeldir)
        
        
    
    def build_data(self):
        
        
        # read edi file to edi object
        edifile = os.path.join(self.epath,self.efile)
        eo = e.Edi()
        eo.readfile(edifile)
        self.lat = eo.lat
        self.lon = eo.lon

        # define z

        zr = np.real(eo.Z.z)
        zi = -np.imag(eo.Z.z)
        ze = eo.Z.zerr
        z = zr + 1j*zi
        if self.errorfloor_type == 'relative':
            zer = ze/np.abs(z)
            zer[(zer<self.errorfloor_z)] = self.errorfloor_z
            ze = np.abs(z)*zer
            
            # set errors in off-diagonals to the minimum error of on-diagonal components
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
    
    
    def make_savepath(self):
        """
        make a savepath that doesn't exist already.
        """        
        
        wkdir = self.wd
        
        # define savepath. need to choose a name that doesn't already exist
        i = 1
        svpath_str = self.efile.split('_')[0]+self.mode
        svpath = svpath_str+'_%02i'%i
        while os.path.exists(os.path.join(wkdir,svpath)):
            i += 1
            svpath = svpath_str+'_%02i'%i
            
        self.savepath = os.path.join(wkdir,svpath)
            
        # make the save path and move into savepath
        os.mkdir(os.path.join(wkdir,svpath))
        os.chdir(os.path.join(wkdir,svpath))
        
    def write_datafile(self):
        """
        write datafile
        """
        
        wkdir = self.wd
        
        self.build_data()
        self.make_savepath()
        svpath = self.savepath
        
        # define format list for writing data file
        fmt = ['%14.5f']+['%12.5e']*16
        
        # define file name and save data file
        fname_bas = os.path.basename(svpath)[:8]
        fname = os.path.join(wkdir,svpath,fname_bas+'.dat')
        self.datafile = os.path.basename(fname)
        np.savetxt(fname,self.data,fmt=fmt,header=self.header,comments='')    
    
    def write_ctlfile(self):
        """
        write control file
        """
        
        # create control file
        ctlfile = open(os.path.join(self.wd,self.savepath,'inregulm.dat'),'wb') # control file name is hardcoded into software!
        
        # define number of weights
        nw_struct = len(self.value_struct)
        nw_aniso = len(self.value_aniso)  
        
        for thing in [(2,self.imax),(nw_struct,nw_aniso),(self.type_struct,self.type_aniso)]:
            ctlfile.write('%1i%6i\n'%thing)
        for thing in [self.value_struct,self.value_aniso]:
            ctlfile.write('  '.join([str(i) for i in thing])+'\n')
        ctlfile.close()

    def build_inmodel(self,vals,modeldir):
        """
        build an inmodel file to be used as a constraint
        need to give it a dictionary containing values (list of rmin,rmax and strike) and bottom depths
        depths are the keys, resistivities are the values
        and a modeldir - needs to have the same steps as
        the model planned to run.
        """
        
        modelf = open(os.path.join(modeldir,'ai1mod.dat'))
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
        
        keys = vals.keys()
        keys.sort()
        for key in keys:
            cond = model[:,1] >= key
            mvals[cond] = np.array(vals[key])
        
        self.inmodel = np.vstack([mi,mthick,mvals.T]).T
    
    
    
    def write_inmodel(self,vals,modeldir):
        """
        """
        
        if not hasattr(self,'inmodel'):
            self.build_inmodel(vals,modeldir)
        
        np.savetxt(os.path.join(self.wd,self.savepath,'inmodel.dat'),self.inmodel,fmt=['%5i','%11.4e','%11.4e','%11.4e','%11.4e'])


class Model():
    """
    deals with outputs from 1d inversions
    """    
    
    def __init__(self,wkdir,**input_parameters):
        self.wd = wkdir
        self.modelfile = 'ai1mod.dat'
        self.respfile = 'ai1dat.dat'
        self.fitfile = 'ai1fit.dat'
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
        f = a*s + (a/s)**2 + (s/a)**2
        
        # define the parameters relating to the best model
        self.params_bestmodel = fit[f==min(f[mis<min(mis)*self.misfit_threshold])][0]
        self.params_fittingmodels = fit[mis<min(mis)*self.misfit_threshold]


    def read_fit(self):
        """
        read fit file to give structure and anisotropy penalties and penalty weights
        """
        # load the file with fit values in it
        fit = np.loadtxt(os.path.join(self.wd,self.fitfile))
        
        # find number of periods
        self.find_nperiods()        
        
        # total misfit
        self.misfit = fit[:,5]/(self.n_periods*8)
        

            
        # structure and anisotropy penalty
        self.penalty_structure = fit[:,6]
        self.penalty_anisotropy = fit[:,7]
        self.weight_structure = fit[:,2]
        self.weight_anisotropy = fit[:,4]
        self.modelno = fit[:,0]
        self.fit = fit
        
        

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
        self.n_periods = n
        
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
        
        
    def plot_model(self,
                   modelno,
                   parameter = None,
                   save = True,
                   c = 'k',
                   smm = (0,180),
                   rmm = (0.1,1000),
                   amm = (0,50),
                   ymm = (10,0),
                   ylabel = True):
        """
        plot the model
        """
        
        # if model file hasn't been read, read it 
        if not hasattr(self,'models'):
            self.read_model()
            
            
        # if modelno is not defined, find the best one according to default parameters
        if modelno is None:
            # find best model
            if not hasattr(self,'params_bestmodel'):
                self.find_bestmodel()
            modelno = self.params_bestmodel[0]
            
        model = self.models[modelno]            
            
        allowed_params = ['tetm_a','a','tetm','strike']#'te','tm',
        
        if parameter not in allowed_params:
            parameter = 'tetm_a'
        
        if self.station is None:
            self.station = ''
              
        while np.any(model[:,4]>180):
            model[:,4][model[:,4]>180] -= 180
        while np.any(model[:,4]<0):
            model[:,4][model[:,4]<0] += 180
        
        
        if 'tetm' in parameter:
            if parameter != 'tetm':
                ax1 = plt.subplot(121)
            plt.plot(model[:,2],model[:,1],color=c,label = 'model{}'.format(modelno))
            plt.plot(model[:,3],model[:,1],'--',color=c,label = 'model{}'.format(modelno))
            plt.title('tetm')
            plt.xscale('log')

#            plt.legend()
            plt.title(self.station)
            if ylabel:
                plt.ylabel('depth')
            plt.ylim(ymm)
            plt.xlim(rmm)
        if 'a' in parameter:
            if parameter != 'a':
                ax2 = plt.subplot(122)
            plt.plot(model[:,3]/model[:,2],model[:,1],c=c,label = 'model{}'.format(modelno))
            plt.xlabel('aniso')
            plt.ylim(ymm)
            plt.xlim(amm)
            if ylabel:
                plt.ylabel('depth')
            plt.title(self.station)
        if parameter == 'strike':
            plt.plot(model[:,4],model[:,1],c=c,label = 'model{}'.format(modelno))
            plt.xlabel('strike')
            plt.ylim(ymm)
            plt.xlim(smm)
            if ylabel:
                plt.ylabel('depth')
            plt.title(self.station)
            plt.grid()
        if save:
            plt.savefig(os.path.join(self.wd,parameter+'.png'))
            plt.close()   
            
    def plot_lcurve_contourmap(self, 
                               levels = None,
                               imethod = 'linear',
                               symbol = 'o',
                               fontsize = 8,
                               draw_threshold = None,
                               lim = 100,
                               cmap = 'rainbow'
                               ):
        """
        plot 'lcurve' contour map
        N  = number of levels to contour
        imethod = type of interpolation for numpy.griddata
                  options are 'nearest','linear', and 'cubic'
        symbol = symbol for plotting
        fontsize = numeric value
        draw_threshold = draw a black contour at a given percentage of the minimum rms error
                         e.g. 110 will draw a contour at 1.1 * minimum error value achieved
        
        """
        
        self.read_fit()
        a = self.penalty_anisotropy
        s = self.penalty_structure
        mis = self.misfit
        aa = [str(round(i,1)) for i in self.weight_anisotropy]
        ss = [str(round(i,1)) for i in self.weight_structure]
        
        # grid the data to make contour plot
        # first, define points to grid
        points = np.vstack([s,a]).T
        # define points xi to grid onto
        xi = np.array(np.meshgrid(np.linspace(0.,max(s)),np.linspace(0.,max(a)))).T
#        print xi
        f1 = si.griddata(points,mis,xi,method = imethod)
        cmap = plt.get_cmap(cmap)
        if levels is not None:
            plt.contour(xi[:,:,0],xi[:,:,1],f1,cmap=cmap,levels=levels)
        else:
            plt.contour(xi[:,:,0],xi[:,:,1],f1,cmap=cmap)
        plt.colorbar()
        if draw_threshold is not None:
            plt.contour(xi[:,:,0],xi[:,:,1],f1,
                        levels=[round(min(mis)*draw_threshold/100,2)],
                                linewidths=[1.5],
                                colors='k')
        plt.gca().set_aspect('equal')
        plt.scatter(s,a,c='k',marker=symbol,lw=0)
        
        for i in range(len(aa)):
            if (s[i]<lim) and (a[i]<lim):
                plt.text(s[i],a[i],','.join([ss[i],aa[i]]),fontsize=fontsize)
        
        plt.xlim(0,lim)
        plt.ylim(0,lim)
        plt.xlabel('structure penalty')
        plt.ylabel('anisotropy penalty')            
            
class Data():
    """
    deals with responses from 1d inversions
    """    
    def __init__(self, wkdir, **input_parameters):
        self.wd = wkdir
        self.respfile = 'ai1dat.dat'
        self.datafile = None
        
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

        self.resistivity_mod = resmod.T.reshape(n,len(resp)/n,2,2)
        self.phase_mod = phsmod.T.reshape(n,len(resp)/n,2,2)
        self.freq = 1./period
        zabs = np.zeros((n,len(resp)/n,2,2))
        for m in range(n):
            for f in range(len(self.freq)):
                zabs[m,f] = (self.resistivity_mod[m,f]/(0.2*period[f]))**0.5
        zr = zabs*np.cos(np.deg2rad(self.phase_mod))
        zi = zabs*np.sin(np.deg2rad(self.phase_mod))
        self.z_mod = zr + 1j*zi
            
      
        
    def read_datafile(self,datafile):
        """
        read data file into the data object.
        calculate resistivity and phase
        """
        self.datafile = datafile        
        
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
        

    def plot_responses(self,modelno=45,adjust_phase = True):
        """
        
        """
        
        if not hasattr(self,'resistivity'):
            self.read_datafile()
            
        T = 1./self.freq
        r = self.resistivity
        re = self.resistivity_err
        p = self.phase
        pe = self.phase_err
        rmod = self.resistivity_mod[modelno]
        pmod = self.phase_mod[modelno]
        
        if adjust_phase:
            p[p<0] += 180
            pmod[pmod<0] += 180
        
        ii = 1
        
        for mode,err,model in [[r,re,rmod],[p,pe,pmod]]:
            for sequence in [[[0,1],[1,0]],[[0,0],[1,1]]]:
                ax = plt.subplot(2,2,ii)
                ax.set_color_cycle(['b','r'])
                for i,j in sequence:
                    plt.errorbar(T,mode[:,i,j],err[:,i,j],fmt='--')
                    plt.plot(T,model[:,i,j],'k--')
#            for i,j in [[0,0],[1,1]]:
#                plt.errorbar(T,mode[:,i,j],err[:,i,j],fmt='.')
#                plt.plot(T,model[:,i,j],'k.')
                    plt.xscale('log')
                if ii <= 2:
                    plt.yscale('log')
                plt.grid()
                ii += 1


def sort_folder_list(wkdir,order_file,indices=[0,9999]):
    """
    sort subfolders in wkdir according to order in order_file
    
    wkdir = working directory containing subfolders
    order = full path to text file containing order.
            needs to contain a string to search on that is the same length
            for each item in the list
    indices = indices to search on; default take the whole string
    
    """
    order = open(os.path.join(order_file)).readlines()
    plst = []
    flst = os.listdir(wkdir)
    for o in order:
        for f in flst:
            if str.lower(o.strip().split('_')[0][indices[0]:indices[1]]) == str.lower(f)[indices[0]:indices[1]]:
                plst.append(f)
    return plst
       
def get_elevation(x,y,elevfn,skiprows = 1):
    """
    get elevation at an arbitrary point, interpolated from an xyz file containing elevations
    x = x location of point, can be a numpy array to return multiple z values
    y = y location of point, can be a numpy array to return multiple z values
    elevfn = full path to elevation filename
    """
    elev = np.loadtxt(elevfn)
    f = si.interp2d(elev[:,0],elev[:,1],elev[:,2])
    return f(x,y)
    
    
