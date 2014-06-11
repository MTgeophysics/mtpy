# -*- coding: utf-8 -*-
"""
Created on Thu Jun 05 14:59:04 2014

@author: a1655681
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

class Plot_model():
    """
    plot model results. Inherits a mtpy.modeling.pek1dclasses.Model object
    
    
    parameter = list of parameters to plot. allowed values are
                ['minmax','aniso','strike']
    xlim = dictionary containing limits for parameters
                if not given then defaults are used.
                keys are the parameter values.
                format {parameter:[xmin,xmax]}
    ylim = list containing depth limits, default [6,0]
    titles = dictionary containing titles for plots. Keys are the
             values in parameter (e.g. 'minmax')
    
    """   
    def __init__(self, Model, **input_parameters):

        self.Model = Model
        self.wd = self.Model.wd
        self.parameters = ['minmax','aniso','strike']
        self.titles = {'minmax':'Minimum and maximum\nresistivity, ohm-m',
                       'aniso':'Anisotropy in resistivity',
                       'strike':'Strike angle of\nminimum resistivity'}
        self.xlim = {'minmax':[0.1,1000],
                     'aniso':[0,20],
                     'strike':[0,180]}
        self.ylim = [6,0]
        self.modelno = 0
        self.modeltype = 'model'
        self.save = True
        self.output_filename = None

        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key]) 
    
    
    def plot_parameter(self,parameter):
        """
        base function for plotting a single model
        
        **parameter** string or list containing 'model', 'inmodel' or both
        tells the function whether to plot the model, inmodel (a priori), or both
        
        """
        data_list = []
        
        if 'model' in self.modeltype:
            if self.modeltype != 'inmodel':
                self.Model.read_model()
                data_list.append(self.models[self.modelno-1])
        
        if 'inmodel' in self.modeltype:
            self.read_inmodel()
            data_list.append(self.inmodel)
        
        
        # assess how many parameters
        allowed_params = ['minmax','aniso','strike']
        symbol = 'k-'
        
        if parameter not in allowed_params:
            print "invalid parameter"
            return
        
        for data in data_list:
            
            dep = data[:,1]
            resmin = data[:,2]
            resmax = data[:,3]
            strike = data[:,4]%180
               
            axes_count = 1
            
            if 'minmax' in parameter:
                plt.plot(resmax,dep,symbol)
                symbol += '-'
                plt.plot(resmin,dep,symbol)
                symbol = symbol[:-1]
                axes_count += 1
                plt.xscale('log')
            elif 'aniso' in parameter:
                plt.plot(resmax/resmin,dep,symbol)
                axes_count += 1
            elif 'strike' in parameter:
                plt.plot(strike,dep,symbol)    
                axes_count += 1      
            
            symbol = 'b-'
            
        plt.set_xlim(self.xlim[parameter][0],self.xlim[parameter][1])
        plt.set_ylim(self.ylim[0],self.ylim[1])
        plt.set_title(self.titles[parameter],fontsize=10)
        plt.grid()
        
        
    def plot(self):
        """
        plot 1 or more parameters in a model
        
        """
        n_axes = len(self.parameters)   
        
        sp = 1
        for parameter in self.parameters:
            plt.subplot(n_axes,1,sp)
            self.plot_parameter(parameter)
            sp += 1
        
        if self.save:
            if self.output_filename is None:
                self.construct_filename()
            plt.savefig(os.path.join(self.wd,self.output_filename))
                
            
    def plot_section(self):
        """
        plot 1d models along a profile
        
        """
        return

    
    def construct_filename(self):
        """
        create a filename to save file to
        
        """
        filename = os.path.basename(self.wd)
        filename += '_'.join(self.parameters)     
        self.output_filename = filename
        
        
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
                        levels=[round(min(mis)*draw_threshold/100.,2)],
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


class Plot_responses():
    """
    plot model responses and/or input data
    inherits a Response and Data class    
    
    """
    
    def __init__(self, Data, Response, **input_parameters):

        self.Data = Data
        self.Response = Response
        self.wd = self.Response.wd
        self.modelno = 0
        self.save = True
        self.output_filename = None

        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])
                     

    def plot_responses(self,datafile,modelno=0,adjust_phase = True):
        """
                
        
        """
        
        if not hasattr(self.Data,'resistivity'):
            self.read_datafile(datafile)
        if not hasattr(self.Response,'resistivity'):
            self.read_respfile()
            
        T = 1./self.freq
        r = self.Data.resistivity
        re = self.Data.resistivity_err
        p = self.Data.phase
        pe = self.Data.phase_err
        rmod = self.Response.resistivity[self.modelno]
        pmod = self.Response.phase[self.modelno]
        
        if adjust_phase:
            p[p<0] += 180
            pmod[pmod<0] += 180
        
        ii = 1
        fig = plt.figure()
        fig.text(0.5,0.95,self.Data.datafile[:-4])
        
        for mode,err,model in [[r,re,rmod],[p,pe,pmod]]:
            for sequence in [[[0,1],[1,0]],[[0,0],[1,1]]]:
                ax = plt.subplot(2,2,ii)
                ax.set_color_cycle(['b','r'])
                for i,j in sequence:
                    plt.errorbar(T,mode[:,i,j],err[:,i,j],fmt='--')
                    plt.plot(T,model[:,i,j],'k--')
                    plt.xscale('log')
                if ii <= 2:
                    plt.yscale('log')
                plt.grid()
                ii += 1