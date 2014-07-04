# -*- coding: utf-8 -*-
"""
Created on Thu Jun 05 14:59:04 2014

@author: a1655681
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import mtpy.modeling.pek1dclasses as pek1dc
import mtpy.utils.elevation_data as ed
from matplotlib.font_manager import FontProperties

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
        self.Inmodel = None
        self.working_directory = self.Model.working_directory
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
                data_list.append(self.Model.models[self.modelno-1])
        
        if 'inmodel' in self.modeltype:
            self.Inmodel.read_inmodel()
            data_list.append(self.Inmodel.inmodel)
        
        
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
            
        plt.xlim(self.xlim[parameter][0],self.xlim[parameter][1])
        plt.ylim(self.ylim[0],self.ylim[1])
        plt.title(self.titles[parameter],fontsize=10)
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
            plt.savefig(os.path.join(self.working_directory,self.output_filename))
                
            
    def plot_section(self):
        """
        plot 1d models along a profile
        
        """
        return

    
    def construct_filename(self):
        """
        create a filename to save file to
        
        """
        filename = os.path.basename(self.working_directory)
        filename += '_'.join(self.parameters)     
        self.output_filename = filename
        

class Plot_fit():
    
    def __init__(self, Fit, **input_parameters):
        self.Fit = Fit
        self.imethod = 'linear'
        self.symbol = 'o'
        self.fontsize = 8
        self.cmap = 'rainbow'
        
        
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key]) 

        
    def plot_lcurve_contourmap(self, draw_threshold = False,
                               xlim = None, ylim = None,
                               contour_step = None
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
        
        self.Fit.read_fit()
        a = self.Fit.penalty_anisotropy
        s = self.Fit.penalty_structure
        mis = self.Fit.misfit
        aa = [str(round(i,1)) for i in self.Fit.weight_anisotropy]
        ss = [str(round(i,1)) for i in self.Fit.weight_structure]
           
        
        # grid the data to make contour plot
        # first, define points to grid
        points = np.vstack([s,a]).T
        # define points xi to grid onto
        xi = np.array(np.meshgrid(np.linspace(0.,max(s)),np.linspace(0.,max(a)))).T
#        print xi
        f1 = si.griddata(points,mis,xi,method = self.imethod)
        cmap = plt.get_cmap(self.cmap)
        
        if contour_step is not None:
            if contour_step < 1.:
                rounding = int(np.ceil(np.abs(np.log10(contour_step))))
            else:
                rounding = 0
            levels = np.arange(round(np.amin(self.Fit.misfit),rounding),
                   round(np.amax(self.Fit.misfit),rounding),
                   contour_step)
            plt.contour(xi[:,:,0],
                        xi[:,:,1],
                        f1,cmap=cmap,
                        levels=levels)
        else:
            plt.contour(xi[:,:,0],xi[:,:,1],f1,cmap=cmap)
        plt.colorbar()
        
        if draw_threshold:
            plt.contour(xi[:,:,0],xi[:,:,1],f1,
                        levels=[round(min(mis)*self.Fit.misfit_threshold,2)],
                                linewidths=[1.5],
                                colors='k')
        plt.gca().set_aspect('equal')
        plt.scatter(s,a,c='k',marker=self.symbol,lw=0)
        

                
        if xlim is None:
            xlim = plt.xlim()
        if ylim is None:
            ylim = plt.ylim()
        plt.xlim(xlim)
        plt.ylim(ylim)
            
        for i in range(len(aa)):
            if (s[i] < xlim[-1]) and (a[i] < ylim[-1]):
                plt.text(s[i],a[i],','.join([ss[i],aa[i]]),fontsize=self.fontsize)
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
        self.working_directory = self.Response.working_directory
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

      
         
         
class Plot_map():
    """
    class dealing with plotting results in map format, from multiple 
    1d models
    
    """
    def __init__(self,**input_parameters):
        self.levels = None
        self.n_levels = 10
        self.escale = 0.001
        self.anisotropy_threshold = 2.
        self.cmap = 'jet_r'
        self.scaleby = 'resmin'
        self.anisotropy_display_factor = 0.75
        self.xlim = None
        self.ylim = None
        self.cbar=True

        
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key]) 

        
    def plot_aniso_depth_map(self,aniso_depth_file,
                             scale='km',
                             header_rows=1):
        """
        """
        import matplotlib.patches as mpatches
        
        surface = np.loadtxt(aniso_depth_file,skiprows=header_rows)
        surface = surface[surface[:,4]/surface[:,3]>self.anisotropy_threshold]
        x,y,z,resmin,resmax,strike = [surface[:,i] for i in range(len(surface[0]))]
        aniso = resmax/resmin
                
        self.plot_interface(x,y,z,scale=scale)
        
        if self.scaleby == 'resmin':
            scale = 1./resmin
        else:
            scale = aniso**self.anisotropy_display_factor
                              
        # make rectangles
        recs = [mpatches.Rectangle(xy=np.array([x[i],y[i]]), 
                                   width = self.escale*scale[i],
                                   height = self.escale,
                                   angle=90-strike[i],
                                   lw=0) for i in range(len(x))]       
        ax1 = plt.gca()        
        for i,e in enumerate(recs):
            ax1.add_artist(e)
            e.set_facecolor('k')
        
    
    def plot_interface(self,x,y,z,scale='km'):

        import pek1dplotting as p1dp
        z = p1dp.update_scale(z,scale)

        xi = np.array(np.meshgrid(np.linspace(min(x),max(x),20),np.linspace(min(y),max(y),50))).T
        
        zi = si.griddata(np.vstack([x,y]).T,z,xi)
        
        cmap = plt.get_cmap(self.cmap)
        
        if self.levels is None:
            plt1 = plt.contourf(xi[:,:,0],xi[:,:,1],zi,cmap=cmap)
            self.levels = plt1.levels
        else:
            plt1 = plt.contourf(xi[:,:,0],xi[:,:,1],zi,
                                levels=self.levels,cmap=cmap)
        
        ax = plt.gca()
        ax.set_aspect('equal')
        
        if self.cbar:
            plt.colorbar()
        
        if self.xlim is not None:
            plt.xlim(self.xlim)
        if self.ylim is not None:
            plt.ylim(self.ylim)

        

        
    def plot_aniso_and_interfaces(self,aniso_depth_file,
                                 xyzfiles, header_rows = [1,1],
                                 scale=['km','km']):
        """
        plot a set of models and up to three interfaces for comparison.        
        
        """        
        import pek1dplotting as p1dp        
        
        if type(xyzfiles) == str:
            xyzfiles = [xyzfiles]
        if len(xyzfiles) == 1:
            s1,s2 = 2,1
            sp = [2]
            sp_labels = ['xy']
        elif len(xyzfiles) == 2:
            s1,s2 = 3,1
            sp = [2,3]
            sp_labels = ['y','xy']
        # if more than 3 interfaces provided, plot the first 3
        elif len(xyzfiles) >= 3:
            s1,s2 = 2,2
            sp = [2,3,4]
            sp_labels = ['','xy','x']
            
        xyz = np.loadtxt(aniso_depth_file,skiprows=header_rows[0])
        z = p1dp.update_scale(xyz[:,2],scale[0])
        zmin,zmax = np.amin(z),np.amax(z)
        
        for f in xyzfiles:
            z = p1dp.update_scale(np.loadtxt(f,skiprows=header_rows[1])[:,2],
                                  scale[1])
            if np.amin(z) < zmin:
                zmin = np.amin(z)
            if np.amax(z) > zmax:
                zmax = np.amax(z)
            
        self.levels = np.linspace(zmin,zmax,self.n_levels)
        
        x,y = xyz[:,0],xyz[:,1]
        ar = ((float(s1)/float(s2))*((np.amax(y) - np.amin(y))/(np.amax(x) - np.amin(x))))**0.9
        plt.figure(figsize=(10,10*ar))
        
        plt.subplot(s1,s2,1)
        self.plot_aniso_depth_map(aniso_depth_file,
                                  scale=scale[0],
                                  header_rows=header_rows[0])     
        plt.gca().set_xticklabels([])
        for s,ss in enumerate(sp):
            ax = plt.subplot(s1,s2,ss)
            print xyzfiles[s]
            xyz = np.loadtxt(xyzfiles[s],skiprows=header_rows[1])
            x,y,z = [xyz[:,i] for i in range(3)]
            self.plot_interface(x,y,z,
                                scale=scale[1],
                                cbar=False)
            if 'x' not in sp_labels[s]:
                ax.set_xticklabels([])
            if 'y' not in sp_labels[s]:
                ax.set_yticklabels([])

        
        ax = plt.axes([0.85,0.1,0.1,0.8])
        ax.set_visible(False)
        if self.cbar:
            plt.colorbar(fraction=0.8)


    def plot_location_map(self,
                          plot_names = True):
        """
        plot location map of all stations.        
        
        """


def update_scale(z,scale):
    
    if 'k' not in scale:
        z = z/1000.
    if '-' in scale:
        z = -1.*z
        
    return z                        
    
        

class Plot_profile():
    """

    """    
    
    
    def __init__(self,Model_suite,**input_parameters):
        
        self.working_directory = '.'
        self.station_listfile = None
        self.station_xyfile = None
        self.Model_suite = Model_suite
        self.modeltype = 'model'
        self.fig_width = 1.
        self.ax_width = 0.03
        self.ax_height = 0.8
        self.ax_bottom = 0.1
        self.plot_spacing = 0.02
        self.ylim = [6,0]
        self.titles = {'minmax':'Minimum and maximum resistivity, $\Omega m$',
                       'aniso':'Anisotropy in resistivity (maximum/minimum resistivity)',
                       'strike':'Strike angle of minimum resistivity, $^\circ$'}
        self.xlim = {'minmax':[0.1,1000],
                     'aniso':[0,20],
                     'strike':[0,180]}
                     
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])
            
        self.working_directory = os.path.abspath(self.working_directory)

            
    def get_profile(self):
        """
        get location of profile by linear regression
        """
        if self.Model_suite.x is None:
            print "can't get profile, no x y locations"
            return
        
        x = self.Model_suite.x
        y = self.Model_suite.y
        
        self.profile = np.polyfit(x,y,1)
        


    def get_profile_origin(self):
        """
        get the origin of the profile in real world coordinates
        
        Author: Alison Kirkby (2013)
        """
        if not hasattr(self,'profile'):
            self.get_profile()
            
        x,y = self.Model_suite.x,self.Model_suite.y
        x1,y1 = x[0],y[0]
        [m,c1] = self.profile
        x0 = (y1+(1.0/m)*x1-c1)/(m+(1.0/m))
        y0 = m*x0+c1
        
        self.profile_origin = [x0,y0]        
 
       
    def get_station_distance(self):
        """
        project a well location onto the profile.
        stores the result in the setup variable well_locations.
        x,y = x,y locations of well
        
        Author: Alison Kirkby
        """
        distances = []
        
        x,y = self.Model_suite.x,self.Model_suite.y
        
        if not hasattr(self,'profile_origin'):
            self.get_profile_origin()

        x0,y0 = self.profile_origin
        
            
        # project drillhole location onto profile and get distance along profile
        [m,c1] = self.profile
        [x0,y0] = self.profile_origin
        xp = (y+(1.0/m)*x-c1)/(m+(1.0/m))
        yp = m*x+c1
        xp -= x0
        yp -= y0
        distances = (xp**2.+yp**2.)**(0.5)
            
        self.station_distances = distances
#        for 
                
        
    def plot_parameter(self,parameter,
                       ylim=[6,0],
                       horizon_list = None,
                       horizon_zscale = 'km',
                       new_figure = True):
        """
        parameter = 'anisotropy', 'minmax', or 'strike' or list containing 
        several of these
        
        """
        import pek1dplotting as p1dp



        
        self.get_station_distance()
        xlim = self.xlim[parameter]        
        
        # define some initial locations to put the plots corresponding to distance along profile
        profile_x = (self.station_distances - np.amin(self.station_distances))
        
        # normalise so max distance is at 1
        profile_x /= (np.amax(self.station_distances)-np.amin(self.station_distances))
        
        # make an empty array to put buffered distances
        profile_x_buf = np.zeros_like(profile_x)      
        
        modelno = self.Model_suite.modelno

        px = self.ax_width    
        
        # shift each station along the profile so that they don't overlap each other
        for i in range(len(profile_x)):
            if i == 0:
                profile_x_buf[i] = profile_x[i]
            else:
                profile_x_buf[i] = max(profile_x[i],profile_x_buf[i-1]+px+self.plot_spacing)

        # renormalise so that end station is still within the plot bounds
        profile_x_buf /= np.amax(profile_x_buf)/(self.fig_width-3.*px)
        profile_x_buf += px

        font0 = FontProperties()
        font = font0.copy()
        font.set_family('serif')        
        
        if new_figure:        
            plt.figure(figsize=(len(profile_x),5*self.ax_height))

        for i in range(len(profile_x_buf)):

            data_list = []
            try:
                Model = self.Model_suite.model_list[i]
    
                data_list.append(Model.models[modelno-1])
    
                if len(self.Model_suite.inmodel_list) > 0:
                    Inmodel = self.Model_suite.inmodel_list[i]
                    data_list.append(Inmodel.inmodel)
                ax = plt.axes([profile_x_buf[i],self.ax_bottom,px,self.ax_height])
                for modelvals in data_list:
                    
                    if 'minmax' in parameter:
                        plt.plot(modelvals[:,3],modelvals[:,1],'0.5')
                        plt.plot(modelvals[:,2],modelvals[:,1],'k-')
                        plt.xscale('log')
                    if 'aniso' in parameter:
                        plt.plot(modelvals[:,3]/modelvals[:,2],modelvals[:,1],'k-')
                    if 'strike' in parameter:
                        plt.plot(modelvals[:,4]%180,modelvals[:,1],'k-')
                if horizon_list is not None:
                    for h in horizon_list:
                        elev = ed.get_elevation(Model.x,Model.y,h)
                        elev = p1dp.update_scale(elev,horizon_zscale)
                        plt.plot(xlim,[elev]*2)
    
                plt.ylim(self.ylim)
                plt.xlim(xlim)
                if i != 0:
                    ax.set_yticklabels([])
                else:
                    for label in ax.get_yticklabels():
                        label.set_fontproperties(font)
                ax.set_xticklabels([])
                plt.title(Model.station)
                if i == 0:
#                    plt.title(self.titles[parameter],
#                              fontproperties=font,
#                              loc = 'left')
                    plt.ylabel('$Depth, km$')

            except IndexError:
                print "station {} omitted".format(self.Model_suite.model_list[i].station)


            
        self.profile_x = profile_x_buf
        
    def plot_location_map(self):
        """
        plot location map of all stations with profile shown on map.        
        
        """
        
        if self.Model_suite.station_xyfile is None:
            print "can't get locations, no x y file"
            return
        
        if not hasattr(self,'profile_origin'):
            self.get_profile_origin()

        font0 = FontProperties()
        font = font0.copy()
        font.set_family('serif')        
            
        xy_all = np.genfromtxt(self.Model_suite.station_xyfile,invalid_raise=False)[:,1:]
        plt.plot(xy_all[:,0],xy_all[:,1],'.',c='0.5')
        
        m,c = self.profile
        x0,y0 = self.profile_origin
        
        if m > 1:
            y1 = max(self.Model_suite.y)
            x1 = (y1-c)/m
        else:
            x1 = max(self.Model_suite.x)
            y1 = m*x1 + c
        
        plt.plot([x0,x1],[y0,y1],'k')
        plt.plot(self.Model_suite.x,self.Model_suite.y,'k.')
        
        ax=plt.gca()
        for label in ax.get_yticklabels():
            label.set_fontproperties(font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(font)        
        
