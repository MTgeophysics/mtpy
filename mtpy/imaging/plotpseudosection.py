# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:39:58 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import mtpy.imaging.mtplottools as mtpl
import mtpy.imaging.mtcolors as mtcl
import matplotlib.gridspec as gridspec

#==============================================================================

class PlotResPhasePseudoSection(object):
    """
    plot a resistivity and phase pseudo section for different components
    

    """
    
    def __init__(self, fn_lst=None, res_object_lst=None, z_object_lst=None, 
                 mt_object_lst=None, fignum=1, dpi=300, rot_z=0, plot_xy='y', 
                 plot_yx='y', plot_xx='n', plot_yy='n', plot_yn='y', 
                 plotnum=1, title=None, ftol=0.1, stationid=(0,4), 
                 linedir='ew', plot_period=None):
        """
        Initialize parameters
        """
        
        #--> get the inputs into a list of mt objects
        self.mt_lst = mtpl.get_mtlst(fn_lst=fn_lst, 
                                     res_object_lst=res_object_lst,
                                     z_object_lst=z_object_lst, 
                                     mt_object_lst=mt_object_lst)
        
        #set some of the properties as attributes much to Lars' discontent
        self.fignum = fignum
        self.fig_size = [8, 4]
        self.title = title
        self.dpi = dpi
        self.font_size = 7
        self.ftol = ftol
        self.stationid = stationid
        self.linedir = linedir
        self.aspect = 'auto'
        self.xtickspace = 1
        self.plot_yn = plot_yn
        self.plot_xx = plot_xx
        self.plot_xy = plot_xy
        self.plot_yx = plot_yx
        self.plot_yy = plot_yy
        self.plot_style = 'imshow'
        self.res_limits = (0, 3)
        self.phase_limits = (0, 90)
        self.pstep = 3
        self.imshow_interp = 'bicubic'
        self.plot_ylimits = None
        
        self.cb_pad = .02
        self.cb_orientation = 'horizontal'
        self.cb_shrink = .7
        
        #if rotation angle is an int or float make an array the length of 
        #mt_lst for plotting purposes
        if type(rot_z) is float or type(rot_z) is int:
            self.rot_z = np.array([rot_z]*len(self.mt_lst))
        
        #if the rotation angle is an array for rotation of different 
        #freq than repeat that rotation array to the len(mt_lst)
        elif type(rot_z) is np.ndarray:
            if rot_z.shape[0]!=len(self.mt_lst):
                self.rot_z = np.repeat(rot_z, len(self.mt_lst))
                
        else:
            self.rot_z = rot_z
            
        self.res_cmap = mtcl.cmapdict['mt_rd2gr2bl']
        self.phase_cmap = mtcl.cmapdict['mt_bl2gr2rd']
        
        #create empty lists to put things into
        self.stationlst = []
        self.offsetlst = []
        
        period_lst = np.array([len(mt.period) for mt in self.mt_lst])
        
        #find index where the longest period is if multiple pick the first one
        max_find = np.where(period_lst==period_lst.max())[0]
        if len(max_find)>0:
            max_find = max_find[0]
    
        if plot_period is None:
            self.plot_period = \
                  self.mt_lst[max_find].period
        
        else:
            self.plot_period = plot_period
            
        #create empty arrays to put data into
        ns = len(self.mt_lst)
        nt = len(self.plot_period)
        
        self.resxx = np.zeros((nt, ns))
        self.resxy = np.zeros((nt, ns))
        self.resyx = np.zeros((nt, ns))
        self.resyy = np.zeros((nt, ns))
        
        self.phasexx = np.zeros((nt, ns))
        self.phasexy = np.zeros((nt, ns))
        self.phaseyx = np.zeros((nt, ns))
        self.phaseyy = np.zeros((nt, ns))
        
        if self.plot_yn == 'y':
            self.plot()
        
    #---need to rotate data on setting rotz
    def _set_rot_z(self, rot_z):
        """
        need to rotate data when setting z
        """
        
        #if rotation angle is an int or float make an array the length of 
        #mt_lst for plotting purposes
        if type(rot_z) is float or type(rot_z) is int:
            rot_z = np.array([rot_z]*len(self.mt_lst))
        
        #if the rotation angle is an array for rotation of different 
        #freq than repeat that rotation array to the len(mt_lst)
        elif type(rot_z) is np.ndarray:
            if rot_z.shape[0]!=len(self.mt_lst):
                rot_z = np.repeat(rot_z, len(self.mt_lst))
                
        else:
            pass
        
        #rotate the data
        for ii,mt in enumerate(self.mt_lst):
            mt.rot_z = rot_z[ii]
            
    rot_z = property(fset=_set_rot_z, doc="rotation angle(s)")

    
    def sort_by_offsets(self):
        """
        get list of offsets to sort the mt list
        
        """
        
        dtype = [('station', 'S10'), ('offset', float), ('spot', int)]
        slst = []
        #get offsets
        for ii, mt in enumerate(self.mt_lst):
            #get offsets between stations
            if ii == 0:
                east0 = mt.lon
                north0 = mt.lat
                offset = 0.0
            else:
                east = mt.lon
                north = mt.lat
                #if line is predominantly e-w
                if self.linedir=='ew': 
                    if east0 < east:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif east0 > east:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset = 0
                #if line is predominantly n-s
                elif self.linedir == 'ns':
                    if north0 < north:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif north0 > north:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset=0
            #append values to list for sorting            
            slst.append((mt.station, offset, ii))
        
        #create a structured array according to the data type and values
        v_array = np.array(slst, dtype=dtype)
        
        #sort the structured array by offsets
        sorted_array = np.sort(v_array, order=['offset'])
        
        #create an offset list as an attribute
        self.offset_lst = np.array([ss[1] for ss in sorted_array])
        
        #create a station lst as an attribute
        self.station_lst = np.array([ss[0][self.stationid[0]:self.stationid[1]]
                                     for ss in sorted_array])
        
        #create an index list of the sorted index values 
        index_lst = [ss[2] for ss in sorted_array]
        
        #create a new mt_lst according to the offsets from the new index_lst
        new_mt_lst = [self.mt_lst[ii] for ii in index_lst]
        
        #set the mt_lst attribute as the new sorted mt_lst
        self.mt_lst_sort = new_mt_lst
        
    def get_rp_arrays(self):
        """
        get resistivity and phase values in the correct order according to 
        offsets and periods.

        """        
        
        self.sort_by_offsets()
        
        #create empty arrays to put data into need to reset to zero in case 
        #something has changed
        ns = len(self.mt_lst)
        nt = len(self.plot_period)
        
        self.resxx = np.zeros((nt, ns))
        self.resxy = np.zeros((nt, ns))
        self.resyx = np.zeros((nt, ns))
        self.resyy = np.zeros((nt, ns))
        
        self.phasexx = np.zeros((nt, ns))
        self.phasexy = np.zeros((nt, ns))
        self.phaseyx = np.zeros((nt, ns))
        self.phaseyy = np.zeros((nt, ns))
        
        #make a dictionary of the periods to plot for a reference
        period_dict = dict([(key, vv) 
                             for vv, key in enumerate(self.plot_period)])
                             
        for ii, mt in enumerate(self.mt_lst_sort):
            #get resisitivity and phase in a dictionary and append to a list
            rp = mt.get_ResPhase()
            
            for rr, rper in enumerate(self.plot_period):
                jj = None
                for kk, iper in enumerate(mt.period):
                    if iper == rper:
                        jj = period_dict[rper]
                        self.resxx[jj, ii] = np.log10(rp.resxx[kk])
                        self.resxy[jj, ii] = np.log10(rp.resxy[kk])
                        self.resyx[jj, ii] = np.log10(rp.resyx[kk])
                        self.resyy[jj, ii] = np.log10(rp.resyy[kk])
                        
                        self.phasexx[jj, ii] = rp.phasexx[kk]
                        self.phasexy[jj, ii] = rp.phasexy[kk]
                        self.phaseyx[jj, ii] = rp.phaseyx[kk]
                        self.phaseyy[jj, ii] = rp.phaseyy[kk]
                        
                        break
                        
                    elif rper*(1-self.ftol) <= iper and \
                         iper <= rper*(1+self.ftol):
                             jj = period_dict[rper]
                             self.resxx[jj, ii] = np.log10(rp.resxx[kk])
                             self.resxy[jj, ii] = np.log10(rp.resxy[kk])
                             self.resyx[jj, ii] = np.log10(rp.resyx[kk])
                             self.resyy[jj, ii] = np.log10(rp.resyy[kk])
                            
                             self.phasexx[jj, ii] = rp.phasexx[kk]
                             self.phasexy[jj, ii] = rp.phasexy[kk]
                             self.phaseyx[jj, ii] = rp.phaseyx[kk]
                             self.phaseyy[jj, ii] = rp.phaseyy[kk]
                             
                             break
                    else:
                        pass
                        
                if jj is None:
                    print 'did not find period {0:.6g} (s) for {1}'.format(
                               rper, self.station_lst[ii])
        
    def plot(self):
        
        #--> set subplot spacing
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['figure.subplot.right'] = .90
        plt.rcParams['figure.subplot.bottom'] = .09
        plt.rcParams['figure.subplot.top'] = .98

        #get apparent resistivity and phase
        self.get_rp_arrays()
        
        #make a list of tuples to see how many subplots are needed
        ynlst = [self.plot_xx+'xx', self.plot_xy+'xy', self.plot_yx+'yx', 
                 self.plot_yy+'yy']
        reslst = [self.resxx, self.resxy, self.resyx, self.resyy]
        phaselst = [self.phasexx, self.phasexy, self.phaseyx, self.phaseyy]
        plst = [(yn[1:], res, phase) for yn, res, phase in zip(ynlst, 
                                                           reslst, 
                                                           phaselst) 
                                                           if yn[0]=='y']
            
        #make a general subplot array
        gs = gridspec.GridSpec(2, len(plst), 
                               height_ratios=[1, 1], 
                               hspace=.00, 
                               wspace=.025)
                               
        #get ylimits for plot
        if self.plot_ylimits == None:
            self.plot_ylimits = (self.plot_period.min(), 
                                 self.plot_period.max())
                               
        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        ns = len(self.station_lst)
        #--> plot data
        self.fig = plt.figure(self.fignum, self.fig_size, dpi=self.dpi)
         
        #plot as a mesh where the data are left as blocks          
        if self.plot_style == 'pcolormesh':
            #need to add another element at the end of the array so pcolor 
            #will plot the full array
            xgrid, ygrid = np.meshgrid(np.append(self.offset_lst,
                                                 self.offset_lst[-1]*1.1), 
                                       np.append(self.plot_period,
                                                 self.plot_period[-1]*1.1))
            
            for ii, tt in enumerate(plst):
                axr = self.fig.add_subplot(gs[0, ii])
                axp = self.fig.add_subplot(gs[1, ii])
                
                #plot apparent resistivity
                axr.pcolormesh(xgrid, ygrid, np.flipud(tt[1]), 
                               cmap=self.res_cmap, 
                               vmin=self.res_limits[0], 
                               vmax=self.res_limits[1])
                               
                axr.set_xticks(self.offset_lst)
                plt.setp(axr.get_xticklabels(), visible=False)
                axr.grid(which='major', alpha=.25)
                axr.set_yscale('log')
                axr.set_xlim(self.offset_lst.min(), self.offset_lst.max()*1.1)
                axr.set_ylim(self.plot_ylimits)
                
                #label the plot with a text box
                axr.text(self.offset_lst.min()*.95,
                         self.plot_ylimits[1]*.8,
                         '$Z_{'+tt[0]+'}$',
                         fontdict=font_dict,
                         verticalalignment='top',
                         horizontalalignment='left',
                         bbox={'facecolor':'white', 'alpha':.5})
                
                #plot phase
                axp.pcolormesh(xgrid, ygrid, np.flipud(tt[2]), 
                               cmap=self.phase_cmap, 
                               vmin=self.phase_limits[0],
                               vmax=self.phase_limits[1])
                axp.grid(which='major', alpha=.25)
                axp.set_xticks(self.offset_lst)
                axp.set_xticklabels([self.station_lst[st] 
                                    for st in range(0,ns,self.xtickspace)])
                axp.set_yscale('log')
                axp.set_xlim(self.offset_lst.min(), self.offset_lst.max()*1.1)
                axp.set_ylim(self.plot_ylimits)
                if ii == 0:
                    axp.set_ylabel('Period (s)', font_dict)
                    axr.set_ylabel('Period (s)', font_dict)
                
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)
                
                #add colorbars                    
                if ii == len(plst)-1:
                    cminr = self.res_limits[0]
                    cmaxr = self.res_limits[1]
                    #add colorbar for res
                    axrpos = axr.get_position()
                    
                    #set position just to the right of the figure
                    cbr_position = (axrpos.bounds[0]+axrpos.bounds[2]+.0375,
                                    axrpos.bounds[1]+.05,
                                    .015,
                                    axrpos.bounds[3]*.75)
                                    
                    self.cbaxr = self.fig.add_axes(cbr_position)
                    self.cbr = mcb.ColorbarBase(self.cbaxr,
                                            cmap=self.res_cmap,
                                            norm=colors.Normalize(vmin=cminr,
                                                                  vmax=cmaxr),
                                            orientation='vertical')
                    tkrmin = np.ceil(cminr)
                    tkrmax = np.floor(cmaxr)
                        
                    self.cbr.set_ticks(np.arange(tkrmin, tkrmax+1))
                    cbr_ticklabels = [mtpl.labeldict[ll] 
                                      for ll in np.arange(tkrmin, tkrmax+1)]
                            
                    self.cbr.set_ticklabels(cbr_ticklabels)
                    self.cbr.ax.yaxis.set_label_position('right')
                    self.cbr.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbr.ax.yaxis.tick_left()
                    self.cbr.ax.tick_params(axis='y', direction='in', pad=1)
                    self.cbr.set_label('App. Res ($\Omega \cdot$m)', 
                                        fontdict={'size':self.font_size})
                    
                    #--> add colorbar for phase                    
                    cminp = self.phase_limits[0]
                    cmaxp = self.phase_limits[1]
                    
                    axppos = axp.get_position()
                    
                    #set position just to the right of the figure
                    cbp_position = (axppos.bounds[0]+axppos.bounds[2]+.0375,
                                    axppos.bounds[1]+.05,
                                    .015,
                                    axppos.bounds[3]*.75)
                                    
                    self.cbaxp = self.fig.add_axes(cbp_position)
                    self.cbp = mcb.ColorbarBase(self.cbaxp,
                                            cmap=self.phase_cmap,
                                            norm=colors.Normalize(vmin=cminp,
                                                                  vmax=cmaxp),
                                            orientation='vertical')
                    self.cbp.set_ticks([cminp, (cmaxp-cminp)/2, cmaxp])
                    self.cbp.set_ticklabels(['{0:.0f}'.format(cminp),
                                             '{0:.0f}'.format((cmaxp-cminp)/2),
                                             '{0:.0f}'.format(cmaxp)])
                    self.cbp.ax.yaxis.set_label_position('right')
                    self.cbp.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbp.ax.yaxis.tick_left()
                    self.cbp.ax.tick_params(axis='y', direction='in', pad=.5)
                    self.cbp.set_label('Phase (deg)', 
                                       fontdict={'size':self.font_size})
                
                #make axes attributes for user editing
                if tt == 'xx':
                    self.ax_rxx = axr
                    self.ax_pxx = axp
                elif tt == 'xy':
                    self.ax_rxy = axr
                    self.ax_pxy = axp
                elif tt == 'yx':
                    self.ax_ryx = axr
                    self.ax_pyx = axp
                elif tt == 'yy':
                    self.ax_ryy = axr
                    self.ax_pyy = axp
                
                    
                    
                
            plt.show()
        
        #plot data as an image which can have interpolation
        elif self.plot_style == 'imshow':   
            #make ticks simulate a log scale in the y direction
            #--> set major and minor ticks with appropriate labels 
            major_yticks = np.arange(np.ceil(np.log10(self.plot_ylimits[0])),
                                     np.floor(np.log10(self.plot_ylimits[1]))+1)
            
            #make minor ticks look like they are on a log scale
            minor_yticks = []
            for ll in major_yticks:
                minor_yticks += [np.arange(1,10)*10**ll]
            minor_yticks = np.array(minor_yticks)
            minor_yticks = np.log10(minor_yticks.flatten())
            
            #set ticklabels as 10** 
            yticklabels = [mtpl.labeldict[ll] for ll in major_yticks]
            
            for ii, tt in enumerate(plst):
                axr = self.fig.add_subplot(gs[0, ii])
                axp = self.fig.add_subplot(gs[1, ii])

                #plot apparent resistivity
                axr.imshow(tt[1], 
                           cmap=self.res_cmap, 
                           vmin=self.res_limits[0], 
                           vmax=self.res_limits[1],
                           aspect=self.aspect,
                           interpolation=self.imshow_interp,
                           extent=(self.offset_lst.min(), 
                                   self.offset_lst.max(),
                                   np.log10(self.plot_period.min()),
                                   np.log10(self.plot_period.max())))
                
                #set x ticks but remove labels
                axr.set_xticks(self.offset_lst)
                plt.setp(axr.get_xticklabels(), visible=False)
                
                #set y-axis ticks
                axr.yaxis.set_ticks(major_yticks)                    
                axr.yaxis.set_ticks(minor_yticks, minor=True)
                axr.set_yticklabels(yticklabels[::-1])
                
                axr.grid(which='major', alpha=.25)
                axr.set_xlim(self.offset_lst.min(), self.offset_lst.max())
                axr.set_ylim(np.log10(self.plot_ylimits[0]),
                             np.log10(self.plot_ylimits[1]))
                
                #label the plot with a text box
                axr.text(self.offset_lst.min()*.95,
                         np.log10(self.plot_ylimits[1])*.95,
                         '$Z_{'+tt[0]+'}$',
                         fontdict=font_dict,
                         verticalalignment='top',
                         horizontalalignment='left',
                         bbox={'facecolor':'white', 'alpha':.5})
                
                if ii == 0:
                    axr.set_ylabel('Period (s)', font_dict)
                
                #plot phase
                axp.imshow(tt[2], 
                           cmap=self.phase_cmap, 
                           vmin=self.phase_limits[0],
                           vmax=self.phase_limits[1],
                           aspect=self.aspect,
                           interpolation=self.imshow_interp,
                           extent=(self.offset_lst.min(), 
                                   self.offset_lst.max(),
                                   np.log10(self.plot_period.min()),
                                   np.log10(self.plot_period.max())))
                                       
                axp.grid(which='major', alpha=.25)
                axp.set_xticks(self.offset_lst)
                axp.set_xticklabels([self.station_lst[st] 
                                    for st in range(0,ns,self.xtickspace)])
                
                #remove tick labels if not the first subplot
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)
                
                #set y-axis ticks
                axp.yaxis.set_ticks(major_yticks)                    
                axp.yaxis.set_ticks(minor_yticks, minor=True)
                axp.set_yticklabels(yticklabels[::-1])
                
                axp.set_xlim(self.offset_lst.min(), self.offset_lst.max())
                axp.set_ylim(np.log10(self.plot_ylimits[0]),
                             np.log10(self.plot_ylimits[1]))
                
                if ii == 0:
                    axp.set_ylabel('Period (s)', font_dict)
                             
                #add colorbars                    
                if ii == len(plst)-1:
                    cminr = self.res_limits[0]
                    cmaxr = self.res_limits[1]
                    #add colorbar for res
                    axrpos = axr.get_position()
                    
                    #set position just to the right of the figure
                    cbr_position = (axrpos.bounds[0]+axrpos.bounds[2]+.0375,
                                    axrpos.bounds[1]+.05,
                                    .015,
                                    axrpos.bounds[3]*.75)
                                    
                    self.cbaxr = self.fig.add_axes(cbr_position)
                    self.cbr = mcb.ColorbarBase(self.cbaxr,
                                            cmap=self.res_cmap,
                                            norm=colors.Normalize(vmin=cminr,
                                                                  vmax=cmaxr),
                                            orientation='vertical')
                    tkrmin = np.ceil(cminr)
                    tkrmax = np.floor(cmaxr)
                        
                    self.cbr.set_ticks(np.arange(tkrmin, tkrmax+1))
                    cbr_ticklabels = [mtpl.labeldict[ll] 
                                      for ll in np.arange(tkrmin, tkrmax+1)]
                            
                    self.cbr.set_ticklabels(cbr_ticklabels)
                    self.cbr.ax.yaxis.set_label_position('right')
                    self.cbr.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbr.ax.yaxis.tick_left()
                    self.cbr.ax.tick_params(axis='y', direction='in', pad=1)
                    self.cbr.set_label('App. Res ($\Omega \cdot$m)', 
                                        fontdict={'size':self.font_size})
                    
                    #--> add colorbar for phase                    
                    cminp = self.phase_limits[0]
                    cmaxp = self.phase_limits[1]
                    
                    axppos = axp.get_position()
                    
                    #set position just to the right of the figure
                    cbp_position = (axppos.bounds[0]+axppos.bounds[2]+.0375,
                                    axppos.bounds[1]+.05,
                                    .015,
                                    axppos.bounds[3]*.75)
                                    
                    self.cbaxp = self.fig.add_axes(cbp_position)
                    self.cbp = mcb.ColorbarBase(self.cbaxp,
                                            cmap=self.phase_cmap,
                                            norm=colors.Normalize(vmin=cminp,
                                                                  vmax=cmaxp),
                                            orientation='vertical')
                    self.cbp.set_ticks([cminp, (cmaxp-cminp)/2+cminp, cmaxp])
                    self.cbp.set_ticklabels(['{0:.0f}'.format(cminp),
                                             '{0:.0f}'.format((cmaxp-cminp)/2+cminp),
                                             '{0:.0f}'.format(cmaxp)])
                    self.cbp.ax.yaxis.set_label_position('right')
                    self.cbp.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbp.ax.yaxis.tick_left()
                    self.cbp.ax.tick_params(axis='y', direction='in', pad=.5)
                    self.cbp.set_label('Phase (deg)', 
                                       fontdict={'size':self.font_size})
                                       
                if tt[0] == 'xx':
                    self.ax_rxx = axr
                    self.ax_pxx = axp
                elif tt[0] == 'xy':
                    self.ax_rxy = axr
                    self.ax_pxy = axp
                elif tt[0] == 'yx':
                    self.ax_ryx = axr
                    self.ax_pyx = axp
                elif tt[0] == 'yy':
                    self.ax_ryy = axr
                    self.ax_pyy = axp
                    
                
            plt.show()

    def save_plot(self, save_fn, file_format='pdf', orientation='portrait', 
                  fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_ResPhase.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            plt.clf()
            plt.close(self.fig)
            
        else:
            save_fn = os.path.join(save_fn, 
                                   self._mt.station+'_ResPhasePseudoSection.'+
                                    file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation)
        
        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)
        
        else:
            pass
        
        self.fig_fn = save_fn
        print 'Saved figure to: '+self.fig_fn

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> [ax.grid(True, which='major') for ax in [p1.axr,p1.axp]]
            >>> p1.update_plot()
        
        """

        self.fig.canvas.draw()
        
    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.xy_color = (.5,.5,.9)
            >>> p1.xy_marker = '*'
            >>> p1.redraw_plot()
        """
        
        plt.close(self.fig)
        self.plot()
        
    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """
        
        return "Plots Resistivity and phase as a pseudo section."
