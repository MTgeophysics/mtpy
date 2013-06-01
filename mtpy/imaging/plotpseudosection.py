# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:39:58 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
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
        
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .07
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .06
        plt.rcParams['figure.subplot.top'] = .94
        plt.rcParams['figure.subplot.wspace'] = .01
        plt.rcParams['figure.subplot.hspace'] = .20

        self.get_rp_arrays()
        
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
                               
        
        #--> plot data
        self.fig = plt.figure(self.fignum, self.fig_size, dpi=self.dpi)
                   
        cbarkwargs = {'shrink':self.cb_shrink,
                      'orientation':self.cb_orientation,
                      'pad':self.cb_pad}
        
        ns = len(self.station_lst)
        
        if self.plot_style == 'pcolormesh':
            xgrid, ygrid = np.meshgrid(self.offset_lst, self.plot_period)
            
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
                axr.set_xlim(self.offset_lst.min(), self.offset_lst.max())
                axr.set_ylim(self.plot_period.min(), self.plot_period.max())
                
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
                axp.set_xlim(self.offset_lst.min(), self.offset_lst.max())
                axp.set_ylim(self.plot_period.min(), self.plot_period.max())
                
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)
                    
                
            plt.show()
        elif self.plot_style == 'imshow':
            nt = len(self.plot_period)
            major_yticks = [np.log10(self.plot_period[ll]) 
                            for ll in np.arange(0, nt, self.pstep)][::-1]
            minor_yticks = [np.log10(self.plot_period[ll]) 
                            for ll in np.arange(0, nt, 1)][::-1]
            yticklabels = ['{0:>4}'.format('{0:.1e}'.format(
                                                       self.plot_period[ll])) 
                            for ll in np.arange(0, nt, self.pstep)]+\
                          ['{0:>4}'.format('{0:.1e}'.format(
                                                       self.plot_period[-1]))]
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
                            
                axr.set_xticks(self.offset_lst)
                plt.setp(axr.get_xticklabels(), visible=False)
                axr.grid(which='major', alpha=.25)
                axr.set_xlim(self.offset_lst.min(), self.offset_lst.max())
                axr.set_ylim(np.log10(self.plot_period.min()), 
                             np.log10(self.plot_period.max()))
                
                #set y-axis major ticks
                axr.yaxis.set_ticks(major_yticks)
                #set y-axis minor ticks                     
                axr.yaxis.set_ticks(minor_yticks, minor=True)
                #set y-axis tick labels
                axr.set_yticklabels(yticklabels[1:-1]+['',''])
                
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
                axp.set_xlim(self.offset_lst.min(), self.offset_lst.max())
                axp.set_ylim(np.log10(self.plot_period.min()),
                             np.log10(self.plot_period.max()))
                
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)
                
                #set y-axis major ticks
                axp.yaxis.set_ticks(major_yticks)
                #set y-axis minor ticks                     
                axp.yaxis.set_ticks(minor_yticks, minor=True)
                #set y-axis tick labels
                axp.set_yticklabels(yticklabels[1:]+[''])
                    
                
            plt.show()
#            plt.title('Log$_{10}\mathbf{(1/\sigma_{xy})}$',fontsize=fs+4,fontweight='bold')
#            ax2.yaxis.set_minor_locator(MultipleLocator(.2))
#            plt.grid()       
#            if tt[0] == 'xx':
#                self.ax_rxx = axr
#                self.ax_pxx = axp
#    
#        ax2=plt.subplot(2,4,2)
#        plt.imshow(resxy[0:nperiod,:],**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#    #    plt.ylabel('Log$_{10}$ Period',fontsize=10,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('Log$_{10}\mathbf{(1/\sigma_{xy})}$',fontsize=fs+4,fontweight='bold')
#        ax2.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#        ax3=plt.subplot(2,4,3)
#        plt.imshow(resyx[0:nperiod,:],**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('Log$_{10}\mathbf{(1/\sigma_{yx})}$',fontsize=fs+4,fontweight='bold')
#        ax3.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#        ax6=plt.subplot(2,4,6)
#        plt.imshow(phasexy[0:nperiod,:],vmin=0,vmax=90,**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('$\mathbf{\phi_{xy}}$',fontsize=fs+4)
#        ax6.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#        ax7=plt.subplot(2,4,7)
#        plt.imshow(phaseyx[0:nperiod,:],vmin=0,vmax=90,**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('$\mathbf{\phi_{yx}}$',fontsize=fs+4)
#        ax7.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#    #    fig2=plt.figure(2,dpi=150)    
#    
#        ax1=plt.subplot(2,4,1)
#        plt.imshow(resxx[0:nperiod,:],**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#        plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('Log$_{10}\mathbf{(1/\sigma_{xx})}$',fontsize=fs+4,fontweight='bold')
#        ax1.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#        ax4=plt.subplot(2,4,4)
#        plt.imshow(resyy[0:nperiod,:],**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('Log$_{10}\mathbf{(1/\sigma_{yy})}$',fontsize=fs+4,fontweight='bold')
#        ax4.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#        ax5=plt.subplot(2,4,5)
#        plt.imshow(phasexx[0:nperiod,:],**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#        plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('$\mathbf{\phi_{xx}}$',fontsize=fs+4)
#        ax5.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        
#        ax8=plt.subplot(2,4,8)
#        plt.imshow(phaseyy[0:nperiod,:],**kwargs)
#        plt.colorbar(**cbarkwargs)
#        plt.xticks(np.arange(0,n,xtickspace),
#                   [stationlst[st] for st in range(0,n,xtickspace)])
#    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#        plt.title('$\mathbf{\phi_{yy}}$',fontsize=fs+4)
#        ax8.yaxis.set_minor_locator(MultipleLocator(.2))
#        plt.grid()
#        plt.show()
