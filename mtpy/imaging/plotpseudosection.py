# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:39:58 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MultipleLocator
import mtpy.imaging.mtplottools as mtpl

#==============================================================================

class PlotResPhasePseudoSection(object):
    """
    plot a resistivity and phase pseudo section for different components
    

    """
    
    def __init__(self, fn_lst=None, res_object_lst=None, z_object_lst=None, 
                 mt_object_lst=None, fignum=1, dpi=300, rot_z=0, plot_xy='y', 
                 plot_yx='y', plot_xx='n', plot_yy='n', plot_yn='y', 
                 plotnum=1, title=None, ftol=0.1):
        """
        Initialize parameters
        """
        
        #--> get the inputs into a list of mt objects
        self.mt_lst = get_mtlst(fn_lst=fn_lst, 
                                 res_object_lst=res_object_lst,
                                 z_object_lst=z_object_lst, 
                                 mt_object_lst=mt_object_lst)
        
        #set some of the properties as attributes much to Lars' discontent
        self.fignum = fignum
        self.title = title
        self.dpi = dpi
        self.font_size = 7
        self.ftol = ftol
        
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
            
        self.res_cmap = 'mt_rd2gr2bl'
        self.phase_cmap = 'mt_bl2gr2rd'
        
        
        
    def plot(self):
        
    
        plt.rcParams['font.size']=self.font_size
        plt.rcParams['figure.subplot.left']=.07
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.06
        plt.rcParams['figure.subplot.top']=.94
        plt.rcParams['figure.subplot.wspace']=.01
        plt.rcParams['figure.subplot.hspace']=.20

        
        #create empty lists to put things into
        stationlst = []
        offsetlst = []
        minperiodlst = []
        maxperiodlst = []
        periodlst = []
        nlst = []
        
        period_lst = [len(mt.period) for mt in self.mt_lst]
        
        self.plot_period = np.where
        
        #create empty arrays to put data into
        n=len(self.mt_lst)
        resxx=np.zeros((maxperiod,n))
        resxy=np.zeros((maxperiod,n))
        resyx=np.zeros((maxperiod,n))
        resyy=np.zeros((maxperiod,n))
        
        phasexx=np.zeros((maxperiod,n))
        phasexy=np.zeros((maxperiod,n))
        phaseyx=np.zeros((maxperiod,n))
        phaseyy=np.zeros((maxperiod,n))
        
        #get data into arrays
        for ii,fn in enumerate(edifilelst):
            #get offsets between stations
            imp=Z.Z(fn)
            stationlst.append(imp.station[stationid[0]:stationid[1]])
            zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
            
            if ii==0:
                east0=east
                north0=north
                offset=0.0
            else:
                if linedir=='ew': 
                    if east0<east:
                        offset=np.sqrt((east0-east)**2+(north0-north)**2)
                    elif east0>east:
                        offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset=0
                elif linedir=='ns':
                    if north0<north:
                        offset=np.sqrt((east0-east)**2+(north0-north)**2)
                    elif north0>north:
                        offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset=0
            offsetlst.append(offset)
            #get resisitivity and phase in a dictionary and append to a list
            rp=imp.getResPhase(ffactor=ffactor,thetar=rotz)
            m=len(imp.period)
            resxx[0:m,ii]=np.log10(rp.resxx)#[::-1]
            resxy[0:m,ii]=np.log10(rp.resxy)#[::-1]
            resyx[0:m,ii]=np.log10(rp.resyx)#[::-1]
            resyy[0:m,ii]=np.log10(rp.resyy)#[::-1]
            
            phasexx[0:m,ii]=rp.phasexx#[::-1]
            phasexy[0:m,ii]=rp.phasexy#[::-1]
            phaseyx[0:m,ii]=rp.phaseyx+180#[::-1]
            phaseyy[0:m,ii]=rp.phaseyy#[::-1]
            
            #get min and max of the period, just in case some edi files have 
            #different period information
            minperiodlst.append(min(imp.period))
            maxperiodlst.append(max(imp.period))
            nlst.append(len(imp.period))
            periodlst.append(imp.period)
            
        nperiod=max(nlst)
        for pp in range(len(periodlst)):
            if nlst[pp]==nperiod:
                period=np.log10(periodlst[pp])#[::-1]
    
        #plot data
        fig=plt.figure(1,[12,4],dpi=150)
        extent=(0,n,period[-1],period[0])
        aspect=aspect
        cmap=cmap
        pad=.2
        kwargs={'aspect':aspect,'cmap':cmap,'extent':extent}
        cbarkwargs={'shrink':.7,'orientation':'horizontal','pad':pad}
    
        
        ax2=plt.subplot(2,4,2)
        plt.imshow(resxy[0:nperiod,:],**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
    #    plt.ylabel('Log$_{10}$ Period',fontsize=10,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('Log$_{10}\mathbf{(1/\sigma_{xy})}$',fontsize=fs+4,fontweight='bold')
        ax2.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
        ax3=plt.subplot(2,4,3)
        plt.imshow(resyx[0:nperiod,:],**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('Log$_{10}\mathbf{(1/\sigma_{yx})}$',fontsize=fs+4,fontweight='bold')
        ax3.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
        ax6=plt.subplot(2,4,6)
        plt.imshow(phasexy[0:nperiod,:],vmin=0,vmax=90,**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('$\mathbf{\phi_{xy}}$',fontsize=fs+4)
        ax6.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
        ax7=plt.subplot(2,4,7)
        plt.imshow(phaseyx[0:nperiod,:],vmin=0,vmax=90,**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('$\mathbf{\phi_{yx}}$',fontsize=fs+4)
        ax7.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
    #    fig2=plt.figure(2,dpi=150)    
    
        ax1=plt.subplot(2,4,1)
        plt.imshow(resxx[0:nperiod,:],**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
        plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('Log$_{10}\mathbf{(1/\sigma_{xx})}$',fontsize=fs+4,fontweight='bold')
        ax1.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
        ax4=plt.subplot(2,4,4)
        plt.imshow(resyy[0:nperiod,:],**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('Log$_{10}\mathbf{(1/\sigma_{yy})}$',fontsize=fs+4,fontweight='bold')
        ax4.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
        ax5=plt.subplot(2,4,5)
        plt.imshow(phasexx[0:nperiod,:],**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
        plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('$\mathbf{\phi_{xx}}$',fontsize=fs+4)
        ax5.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        
        ax8=plt.subplot(2,4,8)
        plt.imshow(phaseyy[0:nperiod,:],**kwargs)
        plt.colorbar(**cbarkwargs)
        plt.xticks(np.arange(0,n,xtickspace),
                   [stationlst[st] for st in range(0,n,xtickspace)])
    #    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
        plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
        plt.title('$\mathbf{\phi_{yy}}$',fontsize=fs+4)
        ax8.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.grid()
        plt.show()
