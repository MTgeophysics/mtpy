# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 10:05:48 2012

@author: Jared Peacock
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec


class Strike:
    """
    Tools to deal with Strike (version 5.0) by:
    ===========================================
    
        McNeice, G. W.and Jones, A. G., 2001, Multisite, multifrequency tensor 
            decomposition of magnetotelluric data: *Geophysics*, **66**, 
            pg. 158--173.
    
    For Further Information:
    ------------------------
        Groom, R. W. and Bailey, R. C., 1989, Decomposition of magnetotelluric
            impedance tensors in the presence of local three-dimensional 
            galvanic distortion: *Journal of Geophysical Research*, *94*, 
            pg. 1913--1925 
        
    """
    
    def __init__(self):
        pass
    
    def readDCMP(self,dcmpfn):
        """
        read in the output file .dcmp by strike
        
        Arguments:
        ----------
        
            **dcmpfn** : string
                         full path to the output file from Strike (.dcmp)
        
        Returns:
        --------
            **Strike** data type
            
        Attributes:
        ----------
            * Strike.station
            * Strike.lat
            * Strike.lon
            * Strike.regional_strike 
            * Strike.shear  
            * Strike.chanel
            * Strike.twist
            * Strike.rho_a
            * Strike.rho_b
            * Strike.phase_a
            * Strike.phase_b
            * Strike.rms
            * Strike.skew
            * Strike.anisotropy
            * Strike.phasediff
        
        :Example: ::
            
            >>> st = Strike()
            >>> st.readDCMP(r"/test.dcmp")
            
        """
        
        dfid=file(dcmpfn,'r')
        dlines=dfid.readlines()
        
        kdict=dict([(key,[]) for key in ['regional','shear','chanel','twist',
                    'rho_a','rho_b','phase_a','phase_b','rms','skew',
                    'anisotropy','phasediff']])
                    
        ldict={'regional azimuth (cf AZIMUTH coordinate frame)':'regional',
               'shear angle':'shear',
               'channelling angle':'chanel',
               'twist angle':'twist',
               'app rho a':'rho_a',
               'app rho b':'rho_b',
               'imped phase a':'phase_a',
               'imped phase b':'phase_b',
               'av. rms error':'rms',
               'skew':'skew',
               'anis':'anisotropy',
               'phadif':'phasediff'}
        
        for dline in dlines:
            #get station name
            if dline.find('input file')>0:
                self.station=dline.strip().split('>')[1][:-4]
            #get location
            elif dline.find('>LATITUDE')==0:
                self.lat=float(dline.strip().split('=')[1])
                
            elif dline.find('>LONGITUDE')==0:
                self.lon=float(dline.strip().split('=')[1])
                
            elif dline.find('>ELEVATION')==0:
                self.elev=float(dline.strip().split('=')[1])
                
            elif dline.find('>AZIMUTH')==0:
                self.thetar=float(dline.strip().split('=')[1])
                
            elif dline.find('#')==0:
                pass
            
            #get data blocks
            else:
                dstr=dline.strip().split()
                #if the string is just a number pass
                if len(dstr)==1:
                    pass
            
                else:
                    #try to get data block
                    try:
                        test=float(dstr[1])
                        tlst=kdict[tkey].append([float(ii) for ii in dstr])
                    
                    #if not the data block figure out the dictionary key
                    except ValueError:
                        dstr=dline.strip().split(None,1)[1]
                        try:
                            tkey=ldict[dstr]
                        
                        #the rms key has an extra float at the end, remove it
                        except KeyError:
                            if dstr.find('rms')>0:
                                try:
                                    dstr=dstr.rsplit(None,1)[0]
                                    tkey=ldict[dstr]
                                except KeyError:
                                    print 'Did not find the key for ',dstr
            
            
        #make the data attributes
        self.regional_strike=np.array(kdict['regional'])
        self.shear=np.array(kdict['shear'])
        self.chanel=np.array(kdict['chanel'])
        self.twist=np.array(kdict['twist'])
        self.rho_a=np.array(kdict['rho_a'])
        self.rho_b=np.array(kdict['rho_b'])
        self.phase_a=np.array(kdict['phase_a'])
        self.phase_b=np.array(kdict['phase_b'])
        self.rms=np.array(kdict['rms'])
        self.skew=np.array(kdict['skew'])
        self.anisotropy=np.array(kdict['anisotropy'])
        self.phasediff=np.array(kdict['phasediff'])
            
    def plotHistogram(self,fignum=1,dpi=300,fs=8,
                      plotlst=['regional','skew','shear','twist','chanel'],):
                                        
        """
        plot the histogram of the different angles
        
        Arguments:
        ----------
            **fignum** : int (figure number label)
            
            **dpi** : int
                      Dots-per-inch resolution of figure
                      
            **fs** : float
                     font size
                     
            **plotlst** : list of angles to plot.  Options are:
                          * 'regional' for regional strike
                          * 'skew' for skew angle
                          * 'shear' for shear angle
                          * 'twist' for twist angle
                          * 'chanel' for channeling angle
        
        :Example: ::
            
            >>> st=Strike()
            >>> st.readDCMP(r"/test.dcmp")
            >>> st.plotHistogram(['regional'])
            
        """
        
        try:
            self.anisotropy
        except AttributeError:
            print 'Need to read in file first'
            
        plt.rcParams['font.size']=fs-2
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.left']=.1
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.hspace']=.05
        
        
        nl=len(plotlst)
        
        fig=plt.figure(fignum,dpi=dpi)
        plt.clf()
        for ii,ll in enumerate(plotlst,1):
            ax=fig.add_subplot(nl,1,ii)
            if ll=='regional':
                l1=ax.hist(self.regional_strike[:,1],color=(.2,0,.8),bins=90,
                           range=(-180,180),rwidth=5)
                label='Regional Strike'
            
            elif ll=='skew':
                l1=ax.hist(self.skew[:,1],color=(.9,0,.1),bins=90,
                           range=(-180,180),rwidth=5)
                label='Skew Angle'
            
            elif ll=='twist':
                l1=ax.hist(self.twist[:,1],color=(.3,.3,0),bins=90,
                           range=(-180,180),rwidth=5)
                label='Twist Angle'
                
            elif ll=='shear':
                l1=ax.hist(self.shear[:,1],color=(.2,.8,0),bins=90,
                           range=(-180,180),rwidth=5)
                label='Shear Angle'
            
            elif ll=='chanel':
                l1=ax.hist(self.chanel[:,1],color=(.5,.2,.5),bins=90,
                           range=(-180,180),rwidth=5)
                label='Channeling Angle'
                

            ax.set_ylim(0,len(self.regional_strike[:,0]))
            ax.set_xlim(-180,180)
            
            ax.xaxis.set_major_locator(MultipleLocator(30))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            
            if ii<nl:
                plt.setp(ax.xaxis.get_ticklabels(),visible=False)
            else:
               ax.set_xlabel('Angle',fontdict={'size':fs,'weight':'bold'}) 

            ax.set_ylabel('Counts',fontdict={'size':fs,'weight':'bold'})
            ax.legend([l1[2][0]],[label],loc='upper left',prop={'size':fs-2},
                      borderaxespad=.05)
            ax.grid(which='both',alpha=.25)
            plt.show()
            
    def plotAngles(self,dpi=300,fs=8,ms=4,fignum=1,
                   plotlst=['regional','skew','shear','twist','chanel','rms']):
        """
        plot the angles vs log period
        
        Arguments:
        ----------
            **fignum** : int (figure number label)
            
            **dpi** : int
                      Dots-per-inch resolution of figure
                      
            **fs** : float
                     font size
                     
            **ms** : float
                     marker size
                     
            **plotlst** : list of angles to plot.  Options are:
                          * 'regional' for regional strike
                          * 'skew' for skew angle
                          * 'shear' for shear angle
                          * 'twist' for twist angle
                          * 'chanel' for channeling angle
        
        :Example: ::
            
            >>> st=Strike()
            >>> st.readDCMP(r"/test.dcmp")
            >>> st.plotAngles(['regional','skew','rms'])
        """
        
        try:
            self.anisotropy
        except AttributeError:
            print 'Need to read in file first'
            
        plt.rcParams['font.size']=fs-2
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.left']=.1
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.hspace']=.05
        
        #set some empty lists for the legend
        labellst=[]
        linelst=[]
        
        fig=plt.figure(fignum,dpi=dpi)
        plt.clf()
        
        #set a grid of subplot
        gs=gridspec.GridSpec(2,2,height_ratios=(5,1))
        ax=fig.add_subplot(gs[0,:])
        
        #------plot angles----------------------------- 
        for ii,ll in enumerate(plotlst,1):
            
            if ll=='regional':
                l1=ax.semilogx(self.regional_strike[:,0],
                               self.regional_strike[:,1],
                               color=(.2,0,.8),ls='None',marker='s',ms=ms)
                labellst.append('Regional Strike')
                linelst.append(l1[0])
            
            elif ll=='skew':
                l1=ax.semilogx(self.skew[:,0],self.skew[:,1],color=(.9,0,.1),
                               ls='None',marker='d',ms=ms)
                labellst.append('Skew Angle')
                linelst.append(l1[0])
            
            elif ll=='twist':
                l1=ax.semilogx(self.twist[:,0],self.twist[:,1],
                               color=(.3,.3,0),ls='None',marker='v',ms=ms)
                labellst.append('Twist Angle')
                linelst.append(l1[0])
                
            elif ll=='shear':
                l1=ax.semilogx(self.shear[:,0],self.shear[:,1],
                               color=(.2,.8,0),ls='None',marker='h',ms=ms)
                labellst.append('Shear Angle')
                linelst.append(l1[0])
            
            elif ll=='chanel':
                l1=ax.semilogx(self.chanel[:,0],self.chanel[:,1],
                               color=(.5,.2,.5),ls='None',marker='p',ms=ms)
                labellst.append('Channeling Angle')
                linelst.append(l1[0])
        
        #get ylimits        
        ymax=np.max([self.regional_strike[:,1].max(),self.skew[:,1].max(),
                     self.shear[:,1].max(),self.twist[:,1].max(),
                     self.chanel[:,1].max()])
        ymin=np.min([self.regional_strike[:,1].min(),self.skew[:,1].min(),
                     self.shear[:,1].min(),self.twist[:,1].min(),
                     self.chanel[:,1].min()])
        ax.set_ylim(ymin-5,ymax+5)
        
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        
        plt.setp(ax.xaxis.get_ticklabels(),visible=False)
        
        ax.set_ylabel('Angle (deg)',fontdict={'size':fs,'weight':'bold'}) 

        ax.legend(linelst,labellst,loc='upper left',prop={'size':fs-2},
                  borderaxespad=.05)
        ax.grid(which='both',alpha=.25)
        
        #----plot the rms---------------------------
        ax2=fig.add_subplot(gs[1,:])
        l1=ax2.semilogx(st.rms[:,0],st.rms[:,1],color=(.2,.6,.8),ls='None',
                       marker='o',ms=ms)
        ax2.set_ylim(0,st.rms[:,1].max()+.5)
        
        ax2.yaxis.set_major_locator(MultipleLocator(1))
        ax2.yaxis.set_minor_locator(MultipleLocator(.2))
        
        ax2.set_ylabel('RMS',fontdict={'size':fs,'weight':'bold'}) 

        ax2.set_xlabel('Period (s)',fontdict={'size':fs,'weight':'bold'})
        ax2.grid(which='both',alpha=.25)
        
        plt.show()
        
    def plotResPhase(self,ms=4,fs=8,fignum=1,dpi=300):
        """
        plot the regional resistivity and phase
        
        Arguments:
        ----------
            **fignum** : int (figure number label)
            
            **dpi** : int
                      Dots-per-inch resolution of figure
                      
            **fs** : float
                     font size
                     
            **ms** : float
                     marker size
          
        :Example: ::
            
            >>> st=Strike()
            >>> st.readDCMP(r"/test.dcmp")
            >>> st.plotResPhase(fignum=2)          
        
        """
        
        try:
            self.anisotropy
        except AttributeError:
            print 'Need to read in file first'
            
        plt.rcParams['font.size']=fs-2
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.left']=.1
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.hspace']=.05

        gs=gridspec.GridSpec(2,2,height_ratios=[2,1])
        
        fig=plt.figure(fignum,dpi=dpi)
        plt.clf()
        
        #plot apparent resistivity
        axr=fig.add_subplot(gs[0,:])
        r1=axr.loglog(self.rho_a[:,0],self.rho_a[:,1],ls='none',marker='s',
                      color=(.1,0,.9),ms=ms)
        r2=axr.loglog(self.rho_b[:,0],self.rho_b[:,1],ls='none',marker='o',
                      color=(.9,0,.1),ms=ms)
        
        plt.setp(axr.xaxis.get_ticklabels(),visible=False)
        axr.set_ylabel('App. Res. ($\Omega \cdot$m)',
                       fontdict={'size':fs,'weight':'bold'})
        axr.legend([r1[0],r2[0]],['Regional_a','Regional_b'],prop={'size':fs-2},
                   borderaxespad=.05,loc='upper left')
        axr.grid(which='both',alpha=.25)    
        
        #plot phase
        axp=fig.add_subplot(gs[1,:])
        axp.semilogx(self.phase_a[:,0],self.phase_a[:,1],ls='none',marker='s',
                     color=(.1,0,.9),ms=ms)
        axp.semilogx(self.phase_b[:,0],self.phase_b[:,1],ls='none',marker='o',
                     color=(.9,0,.1),ms=ms)
                     
        axp.set_xlabel('Period (s)',fontdict={'size':fs,'weight':'bold'})
        axp.set_ylabel('Phase (deg)',
                       fontdict={'size':fs,'weight':'bold'})
        axp.grid(which='both',alpha=.25) 
        
        plt.show()
        








            
        
        
        
    