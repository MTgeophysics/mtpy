# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 17:27:06 2010

Author: Jared Peacock

Based on MTAnalysis.m and LoadEDI.m written by Stephan Thiel 
"""

import numpy as np
import matplotlib.pyplot as plt
from mtpy.core.z import Z
import os
from matplotlib.ticker import FormatStrFormatter,MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.colorbar import *

#make a custom colormap to use for plotting
ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

class ptplots(Z):
    """ptplots will plot different parameters of the phase tensor and invariants
    given a list of plot types for an input of impvar which can be an .imp file,
    .edi file or a dictionary with keys:
     impvar['station'] -> station name 
     impvar['lat'] -> latitude in deg 
     impvar['lon'] -> longitude in deg 
     impvar['periodlst'] -> period list 
     impvar['zlst'] -> impedance tensor list at [[zxx,zxy],[zyx,zyy]]
     impvar['zvar'] -> impedance tensor variance same fmt as z
     impvar['tip'] -> tipper array as [txx,txy]
     impvar['tipvar'] -> tipper variance as [txxvar,tyyvar]
    
    Plots are: 
    
     plotPhaseTensor(save='n',fmt='pdf',fignum=1,thetar=0) -> plots the phase tensor as 
     ellipses with long axis directed towards electrical strike and short axis
     directed across magnetic strike.  Spaces the ellipses across the period axis.
    
     plotStrikeAngle(save='n',fmt='pdf',fignum=1,thetar=0) -> plots the electric 
     strike angle as a function of period with error bars. 
    
     plotMinMaxPhase(save='n',fmt='pdf',fignum=1,thetar=0) -> plots the min and max 
     phase angle as a function of period with error bars.  
    
     plotAzimuth(save='n',fmt='pdf',fignum=1,thetar=0) -> plots the azimuth angle as
     a function of period with error bars. 
    
     plotSkew(save='n',fmt='pdf',fignum=1,thetar=0) -> plots the skew angle as a 
     function of period with error bars.  
    
     plotElliticity(save='n',fmt='pdf',fignum=1,thetar=0) -> plots the ellipticity 
     as a function of period with error bars. 
    
     plotAll(save='n',fmt='pdf',fignum=1,thetar=0) -> plots phase tensor, strike angle, 
     min and max phase angle, azimuth, skew, and ellipticity as subplots on one
     plot. 
    
    You can save the plots by making save='y' or a path you input. fmt is the
    fmt you want to save the figure as. Can be: pdf,eps,ps,png or svg. Fignum
    is the figure number in the event you want to plot multipl plots."""

    colorlst=['b','r','k','g','c','m','y']
    markerlst=['s','o','d','h','v','^','>','<','p','*','x']
    pcmlst=[(mm,cc) for mm in markerlst for cc in colorlst]
        
    
    def __init__(self,impvar):
        self.impvar=impvar
        if type(impvar) is list:
            self.z=[]
            for fn in impvar:
                self.z.append(Z(fn))
            self.savepath=os.path.join(os.path.dirname(self.impvar[0]),'PTPlots')
        else:
            self.z=[Z(impvar)]
        if type(impvar) is str:
            self.savepath=os.path.join(os.path.dirname(self.impvar),'PTPlots')
        else:
            self.savepath=os.path.join(os.getcwd(),'PTPlots')
        #make a list of colors and makers to use for plotting
        self.mclst=[]
        
        for cc in ptplots.colorlst:
            for mm in ptplots.markerlst:
                self.mclst.append([cc,mm])
    def plotPhaseTensor(self,spacing=6,esize=5,colorkey='beta',save='n',fmt='pdf',
                        fignum=1,thetar=0):
        """plotPhaseTensor will plot phase tensor ellipses. Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. """
        stationstr=self.z[0].station
        stationlst=[]
        
        plt.rcParams['font.size']=7
        plt.rcParams['figure.subplot.left']=.07
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.08
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.wspace']=.2
        plt.rcParams['figure.subplot.hspace']=.4
          
        
        for dd in range(len(self.z)):
            
            #get phase tensor elements 
            pt=self.z[dd].getPhaseTensor(thetar=thetar)
            periodlst=self.z[dd].period
            n=len(periodlst)
            station=self.z[dd].station
            stationlst.append(station)
            if dd!=0:
                stationstr+=','+station
            else:
                pass
            
            #create a plot instance
            fig=plt.figure(fignum,[10,2],dpi=150)
            ax=fig.add_subplot(1,1,1,aspect='equal')
            for ii in range(n):
                #make sure the ellipses will be visable
                scaling=esize/pt.phimax[ii]
                eheight=pt.phimin[ii]*scaling
                ewidth=pt.phimax[ii]*scaling
                    
                #create an ellipse scaled by phimin and phimax and oriented along
                #the azimuth    
                ellip=Ellipse((spacing*ii,0),width=ewidth,
                              height=eheight,
                              angle=pt.azimuth[ii])
                ax.add_artist(ellip)
                
                if pt.phimin[ii]<0 or pt.phimin[ii]=='nan':
                    cvars=0
#                    print 'less than 0 or nan',cvars
                else:
                    cvars=(pt.phimin[ii]/(np.pi/2))
#                    print 'calculated ',cvars
                    if cvars>1.0:
                        cvars=0
#                        print 'greater than 1 ',cvars
#                print cvars
                    
                ellip.set_facecolor((1,1-cvars,.1))
            
        xticklabels=['%2.2g' % periodlst[ii] for ii in np.arange(start=0,stop=n,
                     step=3)]
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.title('Phase Tensor Ellipses for '+stationstr,fontsize=14)
        plt.xticks(np.arange(start=0,stop=spacing*n,step=3*spacing),xticklabels)
        ax.set_ylim(-1*(spacing+3),spacing+3)
        ax.set_xlim(-spacing,n*spacing+3)
        plt.grid()
        ax2=make_axes(ax,shrink=.4)
        cb=ColorbarBase(ax2[0],cmap=ptcmap,norm=Normalize(vmin=0,vmax=90))
        cb.set_label('Min Phase')
        
        if len(stationlst)>1:
            plt.legend(stationlst,loc=0,markerscale=.4,borderaxespad=.05,
                   labelspacing=.1,handletextpad=.2)
            leg=plt.gca().get_legend()
            ltext=leg.get_texts()  # all the text.Text instance in the legend
            plt.setp(ltext, fontsize=10)    # the legend text fontsize
        
        plt.show()
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                     self.z.station+'PT.'+fmt),fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                                   self.z.station+'PT.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'PT.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                                   self.z.station+'PT.'+fmt)
        elif save=='n':
            pass
        
    def plotStrikeangle(self,save='n',fmt='pdf',fignum=1,thetar=0):
        """Plot the strike angle as calculated from the invariants. Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. """
        
        stationstr=self.z[0].station
        stationlst=[]
        mlst=[]
        
        if len(self.z)>0:
            medlst=[]
        for dd in range(len(self.z)):
            mm=np.remainder(dd,4)
            zinv=self.z[dd].getInvariants(thetar=thetar)
            period=self.z[dd].period
            medlst.append(zinv.strike)
            
            fig=plt.figure(fignum,[6,8],dpi=100)
            ax=plt.subplot(1,1,1)
            erxy=plt.errorbar(period,zinv.strike,
                        marker=self.mclst[mm][1],ms=4,mfc='None',
                        mec=self.mclst[mm][0],mew=1,ls='None',
                        yerr=zinv.strikeerr,ecolor=self.mclst[mm][0])
            mlst.append(erxy[0])
            stationlst.append(self.z[dd].station)
            if dd!=0:
                stationstr+=','+self.z[dd].station
            else:
                pass
        
        if medlst:
            medlst=np.array(medlst)
            medp=plt.plot(period,np.median(medlst,axis=0),lw=4,marker='x',ms=5,
                          color='k',mec='k')
            mlst.append(medp[0])
            stationlst.append('Median')
            meanp=plt.plot(period,np.mean(medlst,axis=0),lw=4,marker='x',ms=5,
                          color='g',mec='g')
            mlst.append(meanp[0])
            stationlst.append('Mean')
        ax.set_yscale('linear')
        ax.set_xscale('log')
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=-180,ymax=180)
        plt.grid(True)
        if len(stationlst)>1:
            plt.legend(tuple(mlst),tuple(stationlst),loc=0,markerscale=.4,borderaxespad=.05,
                   labelspacing=.1,handletextpad=.2)
            leg = plt.gca().get_legend()
            ltext=leg.get_texts()  # all the text.Text instance in the legend
            plt.setp(ltext, fontsize=10)    # the legend text fontsize
    
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Strike Angle (deg)',fontsize=10,fontweight='bold')
        plt.title('Strike Angle for '+stationstr,fontsize=12,fontweight='bold')
        plt.show() 
      
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                 self.z.station+'Strike.'+fmt),fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z.station+'Strike.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'Strike.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z.station+'Strike.'+fmt)
        elif save=='n':
            pass 
        
    def plotMinMaxPhase(self,save='n',fmt='pdf',fignum=1,thetar=0):
        """Plot the minimum and maximum phase of phase tensor. Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. """
        
        stationstr=self.z[0].station
        stationlst=[]
        for dd in range(len(self.z)):
            pt=self.z[dd].getPhaseTensor(thetar=thetar)
            period=self.z[dd].period
            minphi=pt.phiminang
            minphivar=pt.phiminangvar
            maxphi=pt.phimaxang
            maxphivar=pt.phimaxangvar
            
            fig=plt.figure(fignum,[6,8],dpi=100)
            ax=plt.subplot(1,1,1)
            ermin=plt.errorbar(period,minphi,marker='s',ms=4,mfc='None',
                               mec=self.mclst[dd][0],mew=1,ls='None',
                               yerr=minphivar,ecolor='k')
            ermax=plt.errorbar(period,maxphi,marker='o',ms=4,mfc='None',
                               mec=self.mclst[dd][0],mew=1,ls='None',
                               yerr=maxphivar,ecolor='k')
            stationlst.append(self.z[dd].station)
            if dd!=0:
                stationstr+=','+self.z[dd].station
            else:
                pass
            
        ax.set_xscale('log')
        ax.set_yscale('linear')
        plt.legend((ermin[0],ermax[0]),('Minimum','Maximum'),loc='udder left')
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=0,ymax=90)
        plt.grid(True)
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Phase (deg)',fontsize=10,fontweight='bold')
        plt.title('Min and Max Phase for '+stationstr,fontsize=12,
                  fontweight='bold')
        plt.show()
        
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                 self.z.station+'MMPhase.'+fmt),fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z.station+'MMPhase.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'MMPhase.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z.station+'MMPhase.'+fmt)
        elif save=='n':
            pass
    
    def plotAzimuth(self,save='n',fmt='pdf',fignum=1,thetar=0,dpi=300):
        """plot the azimuth of maximum phase. Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. """
        
        stationstr=self.z[0].station
        stationlst=[]
        mlst=[]
        if len(self.z)>0:
            medazm=[]
        plt.rcParams['font.size']=7
        fig=plt.figure(fignum,[6,8],dpi=dpi)
        ax=plt.subplot(1,1,1)
        for dd in range(len(self.z)):
            
            pt=self.z[dd].getPhaseTensor(thetar=thetar)
            #tipdict=self.z[dd].getTipper()
            az=90-np.array(pt.azimuth)
            azvar=np.array(pt.azimuthvar)
            medazm.append(az)
            #realarrow=tipdict['magreal']
            #realarrowvar=np.zeros(len(realarrow))+.00000001
            period=self.z[dd].period
            

            eraz=plt.errorbar(period,az,marker=self.mclst[dd][1],ms=.5,
                              mfc='None',mec=self.mclst[dd][0],mew=.5,
                              ls='None',yerr=azvar,ecolor=self.mclst[dd][0])
#            ertip=plt.errorbar(period,realarrow,marker='>',ms=4,mfc='None',mec='k',
#                               mew=1,ls='None',yerr=realarrowvar,ecolor='k')
            mlst.append(eraz[0])
            stationlst.append(self.z[dd].station)
            if dd!=0:
                stationstr=''
#                stationstr+=','+self.z[dd].station
            else:
                pass
        
        if medazm:
            medazm=np.array(medazm)
            medp=plt.plot(period,np.median(medazm,axis=0),lw=2,marker='x',ms=.5,
                          color='k',mec='k')
            mlst.append(medp[0])
            stationlst.append('Median')
            meanp=plt.plot(period,np.mean(medazm,axis=0),lw=2,marker='x',ms=.5,
                          color='g',mec='g')
            mlst.append(meanp[0])
            stationlst.append('Mean')
        ax.set_xscale('log')
        ax.set_yscale('linear')
        if len(stationlst)>1:
            plt.legend(tuple(mlst),tuple(stationlst),loc='lower left',
                       markerscale=2,borderaxespad=.05,
                       labelspacing=.1,handletextpad=.2,prop={'size':7},
                       ncol=int(np.ceil(dd/8)))
#            leg=plt.gca().get_legend()
#            ltext=leg.get_texts()  # all the text.Text instance in the legend
#            plt.setp(ltext, fontsize=8)    # the legend text fontsize
        
        ax.yaxis.set_major_locator(MultipleLocator(15))
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=-200,ymax=200)

        plt.grid(alpha=.25)
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Azimuth (deg. clockwise from N)',fontsize=10,
                   fontweight='bold')
        plt.title('Azimuth of Max Phase for '+ \
                                stationstr,fontsize=10,fontweight='bold')
        plt.show()
        
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,self.z.station+'Az.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                                   self.z.station+'Az.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'Az.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                                   self.z.station+'Az.'+fmt)
        elif save=='n':
            pass
        
    def plotSkew(self,save='n',fmt='pdf',fignum=1,thetar=0):
        """Plot the skew angle. Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. """
        
        stationstr=self.z[0].station
        stationlst=[]
        mlst=[]
        for dd in range(len(self.z)):
        
            pt=self.z[dd].getPhaseTensor(thetar=thetar)
            skew=pt.beta
            skewvar=pt.betavar
            period=self.z[dd].period
            
            fig=plt.figure(fignum,[6,8],dpi=100)
            ax=plt.subplot(1,1,1)
            erskew=plt.errorbar(period,skew,marker=self.mclst[dd][1],ms=4,
                                mfc='None',mec=self.mclst[dd][0],mew=1,
                                ls='None',yerr=skewvar,ecolor=self.mclst[dd][0])
            mlst.append(erskew[0])
            stationlst.append(self.z[dd].station)
            if dd!=0:
                stationstr+=','+self.z[dd].station
            else:
                pass
              
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.yaxis.set_major_locator(MultipleLocator(10))
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=-45,ymax=45)
        plt.grid(True)
        if len(stationlst)>1:
            plt.legend(tuple(mlst),tuple(stationlst),loc=0,markerscale=.4,
                       borderaxespad=.05,labelspacing=.1,handletextpad=.2)
            leg=plt.gca().get_legend()
            ltext=leg.get_texts()  # all the text.Text instance in the legend
            plt.setp(ltext, fontsize=10)    # the legend text fontsize
        
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Skew Angle (deg)',fontsize=10,fontweight='bold')
        plt.title('Skew Angle for '+stationstr,fontsize=12,fontweight='bold')
        plt.show()
        
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                 self.z.station+'Skew.'+fmt),fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z.station+'Skew.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'Skew.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z.station+'Skew.'+fmt)
        elif save=='n':
            pass
        
    def plotEllipticity(self,save='n',fmt='pdf',fignum=1,thetar=0):
        """Plot the ellipticity as (phimax-phimin)/phimax+phimin). Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. """
        
        plt.rcParams['font.size']=8
        plt.rcParams['figure.subplot.left']=.07
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.08
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.wspace']=.2
        plt.rcParams['figure.subplot.hspace']=.4
           
    
        stationstr=self.z[0].station
        stationlst=[]
        mlst=[]
        for dd in range(len(self.z)):    
       
            period=self.z[dd].period
            pt=self.z[dd].getPhaseTensor(thetar=thetar)
            ellipticity=pt.ellipticity
            ellipticityvar=pt.ellipticityvar
            
            fig=plt.figure(fignum,[6,8],dpi=100)
            ax=plt.subplot(1,1,1)
            erskew=plt.errorbar(period,ellipticity,marker=self.mclst[dd][1],
                                ms=4,mfc='None',mec=self.mclst[dd][0],mew=1,
                                ls='None',yerr=ellipticityvar,
                                ecolor=self.mclst[dd][0])
            stationlst.append(self.z[dd].station)
            mlst.append(erskew[0])
            if dd!=0:
                stationstr+=','+self.z[dd].station
            else:
                pass
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.yaxis.set_major_locator(MultipleLocator(.1))
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=0,ymax=1)
        #plt.yticks(range(10),np.arange(start=0,stop=1,step=.1))
        plt.grid(True)
        if len(stationlst)>1:
            plt.legend(tuple(mlst),tuple(stationlst),loc=0,markerscale=.4,borderaxespad=.05,
                   labelspacing=.1,handletextpad=.2)
            leg=plt.gca().get_legend()
            ltext=leg.get_texts()  # all the text.Text instance in the legend
            plt.setp(ltext, fontsize=10)    # the legend text fontsize
        
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('$\mathbf{\phi_{max}-\phi_{min}/\phi_{max}+\phi_{min}}$',
                   fontsize=12,fontweight='bold')
        plt.title('Ellipticity for '+stationstr,fontsize=12,
                  fontweight='bold')
        plt.show()
        
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                 self.z.station+'Ellip.'+fmt),fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z.station+'Ellip.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'Ellip.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z.station+'Ellip.'+fmt)
        elif save=='n':
            pass
        
    def plotAll(self,xspacing=6,esize=5,save='n',fmt='pdf',
                fignum=1,thetar=0):
        """plotAll will plot phase tensor, strike angle, min and max phase angle, 
        azimuth, skew, and ellipticity as subplots on one plot.  Save='y' if you 
        want to save the figure with path similar to input file or Save=savepath
        if you want to define the path yourself.  fmt can be pdf,eps,ps,png,
        svg. Fignum is the number of the figure."""
        
        stationstr=self.z[0].station
        stationlst=[]
        #Set plot parameters
        plt.rcParams['font.size']=8
        plt.rcParams['figure.subplot.left']=.07
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.08
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.wspace']=.2
        plt.rcParams['figure.subplot.hspace']=.4
        
        fs=8
        tfs=10
        #begin plotting
        fig=plt.figure(fignum,[8,10],dpi=150)
        for dd in range(len(self.z)):
            #get phase tensor infmtion
            pt=self.z[dd].getPhaseTensor(thetar=thetar)
            zinv=self.z[dd].getInvariants(thetar=thetar)
            tip=self.z[dd].getTipper(thetar=thetar)
            period=self.z[dd].period
            n=len(period)
            stationlst.append(self.z[dd].station)
            if dd!=0:
                stationstr+=','+self.z[dd].station
            else:
                pass

            
            #plotPhaseTensor
            ax1=plt.subplot(3,1,1,aspect='equal')
            for ii in range(n):
                #make sure the ellipses will be visable
                scaling=esize/pt.phimax[ii]
                eheight=pt.phimin[ii]*scaling
                ewidth=pt.phimax[ii]*scaling
                    
                #create an ellipse scaled by phimin and phimax and oriented along
                #the azimuth    
                ellip=Ellipse((xspacing*ii,0),width=ewidth,
                              height=eheight,
                              angle=pt.azimuth[ii])
                ax1.add_artist(ellip)
                
                if pt.phimin[ii]<0 or pt.phimin[ii]=='nan':
                    cvars=0
#                    print 'less than 0 or nan',cvars
                else:
                    cvars=(pt.phimin[ii]/(np.pi/2))
#                    print 'calculated ',cvars
                    if cvars>1.0:
                        cvars=0
#                        print 'greater than 1 ',cvars
#                print cvars
                    
                ellip.set_facecolor((1,1-cvars,.1))
            
            xticklabels=['%2.2g' % period[ii] for ii in np.arange(start=0,stop=n,
                         step=3)]
            plt.xlabel('Period (s)',fontsize=8,fontweight='bold')
            #plt.title('Phase Tensor Ellipses for '+stationstr,fontsize=14)
            plt.xticks(np.arange(start=0,stop=xspacing*n,step=3*xspacing),
                       xticklabels)
            ax1.set_ylim(-1*(xspacing+3),xspacing+3)
            ax1.set_xlim(-xspacing,n*xspacing+3)
            plt.grid()
            if dd==0:
                ax1cb=make_axes(ax1,shrink=.3,orientation='horizontal',pad=.30)
                cb=ColorbarBase(ax1cb[0],cmap=ptcmap,
                                norm=Normalize(vmin=min(pt.phiminang),
                                               vmax=max(pt.phiminang)),
                                orientation='horizontal')
                cb.set_label('Min Phase')
        
            if len(stationlst)>1:
                plt.legend(stationlst,loc=0,markerscale=.4,borderaxespad=.05,
                       labelspacing=.1,handletextpad=.2)
                leg=plt.gca().get_legend()
                ltext=leg.get_texts()  # all the text.Text instance in the legend
                plt.setp(ltext, fontsize=10)    # the legend text fontsize
        
            
            #plotStrikeAngle
            
            az=90-np.array(pt.azimuth)
            azvar=np.array(pt.azimuthvar)
            realarrow=tip.magreal
            realarrowvar=np.zeros(len(realarrow))+.00000001
            
            ax2=plt.subplot(3,2,3)
            erxy=plt.errorbar(period,zinv.strike,
                              marker=self.pcmlst[dd][0],ms=4,mfc='None',
                              mec=self.pcmlst[dd][1],mew=1,ls='None',
                              yerr=zinv.strikeerr,
                              ecolor=self.pcmlst[dd][1])
            eraz=plt.errorbar(period,az,marker=self.pcmlst[dd+1][0],ms=4,
                              mfc='None',mec=self.pcmlst[dd+1][1],mew=1,
                              ls='None',yerr=azvar,ecolor=self.pcmlst[dd+1][1])
            #ertip=plt.errorbar(period,realarrow,marker='>',ms=4,mfc='None',mec='k',
            #                   mew=1,ls='None',yerr=realarrowvar,ecolor='k')
            plt.legend((erxy[0],eraz[0]),('Strike','Azimuth'),loc=0,
                       markerscale=.2,borderaxespad=.01,labelspacing=.1,
                       handletextpad=.2)
            leg = plt.gca().get_legend()
            ltext  = leg.get_texts()  # all the text.Text instance in the legend
            plt.setp(ltext, fontsize=6)    # the legend text fontsize
    
            
            ax2.set_yscale('linear')
            ax2.set_xscale('log')
            plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                     xmin=10**(np.floor(np.log10(period[0]))))
            plt.ylim(ymin=-200,ymax=200)
            plt.grid(True)
            #plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
            plt.ylabel('Angle (deg)',fontsize=fs,fontweight='bold')
            plt.title('Strike Angle, Azimuth',fontsize=tfs,
                      fontweight='bold')
            
            #plotMinMaxPhase
            
            minphi=pt.phiminang
            minphivar=pt.phiminangvar
            maxphi=pt.phimaxang
            maxphivar=pt.phimaxangvar
    
            ax3=plt.subplot(3,2,4,sharex=ax2)
            ermin=plt.errorbar(period,minphi,marker=self.pcmlst[dd][0],ms=4,
                               mfc='None',mec=self.pcmlst[dd][1],mew=1,ls='None',
                               yerr=minphivar,ecolor=self.pcmlst[dd][1])
            ermax=plt.errorbar(period,maxphi,marker=self.pcmlst[dd+1][0],ms=4,
                               mfc='None',mec=self.pcmlst[dd+1][1],mew=1,
                               ls='None',yerr=maxphivar,
                               ecolor=self.pcmlst[dd+1][1])
            ax3.set_xscale('log')
            ax3.set_yscale('linear')
            plt.legend((ermin[0],ermax[0]),('$\phi_{min}$','$\phi_{max}$'),
                       loc='upper left',markerscale=.2,borderaxespad=.01,
                       labelspacing=.1,handletextpad=.2)
            leg = plt.gca().get_legend()
            ltext  = leg.get_texts()  # all the text.Text instance in the legend
            plt.setp(ltext, fontsize=6.5)    # the legend text fontsize
            plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                     xmin=10**(np.floor(np.log10(period[0]))))
            plt.ylim(ymin=0,ymax=90)
            plt.grid(True)
            #plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
            plt.ylabel('Phase (deg)',fontsize=fs,fontweight='bold')
            plt.title('$\mathbf{\phi_{min}}$ and $\mathbf{\phi_{max}}$',fontsize=tfs,
                      fontweight='bold')

            
            #plotSkew
            
            skew=pt.beta
            skewvar=pt.betavar
    
            ax5=plt.subplot(3,2,5,sharex=ax2)
            erskew=plt.errorbar(period,skew,marker=self.pcmlst[dd][0],ms=4,
                                mfc='None',mec=self.pcmlst[dd][1],mew=1,
                                ls='None',yerr=skewvar,
                                ecolor=self.pcmlst[dd][1])
            ax5.set_xscale('log')
            ax5.set_yscale('linear')
            ax5.yaxis.set_major_locator(MultipleLocator(10))
            plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),xmin=10**(
                                np.floor(np.log10(period[0]))))
            plt.ylim(ymin=-45,ymax=45)
            plt.grid(True)
            plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
            plt.ylabel('Skew Angle (deg)',fontsize=fs,fontweight='bold')
            plt.title('Skew Angle',fontsize=tfs,fontweight='bold')
            
            #plotEllipticity
            
            ellipticity=pt.ellipticity
            ellipticityvar=pt.ellipticityvar
    
            ax6=plt.subplot(3,2,6,sharex=ax2)
            erskew=plt.errorbar(period,ellipticity,marker=self.pcmlst[dd][0],
                                ms=4,mfc='None',mec=self.pcmlst[dd][1],mew=1,
                                ls='None',yerr=ellipticityvar,
                                ecolor=self.pcmlst[dd][1])
            ax6.set_xscale('log')
            ax6.set_yscale('linear')
            ax6.yaxis.set_major_locator(MultipleLocator(.1))
            plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                     xmin=10**(np.floor(np.log10(period[0]))))
            plt.ylim(ymin=0,ymax=1)
            #plt.yticks(range(10),np.arange(start=0,stop=1,step=.1))
            plt.grid(True)
            plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
            plt.ylabel('$\mathbf{\phi_{max}-\phi_{min}/\phi_{max}+\phi_{min}}$',
                       fontsize=fs,fontweight='bold')
            plt.title('Ellipticity',fontsize=tfs,fontweight='bold')
            #plt.suptitle(self.z.station,fontsize=tfs,fontweight='bold')
        plt.suptitle('Phase Tensor Elements for: '+stationstr,fontsize=12,
                     fontweight='bold')
            
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                     self.z[0].station+'All.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z[0].station+'All.'+fmt)
            plt.close()
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z[0].station+'All.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z[0].station+'All.'+fmt)
            plt.close()
        elif save=='n':
            pass

    def plotResPhase(self,ffactor=1,save='n',fmt='pdf',fignum=1,thetar=0):
        """plotResPhase(self,df=100.) will plot the resistivity and phase for
        all impedance tensor polarizations.  2 plots, one containing xy,yx 
        polarizations, the other xx,yy.  Save='y' if you want to save the figure
        with path similar to input file or Save=savepathif you want to define 
        the path yourself.  fmt can be pdf,eps,ps,png, svg. Fignum is the 
        number of the figure."""
        
        rp=self.z[0].getResPhase(ffactor=ffactor)
        period=self.z[0].period
        #set some plot parameters
        plt.rcParams['font.size']=10
        plt.rcParams['figure.subplot.left']=.13
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.07
        plt.rcParams['figure.subplot.top']=.96
        plt.rcParams['figure.subplot.wspace']=.25
        plt.rcParams['figure.subplot.hspace']=.20
        
            
        #make figure for xy,yx components
        fig=plt.figure(fignum,[6,8],dpi=100)
        ax=plt.subplot(2,1,1)
        ax.yaxis.set_label_coords(-.1, 0.5)
        erxy=plt.errorbar(period,rp.resxy,marker='s',ms=4,mfc='None',
                          mec='b',mew=1,ls='None',yerr=rp.resxyerr,
                          ecolor='b')
        eryx=plt.errorbar(period,rp.resyx,marker='s',ms=4,mfc='None',
                          mec='r',mew=1,ls='None',yerr=rp.resyxerr,
                          ecolor='r')
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Addarent Resistivity ($\mathbf{\Omega \cdot m}$)',
                   fontsize=10,fontweight='bold')
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlim(xmin=10**(np.ceil(np.log10(period[0]))),
                 xmax=10**(np.floor(np.log10(period[-1]))))
        plt.grid(True)
        plt.legend((erxy[0],eryx[0]),('$E_x/B_y$','$E_y/B_x$'),loc=0,
                markerscale=1,borderaxespad=.01,labelspacing=.1,
                handletextpad=.2)
        
        ax2=plt.subplot(2,1,2,sharex=ax)
        ax2.yaxis.set_label_coords(-.1, 0.5)
        plt.errorbar(period,rp.phasexy,marker='s',ms=4,mfc='None',
                     mec='b',mew=1,ls='None',yerr=rp.phasexyerr,
                     ecolor='b')
        plt.errorbar(period,rp.phaseyx,marker='s',ms=4,mfc='None',
                     mec='r',mew=1,ls='None',yerr=rp.phaseyxerr,
                     ecolor='r')
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Imepdance Phase (deg)',fontsize=10,fontweight='bold')
        ax2.set_xscale('log')
        plt.xlim(xmin=10**(np.ceil(np.log10(period[0]))),
                 xmax=10**(np.floor(np.log10(period[-1]))))
        if min(rp.phasexy)<0 or min(rp.phaseyx)<0:
            ymin=-90
            ax2.yaxis.set_major_locator(MultipleLocator(30))
        else:
            ymin=0
        plt.ylim(ymin=ymin,ymax=90)
        plt.grid(True)
        plt.show()
        
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig.savefig(os.path.join(self.savepath,
                                self.z.station+'Resxy.'+fmt), fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z.station+'Resxy.'+fmt)
        elif len(save)>1:
            fig.savefig(os.path.join(save,self.z.station+'Resxy.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z.station+'Resxy.'+fmt)
        elif save=='n':
            pass
        
        #plot figure for xx,yy components
        fig2=plt.figure(fignum+1,[6,8],dpi=100)
        ax3=plt.subplot(2,1,1)
        ax3.yaxis.set_label_coords(-.1, 0.5)
        erxy2=plt.errorbar(period,rp.resxx,marker='s',ms=4,mfc='None',
                           mec='b',mew=1,ls='None',yerr=rp.resxxerr,
                           ecolor='b')
        eryx2=plt.errorbar(period,rp.resyy,marker='s',ms=4,mfc='None',
                           mec='r',mew=1,ls='None',yerr=rp.resyyerr,
                           ecolor='r')
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Addarent Resistivity ($\mathbf{\Omega \cdot m}$)',
                   fontsize=10,fontweight='bold')
        ax3.set_yscale('log')
        ax3.set_xscale('log')
        plt.xlim(xmin=10**(np.ceil(np.log10(period[0]))),
                 xmax=10**(np.floor(np.log10(period[-1]))))
        plt.grid(True)
        plt.legend((erxy2[0],eryx2[0]),('$E_x/B_x$','$E_y/B_y$'),
                    loc=0, markerscale=1,borderaxespad=.01,labelspacing=.1,
                    handletextpad=.2)
        
        ax3=plt.subplot(2,1,2,sharex=ax3)
        ax3.yaxis.set_label_coords(-.1, 0.5)
        plt.errorbar(period,rp.phasexx,marker='s',ms=4,mfc='None',
                     mec='b',mew=1,ls='None',yerr=rp.phasexxerr,
                     ecolor='b')
        plt.errorbar(period,rp.phaseyy,marker='s',ms=4,mfc='None',
                     mec='r',mew=1,ls='None',yerr=rp.phaseyyerr,
                     ecolor='r')
        plt.xlabel('Period (s)',fontsize=10,fontweight='bold')
        plt.ylabel('Imepdance Phase (deg)',fontsize=10,fontweight='bold')
        ax3.set_xscale('log')
        plt.xlim(xmin=10**(np.ceil(np.log10(period[0]))),
                 xmax=10**(np.floor(np.log10(period[-1]))))
        if min(rp.phasexx)<0 or min(rp.phaseyy)<0:
            ymin=-90
            ax3.yaxis.set_major_locator(MultipleLocator(30))
        else:
            ymin=0
        plt.ylim(ymin=ymin,ymax=90)
        plt.grid(True)
        plt.show()
        
        if save=='y':
            if not os.path.exists(self.savepath):
                os.mkdir(self.savepath)
                print 'Made Directory: '+ self.savepath
            fig2.savefig(os.path.join(self.savepath,
                                self.z.station+'Resxx.'+fmt),fmt=fmt)
            print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z.station+'Resxx.'+fmt)
        elif len(save)>1:
            fig2.savefig(os.path.join(save,self.z.station+'Resxx.'+fmt),
                        fmt=fmt)
            print 'Saved figure to: '+os.path.join(save,
                                               self.z.station+'Resxx.'+fmt)
        elif save=='n':
            pass
