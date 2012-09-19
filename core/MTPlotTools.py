'''
This module contains functions to plot different MT things.
'''

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import BIRRPTools as brp
import os
import MTtools as mt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
import Z
import LatLongUTMconversion as utm2ll
from matplotlib.patches import Ellipse,Rectangle,Arrow
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.colorbar import *
import matplotlib.gridspec as gridspec

#==============================================================================
# Make some color maps for plotting
#==============================================================================
#yellow to red
ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

#blue to yellow to red
skcmapdict={'red':((0.0,0.0,0.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,1.0,0.0),(.5,1.0,0.0),(.5,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,1.0),(.5,0.0,1.0),(0.5,0.1,0.1),(1.0,0.1,0.1))}
skcmap=LinearSegmentedColormap('skcmap',skcmapdict,256)

#blue to white to red
skcmapdict2={'red':((0.0,0.0,0.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,1.0,0.0),(.5,1.0,0.0),(.5,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,1.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,0.0,0.0))}
skcmap2=LinearSegmentedColormap('skcmap2',skcmapdict2,256)


#white to blue
ptcmapdict3={'red':((0.0,1.0,1.0),(1.0,0.0,0.0)),
            'green':((0.0,1.0,1.0),(1.0,0.0,0.0)),
            'blue':((0.0,1.0,1.0),(1.0,1.0,1.0))}
ptcmap3=LinearSegmentedColormap('ptcmap3',ptcmapdict3,256)

#red to blue
rtcmapdict={'red':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'green':((0.0,0.0,0.0),(1.0,0.0,0.0)),
            'blue':((0.0,1.0,0.0),(1.0,1.0,0.0))}
rtcmap=LinearSegmentedColormap('rtcmap',rtcmapdict,256)

ckdict={'phiminang':'$\Phi_{min}$ (deg)','phimin':'$\Phi_{min}$ (deg)',
        'phimaxang':'$\Phi_{max}$ (deg)','phimax':'$\Phi_{max}$ (deg)',
        'phidet':'Det{$\Phi$} (deg)','beta':'Skew (deg)',
        'ellipticity':'Ellipticity'}

zonedict=dict([(a,ii) for ii,a in enumerate(['a','b','c','d','e','f','g','h',
               'i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x',
               'y','z'])])

#==============================================================================
# Plotting tools
#==============================================================================


def plotcoh(filename,fignum=1,savefigfilename=None,dpi=None,format=None,
            orientation=None):
    """Will plot coherence output from birrp_bbconvert.  If you want to save the
    plot do so using the interactive gui.  If coh3=0 only two plots.
    
     Inputs:
        filename = filename containing coherence (.coh)
        fignum = figure number
        savefigfilenames = supply filenames to save figures to if desired
        dpi = figure resolution
        format = file type of saved figure pdf,svg,eps...
        orientation = orientation of figure on A4 paper
    Outputs:
        none"""
    
  
    station,period,freq,coh1,zcoh1,coh2,zcoh2,coh3,zcoh3=brp.readcoh(filename)
    nd=len(period)
    
    #===============================================================================
    # Plot coherence
    #===============================================================================
    if coh3[0]==0 and coh3[len(coh3)-1]==0 and zcoh3[0]==0 and zcoh3[len(coh3)-1]==0:
        plt.rcParams['font.size']=10
        plt.rcParams['figure.subplot.left']=.09
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.075
        plt.rcParams['figure.subplot.top']=.92
        plt.rcParams['figure.subplot.wspace']=.25
        plt.rcParams['figure.subplot.hspace']=.35
        #plt.rcParams['font.family']='helvetica'
        
        fig1=plt.figure(fignum, [8,6], 100)
        ax1=plt.subplot(1,2,1)
        ax1.semilogx(period,coh1,'k',linewidth=2)
        ax1.semilogx(period,zcoh1,'k--', linewidth=2,dashes=[4,1,4,1])
        ax1.grid(True)
        ax1.legend(('Normal','Zero'),'lower left',shadow=False,borderaxespad=.50)
        ax1.set_title('$\mathbf{E_x/B_y}$',fontsize=12,fontweight='bold')
        ax1.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax1.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax1.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax1.set_ylim(ymin=0,ymax=1)
        
        ax2=plt.subplot(1,2,2)
        ax2.semilogx(period,coh2,'k',linewidth=2)
        ax2.semilogx(period,zcoh2,'k--', linewidth=2,dashes=[4,1,4,1])
        ax2.grid(True)
        ax2.legend(('Normal','Zero'),'lower left',shadow=False,borderaxespad=.50)
        ax2.set_title('$\mathbf{E_y/B_x}$',fontsize=12,fontweight='bold')
        ax2.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax2.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax2.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax2.set_ylim(ymin=0,ymax=1)
        
        plt.suptitle('Coherence for Station: '+station,fontsize=12,fontweight='bold')
        plt.show()
        if savefigfilename!=None:
            if dpi==None:
                dpi=100
            if format==None:
                format='pdf'
            if orientation==None:
                orientation='landscape'
            plt.savefig(savefigfilename,dpi=dpi,format=format,
                        orientation=orientation)
            plt.close(fig1)
    
    else:
        plt.rcParams['font.size']=10
        plt.rcParams['figure.subplot.left']=.07
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.075
        plt.rcParams['figure.subplot.top']=.90
        plt.rcParams['figure.subplot.wspace']=.25
        plt.rcParams['figure.subplot.hspace']=.35
        #plt.rcParams['font.family']='helvetica'
        
        fig1=plt.figure(fignum, [10,6], 100)
        ax1=plt.subplot(1,3,1)
        ax1.semilogx(period,coh1,'k',linewidth=2)
        ax1.semilogx(period,zcoh1,'k--', linewidth=2,dashes=[4,1,4,1])
        ax1.grid(True)
        ax1.legend(('Normal','Zero'),'best',shadow=False,borderaxespad=.50)
        ax1.set_title('$\mathbf{E_x/B_y}$',fontsize=12,fontweight='bold')
        ax1.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax1.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax1.set_xscale('log')
        ax1.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax1.set_ylim(ymin=0,ymax=1)
        
        ax2=plt.subplot(1,3,2)
        ax2.semilogx(period,coh2,'k',linewidth=2)
        ax2.semilogx(period,zcoh2,'k--', linewidth=2,dashes=[4,1,4,1])
        ax2.grid(True)
        ax2.legend(('Normal','Zero'),'best',shadow=False,borderaxespad=.50)
        ax2.set_title('$\mathbf{E_y/B_x}$',fontsize=12,fontweight='bold')
        ax2.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax2.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax2.set_xscale('log')
        ax2.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax2.set_ylim(ymin=0,ymax=1)
        
        ax3=plt.subplot(1,3,3)
        ax3.semilogx(period,coh3,'k',linewidth=2)
        ax3.semilogx(period,zcoh3,'k--', linewidth=2,dashes=[4,1,4,1])
        ax3.grid(True)
        ax3.legend(('Normal','Zero'),'best',shadow=False,borderaxespad=.50)
        ax3.set_title('$\mathbf{E_x/B_z}$',fontsize=12,fontweight='bold')
        ax3.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax3.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax3.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax3.set_ylim(ymin=0,ymax=1)
        
        plt.suptitle('Coherence for Station: '+station,fontsize=12,fontweight='bold')
        plt.show()
        if savefigfilename!=None:
            if dpi==None:
                dpi=100
            if format==None:
                format='pdf'
            if orientation==None:
                orientation='landscape'
            plt.savefig(savefigfilename,dpi=dpi,format=format,
                        orientation=orientation)
            plt.clf()
            plt.close(fig1)

    
def plotResPhase(filename,fignum=1,df=100,ffactor=1,plotnum=1,title=None,
                 savefigfilename=None,dpi=None,format=None,orientation=None,
                 rotz=0):
    """
    plotResPhase(filename,fignum) will plot the apparent resistivity and 
    phase for TE and TM modes 
    from a .dat file produced by writedat.  If you want to save the plot 
    use the save button on top left.
    
    Inputs:
        filename = filename containing impedance (.edi) or resistivity and phase
                   information (.dat)
        fignum = figure number
        df = sampling frequency (Hz)
        ffactor = fudge factor for computing resistivity from impedances
        thetar = rotation angle of impedance tensor (deg or radians)
        plotnum = 1 for just Ex/By and Ey/Bx
                  2 for all 4 components
                  3 for off diagonal plus the determinant
        title = title of plot
        savefigfilename = supply filename to save figure to if desired
        dpi = figure resolution
        format = file type of saved figure pdf,svg,eps...
        orientation = orientation of figure on A4 paper
        
    Outputs:
        none
        
    """
    
    if dpi==None:
        dpi=100
    #read in .dat file
    if filename.find('dat',-4,len(filename))>=0:
        print 'Reading .dat file'
        station,period,resdat,phasedat=brp.readdat(filename)
#        if period[0]>period[-1]:
#            res=res.reverse()
        resxy=resdat[0][1][:]
        resxyerr=resdat[1][1][:]
        resyx=resdat[0][2][:]
        resyxerr=resdat[1][2][:]
        resxx=resdat[0][0][:]
        resxxerr=resdat[1][0][:]
        resyy=resdat[0][3][:]
        resyyerr=resdat[1][3][:]
        
        phasexy=phasedat[0][1][:]
        phasexyerr=phasedat[1][1][:]
        phaseyx=phasedat[0][2][:]
        phaseyxerr=phasedat[1][2][:]
        phasexx=phasedat[0][0][:]
        phasexxerr=phasedat[1][0][:]
        phaseyy=phasedat[0][3][:]
        phaseyyerr=phasedat[1][3][:] 
        
        rpdict={}
        rpdict['resxx']=np.array(resxx)
        rpdict['resxy']=np.array(resxy)
        rpdict['resyx']=np.array(resyx)
        rpdict['resyy']=np.array(resyy)
        rpdict['resxxerr']=np.array(resxxerr)
        rpdict['resxyerr']=np.array(resxyerr)
        rpdict['resyxerr']=np.array(resyxerr)
        rpdict['resyyerr']=np.array(resyyerr)
        rpdict['phasexx']=np.array(phasexx)
        rpdict['phasexy']=np.array(phasexy)
        rpdict['phaseyx']=np.array(phaseyx)
        rpdict['phaseyy']=np.array(phaseyy)
        rpdict['phasexxerr']=np.array(phasexxerr)
        rpdict['phasexyerr']=np.array(phasexyerr)
        rpdict['phaseyxerr']=np.array(phaseyxerr)
        rpdict['phaseyyerr']=np.array(phaseyyerr)

            
        
    elif filename.find('edi',-4,len(filename))>=0:
        print 'Reading .edi file'
        impz=Z.Z(filename)
        station=impz.station
        period=impz.period
        rp=impz.getResPhase(thetar=rotz)
    else:
        raise ValueError('Could not read file: '+filename)
    
 
    plt.rcParams['font.size']=10
    plt.rcParams['figure.subplot.left']=.13
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.1
    plt.rcParams['figure.subplot.top']=.95
    plt.rcParams['figure.subplot.wspace']=.25
    plt.rcParams['figure.subplot.hspace']=.05
    #plt.rcParams['font.family']='helvetica'
    fontdict={'size':14,'weight':'bold'}
    
    gs=gridspec.GridSpec(2,2,height_ratios=[2,1.5],hspace=.05)
    
        
    #make figure for xy,yx components
    if plotnum==1 or plotnum==3: 
        fig=plt.figure(fignum,[8,10],dpi=dpi)
        gs.update(hspace=.05,wspace=.15,left=.1)
    elif plotnum==2:
        fig=plt.figure(fignum,[10,10],dpi=dpi)
        gs.update(hspace=.05,wspace=.15,left=.07)
    
    #---------plot the apparent resistivity-----------------------------------
    if plotnum==1  or plotnum==3:
        ax=plt.subplot(gs[0,:])
        ax2=plt.subplot(gs[1,:],sharex=ax)
        ax.yaxis.set_label_coords(-.055, 0.5)
        ax2.yaxis.set_label_coords(-.055, 0.5)
    elif plotnum==2:
        ax=plt.subplot(gs[0,0])
        ax2=plt.subplot(gs[1,0],sharex=ax)
        ax.yaxis.set_label_coords(-.075, 0.5)
        ax2.yaxis.set_label_coords(-.075, 0.5)
    
    
    erxy=ax.errorbar(period,rp.resxy,marker='s',ms=4,mfc='None',mec='b',
                      mew=1,ls='None',yerr=rp.resxyerr,ecolor='b')
    eryx=ax.errorbar(period,rp.resyx,marker='o',ms=4,mfc='None',mec='r',
                      mew=1,ls='None',yerr=rp.resyxerr,ecolor='r')
    #ax.set_xlabel('Period (s)',fontdict=fontdict)
    pylab.setp( ax.get_xticklabels(), visible=False)
    ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
               fontdict=fontdict)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
             xmax=10**(np.ceil(np.log10(period[-1]))))
    ax.grid(True)
    ax.legend((erxy[0],eryx[0]),('$E_x/B_y$','$E_y/B_x$'),loc=3,
                markerscale=1,borderaxespad=.01,labelspacing=.07,
                handletextpad=.2,borderpad=.02)
    
    #-----Plot the phase----------------------------------------------------
    
    ax2.errorbar(period,rp.phasexy,marker='s',ms=4,mfc='None',mec='b',
                 mew=1,ls='None',yerr=rp.phasexyerr,ecolor='b')
    ax2.errorbar(period,np.array(rp.phaseyx)+180,marker='o',ms=4,
                 mfc='None',mec='r',mew=1,ls='None',yerr=rp.phaseyxerr,
                 ecolor='r')
    ax2.set_xlabel('Period (s)',fontdict)
    ax2.set_ylabel('Phase (deg)',fontdict)
    ax2.set_xscale('log')
    #ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
    #         xmax=10**(np.ceil(np.log10(period[-1]))))
    #check the phase to see if any point are outside of [0:90]    
    if min(rp.phasexy)<0 or min(rp.phaseyx+180)<0:
        pymin=min([min(rp.phasexy),min(rp.phaseyx+180)])
        if pymin>0:
            pymin=0
    else:
        pymin=0
    
    if max(rp.phasexy)>90 or max(rp.phaseyx+180)>90:
        pymax=min([max(rp.phasexy),max(rp.phaseyx+180)])
        if pymax<91:
            pymax=90
    else:
        pymax=90
    
    ax2.set_ylim(ymin=pymin,ymax=pymax)        
    ax2.yaxis.set_major_locator(MultipleLocator(30))
    ax2.yaxis.set_minor_locator(MultipleLocator(1))
    ax2.grid(True)
    
    if plotnum==2:
        #---------plot the apparent resistivity-----------------------------------
        ax3=plt.subplot(gs[0,1])
        ax3.yaxis.set_label_coords(-.1, 0.5)
        erxx=ax3.errorbar(period,rp.resxx,marker='s',ms=4,mfc='None',mec='b',
                          mew=1,ls='None',yerr=rp.resxxerr,ecolor='b')
        eryy=ax3.errorbar(period,rp.resyy,marker='o',ms=4,mfc='None',mec='r',
                          mew=1,ls='None',yerr=rp.resyyerr,ecolor='r')
        #ax.set_xlabel('Period (s)',fontdict=fontdict)
    #    ax3.set_ylabel('Apparent Resistivity ($\mathbf{\Omega \cdot m}$)',
    #               fontdict=fontdict)
        ax3.set_yscale('log')
        ax3.set_xscale('log')
        pylab.setp( ax3.get_xticklabels(), visible=False)
        ax3.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
                 xmax=10**(np.ceil(np.log10(period[-1]))))
        ax3.grid(True)
        ax3.legend((erxx[0],eryy[0]),('$E_x/B_x$','$E_y/B_y$'),loc=3,
                    markerscale=1,borderaxespad=.01,labelspacing=.07,
                    handletextpad=.2,borderpad=.02)
        
        #-----Plot the phase----------------------------------------------------
        ax4=plt.subplot(gs[1,1],sharex=ax3)
        
        ax4.yaxis.set_label_coords(-.1, 0.5)
        ax4.errorbar(period,rp.phasexx,marker='s',ms=4,mfc='None',mec='b',
                     mew=1,ls='None',yerr=rp.phasexxerr,ecolor='b')
        ax4.errorbar(period,np.array(rp.phaseyy),marker='o',ms=4,
                     mfc='None',mec='r',mew=1,ls='None',yerr=rp.phaseyyerr,
                     ecolor='r')
        ax4.set_xlabel('Period (s)',fontdict)
        #ax4.set_ylabel('Imepdance Phase (deg)',fontdict)
        ax4.set_xscale('log')
        #ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
        #         xmax=10**(np.ceil(np.log10(period[-1]))))
        ax4.set_ylim(ymin=-180,ymax=180)        
        ax4.yaxis.set_major_locator(MultipleLocator(30))
        ax4.yaxis.set_minor_locator(MultipleLocator(5))
        ax4.grid(True)
        
    if plotnum==3:
            
        
        erdet=ax.errorbar(period,rp.resdet,marker='d',ms=4,mfc='None',mec='g',
                      mew=1,ls='None',yerr=rp.resdeterr,ecolor='g')
    
        ax.legend((erxy[0],eryx[0],erdet[0]),('$E_x/B_y$','$E_y/B_x$',
                  '$\det(\mathbf{\hat{Z}})$'),loc=3,markerscale=1,borderaxespad=.01,
                  labelspacing=.07,handletextpad=.2,borderpad=.02)
        
        #-----Plot the phase----------------------------------------------------
        
        ax2.errorbar(period,rp.phasedet,marker='d',ms=4,mfc='None',mec='g',
                     mew=1,ls='None',yerr=rp.phasedeterr,ecolor='g')
    
    
    #make title and show
    if title!=None:
        plt.suptitle(title,fontsize=14,fontweight='bold')
    else:
        plt.suptitle(station,fontsize=14,fontweight='bold')
    plt.show()
    if savefigfilename!=None:
            if dpi==None:
                dpi=100
            if format==None:
                format='pdf'
            if orientation==None:
                orientation='landscape'
            fig.savefig(savefigfilename,dpi=dpi,format=format,
                        orientation=orientation)
            plt.clf()
            plt.close(fig)

def resPhasePlots(filenamelst,plottype=1,ffactor=1,kwdict=None,plotnum=1,
                  ylim=[0,3],fignum=1,rotz=0):
    """resPhasePlots will plot multiple responses given full path filenames.  
    Can input key word list dictionary if parameters for each plot are different
    or just single values.
    
    Input:
        for plottype==1 (plots a new figure for each station with resitivity and
                         phase):
            plotnum = 1,2,3
                        1 -> will plot just the off diagonal components
                        2 -> will plot all 4 components one subplot for each 
                             diagonal pair
                        3 -> will plot the determinant of Z only
            kwdict = keywords for plotResPhase, can input as a list of 
                     keywords that is the same length as filenamelst
            
            
            ylim = min and max powers of 10 to plot apparent resistivity
        \t 1 => each filename has its own figure window \n
        \t 2 => plotted in one figure with subplots for each filename \n
        \t 3 => plotted onto one figure and one plot"""
    
    plt.rcParams['font.size']=12    
    
    if plottype==1:
        if kwdict==None:
            for ii,filename in enumerate(filenamelst):
                station=os.path.basename(filename)[0:-4]
                print filename
                plotResPhase(filename,fignum=ii+1,ffactor=ffactor,
                             plotnum=plotnum,title=station,rotz=rotz)
        else:
            for ii,filename in enumerate(filenamelst):
                station=os.path.basename(filename)[0:-4]
                if type(kwdict) is list:
                    plotResPhase(filename,fignum=ii+1,plotnum=plotnum,
                                 title=station,rotz=rotz,**kwdict[ii])
                elif type(kwdict) is dict:
                     plotResPhase(filename,fignum=ii+1,plotnum=plotnum,
                                 title=station,rotz=rotz,**kwdict)
    
    if plottype==2:
    
        nsp=len(filenamelst)
        #get the number of rows with maximum of 8 columns        
        nrows=int(np.ceil(nsp/8.))
        if nsp<8:
            ncols=nsp
        else:
            ncols=8
        
        fig=plt.figure(fignum,[24,14])
        gs2=gridspec.GridSpec(nrows,1,hspace=.2,left=.06,right=.99)
        
        fontdict={'size':12,'weight':'bold'}
        for rr in range(nrows):
            g1=gridspec.GridSpecFromSubplotSpec(2,ncols,subplot_spec=gs2[rr])
            g1.set_height_ratios([2,1])
            g1._wspace=.08
            g1._hspace=.02
            for cc in range(ncols):
                try:
                    impz=Z.Z(filenamelst[(rr*ncols)+cc])
                    station=impz.station
                    period=impz.period
                    rp=impz.getResPhase(thetar=rotz)
                    
                    axr=plt.Subplot(fig,g1[0,cc])
                    ax=fig.add_subplot(axr)
                    ax.set_title(station,fontdict=fontdict)
                    ax.yaxis.set_label_coords(-.15, 0.5)
                    
                    #--Plot Res------------------------------------------------
                    erxy=ax.errorbar(period,rp.resxy,marker='s',ms=4,
                                     mfc='None',mec='b',mew=1,ls='None',
                                     yerr=rp.resxyerr,ecolor='b')
                    eryx=ax.errorbar(period,rp.resyx,marker='o',ms=4,
                                     mfc='None',mec='r',mew=1,ls='None',
                                     yerr=rp.resyxerr,ecolor='r')
                    if cc==0:
                        ax.set_ylabel('App. Res.'+
                                      '($\mathbf{\Omega \cdot m}$)',
                                       fontdict=fontdict)
                    ax.set_yscale('log')
                    ax.set_xscale('log')
                    ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
                             xmax=10**(np.ceil(np.log10(period[-1]))))
                    ax.grid(True)
                    ax.set_xticklabels(['' for tt in range(6)])
                    ax.set_ylim(ymin=10**ylim[0],ymax=10**ylim[1])
                    
                    #--Plot Phase----------------------------------------------
                    axp=plt.Subplot(fig,g1[1,cc])
                    ax2=fig.add_subplot(axp,sharex=ax)
            
                    ax2.yaxis.set_label_coords(-.15, 0.5)
                    ax2.errorbar(period,rp.phasexy,marker='s',ms=4,
                                 mfc='None',mec='b',mew=1,ls='None',
                                 yerr=rp.phasexyerr,ecolor='b')
                    ax2.errorbar(period,np.array(rp.phaseyx)+180,
                                 marker='o',ms=4,mfc='None',mec='r',mew=1,
                                 ls='None',yerr=rp.phaseyxerr,ecolor='r')
                    ax2.set_xlabel('Period (s)',fontdict)
                    xtl=[str(int(ss)) for ss in np.arange(np.floor(np.log10(period[0])),
                                  np.ceil(np.log10(period[-1])))]
                    ax2.set_xticklabels(xtl)
                    if cc==0:
                        ax2.set_ylabel('Phase (deg)',
                                       fontdict=fontdict)
                    if rr==nrows-1:
                        ax2.set_xlabel('Period (s)',fontdict=fontdict)
                    ax2.set_xscale('log')
                    #ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
                    #         xmax=10**(np.ceil(np.log10(period[-1]))))
                    ax2.set_ylim(ymin=0,ymax=90)        
                    ax2.yaxis.set_major_locator(MultipleLocator(15))
                    ax2.yaxis.set_minor_locator(MultipleLocator(5))
                    ax2.grid(True)
                    
                    if cc>0:
                        ax.set_yticklabels(['' for tt in range(6)])
                        ax2.set_yticklabels(['' for tt in range(6)])
                except IndexError:
                    break
            
        
    #plot on same plot        
    elif plottype==3:
        
        colorlst=[['b','r'],['c','m'],['g','y'],['b','r'],['c','m'],['g','y']]*5
        markerlst=[['s','s'],['D','D'],['o','o'],['x','x'],['v','v'],['^','^'],
                   ['p','p'],['+','+']]*5
        
        plt.rcParams['font.size']=12
        plt.rcParams['figure.subplot.left']=.08
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.92
        plt.rcParams['figure.subplot.wspace']=.25
        plt.rcParams['figure.subplot.hspace']=.20
        #plt.rcParams['font.family']='helvetica'
        
        #make figure
        fig=plt.figure(fignum,[10,6],dpi=100)
        if plotnum==2:
            fig2=plt.figure(fignum+1,[10,6],dpi=100)
            legendlst2=[]
            legendlabellst2=[]
        elif plotnum==3:
            fig3=plt.figure(fignum+2,[10,6],dpi=100)
            legendlst3=[]
            legendlabellst3=[]
        #read in filenames and plot on subplots
        mspot=0
        legendlst=[]
        legendlabellst=[]
        periodmin=[]
        periodmax=[]
        for ii,filename in enumerate(filenamelst):
            
            xycolor=colorlst[mspot][0]
            yxcolor=colorlst[mspot][1]
            
            xymarker=markerlst[mspot][0]
            yxmarker=markerlst[mspot][1]
            mspot+=1
            
            station=os.path.basename(filename)[0:-4]
            #read in .dat file
            if filename.find('dat',-4,len(filename))>=0:
                print 'Reading .dat file'
                ofil,period,resdat,phasedat=mt.readdat(filename)
                resxy=resdat[0][1][:]
                resxyerr=resdat[1][1][:]
                resyx=resdat[0][2][:]
                resyxerr=resdat[1][2][:]
                resxx=resdat[0][0][:]
                resxxerr=resdat[1][0][:]
                resyy=resdat[0][3][:]
                resyyerr=resdat[1][3][:]
                
                phasexy=phasedat[0][1][:]
                phasexyerr=phasedat[1][1][:]
                phaseyx=phasedat[0][2][:]
                phaseyxerr=phasedat[1][2][:]
                phasexx=phasedat[0][0][:]
                phasexxerr=phasedat[1][0][:]
                phaseyy=phasedat[0][3][:]
                phaseyyerr=phasedat[1][3][:] 
                
                res=[[resxx,resxxerr],[resxy,resxyerr],[resyx,resyxerr],
                     [resyy,resyyerr]]
                phase=[[phasexx,phasexxerr],[phasexy,phasexyerr],[phaseyx,phaseyxerr],
                       [phaseyy,phaseyyerr]]
                    
            #read in .edi file
            elif filename.find('edi',-4,len(filename))>=0:
                print 'Reading .edi file'
                impz=Z.Z(filename)
                station=impz.station
                rp=impz.getResPhase(thetar=rotz)
                period=impz.period
         
            periodmin.append(min(period))
            periodmax.append(max(period))
            #make figure for xy,yx components
            #plot resisitivity 
            ax=fig.add_subplot(2,1,1)

            erxy=ax.errorbar(period,rp.resxy,marker=xymarker,
                              ms=4,mfc='None',
                              mec=xycolor,mew=1,ls='None',
                              yerr=rp.resxyerr,
                              ecolor=xycolor)
            eryx=ax.errorbar(period,rp.resyx,marker=yxmarker,
                              ms=4,mfc='None',
                              mec=yxcolor,mew=1,ls='None',
                              yerr=rp.resyxerr,
                              ecolor=yxcolor)
            legendlst.append(erxy[0])
            legendlst.append(eryx[0])
            legendlabellst.append('{0} $E_x/B_y$'.format(ii))
            legendlabellst.append('{0} $E_y/B_x$'.format(ii))
            
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.grid(True)
            
            if plotnum==2:
                ax3=fig2.add_subplot(2,1,1)

                erxx=ax3.errorbar(period,rp.resxx,marker=xymarker,
                                  ms=4,mfc='None',
                                  mec=xycolor,mew=1,ls='None',
                                  yerr=rp.resxxerr,
                                  ecolor=xycolor)
                eryy=ax3.errorbar(period,rp.resyy,marker=yxmarker,
                                  ms=4,mfc='None',
                                  mec=yxcolor,mew=1,ls='None',
                                  yerr=rp.resyyerr,
                                  ecolor=yxcolor)
                legendlst2.append(erxx[0])
                legendlst2.append(eryy[0])
                legendlabellst2.append('$E_x/B_x$')
                legendlabellst2.append('$E_y/B_y$')
                
                ax3.set_yscale('log')
                ax3.set_xscale('log')
                ax3.grid(True)

            #plot phase
            ax2=fig.add_subplot(2,1,2,sharex=ax)
            ax2.errorbar(period,rp.phasexy,marker=xymarker,
                         ms=4,mfc='None',mec=xycolor,
                         mew=1,ls='None',yerr=rp.phasexyerr,
                         ecolor=xycolor)
            ax2.errorbar(period,np.array(rp.phaseyx)+180,
                         marker=yxmarker,
                         ms=4,mfc='None',mec=yxcolor,
                         mew=1,ls='None',yerr=rp.phaseyxerr,
                         ecolor=yxcolor)
            
            if plotnum==2:
                ax4=fig2.add_subplot(2,1,2,sharex=ax3)
                ax4.errorbar(period,rp.phasexx,marker=xymarker,
                         ms=4,mfc='None',mec=xycolor,
                         mew=1,ls='None',yerr=rp.phasexxerr,
                         ecolor=xycolor)
                ax4.errorbar(period,np.array(rp.phaseyy)+180,
                         marker=yxmarker,
                         ms=4,mfc='None',mec=yxcolor,
                         mew=1,ls='None',yerr=rp.phaseyyerr,
                         ecolor=yxcolor)
                         
        ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                         fontdict={'size':14,'weight':'bold'})
        ax2.set_ylabel('Phase(rad)', fontdict={'size':14,'weight':'bold'})
        ax2.set_xlabel('Period (s)',fontdict={'size':14,'weight':'bold'})
        ax2.set_xscale('log')
        ax.set_xlim(xmin=10**(np.floor(np.log10(min(periodmin)))),
                 xmax=10**(np.ceil(np.log10(max(periodmax)))))
        ax2.set_ylim(ymin=0,ymax=90)
        ax2.yaxis.set_major_locator(MultipleLocator(10))
        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.grid(True)
        
        ax2.legend(legendlst,legendlabellst,loc=2,
                   markerscale=1,borderaxespad=.05,labelspacing=.08,
                   handletextpad=.15,borderpad=.05,ncol=int(len(legendlst)/4.))
        plt.suptitle(station,fontsize=13,fontweight='bold')
        plt.show()
        if plotnum==2:
            ax3.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                         fontdict={'size':14,'weight':'bold'})
            ax4.set_ylabel('Phase(rad)', fontdict={'size':14,'weight':'bold'})
            ax4.set_xlabel('Period (s)',fontdict={'size':14,'weight':'bold'})
            ax4.set_xscale('log')
            ax3.set_xlim(xmin=10**(np.floor(np.log10(min(periodmin)))),
                     xmax=10**(np.ceil(np.log10(max(periodmax)))))
            ax4.set_ylim(ymin=-180,ymax=180)
            ax4.yaxis.set_major_locator(MultipleLocator(30))
            ax4.yaxis.set_minor_locator(MultipleLocator(5))
            ax4.grid(True)
            
            ax4.legend(legendlst,legendlabellst,loc=2,
                       markerscale=1,borderaxespad=.05,labelspacing=.08,
                       handletextpad=.15,borderpad=.05)
            plt.suptitle(station,fontsize=13,fontweight='bold')
            plt.show()
            
def plotTipper(tipper):
    """
    plotTipper will plot the real and imaginary induction arrow
    
    Input: 
        edifile = full path to edifile
        
    Output:
        None
    """
    
    impz=Z.Z(edifile)
    
    
    
    pass
    
def plotTS(combinefilelst,df=1.0,fignum=1,start=None,stop=None,
           savefigname=None,dpi=None,format=None,orientation=None):
    """
    plotTS will plot timeseries that have been combined using combineFewFiles.
    combinefilelst is the output list of filenames
    
    Inputs:
        combinefilelst = list of combined file paths
        fignum = figure number
        start = start time string in hhmmss
        stop = stop time string in hhmmss
        savefig = 'y' to save, 'n' to not save
        savefigname = path to saved figure name
        dpi = resolution of figure
        format = type of saved file, pdf,svg,eps...
        orientation = orientation of figure on A4 paper
    
    Outputs:
        none
        
    """
    

    
    if start==None:
        combstr=os.path.basename(combinefilelst[0])
        tspot=combstr.find('to')
        print tspot
        #if combined files are days
        try:
            startf=0.0
            majortick=6
            minortick=1
        #if combined files are for just one day
        except ValueError:
            start=combstr[tspot-2:tspot]+'0000'
            startf=float(start[0:2]+'.'+str(float(start[2:4])/60)[2])
            majortick=1
            minortick=.3
    else:
        startf=int(start[0:2])
        majortick=1*df
        minortick=.3*df
            
    complst=[]
    lc=len(combinefilelst)
    fig1=plt.figure(fignum,[10,7])
    
    plt.rcParams['font.size']=10
    plt.rcParams['figure.subplot.left']=.10
    plt.rcParams['figure.subplot.right']=.96
    plt.rcParams['figure.subplot.bottom']=.07
    plt.rcParams['figure.subplot.top']=.96
    plt.rcParams['figure.subplot.wspace']=.25
    plt.rcParams['figure.subplot.hspace']=.20
    #plt.rcParams['font.family']='helvetica'
    
    for ii in range(lc):
        filename=combinefilelst[ii]
        complst.append(filename[-2:])
        ts=np.loadtxt(filename)
        t=np.arange(len(ts))/(3600.*df)+startf
        t=t[0:len(ts)]
        if ii==0:
            ax=plt.subplot(lc,1,ii+1)
        else:
            ax=plt.subplot(lc,1,ii+1,sharex=ax)
        #plt.plot(t[0:len(ts)],mt.normalizeL2(mt.dctrend(ts)),lw=2)
        plt.plot(t,ts,lw=.8)
        plt.ylabel(complst[ii],fontsize=12,fontweight='bold')
        if ii==len(complst):
    		plt.xlabel('Time (hrs)',fontsize=12,fontweight='bold')      
        ax.xaxis.set_major_locator(MultipleLocator(majortick))
        ax.xaxis.set_minor_locator(MultipleLocator(minortick))        
        plt.axis('tight')
    plt.show()
    if savefigname!=None:
        if dpi==None:
            dpi=100
        if format==None:
            format='pdf'
        if orientation==None:
            orientation='landscape'
        plt.savefig(savefigname,dpi=dpi,format=format,
                    orientation=orientation)
        plt.clf()
        plt.close(fig1)
    
def plotPTpseudoSection(filenamelst,colorkey='phiminang',esize=2,
                        offsetscaling=.005,colorkeymm=[0,90],stationid=[0,4],
                        title=None,cbshrink=.8,linedir='ns',fignum=1,rotz=0,
                        yscale='period',pxy=[8,8],dpi=300,
                        indarrows='n',ascale=5,cmap='ptcmap',tscale='period'):
    
    """
    plotPTpseudoSection(filenamelst,colorkey='beta',esize=2,offsetscaling=
    .005) will plot a pseudo section of phase tensor ellipses given a list of 
    full path filenames. 
     
    colorkey =  the fill color of the ellipses and can be:
            'phimin' for minimum phase
            'beta'  for phase tensor skew angle
            'ellipticity' for phase tensor ellipticity
            'phidet' for the determinant of the phase tensor
    colorkeymm = [min,max] min and max of colorkey to which colorbar is
                    normalized to. 
    esize = size of ellipse, float 
    
    offsetscaling is a factor that scales the distance from one station to the 
    next to make the plot readable. 
    
    stationid is start and stop of station name [start,stop]
    
    title is figure title added to Phase Tensor Elements + title
    
    indarrow = 'yri' to plot induction both real and imaginary induction
                    arrows 
                    'yr' to plot just the real induction arrows
                    'yi' to plot the imaginary induction arrows
                    'n' to not plot them
                    **Note: convention is to point towards a conductor **
    ascale = scaling factor to make induction arrows bigger
    
    cmap = color map of ellipse facecolor.  So far the colormaps are:
            ptcmap = yellow (low phase) to red (high phase)
            ptcmap3 = white (low numbers) to blue (high numbers)
            skcmap = blue to yellow to red
            skcmap2 = blue to white to red
            rtcmap = blue to purple to red
    tscale = period or frequency for the title of the plot

    """
    
    fs=int(dpi/30.)
    plt.rcParams['font.size']=fs
    plt.rcParams['figure.subplot.left']=.08
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.06
    plt.rcParams['figure.subplot.top']=.96
    plt.rcParams['figure.subplot.wspace']=.55
    plt.rcParams['figure.subplot.hspace']=.70
    #plt.rcParams['font.family']='helvetica'
    
    ckmin=colorkeymm[0]
    ckmax=colorkeymm[1]
    
    #create a plot instance
    fig=plt.figure(fignum,pxy,dpi=300)
    ax=fig.add_subplot(1,1,1,aspect='equal')
    stationlst=[]
    offsetlst=[]
    minlst=[]
    maxlst=[]
    #plot phase tensor ellipses
    for ii,fn in enumerate(filenamelst):
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
        #get phase tensor elements and flip so the top is small periods/high 
        #frequencies
        pt=imp.getPhaseTensor(thetar=rotz)
        periodlst=imp.period[::-1]
        phimax=pt.phimax[::-1]
        phimin=pt.phimin[::-1]
        azimuth=pt.azimuth[::-1]
        if colorkey=='phiminang' or colorkey=='phimin':
            colorarray=pt.phiminang[::-1]
        if colorkey=='phidet':
            colorarray=np.sqrt(abs(pt.phidet[::-1]))*(180/np.pi)
        if colorkey=='beta':
            colorarray=pt.beta[::-1]
        if colorkey=='ellipticity':
            colorarray=pt.ellipticity[::-1]
            
        n=len(periodlst)
        minlst.append(min(colorarray))
        maxlst.append(max(colorarray))

        for jj in range(n):
            
            #make sure the ellipses will be visable
            eheight=phimin[jj]/phimax[jj]*esize
            ewidth=phimax[jj]/phimax[jj]*esize
        
                
            #create an ellipse scaled by phimin and phimax and oriented along
            #the azimuth    
            ellipd=Ellipse((offset*offsetscaling,3*jj),width=ewidth,
                          height=eheight,
                          angle=azimuth[jj])
            ax.add_artist(ellipd)
            
            #get face color info
            if colorkey=='phiminang' or  colorkey=='phimin':
                cvar=(pt.phiminang[jj]-ckmin)/(ckmax-ckmin)
            elif colorkey=='phidet':
                cvar=(pt.phidet[jj]-ckmin)/(ckmax-ckmin)
            elif colorkey=='beta':
                cvar=(pt.beta[jj]-abs(ckmin))/(ckmax-ckmin)
            elif colorkey=='ellipticity':
                cvar=(pt.ellipticity[jj]-ckmin)/(ckmax-ckmin)
            else:
                raise NameError('color key '+colorkey+' not supported')
                
            
            #set facecolor depending on the colormap
            #yellow to red
            if cmap=='ptcmap':
                if abs(cvar)>1:
                    ellipd.set_facecolor((1,0,0))
                elif cvar<0:
                    ellipd.set_facecolor((1-abs(cvar),1,abs(cvar)))
                else:
                    ellipd.set_facecolor((1,1-abs(cvar),.1))
            #white to blue
            elif cmap=='ptcmap3':
                if abs(cvar)>1:
                    ellipd.set_facecolor((0,0,0))
                else:
                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
            #blue to yellow to red
            elif cmap=='skcmap2':
                if cvar<0 and cvar>-1:
                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                elif cvar<-1:
                    ellipd.set_facecolor((0,0,1))
                elif cvar>0 and cvar<1:
                    ellipd.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                elif cvar>1:
                    ellipd.set_facecolor((1,0,0))
            #blue to white to red
            elif cmap=='skcmap':
                if cvar<0 and cvar>-1:
                    ellipd.set_facecolor((abs(cvar),abs(cvar),1-abs(cvar)))
                elif cvar<-1:
                    ellipd.set_facecolor((0,0,1))
                elif cvar>0 and cvar<1:
                    ellipd.set_facecolor((1,1-abs(cvar),.01))
                elif cvar>1:
                    ellipd.set_facecolor((1,0,0))
                    
            if indarrows.find('y')==0:
                tip=imp.getTipper(thetar=rotz)
                aheight=.5*esize
                awidth=.2*esize
                #plot real tipper
                if indarrows=='yri' or indarrows=='yr':
                    txr=tip.magreal[jj]*np.cos(tip.anglereal[jj]*np.pi/180)*\
                        esize*5
                    tyr=tip.magreal[jj]*np.sin(tip.anglereal[jj]*np.pi/180)*\
                        esize*5

                    ax.arrow(offset*offsetscaling,3*jj,txr,tyr,lw=.75*awidth,
                         facecolor='k',edgecolor='k',
                         length_includes_head=False,
                         head_width=awidth,head_length=aheight)
                #plot imaginary tipper
                if indarrows=='yri' or indarrows=='yi':
                    txi=tip.magimag[jj]*np.cos(tip.angleimag[jj]*np.pi/180)*\
                        esize*5
                    tyi=tip.magimag[jj]*np.sin(tip.angleimag[jj]*np.pi/180)*\
                        esize*5

                    ax.arrow(offset*offsetscaling,3*jj,txi,tyi,lw=.75*awidth,
                             facecolor='b',edgecolor='b',
                             length_includes_head=False,
                             head_width=awidth,head_length=aheight)
                    
    offsetlst=np.array(offsetlst)
    ax.set_xlim(min(offsetlst)*offsetscaling-4,max(offsetlst)*offsetscaling+4)
    ax.set_ylim(-5,n*3+5)
    if tscale=='period':
        yticklabels=['{0:.3g}'.format(periodlst[ii]) for ii in np.arange(start=0,stop=n,
                     step=2)]
        ax.set_ylabel('Period (s)',fontsize=fs+5,fontweight='bold')
    elif tscale=='frequency':
        yticklabels=['{0:.4g}'.format(1./periodlst[ii]) for ii in np.arange(start=0,stop=n,
                     step=2)]
        ax.set_ylabel('Frequency (Hz)',fontsize=fs+5,fontweight='bold')
    ax.set_xlabel('Station',fontsize=fs+5,fontweight='bold')
    
    if title==None:
        pass
    else:
        ax.set_title(title,fontsize=fs+4)
    plt.yticks(np.arange(start=0,stop=3*n,step=6),yticklabels)
    plt.xticks(np.array(offsetlst)*offsetscaling,stationlst)
    
    if indarrows.find('y')==0:
        if indarrows=='yri':
            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
            ax.legend([treal[0],timag[0]],['Tipper_real','Tipper_imag'],
                      loc='lower right',
                      prop={'size':10,'weight':'bold'},
                      ncol=2,markerscale=.5,borderaxespad=.005,
                      borderpad=.25)
        elif indarrows=='yr':
            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
            ax.legend([treal[0]],['Tipper_real'],
                      loc='lower right',
                      prop={'size':10,'weight':'bold'},
                      ncol=2,markerscale=.5,borderaxespad=.005,
                      borderpad=.25)
        elif indarrows=='yi':
            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
            ax.legend([timag[0]],['Tipper_imag'],
                      loc='lower right',
                      prop={'size':10,'weight':'bold'},
                      ncol=2,markerscale=.5,borderaxespad=.005,
                      borderpad=.25)
    
    ax.grid(alpha=.25,which='both')
    
    print 'Colorkey min = ',min(minlst)
    print 'Colorkey max = ',max(maxlst)
    
    #make a colorbar with appropriate colors             
    ax2=make_axes(ax,shrink=.8)
    if cmap=='ptcmap': 
        cb=ColorbarBase(ax2[0],cmap=ptcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='ptcmap3': 
        cb=ColorbarBase(ax2[0],cmap=ptcmap3,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='skcmap': 
        cb=ColorbarBase(ax2[0],cmap=skcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='skcmap2': 
        cb=ColorbarBase(ax2[0],cmap=skcmap2,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='rtcmap': 
        cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))

    #label the color bar accordingly
    if colorkey=='phimin' or colorkey=='phiminang':
        cb.set_label('$\Phi_{min}$ (deg)',fontdict={'size':10,'weight':'bold'})
    elif colorkey=='beta':
        cb.set_label('Skew (deg)',fontdict={'size':10,'weight':'bold'})
    elif colorkey=='phidet':
        cb.set_label('Det{$\Phi$} (deg)',fontdict={'size':10,'weight':'bold'})
    elif colorkey=='ellipticity':
        cb.set_label('Ellipticity',fontdict={'size':10,'weight':'bold'})
    
    plt.show()

def plotRTpseudoSection(filenamelst,colorkey='rhodet',esize=2,
                        offsetscaling=.005,colorkeymm=[0,90],stationid=[0,4],
                        title=None,cbshrink=.8,linedir='ns',fignum=1,rotz=0,
                        yscale='period',pxy=[8,8],dpi=300):
    
    """
    plotRTpseudoSection(filenamelst,colorkey='beta',esize=2,offsetscaling=
    .005) will plot a pseudo section of resistivity tensor ellipses given a list of 
    full path filenames. (Weckmann et al. 2002)
    
    colorkey is the fill color of the ellipses and can be any of the dictionary 
    keys returned by Z.getPhaseTensor(), note skew is beta:
        'phimin','phi', 'phiminvar', 'azimuthvar', 'azimuth', 'betavar', 
        'phivar', 'alphavar', 'beta', 'ellipticityvar', 'phiminangvar', 
        'ellipticity', 'phimaxangvar', 'alpha', 'phiminang', 'phimaxvar', 
        'phimaxang', 'phimax'
        . 
    esize is the normalized size of the ellipse, meaning the major axis will be
    this number. 
    
    offsetscaling is a factor that scales the distance from one station to the 
    next to make the plot readable. 
    
    colorkeymm is colorkey min and max input as [min,max]
    
    stationid is start and stop of station name [start,stop]
    
    title is figure title added to Phase Tensor Elements + title
    
    """
    
    fs=int(dpi/30.)
    plt.rcParams['font.size']=fs
    plt.rcParams['figure.subplot.left']=.08
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.06
    plt.rcParams['figure.subplot.top']=.96
    plt.rcParams['figure.subplot.wspace']=.55
    plt.rcParams['figure.subplot.hspace']=.70
    #plt.rcParams['font.family']='helvetica'
    
    
    #create a plot instance
    fig=plt.figure(fignum,pxy,dpi=dpi)
    ax=fig.add_subplot(1,1,1,aspect='equal')
    stationlst=[]
    offsetlst=[]
    minlst=[]
    maxlst=[]
    #plot phase tensor ellipses
    for ii,fn in enumerate(filenamelst):
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
        #get phase tensor elements and flip so the top is small periods/high 
        #frequencies
        rt=imp.getResTensor(thetar=rotz)
        periodlst=imp.period[::-1]
        phimax=rt.rhomax[::-1]
        phimin=rt.rhomin[::-1]
        azimuth=rt.rhoazimuth[::-1]
        if colorkey=='rhomin':
            colorarray=phimin
        elif colorkey=='rhomax':
            colorarray=phimax
        elif colorkey=='rhobeta':
            colorarray=rt.rhobeta[::-1]
        elif colorkey=='rhodet':
            colorarray=rt.rhodet[::-1]
        n=len(periodlst)
        minlst.append(min(colorarray))
        maxlst.append(max(colorarray))

        for jj in range(n):
            #make sure the ellipses will be visable
            eheight=phimax[jj]/phimax[jj]*esize
            ewidth=phimin[jj]/phimax[jj]*esize
        
            #create an ellipse scaled by phimin and phimax and oriented along
            #the azimuth 
            
            emax=10*esize
            if eheight>emax or ewidth>emax:
                pass
            else:
                ellip=Ellipse((offset*offsetscaling,3*jj),width=ewidth,
                              height=eheight,
                              angle=azimuth[jj])
                #color the ellipse as red high conductivity, blue low
                cvars=abs(np.log10(colorarray[jj]))/colorkeymm[1]
                if cvars>1:
                    ellip.set_facecolor((0,0,1))
                elif cvars<float(colorkeymm[0])/colorkeymm[1] or cvars<0:
                    ellip.set_facecolor((1,0,0))
                else:
                    ellip.set_facecolor((1-cvars,0,cvars))
                ax.add_artist(ellip)
            
    offsetlst=np.array(offsetlst)
    ax.set_xlim(min(offsetlst)*offsetscaling-4,max(offsetlst)*offsetscaling+4)
    ax.set_ylim(-5,n*3+5)
    if yscale=='period':
        yticklabels=['{0:.3g}'.format(periodlst[ii]) for ii in np.arange(start=0,stop=n,
                     step=2)]
        ax.set_ylabel('Period (s)',fontsize=fs+5,fontweight='bold')
    elif yscale=='frequency':
        yticklabels=['{0:.4g}'.format(1./periodlst[ii]) for ii in np.arange(start=0,stop=n,
                     step=2)]
        ax.set_ylabel('Frequency (Hz)',fontsize=fs+5,fontweight='bold')
    ax.set_xlabel('Station',fontsize=fs+5,fontweight='bold')
    
    if title==None:
        pass
#        plt.title('Phase Tensor Pseudo Section '+title,fontsize=16)
    else:
        ax.set_title(title,fontsize=fs+4)
    plt.yticks(np.arange(start=0,stop=3*n,step=6),yticklabels)
    plt.xticks(np.array(offsetlst)*offsetscaling,stationlst)
    
    ax.grid(alpha=.25,which='both')
    
    print 'Colorkey min = ',min(minlst)
    print 'Colorkey max = ',max(maxlst)
    
    #make colorbar
    ax2=make_axes(ax,shrink=cbshrink)
    cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=colorkeymm[0],
                        vmax=colorkeymm[1]))
    cb.set_label(colorkey,fontdict={'size':14,'weight':'bold'})
    
    plt.show()

def plotPTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='phimin',xpad=.2,
               ypad=.2,tickstrfmt='%2.2f',cborientation='vertical',
               colorkeymm=[0,90.],figsave='y',fmt=['png'],rotz=0,pxy=[10,12],
                galpha=.25,stationid=None,stationpad=.0005,
                sfdict={'size':12,'weight':'bold'},indarrows='n',
                cmap='ptcmap',tscale='period',mapscale='latlon',fignum=1,
                imagefile=None,image_extent=None):
    """ 
    plotPTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='phimin',xpad=.2,
               ypad=.2,tickstrfmt='%2.4f',cborientation='vertical',
               colorkeymax=90.,figsave='y',fmt=['png']) 
    plots phase tensor ellipses in map view from a list of edifiles with full 
    path.  Parameters are:
        
        freqspot = position in frequency list for plotting, an integer
        esize = size of ellipse, float
        colorkey =  the fill color of the ellipses and can be:
            'phimin' for minimum phase
            'beta'  for phase tensor skew angle
            'ellipticity' for phase tensor ellipticity
            'phidet' for the determinant of the phase tensor
        colorkeymm = [min,max] min and max of colorkey to which colorbar is
                    normalized to.
        xpad = pad from xmin and xmax on plot
        ypad = pad from ymin and ymax on plot
        tickstrfmt = format of tick strings needs to be a string format
        cborientation = colorbar orientation horizontal or vertical
        figsave = y or n, if yes figure will be saved to edifilelist path in 
                a folder called PTfigures.
        fmt = ['png'] is format of save figure can be pdf,eps or any other 
            formats supported by matplotlib. Can be a list of multiple formats.
            Note that pdf and eps do not properly yet.
        rotz = rotation angle clockwise from north
        pxy = dimensions of the figure in inches
        galpha = opacity of the grid
        stationid = first and last index of the station name, default is None
                    which is no station names.  If input use 
                    stationid=(0,4) for the 1st through 4th characters
        stationpad = padding for station name in the y direction
        sfdict = dictionary for station name where size is the font size
                and weight is the font weight
        indarrow = 'yri' to plot induction both real and imaginary induction
                    arrows 
                    'yr' to plot just the real induction arrows
                    'yi' to plot the imaginary induction arrows
                    'n' to not plot them
                    **Note: convention is to point towards a conductor **
        cmap = color map of ellipse facecolor.  So far the colormaps are:
            ptcmap = yellow (low phase) to red (high phase)
            ptcmap3 = white (low numbers) to blue (high numbers)
            skcmap = blue to yellow to red
            skcmap2 = blue to white to red
            rtcmap = blue to purple to red
        tscale = period or frequency for the title of the plot
        mapscale = latlon for lats and lons or
                   eastnorth for easting and northing, this is recomended if
                   you want to plot tipper data for small surveys.
        imagefile = path to an image file jpg or png or svg
        image_extent=(xmin,xmax,ymin,ymax) in coordinates accorting to mapscale
        
    """
    jj=freqspot
    fig=plt.figure(fignum,pxy,dpi=200)
    plt.clf()
    ax=fig.add_subplot(1,1,1,aspect='equal')
    
    if imagefile!=None:
        if image_extent==None:
            raise ValueError('Need to put in image extent')
        im=plt.imread(imagefile)
        ax.imshow(im,origin='lower',extent=image_extent,aspect='auto')
        
    elliplst=[]
    latlst=[]
    lonlst=[]
    if not esize is float:
        esize=float(esize)
    
    if mapscale=='latlon':
        tickstrfmt='%.3f'
    elif mapscale=='eastnorth':
        tickstrfmt='%.0f'
    
    ckmin=colorkeymm[0]
    ckmax=colorkeymm[1]
    
    plt.rcParams['font.size']=8
    plt.rcParams['figure.subplot.left']=.1
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.1
    plt.rcParams['figure.subplot.top']=.93
    plt.rcParams['figure.subplot.wspace']=.55
    plt.rcParams['figure.subplot.hspace']=.70
    #plt.rcParams['font.family']='helvetica'
    
    for ii,filename in enumerate(edifilelst):
        #get phase tensor info
        imp=Z.Z(filename)
        #check to see if the period is there
        try:
            freq=1./imp.period[jj]
            if mapscale=='latlon':
                latlst.append(imp.lat)
                lonlst.append(imp.lon)
                plotx=imp.lon
                ploty=imp.lat
            elif mapscale=='eastnorth':
                zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
                if ii==0:
                    zone1=zone
                    plotx=east
                    ploty=north
                else:
                    if zone1!=zone:
                        if zone1[0:2]==zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east=east+500000
                        else:
                            east=east-500000
                        latlst.append(north)
                        lonlst.append(east)
                        plotx=east
                        ploty=north
                    else:
                        latlst.append(north)
                        lonlst.append(east)
                        plotx=east
                        ploty=north
            pt=imp.getPhaseTensor(thetar=rotz)
            phimin=pt.phimin[jj]
            phimax=pt.phimax[jj]
            eangle=pt.azimuth[jj]
            #create an ellipse object
            scaling=esize/phimax
            eheight=phimin*scaling
            ewidth=phimax*scaling
            ellipd=Ellipse((plotx,ploty),width=ewidth,height=eheight,
                          angle=eangle)
            #print imp.lon,imp.lat,scaling,ewidth,eheight,phimin,phimax
            elliplst.append(ellipd)
            ax.add_artist(ellipd)
            
            #get face color info
            if colorkey=='phiminang' or  colorkey=='phimin':
                cvar=(pt.phiminang[jj]-ckmin)/(ckmax-ckmin)
            elif colorkey=='phidet':
                cvar=(pt.phidet[jj]-ckmin)/(ckmax-ckmin)
            elif colorkey=='beta':
                cvar=(pt.beta[jj]-abs(ckmin))/(ckmax-ckmin)
            elif colorkey=='ellipticity':
                cvar=(pt.ellipticity[jj]-ckmin)/(ckmax-ckmin)
            else:
                raise NameError('color key '+colorkey+' not supported')
            
            #set facecolor depending on the colormap
            #yellow to red
            if cmap=='ptcmap':
                if abs(cvar)>1:
                    ellipd.set_facecolor((1,0,0))
                elif cvar<0:
                    ellipd.set_facecolor((1-abs(cvar),1,abs(cvar)))
                else:
                    ellipd.set_facecolor((1,1-abs(cvar),.1))
            #white to blue
            elif cmap=='ptcmap3':
                if abs(cvar)>1:
                    ellipd.set_facecolor((0,0,0))
                else:
                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
            #blue to yellow to red
            elif cmap=='skcmap2':
                if cvar<0 and cvar>-1:
                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                elif cvar<-1:
                    ellipd.set_facecolor((0,0,1))
                elif cvar>0 and cvar<1:
                    ellipd.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                elif cvar>1:
                    ellipd.set_facecolor((1,0,0))
            #blue to white to red
            elif cmap=='skcmap':
                if cvar<0 and cvar>-1:
                    ellipd.set_facecolor((abs(cvar),abs(cvar),1-abs(cvar)))
                elif cvar<-1:
                    ellipd.set_facecolor((0,0,1))
                elif cvar>0 and cvar<1:
                    ellipd.set_facecolor((1,1-abs(cvar),.01))
                elif cvar>1:
                    ellipd.set_facecolor((1,0,0))
                    
            #-----------Plot Induction Arrows---------------------------
            if indarrows.find('y')==0:
                if mapscale=='latlon':
                    print 'Might try mapscale=latlon for better scale of arrows'
                    
                tip=imp.getTipper(thetar=rotz)
                aheight=.5*esize
                awidth=.25*esize
                #plot real tipper
                if indarrows=='yri' or indarrows=='yr':
                    txr=tip.magreal[jj]*np.cos(tip.anglereal[jj]*np.pi/180)*\
                        scaling
                    tyr=tip.magreal[jj]*np.sin(tip.anglereal[jj]*np.pi/180)*\
                        scaling

                    ax.arrow(plotx,ploty,txr,tyr,lw=.0075*awidth,
                         facecolor='k',edgecolor='k',
                         length_includes_head=False,
                         head_width=awidth,head_length=aheight)
                #plot imaginary tipper
                if indarrows=='yri' or indarrows=='yi':
                    txi=tip.magimag[jj]*np.cos(tip.angleimag[jj]*np.pi/180)*\
                        scaling
                    tyi=tip.magimag[jj]*np.sin(tip.angleimag[jj]*np.pi/180)*\
                        scaling

                    ax.arrow(plotx,ploty,txi,tyi,lw=.0075*awidth,
                             facecolor='b',edgecolor='b',
                             length_includes_head=False,
                             head_width=awidth,head_length=aheight)
                         
            
            #------------Plot station name------------------------------
            if stationid!=None:
                ax.text(plotx,ploty+stationpad,
                        imp.station[stationid[0]:stationid[1]],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=sfdict)
                
        #if the period is not there 
        except IndexError:
            print 'Did not find index for station'.format(jj)+imp.station
    
    if mapscale=='latlon':    
        ax.set_xlabel('longitude',fontsize=10,fontweight='bold')
        ax.set_ylabel('latitude',fontsize=10,fontweight='bold')
        ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
        ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
        ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
        ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
    elif mapscale=='eastnorth':
        ax.set_xlabel('Easting (m)',fontsize=10,fontweight='bold')
        ax.set_ylabel('Northing (m)',fontsize=10,fontweight='bold')
        ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
        ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
        ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
        ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
    if tscale=='period':
        titlefreq='{0:.5g} (s)'.format(1./freq)
    else:
        titlefreq='{0:.5g} (Hz)'.format(freq)
    ax.set_title('Phase Tensor Map for '+titlefreq,
                 fontsize=10,fontweight='bold')
    ax.grid(alpha=galpha)
    
    if indarrows.find('y')==0:
        if indarrows=='yri':
            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
            ax.legend([treal[0],timag[0]],['Tipper_real','Tipper_imag'],
                      loc='upper center',
                      prop={'size':10,'weight':'bold'},
                      ncol=2,markerscale=.5,borderaxespad=.005,
                      borderpad=.25)
        elif indarrows=='yr':
            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
            ax.legend([treal[0]],['Tipper_real'],
                      loc='upper center',
                      prop={'size':10,'weight':'bold'},
                      ncol=2,markerscale=.5,borderaxespad=.005,
                      borderpad=.25)
        elif indarrows=='yi':
            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
            ax.legend([timag[0]],['Tipper_imag'],
                      loc='upper center',
                      prop={'size':10,'weight':'bold'},
                      ncol=2,markerscale=.5,borderaxespad=.005,
                      borderpad=.25)
    
    
    
    #make a colorbar with appropriate colors             
    ax2=make_axes(ax,shrink=.8)
    if cmap=='ptcmap': 
        cb=ColorbarBase(ax2[0],cmap=ptcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='ptcmap3': 
        cb=ColorbarBase(ax2[0],cmap=ptcmap3,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='skcmap': 
        cb=ColorbarBase(ax2[0],cmap=skcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='skcmap2': 
        cb=ColorbarBase(ax2[0],cmap=skcmap2,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    elif cmap=='rtcmap': 
        cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))

    #label the color bar accordingly
    if colorkey=='phimin' or colorkey=='phiminang':
        cb.set_label('$\Phi_{min}$ (deg)',fontdict={'size':10,'weight':'bold'})
    elif colorkey=='beta':
        cb.set_label('Skew (deg)',fontdict={'size':10,'weight':'bold'})
    elif colorkey=='phidet':
        cb.set_label('Det{$\Phi$} (deg)',fontdict={'size':10,'weight':'bold'})
    elif colorkey=='ellipticity':
        cb.set_label('Ellipticity',fontdict={'size':10,'weight':'bold'})
        
    plt.show()
    
    #save the figure if desired
    if figsave=='y':
        sf='{0:.5g}'.format(int(np.round(freq)))
        savepath=os.path.join(os.path.dirname(edifilelst[0]),'PTfigures')
        if not os.path.exists(savepath):
            os.mkdir(savepath)
            print 'Made directory: '+savepath
        else:
            pass
        for f in fmt:
            fig.savefig(os.path.join(savepath,
                                 'PTmap'+sf+'Hz.'+f),
                                 format=f)
            print 'Saved file figures to: '+ os.path.join(savepath, \
                                     'PTmap'+sf+'Hz.'+f)
        plt.close()
        

def plotResPhasePseudoSection(edifilelst,stationid=[0,4],ffactor=1,
                              maxperiod=60,aspect=4,cmap='jet_r',xtickspace=2,
                              linedir='ns',rotz=0,dpi=300):
    """
    plotResPhasePseudoSection(edifilelst,stationid=4,ffactor=10E3,df=100.,
                              maxperiod=24,aspect=4,cmap='jet_r') plots a 
    pseudo section from a list of edifiles with full path with descending 
    frequency or ascending period on the y axis and relative position on the x. 
    keyword arguments can be:
        stationid = start and finish of station string for plotting
        ffactor = fudge factor to make sure resistivities plot correctly
        df = sampling frequency (Hz)
        maxperiod = maximum period to plot
        aspect = horizontal to vertical ratio of plot
        cmap = colormap of image
    """
    
    fs=int(dpi/40)
    plt.rcParams['font.size']=fs
    plt.rcParams['figure.subplot.left']=.07
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.06
    plt.rcParams['figure.subplot.top']=.94
    plt.rcParams['figure.subplot.wspace']=.01
    plt.rcParams['figure.subplot.hspace']=.20
    #plt.rcParams['font.family']='helvetica'
    
    
    #create a plot instance
#    ax1=fig.add_subplot(2,2,1)
#    ax2=fig.add_subplot(2,2,2,sharex=ax1)
#    ax3=fig.add_subplot(2,2,3,sharex=ax1)
#    ax4=fig.add_subplot(2,2,4,sharex=ax1)
    
    #create empty lists to put things into
    stationlst=[]
    offsetlst=[]
    minperiodlst=[]
    maxperiodlst=[]
    periodlst=[]
    nlst=[]
    
    #create empty arrays to put data into
    n=len(edifilelst)
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
    
def plotRTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='rhodet',xpad=.2,
               ypad=.2,tickstrfmt='%2.4f',cborientation='vertical',
               colorkeymm=[0,4],figsave='y',fmt=['png'],rotz=0,pxy=[10,12],
                galpha=.25):
    """ 
    plotPTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='phimin',xpad=.2,
               ypad=.2,tickstrfmt='%2.4f',cborientation='vertical',
               colorkeymax=90.,figsave='y',fmt=['png']) 
    plots phase tensor ellipses in map view from a list of edifiles with full 
    path.  Parameters are:
        
        freqspot = position in frequency list for plotting, an integer
        esize = size of ellipse, float
        colorkey =  the fill color of the ellipses and can be any of the 
        dictionary keys returned by Z.getPhaseTensor(), note skew is beta:
        'phimin','phi', 'phiminvar', 'azimuthvar', 'azimuth', 'betavar', 
        'phivar', 'alphavar', 'beta', 'ellipticityvar', 'phiminangvar', 
        'ellipticity', 'phimaxangvar', 'alpha', 'phiminang', 'phimaxvar', 
        'phimaxang', 'phimax'
        colorkeymm = [min,max] min and max of colorkey to which colorbar is
                    normalized to.
        xpad = pad from xmin and xmax on plot
        ypad = pad from ymin and ymax on plot
        tickstrfmt = format of tick strings needs to be a string format
        cborientation = colorbar orientation horizontal or vertical
        figsave = y or n, if yes figure will be saved to edifilelist path in 
                a folder called PTfigures.
        fmt = ['png'] is format of save figure can be pdf,eps or any other 
            formats supported by matplotlib. Can be a list of multiple formats.
            Note that pdf and eps do not properly yet.
        
    """
    jj=freqspot
    fig=plt.figure(1,pxy,dpi=150)
    ax=fig.add_subplot(1,1,1,aspect='equal')
    elliplst=[]
    latlst=[]
    lonlst=[]
    if not esize is float:
        esize=float(esize)
    
    ckmin=colorkeymm[0]
    ckmax=colorkeymm[1]
    
    plt.rcParams['font.size']=8
    plt.rcParams['figure.subplot.left']=.1
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.1
    plt.rcParams['figure.subplot.top']=.93
    plt.rcParams['figure.subplot.wspace']=.55
    plt.rcParams['figure.subplot.hspace']=.70
    #plt.rcParams['font.family']='helvetica'
    
    for ii,filename in enumerate(edifilelst):
        #get phase tensor info
        imp=Z.Z(filename)
        try:
            freq=1./imp.period[jj]
            latlst.append(imp.lat)
            lonlst.append(imp.lon)
            rt=imp.getResTensor(thetar=rotz)
            phimax=rt.rhomax[jj]
            phimin=rt.rhomin[jj]
            eangle=rt.rhoazimuth[jj]
            if colorkey=='rhomin':
                colorarray=phimin
            elif colorkey=='rhomax':
                colorarray=phimax
            elif colorkey=='rhobeta':
                colorarray=rt.rhobeta[jj]
            elif colorkey=='rhodet':
                colorarray=rt.rhodet[jj]
            eangle=rt.rhoazimuth[jj]
            #create an ellipse object
            scaling=esize/phimax
            eheight=phimin*scaling
            ewidth=phimax*scaling
            
#            ellipd=Ellipse((imp.lon,imp.lat),width=ewidth,height=eheight,
#                          angle=eangle)
#            #print imp.lon,imp.lat,scaling,ewidth,eheight,phimin,phimax
#            elliplst.append(ellipd)
#            ax.add_artist(ellipd)
#            #get face color info
        
            #create an ellipse scaled by phimin and phimax and oriented along
            #the azimuth 
            emax=10*esize
#            print eheight,ewidth,emax
            if eheight>emax or ewidth>emax:
                pass
            else:
                ellip=Ellipse((imp.lon,imp.lat),width=ewidth,
                              height=eheight,
                              angle=eangle)
                #color the ellipse as red high conductivity, blue low
                cvars=(np.log10(abs(colorarray))-ckmin)/(ckmax-ckmin)
                #print cvars
                if cvars>1:
                    ellip.set_facecolor((0,0,1))
                elif cvars<float(ckmin)/ckmax or cvars<0:
                    ellip.set_facecolor((1,0,0))
                else:
                    ellip.set_facecolor((1-cvars,0,cvars))
                ax.add_artist(ellip)            
            
        except IndexError:
            pass
        
    ax.set_xlabel('longitude',fontsize=10,fontweight='bold')
    ax.set_ylabel('latitude',fontsize=10,fontweight='bold')
    ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
    ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
    ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
    ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
    titlefreq='%2.3f' % (1/freq)
    ax.set_title('Resistivity Tensor for '+titlefreq+'(s)',
                 fontsize=10,fontweight='bold')
    ax.grid(alpha=galpha)
                 
    ax2=make_axes(ax,shrink=.8)
    if colorkey=='rhodet' or colorkey=='phidet':
        cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label('Log$_{10}$ '+colorkey+' ($\Omega \cdot$m)')
    elif colorkey=='beta' or colorkey=='ellipticity':
        cb=ColorbarBase(ax2[0],cmap=ptcmap3,norm=Normalize(vmin=ckmin,vmax=ckmax))
        cb.set_label(colorkey+' (deg)')
    plt.show()
    
    if figsave=='y':
        sf='%2.3g' % freq
        savepath=os.path.join(os.path.dirname(edifilelst[0]),'RTfigures')
        if not os.path.exists(savepath):
            os.mkdir(savepath)
            print 'Made directory: '+savepath
        else:
            pass
        for f in fmt:
            fig.savefig(os.path.join(savepath,
                                 'RTmap'+sf+'Hz.'+f),
                                 format=f)
            print 'Saved file figures to: '+ os.path.join(savepath, \
                                     'RTmap'+sf+'Hz.'+f)
        plt.close()
    
def comparePT(edilst,esize=5,xspacing=5,yspacing=3,savepath=None,show='y',
              title=None,stationfind=[0,4],xlabel=None,rotz=0):
    """
    comparePT will compare phase tensors given by edilst.  This formulation is
    based on Weise et al. [2008]. A perfect match is equivalent to a unit sphere
    with a fill color of yellow.
    
    Inputs:
        edilst = list of edi files to compare input as [[edi1,edi2],[1n,2n]]
        esize = max size of ellipse, anything larger than this will be set to 
                zero.
        xspacing = spacing along x-axis
        yspaceing = spacing along y-axis
        savepath = full path to save file to, if None file is not saved.
        show = 'y' or 'n', if yes the figure will be displayed and 'n' the 
                figure will not be displayed, useful if producing and saving
                multiple figures.
        title = title of figure, if None, no title is put on figure
        stationfind = index of station in edifile
    
    """
    
    plt.rcParams['font.size']=8
    plt.rcParams['figure.subplot.left']=.1
    plt.rcParams['figure.subplot.right']=.94
    plt.rcParams['figure.subplot.bottom']=.08
    plt.rcParams['figure.subplot.top']=.95
    plt.rcParams['figure.subplot.hspace']=.05
    #create a plot instance
    fig=plt.figure(1,[8,8],dpi=150)
    ax1=fig.add_subplot(1,1,1,aspect='equal')
    
    colorlst=[]
    maxlst=[]
    stationclst=[] 
    sf=stationfind
    
    azimutharr=np.zeros((29,len(edilst)))
    for ss,station in enumerate(edilst):
        #make a data type Z      
        imp1=Z.Z(station[0])
        imp2=Z.Z(station[1])
        stationclst.append(os.path.basename(station[0])[sf[0]:sf[1]])
        
        #get the phase tensor information
        pt1=imp1.getPhaseTensor(thetar=rotz)
        pt2=imp2.getPhaseTensor(thetar=rotz)
        
        #loop over period plotting the difference between phase tensors
        period=imp1.period
        n=len(period)
        for ii in range(n):
            #calculate the difference between the two phase tensor ellipses
            phi=np.eye(2)-(pt1.phi[ii]/pt2.phi[ii]+
                                pt2.phi[ii]/pt1.phi[ii])/2         
            
            #compute the trace        
            tr=phi[0,0]+phi[1,1]
            #Calculate skew of phi and the cooresponding error
            skew=phi[0,1]-phi[1,0]
            #calculate the determinate and determinate error of phi
            phidet=abs(np.linalg.det(phi))
            
            #calculate reverse trace and error
            revtr=phi[0,0]-phi[1,1]
            
            #calculate reverse skew and error
            revskew=phi[1,0]+phi[0,1]
            
            beta=.5*np.arctan2(skew,tr)*(180/np.pi)
            alpha=.5*np.arctan2(revskew,revtr)*(180/np.pi)
            azimuth=alpha-beta  
    #        azimuth=pt1.azimuth[ii]-pt2.azimuth[ii]                  
            azimutharr[ii,ss]=azimuth
            phimax=np.sqrt(abs((.5*tr)**2+(.5*skew)**2))+\
                    np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
                
            #calculate minimum value for phi
            if phidet>=0:
                phimin=np.sqrt(abs((.5*tr)**2+(.5*skew)**2))-\
                np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
            elif phidet<0:
                phimin=-1*np.sqrt(abs((.5*tr)**2+(.5*skew)**2))-np.sqrt(abs(
                            (.5*tr)**2+(.5*skew)**2-(np.sqrt(phidet))**2))            
            #make ellipse
            if abs(phimax)>esize:
                phimax=.1
            if abs(phimin)>esize:
                phimin=.1
            
            #compute width and height of ellipse
            eheightd=phimax
            ewidthd=phimin
            
            #add ellipse to the plot
            ellipd=Ellipse((xspacing*(ss)+1/2,yspacing*(n-ii)),
                           width=ewidthd,
                           height=eheightd,
                           angle=azimuth)
            ax1.add_artist(ellipd)
            
            ecolor=1/((abs(eheightd)+abs(ewidthd))/2)
            if ecolor>1:
                ecolor=.99
            else:
                ecolor=1-ecolor
            colorlst.append(ecolor)
            if ecolor>.98:
                ecolor=.98
            #print 'cvrs3: ',15*cvars3*(np.pi/2)*(180/np.pi) 
            ellipd.set_facecolor((1,1-ecolor,.1)) 
        maxlst.append(np.mean(colorlst))
    
    nellip=Ellipse((xspacing*(ss+1),yspacing*n/2),width=1,height=1,angle=0)
    
    ax1.add_artist(nellip)
    nellip.set_facecolor((1,.1,.1))
    textrect=Rectangle((xspacing*(ss+1)-xspacing*.4,yspacing*n/2-yspacing*.8),
                       5,
                       4,
                       fill=False,
                       clip_on=False)
    ax1.add_artist(textrect)
    ax1.text(xspacing*(ss+1),yspacing*n/2-yspacing/2,'$\Delta$=1',
             fontdict={'size':10},
             horizontalalignment='center',
             verticalalignment='center')
    
    stationclst=[stationclst[kk][-2:] for kk in range(len(stationclst))]                
    
    yticklabels=['%2.3g' % period[n-ii] for ii in np.arange(start=1,stop=n+1,
                 step=2)]
    plt.ylabel('Period (s)',fontsize=10,fontweight='bold')
    plt.yticks(np.arange(start=yspacing,stop=yspacing*n+1,step=2*yspacing),
               yticklabels)
    
    plt.xticks(np.arange(start=0,stop=xspacing*len(stationclst),step=xspacing),
               stationclst)
               
    ax1.set_xlim(-1*(xspacing/2),len(stationclst)*xspacing+5)
    ax1.set_ylim(0,yspacing*n+2)
    if xlabel==None:
        ax1.set_xlabel('Station',fontdict={'size':10,'weight':'bold'})
    else:
        ax1.set_xlabel(xlabel,fontdict={'size':10,'weight':'bold'})
    
    if title!=None:
        ax1.set_title(title,
                  fontdict={'size':12,'weight':'bold'})
    ax1.grid()
    
    
    ax4=make_axes(ax1,shrink=.5,fraction=.1,orientation='vertical',pad=.005)
    cb1=ColorbarBase(ax4[0],cmap=ptcmap,norm=Normalize(vmin=0,vmax=1),
                     orientation='vertical')
    cb1.set_label('(|$\Delta_{max}$|+|$\Delta_{min}$|)/2')
    
    if savepath!=None:
        plt.savefig(savepath,fmt='png')
    
    if show=='n':
        plt.close()
    else:
        plt.show()

def comparePT2(edilst,esize=5,xspacing=5,yspacing=3,savepath=None,show='y',
              title=None,stationfind=[0,4],xlabel=None,colorkeymm=[-5,5],
              fignum=1,rotz=0):
    """
    comparePT will compare phase tensors given by edilst.  This formulation is
    based on Weise et al. [2008]. Instead of plotting the ellipse of the 
    difference, this will plot the first phase tensor with perturbation due
    to the second phase tensor.
    
    Inputs:
        edilst = list of edi files to compare input as [[edi1,edi2],[1n,2n]]
        esize = max size of ellipse, anything larger than this will be set to 
                zero.
        xspacing = spacing along x-axis
        yspaceing = spacing along y-axis
        savepath = full path to save file to, if None file is not saved.
        show = 'y' or 'n', if yes the figure will be displayed and 'n' the 
                figure will not be displayed, useful if producing and saving
                multiple figures.
        title = title of figure, if None, no title is put on figure
        stationfind = index of station in edifile
    
    """
    
    plt.rcParams['font.size']=8
    plt.rcParams['figure.subplot.left']=.1
    plt.rcParams['figure.subplot.right']=.94
    plt.rcParams['figure.subplot.bottom']=.08
    plt.rcParams['figure.subplot.top']=.95
    plt.rcParams['figure.subplot.hspace']=.05
    #create a plot instance
    fig=plt.figure(fignum,[8,8],dpi=150)
    ax1=fig.add_subplot(1,1,1,aspect='equal')
    
    stationclst=[] 
    sf=stationfind
    
    for ss,station in enumerate(edilst):
        #make a data type Z      
        imp1=Z.Z(station[0])
        imp2=Z.Z(station[1])
        stationclst.append(os.path.basename(station[0])[sf[0]:sf[1]])
        
        #get the phase tensor information
        pt1=imp1.getPhaseTensor(thetar=rotz)
        pt2=imp2.getPhaseTensor(thetar=rotz)
        
        #set colorkey
#        colorarray=2*ptdict1[colorkey]-ptdict2[colorkey]
        
        #loop over period plotting the difference between phase tensors
        period=imp1.period
        n=len(period)
        for ii in range(n):
            #calculate the difference between the two phase tensor ellipses         
                              
            azimuth=2*pt1.azimuth[ii]-pt2.azimuth[ii]                                    
            phimax=2*pt1.phimax[ii]-pt2.phimax[ii]  
            phimin=2*pt1.phimin[ii]-pt2.phimin[ii]  
            
            #compute width and height of ellipse
            scaling1=esize/max([phimax,phimin])
            eheight1=pt1.phimin[ii]*scaling1
            ewidth1=pt1.phimax[ii]*scaling1
            
#            scaling2=esize/max([phimax,phimin])
            eheight2=pt2.phimin[ii]*scaling1
            ewidth2=pt2.phimax[ii]*scaling1
            
#            scaling=esize/max([phimax,phimin])
            eheightd=phimin*scaling1
            ewidthd=phimax*scaling1

            
            #add ellipse to the plot            
            ellip1=Ellipse((xspacing*(ss)+1/2,yspacing*(n-ii)),
                           width=ewidth1,
                           height=eheight1,
                           angle=pt1.azimuth[ii])
            ax1.add_artist(ellip1)
            ellip1.set_edgecolor('black')
            ellip1.set_facecolor('none')

            ellip2=Ellipse((xspacing*(ss)+1/2,yspacing*(n-ii)),
                           width=ewidth2,
                           height=eheight2,
                           angle=pt2.azimuth[ii])
            ax1.add_artist(ellip2)
            ellip2.set_edgecolor('blue')
            ellip2.set_facecolor('none')
            
            ellipd=Ellipse((xspacing*(ss)+1/2,yspacing*(n-ii)),
                           width=ewidthd,
                           height=eheightd,
                           angle=azimuth)
            ax1.add_artist(ellipd)
            ellipd.set_edgecolor('red')
            ellipd.set_facecolor('none')
            
#            if abs(azimuthd)<=colorkeymm[1]:
#                cvars=azimuthd/colorkeymm[1]
#                if cvars<0:
#                    ellipd.set_facecolor((.1,1+cvars,abs(cvars)))
#                else:
#                    ellipd.set_facecolor((1,1-cvars,.1))
#            elif abs(azimuthd)>colorkeymm[1]:
#                cvars=azimuthd/colorkeymm[1]
#                if cvars<0:
#                    cvars=-1
#                    ellipd.set_facecolor((.1,1+cvars,abs(cvars)))
#                else:
#                    cvars=1
#                    ellipd.set_facecolor((1,1-cvars,.1))
                    
            #plot arrows of length phimax and angle azimuth
            xarrow=xspacing*(ss)+1/2
            yarrow=yspacing*(n-ii)
            theta1=pt1.azimuth[ii]*(np.pi/180.)
            theta2=pt2.azimuth[ii]*(np.pi/180.)
            arrow1=Arrow(xarrow,yarrow,
                         pt1.phimax[ii]*np.cos(theta1),
                         pt1.phimax[ii]*np.sin(theta1),
                         ec='black',fc='black')
            ax1.add_artist(arrow1)
            
            arrow2=Arrow(xarrow,yarrow,
                         pt2.phimax[ii]*np.cos(theta2),
                         pt1.phimax[ii]*np.sin(theta2),
                         ec='blue',fc='blue')
            ax1.add_artist(arrow2)
    
    stationclst=[stationclst[kk][-2:] for kk in range(len(stationclst))]                
    
    yticklabels=['%2.3g' % period[n-ii] for ii in np.arange(start=1,stop=n+1,
                 step=2)]
    plt.ylabel('Period (s)',fontsize=10,fontweight='bold')
    plt.yticks(np.arange(start=yspacing,stop=yspacing*n+1,step=2*yspacing),
               yticklabels)
    
    plt.xticks(np.arange(start=0,stop=xspacing*len(stationclst),step=xspacing),
               stationclst)
               
    ax1.set_xlim(-1*(xspacing/2),len(stationclst)*xspacing+5)
    ax1.set_ylim(0,yspacing*n+2)
    if xlabel==None:
        ax1.set_xlabel('Station',fontdict={'size':10,'weight':'bold'})
    else:
        ax1.set_xlabel(xlabel,fontdict={'size':10,'weight':'bold'})
    
    if title!=None:
        ax1.set_title(title,
                  fontdict={'size':12,'weight':'bold'})
    ax1.grid()
    
    
    ax4=make_axes(ax1,shrink=.5,fraction=.1,orientation='vertical',pad=.005)
    cb1=ColorbarBase(ax4[0],cmap=ptcmap2,
                     norm=Normalize(vmin=colorkeymm[0],vmax=colorkeymm[1]),
                     orientation='vertical')
    cb1.set_label('$\gamma_1 - \gamma_2$')
    
    if savepath!=None:
        plt.savefig(savepath,fmt='png')
    
    if show=='n':
        plt.close()
    else:
        plt.show()