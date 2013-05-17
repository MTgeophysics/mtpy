'''
This module contains functions to plot different MT things.

J Peacock
'''

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import mtpy.legacy.birrptools as brp
import os
import mtpy.legacy.mttools as mt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
import mtpy.legacy.old_z as Z
reload(Z)
import mtpy.utils.latlongutmconversion as utm2ll
from matplotlib import colors
from matplotlib import patches
import matplotlib.colorbar as mcb
from matplotlib import gridspec

np.random.seed(101)
zvals = np.random.rand(100, 100) * 10

# make a color map of fixed colors
cmap = colors.ListedColormap(['white', 'red'])
bounds=[0,5,10]
norm = colors.BoundaryNorm(bounds, cmap.N)

#==============================================================================
# Make some color maps for plotting
#==============================================================================
#yellow to red
ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=colors.LinearSegmentedColormap('ptcmap',ptcmapdict,256)

#blue to yellow to red
skcmapdict={'red':((0.0,0.0,0.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,1.0,0.0),(.5,1.0,0.0),(.5,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,1.0),(.5,0.0,1.0),(0.5,0.1,0.1),(1.0,0.1,0.1))}
skcmap=colors.LinearSegmentedColormap('skcmap',skcmapdict,256)

#blue to white to red
skcmapdict2={'red':((0.0,0.0,0.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,1.0,0.0),(.5,1.0,0.0),(.5,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,1.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,0.0,0.0))}
skcmap2=colors.LinearSegmentedColormap('skcmap2',skcmapdict2,256)

#color segmented map from blue to red
skmapdict3={'red':((0.0,0.0,0.0),()),
            'green':((0.0,0.0,0.0)),
            'blue':((0.0,0.0,1.0))}


#white to blue
ptcmapdict3={'red':((0.0,1.0,1.0),(1.0,0.0,0.0)),
            'green':((0.0,1.0,1.0),(1.0,0.0,0.0)),
            'blue':((0.0,1.0,1.0),(1.0,1.0,1.0))}
ptcmap3=colors.LinearSegmentedColormap('ptcmap3',ptcmapdict3,256)

#red to blue
rtcmapdict={'red':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'green':((0.0,0.0,0.0),(1.0,0.0,0.0)),
            'blue':((0.0,1.0,0.0),(1.0,1.0,0.0))}
rtcmap=colors.LinearSegmentedColormap('rtcmap',rtcmapdict,256)

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
    
     Arguments:
    -----------
        **filename** : string
                       filename containing coherence (.coh)
                       
        **fignum** : int
                     figure number
                     
        **savefigfilename** : list of strings
                             supply filenames to save figures to if desired
        
        **dpi** : int
                  dots-per-inch of figure resolution
        
        **format** : [ pdf | eps | jpg | png | svg ]
                    file type of saved figure pdf,svg,eps...
            
        **orientation** : [ landscape | portrait ]
                         orientation of figure on A4 paper
    :Example: ::
        
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> cohfile = r"/home/MT01/MT01.coh"
        >>> mtplot.plotcoh(cohfile)
        
        """
    
  
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
        ax1.legend(('Normal','Zero'),'lower left',shadow=False,
                   borderaxespad=.50)
        ax1.set_title('$\mathbf{E_x/B_y}$',fontsize=12,fontweight='bold')
        ax1.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax1.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax1.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax1.set_ylim(ymin=0,ymax=1)
        
        ax2=plt.subplot(1,2,2)
        ax2.semilogx(period,coh2,'k',linewidth=2)
        ax2.semilogx(period,zcoh2,'k--', linewidth=2,dashes=[4,1,4,1])
        ax2.grid(True)
        ax2.legend(('Normal','Zero'),'lower left',shadow=False,
                   borderaxespad=.50)
        ax2.set_title('$\mathbf{E_y/B_x}$',fontsize=12,fontweight='bold')
        ax2.set_xlabel('Period (s)',fontsize=12,fontweight='bold')
        ax2.set_ylabel('Coherence',fontsize=12,fontweight='bold')
        ax2.set_xlim(xmax=period[0],xmin=period[nd-1])
        ax2.set_ylim(ymin=0,ymax=1)
        
        plt.suptitle('Coherence for Station: '+station,fontsize=12,
                     fontweight='bold')
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
        
        plt.suptitle('Coherence for Station: '+station,fontsize=12,
                     fontweight='bold')
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

    
def plotResPhase(filename,fignum=1,ffactor=1,plotnum=1,title=None,
                 savefigfilename=None,dpi=None,format=None,orientation=None,
                 rotz=0):
    """
    plotResPhase(filename,fignum) will plot the apparent resistivity and 
    phase for TE and TM modes 
    from a .dat file produced by writedat.  If you want to save the plot 
    use the save button on top left.
    
    Arguments:
    ----------
        **filename** : string
                       filename containing impedance (.edi) or resistivity and
                       phase information (.dat)
                      
        **fignum** : int
                     figure number
                     
        **ffactor** : float
                      scaling factor for computing resistivity from impedances
        
        **thetar** : float
                     rotation angle of impedance tensor (deg or radians)
        
        **plotnum** : [ 1 | 2 | 3 ]
                        * 1 for just Ex/By and Ey/Bx
                        * 2 for all 4 components
                        * 3 for off diagonal plus the determinant
                        
        **title** : string
                    title of plot
                    
        **savefigfilename** : string
                              supply filename to save figure to if desired
                              
        **dpi** : int
                  dots-per-inch of figure resolution
        
        **format** : [ pdf | eps | jpg | png | svg ]
                    file type of saved figure pdf,svg,eps...
            
        **orientation** : [ landscape | portrait ]
                         orientation of figure on A4 paper
    :Example: ::
        
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> edifile = r"/home/MT01/MT01.edi"
        >>> mtplot.plotResPhase(edifile)
        >>> # plot all 4 components
        >>> mtplot.plotResPhase(edifile,plotnum=2)
        
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
    ax2.errorbar(period,(np.array(rp.phaseyx))%360,marker='o',ms=4,
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
        pymax=max([max(rp.phasexy),max((rp.phaseyx)%360)])
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
    
    Arguments:
    ----------
        **filenamelst** : list of strings
                          full path to .edi files to plot
                          
        **plottype** : [ 1 | 2 | 3 ]
                        * 1 -> each filename has its own figure window 
                        * 2 -> plots in one figure with subplots for each 
                                file in filenamelst
                        * 3 -> plots onto one figure and one plot
                        * *Default* is 1
        
        **plotnum** : [ 1 | 2 | 3 ]
                        * 1 -> will plot just the off diagonal components
                        * 2 -> will plot all 4 components one subplot for each 
                             diagonal pair
                        * 3 -> will plot the determinant of Z only
                        * *Default* is 1
                        
        **ffactor** : float
                      scaling factor for apparent resistivity. *Default* is 1
            
        **ylim** : tuple (lo10(min),log10(max))
                   minimum and maximu of apparent resistivity that scales the
                   plot, on a log scale. *Default* is 0
        
        **rotz** : float
                   angle in degrees to rotate the impedance tensor clockwise
                   positive. *Default* is 0
                    
        **kwdict** : dictionary or list
                     keywords for plotResPhase, can input as a list of 
                     keywords that is the same length as filenamelst
            
    :Example: ::
        
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> import os
        >>> edipath = r"/home/EDIfiles"
        >>> edilst = [os.path.join(edipath,MT{0}.format(ii) 
        >>> ...       for ii in range(1,4))]
        >>> # plot all 3 edi files together for all 4 components
        >>> mtplot.resPhasePlots(edilst,plottype=3,plotnum=2)
        >>> #plot in individual plots with different parameters for each
        >>> klst=[{'thetar':10},{'thetar':45},{'thetar':75}]
        >>> mtplot.resPhasePlots(edilst,plottype=1,plotnum=2,kwdict=klst)
    """
    
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
    
def plotPTpseudoSection(filenamelst,colorkey='phimin',esize=2,
                        offsetscaling=.005,colorkeymm=(0,90),stationid=[0,4],
                        title=None,cbshrink=.8,linedir='ns',fignum=1,rotz=0,
                        pxy=[8,8],dpi=300,indarrows='n',ascale=5,
                        cmap='ptcmap',tscale='period'):
    
    """
    plotPTpseudoSection(filenamelst,colorkey='beta',esize=2,offsetscaling=
    .005) will plot a pseudo section of phase tensor ellipses given a list of 
    full path filenames. 
    
    Arguments:
    ----------
    
        **filenamelst** : list of strings
                          full paths to .edi files to plot
        
        **colorkey** : [ 'phimin' | 'beta' | 'ellipticity' | 'phidet' ]
                       fill color of the ellipses and can be:
                       * 'phimin'      -> minimum phase
                       * 'beta'        -> phase tensor skew angle
                       * 'ellipticity' -> phase tensor ellipticity
                       * 'phidet'      -> determinant of the phase tensor
                       * *Default* is 'phimin'
                       
        **colorkeymm** : tuple (min,max)
                        min and max of colorkey to which colorbar is normalized
                        to.  In degrees except for ellipticity which is already
                        normalized.
                        
        **esize** : float
                    size of ellipse in map units 
        
        **offsetscaling** : float
                            is a factor that scales the distance from one 
                            station to the next to make the plot readable.
                            *Default* is 0.005
        **linedir** : [ 'ns' | 'ew' ]
                      predominant direction of profile line
                      * 'ns' -> North-South Line
                      * 'ew' -> East-West line
                      * *Default* is 'ns'
        
        **stationid** : tuple or list 
                        start and stop of station name indicies.  
                        ex: for MT01dr stationid=(0,4) will be MT01
        
        **rotz** : float
                   angle in degrees to rotate the data clockwise positive.
                   *Default* is 0
        
        **title** : string
                    figure title added to Phase Tensor + title
                    
        **cbshrink** : float
                       percent to shrink the color bar to fit the image better
                       *Default* is 0.8 -> 80 percent
                       
        **fignum** : int
                     figure number.  *Default* is 1
        
        **indarrow** : [ 'yri' | 'yr' | 'yi' | 'n' ]
                        * 'yri' to plot induction both real and imaginary 
                           induction arrows 
                        * 'yr' to plot just the real induction arrows
                        * 'yi' to plot the imaginary induction arrows
                        * 'n' to not plot them
                        * *Default* is 'n'                        
                        **Note: convention is to point towards a conductor **
                         
        **ascale** : float
                     scaling factor to make induction arrows bigger
        
        **cmap** :[ 'ptcmap' | 'ptcmap3' | 'skcmap' | 'skcmap2' | 'rtcmap' ]
                  color map of ellipse facecolor.  So far the colormaps are:
                      * 'ptcmap'  -> yellow (low phase) to red (high phase)
                      * 'ptcmap3' -> white (low numbers) to blue (high numbers)
                      * 'skcmap'  -> blue to yellow to red
                      * 'skcmap2' -> blue to white to red
                      * 'rtcmap'  -> blue to purple to red
                      
        **tscale** : [ 'period' | 'frequency' ]
                     * 'period'    -> plot vertical scale in period
                     * 'frequency' -> plot vertical scale in frequency
                     
    :Example: ::
        
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> import os
        >>> edipath = r"/home/EDIfiles"
        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')]
        >>> # color by phimin with a range of 20-70 deg
        >>> mtplot.plotPTpseudoSection(edilst,colorkeymm=(20,70))
        >>> 
        >>> # add induction arrows to skew angle with range (-5,5)
        >>> mtplot.plotPTpseudoSection(edilst,colorkey='beta',
        >>> ...                        colorkeymm=(-5,5),indarrows='yri')

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
    
    **filenamelst** : list of strings
                          full paths to .edi files to plot
        
        **colorkey** : [ 'rhodet' ]
                       fill color of the ellipses and can be:
                       * 'rhodet'      -> minimum phase
                       * more to come
                       * *Default* is 'rhodet'
                       
        **colorkeymm** : tuple (min,max)
                        min and max of colorkey to which colorbar is normalized
                        to.  In log10 resistivity
                        
        **esize** : float
                    size of ellipse in map units 
        
        **offsetscaling** : float
                            is a factor that scales the distance from one 
                            station to the next to make the plot readable.
                            *Default* is 0.005
        **linedir** : [ 'ns' | 'ew' ]
                      predominant direction of profile line
                      * 'ns' -> North-South Line
                      * 'ew' -> East-West line
                      * *Default* is 'ns'
        
        **stationid** : tuple or list 
                        start and stop of station name indicies.  
                        ex: for MT01dr stationid=(0,4) will be MT01
        
        **rotz** : float
                   angle in degrees to rotate the data clockwise positive.
                   *Default* is 0
        
        **title** : string
                    figure title added to Phase Tensor + title
                    
        **cbshrink** : float
                       percent to shrink the color bar to fit the image better
                       *Default* is 0.8 -> 80 percent
                       
        **fignum** : int
                     figure number.  *Default* is 1
                      
        **yscale** : [ 'period' | 'frequency' ]
                     * 'period'    -> plot vertical scale in period
                     * 'frequency' -> plot vertical scale in frequency
                     
    :Example: ::
        
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> import os
        >>> edipath = r"/home/EDIfiles"
        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')]
        >>> # color by rhodet with a range of 0-4 log10 resistivity
        >>> mtplot.plotRTpseudoSection(edilst,colorkeymm=(0,4))
    
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
               imagefile=None,image_extent=None,refpoint=(0,0),cbshrink=.8,
               arrowprop={'headheight':0.25,'headwidth':0.25,'linewidth':0.5,
                          'arrowscale':1},
               arrowlegend={'placement':'lower right','xborderpad':.2,
                            'yborderpad':.2,'fontpad':.05,
                            'fontdict':{'size':10,'weight':'bold'}}):
    """  
    Plots phase tensor ellipses in map view from a list of edifiles with full 
    path.
    
    Arguments:
    ----------
        **edilst** : list of strings
                          full paths to .edi files to plot
                          
        **freqspot** : position in frequency list for plotting, an integer
        
        **esize** : float
                    size of ellipse in map units 
        
        **colorkey** : [ 'phimin' | 'beta' | 'ellipticity' | 'phidet' ]
                       fill color of the ellipses and can be:
                       * 'phimin'      -> minimum phase
                       * 'beta'        -> phase tensor skew angle
                       * 'ellipticity' -> phase tensor ellipticity
                       * 'phidet'      -> determinant of the phase tensor
                       * *Default* is 'phimin'
                       
        **colorkeymm** : tuple (min,max)
                        min and max of colorkey to which colorbar is normalized
                        to.  In degrees except for ellipticity which is already
                        normalized.
                    
        **xpad** : float
                   padding in x-direction in map units from furthest ellipses
        
        **ypad** : float
                   padding in y-direction in map units from furthest ellipses
        
        **tickstrfmt** : string
                        format of tick strings needs to be a string format.
                        ex: '%.2f' for 2 decimal place in floating number
        
        **cborientation** : [ 'horizontal' | 'vertical' ]
                            colorbar orientation horizontal or vertical
        
        **figsave** :  [ 'y' | 'n' ]
                        * 'y' will save figure to edifilelist path in a folder 
                           called PTfigures.  The figure will also close once
                           saved.
                        * 'n' will not save the figure
                
        **fmt** : list of formats 
                 can be pdf,eps or any other formats supported by matplotlib. 
                 Can be a list of multiple formats.
                 **Note: that pdf and eps sometimes do not work properly.
            
        **rotz** : float
                   angle in degrees to rotate the data clockwise positive.
                   *Default* is 0
        
        **pxy** : tuple (width,height)
                  dimensions of the figure in inches.  *Default* is (10,12)
        
        **galpha** : float [0:1]
                     opacity of the grid.  0 for transparent, 1 for opaque
                     *Default* is 0.25
        
        **stationid** : tuple or list 
                        start and stop of station name indicies.  
                        ex: for MT01dr stationid=(0,4) will be MT01.
                        *Default* is None
                    
        **stationpad** : float
                         padding for station name in the y-direction in map
                         units
        
        **sfdict** : dictionary
                     font dictionary for station name. Keys can be
                     matplotlib.text properties, common ones are:
                         * 'size'   -> for font size
                         * 'weight' -> for font weight
                         * 'color'  -> for color of font
                                     
        **indarrow** : [ 'yri' | 'yr' | 'yi' | 'n' ]
                        * 'yri' to plot induction both real and imaginary 
                           induction arrows 
                        * 'yr' to plot just the real induction arrows
                        * 'yi' to plot the imaginary induction arrows
                        * 'n' to not plot them
                        * *Default* is 'n'                        
                        **Note: convention is to point towards a conductor **
       
       **arrowprop** : dictionary of arrow properties with keys:
                        * 'linewidth'  -> width of the arrow line
                        * 'headheight' -> height of arrow head
                        * 'headwidth'  -> width of the arrow head
                        * 'arrowscale' -> scale size of the arrow
        
        
        **arrowlegend** : dictionary of properties for legend with keys:
                        * 'placement -> placement of arrow legend can be:
                            - 'upper right'
                            - 'lower right'
                            - 'upper left'
                            - 'lower left'
                        * 'xborderpad' -> padding from x axis
                        * yborderpad'  -> padding from y axis
                        * fontpad'     -> padding between arrow and legend text
                        * fontdict'    -> dictionary of font properties
        
        **cmap** :[ 'ptcmap' | 'ptcmap3' | 'skcmap' | 'skcmap2' | 'rtcmap' ]
                  color map of ellipse facecolor.  So far the colormaps are:
                      * 'ptcmap'  -> yellow (low phase) to red (high phase)
                      * 'ptcmap3' -> white (low numbers) to blue (high numbers)
                      * 'skcmap'  -> blue to yellow to red
                      * 'skcmap2' -> blue to white to red
                      * 'rtcmap'  -> blue to purple to red
            
        **tscale** : [ 'period' | 'frequency' ]
                     * 'period'    -> plot vertical scale in period
                     * 'frequency' -> plot vertical scale in frequency
        
        **mapscale** : [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
                      scale of map
                          * 'latlon'     -> lats and longs
                          * 'eastnorth'  -> for easting and northing, this is 
                             recomended if you want to plot tipper data for 
                             small surveys.
                          * 'eastnorthkm' -> for kilometer scaling
                          ** Note: if the ellipses are plotting in the wrong 
                          spot in east-north scale, check to see if your survey
                          crosses UTM grids.  Currently there is no way of 
                          dealing with this.
                   
        **imagefile** : string
                        path to an image file jpg or png or svg
        
        **image_extent** : tuple (xmin,xmax,ymin,ymax)
                           coordinates according to mapscale. Must be input if
                           image file is not None.
        
        **refpoint** : tuple (x,y)
                       reference point estimate relative distance to.  This 
                       point will be (0,0) on the map and everything else is 
                       referenced to this point
         
    :Example: ::
        
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> import os
        >>> edipath = r"/home/EDIfiles"
        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')]
        >>> # color by phimin with a range of 20-70 deg in meter scale
        >>> mtplot.plotPTMaps(edilst,colorkeymm=(20,70),mapscale='eastnorth')
        >>> 
        >>> # add induction arrows to skew angle with range (-5,5)
        >>> mtplot.plotPTMaps(edilst,colorkey='beta',
        >>> ...               colorkeymm=(-5,5),indarrows='yri')
        >>>
        >>> # put an underlying image as a basemap in km scale
        >>> mtplot.plotPTMaps(edilst,imagefile=r"/home/Maps/Basemap.jpg",
        >>> ...               image_extent=(0,20,0,20),mapscale='eastnorthkm')


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
    elif mapscale=='eastnorth' or mapscale=='eastnorthkm':
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
        #get phase tensor
        pt=imp.getPhaseTensor(thetar=rotz)
        pt.phimax=np.nan_to_num(pt.phimax)
        pt.phimin=np.nan_to_num(pt.phimin)
        #check to see if the period is there
        try:
            freq=1./imp.period[jj]
            if mapscale=='latlon':
                latlst.append(imp.lat)
                lonlst.append(imp.lon)
                plotx=imp.lon-refpoint[0]
                ploty=imp.lat-refpoint[1]
            elif mapscale=='eastnorth':
                zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
                if ii==0:
                    zone1=zone
                    plotx=east-refpoint[0]
                    ploty=north-refpoint[1]
                else:
                    if zone1!=zone:
                        if zone1[0:2]==zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east=east+500000
                        else:
                            east=east-500000
                        latlst.append(north-refpoint[1])
                        lonlst.append(east-refpoint[0])
                        plotx=east-refpoint[0]
                        ploty=north-refpoint[1]
                    else:
                        latlst.append(north-refpoint[1])
                        lonlst.append(east-refpoint[0])
                        plotx=east-refpoint[0]
                        ploty=north-refpoint[1]
            elif mapscale=='eastnorthkm':
                zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
                if ii==0:
                    zone1=zone
                    plotx=(east-refpoint[0])/1000.
                    ploty=(north-refpoint[1])/1000.
                else:
                    if zone1!=zone:
                        if zone1[0:2]==zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east=east+500000
                        else:
                            east=east-500000
                        latlst.append((north-refpoint[1])/1000.)
                        lonlst.append((east-refpoint[0])/1000.)
                        plotx=(east-refpoint[0])/1000.
                        ploty=(north-refpoint[1])/1000.
                    else:
                        latlst.append((north-refpoint[1])/1000.)
                        lonlst.append((east-refpoint[0])/1000.)
                        plotx=(east-refpoint[0])/1000.
                        ploty=(north-refpoint[1])/1000.
            else:
                raise NameError('mapscale not recognized')
                    
            
            phimin=pt.phimin[jj]
            phimax=pt.phimax[jj]
            eangle=pt.azimuth[jj]
            #create an ellipse object
            if phimax==0 or phimax>100 or phimin==0 or pt.phiminang[jj]<0 or \
                pt.phiminang[jj]>100:
                eheight=.0000001*esize
                ewidth=.0000001*esize
            else:
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
#                if mapscale=='latlon':
#                    print 'Might try mapscale=latlon for better scale of arrows'
                    
                tip=imp.getTipper(thetar=rotz)
                aheight=arrowprop['headheight']
                awidth=arrowprop['headwidth']
                ascale=arrowprop['arrowscale']
                #plot real tipper
                if indarrows=='yri' or indarrows=='yr':
                    if tip.magreal[jj]<=1.0:
                        txr=tip.magreal[jj]*ascale*\
                            np.sin((tip.anglereal[jj])*np.pi/180)
                        tyr=tip.magreal[jj]*ascale*\
                            np.cos((tip.anglereal[jj])*np.pi/180)
    
                        ax.arrow(plotx,ploty,txr,tyr,lw=arrowprop['linewidth'],
                             facecolor='k',edgecolor='k',
                             length_includes_head=False,
                             head_width=awidth,head_length=aheight)
                    else:
                        pass
                #plot imaginary tipper
                if indarrows=='yri' or indarrows=='yi':
                    if tip.magimag[jj]<=1.0:
                        txi=tip.magimag[jj]*\
                            np.sin((tip.angleimag[jj])*np.pi/180)*scaling
                        tyi=tip.magimag[jj]*\
                            np.cos((tip.angleimag[jj])*np.pi/180)*scaling
    
                        ax.arrow(plotx,ploty,txi,tyi,lw=arrowprop['linewidth'],
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
    elif mapscale=='eastnorthkm':
        ax.set_xlabel('Easting (km)',fontsize=10,fontweight='bold')
        ax.set_ylabel('Northing (km)',fontsize=10,fontweight='bold')
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
    #plot induction arrow scale bar
    if indarrows.find('y')==0:
        parrx=ax.get_xlim()
        parry=ax.get_ylim()
        try:
            axpad=arrowlegend['xborderpad']
        except KeyError:
            axpad=xpad+arrowprop['arrowscale']
        try:
            aypad=arrowlegend['yborderpad']
        except KeyError:
            aypad=ypad
        try:
            txtpad=arrowlegend['fontpad']
        except KeyError:
            txtpad=.25*esize
            
        
        if arrowlegend['placement']=='lower right':
            pax=parrx[1]-axpad
            pay=parry[0]+aypad
            ptx=arrowprop['arrowscale']
            pty=0
            txa=parrx[1]-axpad+arrowprop['arrowscale']/2.
            txy=pay+txtpad
        elif arrowlegend['placement']=='upper right':
            pax=parrx[1]-axpad
            pay=parry[1]-aypad
            ptx=arrowprop['arrowscale']
            pty=0
            txa=parrx[1]-axpad+arrowprop['arrowscale']/2.
            txy=pay+txtpad
        elif arrowlegend['placement']=='lower left':
            pax=parrx[0]+axpad
            pay=parry[0]+aypad
            ptx=arrowprop['arrowscale']
            pty=0
            txa=parrx[0]+axpad+arrowprop['arrowscale']/2.
            txy=pay+txtpad
        elif arrowlegend['placement']=='upper left':
            pax=parrx[0]+axpad
            pay=parry[1]-aypad
            ptx=arrowprop['arrowscale']
            pty=0
            txa=parrx[0]+axpad+arrowprop['arrowscale']/2.
            txy=pay+txtpad
        else:
            raise NameError('arrowlegend not supported.')
            
        ax.arrow(pax,pay,ptx,pty,lw=arrowprop['linewidth'],
                             facecolor='k',edgecolor='k',
                             length_includes_head=False,
                             head_width=arrowprop['headwidth'],
                             head_length=arrowprop['headheight'])
        
        ax.text(txa,txy,'|T|=1',
                horizontalalignment='center',
                verticalalignment='baseline',
                fontdict={'size':10,'weight':'bold'})
                
    ax.grid(alpha=galpha)
    
#    if indarrows.find('y')==0:
#        if indarrows=='yri':
#            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
#            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
#            ax.legend([treal[0],timag[0]],['Tipper_real','Tipper_imag'],
#                      loc='upper center',
#                      prop={'size':10,'weight':'bold'},
#                      ncol=2,markerscale=.5,borderaxespad=.005,
#                      borderpad=.25)
#        elif indarrows=='yr':
#            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
#            ax.legend([treal[0]],['Tipper_real'],
#                      loc='upper center',
#                      prop={'size':10,'weight':'bold'},
#                      ncol=2,markerscale=.5,borderaxespad=.005,
#                      borderpad=.25)
#        elif indarrows=='yi':
#            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
#            ax.legend([timag[0]],['Tipper_imag'],
#                      loc='upper center',
#                      prop={'size':10,'weight':'bold'},
#                      ncol=2,markerscale=.5,borderaxespad=.005,
#                      borderpad=.25)
    
    
    
    #make a colorbar with appropriate colors             
    ax2=make_axes(ax,shrink=cbshrink)
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
        sf='_{0:.5g}'.format(freq)
        savepath=os.path.join(os.path.dirname(edifilelst[0]),'PTfigures')
        if not os.path.exists(savepath):
            os.mkdir(savepath)
            print 'Made directory: '+savepath
        else:
            pass
        for f in fmt:
            fig.savefig(os.path.join(savepath,
                                 'PTmap_'+colorkey+sf+'Hz.'+f),
                                 format=f)
            print 'Saved file figures to: '+ os.path.join(savepath, \
                                     'PTmap_'+colorkey+sf+'Hz.'+f)
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
    
def plotRoseStrikeAngles(edilst,fignum=1,fs=10,dpi=300,thetar=0,ptol=.05,
                         tpad=.60,galpha=.25,prange='data',plottype=1,
                         tipper='n',pterr=None):
    """
    plots the strike angle as determined by phase tensor azimuth (Caldwell et 
    al. [2004]) and invariants of the impedance tensor (Weaver et al. [2003]).
    
    The data is split into decades where the histogram for each is plotted in 
    the form of a rose diagram with a range of 0 to 180 degrees.
    Where 0 is North and 90 is East.   The median angle of the period band is 
    set in polar diagram.  The top row is the strike estimated from
    the invariants of the impedance tensor.  The bottom row is the azimuth
    estimated from the phase tensor.  If tipper is 'y' then the 3rd row is the
    strike determined from the tipper, which is orthogonal to the induction
    arrow direction.  
    
    Arguments:
    ----------
        **edilst** : list of full paths to edifiles to be used
        
        **fignum** : int
                     figure number to be plotted. *Default* is 1
        
        **fs** : float
                 font size for labels of plotting. *Default* is 10
                 
        **dpi** : int
                  dots-per-inch resolution of figure, 300 is needed for 
                  publications. *Default* is 300
                  
        **thetar** : float
                     angle of rotation clockwise positive. *Default* is 0
                     
        **ptol** : float
                   Tolerance level to match periods from different edi files.
                   *Default* is 0.05
                   
        **tpad** : float
                   padding of the angle label at the bottom of each polar 
                   diagram.  *Default* is 1.65
                   
        **galpha** : float (0:1)
                     transparency of the font boxes and grid lines.  0 is fully
                     tranparent and 1 is opaque.  *Default* is 0.25
                     
        **prange** : [ 'data' | (period_min,period_max) ]
                    period range to estimate the strike angle. Options are:
                        * *'data'* for estimating the strike for all periods
                            in the data.
                        * (pmin,pmax) for period min and period max, input as
                          (log10(pmin),log10(pmax))
        
        **plottype** : [ 1 | 2 ]
                        * *1* to plot individual decades in one plot
                        * *2* to plot all period ranges into one polar diagram
                              for each strike angle estimation
        
        **tipper** : [ 'y' | 'n' ]
                      * *'y'* to plot the tipper strike
                      * *'n'* to not plot tipper strike
                      
        **pterr** : float   
                    Maximum error in degrees that is allowed to estimate strike.
                    *Default* is None allowing all estimates to be used.
                

    :Example: ::
        
        >>> import os
        >>> import mtpy.imaging.mtplottools as mtplot
        >>> edipath = r"/home/EDIFiles"
        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')>0]
        >>> # plot rose plots in decades with tipper and an error floor on pt
        >>> mtplot.plotRoseStrikeAngles(edilst,plottype=1,pterr=5)
        >>> # plot all decades into one rose plot for each estimation
        >>> mtplot.plotRoseStrikeAngles(edilst,plottype=2,pterr=5)
    """
    plt.rcParams['font.size']=fs-2
    plt.rcParams['figure.subplot.left']=.07
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.09
    plt.rcParams['figure.subplot.top']=.90
    plt.rcParams['figure.subplot.wspace']=.2
    plt.rcParams['figure.subplot.hspace']=.4   
    
    invlst=[]
    ptlst=[]
    
    if tipper=='y':
        tiprlst=[]
    
    nc=len(edilst)
    nt=0
    kk=0
    
    for dd,edi in enumerate(edilst):
        z1=Z.Z(edi)
        mm=np.remainder(dd,4)
        period=z1.period
    
        #get maximum length of periods and which index that is
        if len(period)>nt:
            nt=len(period)
            kk=dd
        #-----------get strike angle from invariants--------------------------
        zinv=z1.getInvariants(thetar=thetar)
        
        #add 90 degrees because invariants assume 0 is north, but plotting assumes
        #that 90 is north and measures clockwise, thus the negative
        zs=90-zinv.strike
        
        #for plotting put the NW angles into the SE quadrant 
        zs[np.where(zs>90)]=zs[np.where(zs>90)]-180
        zs[np.where(zs<-90)]=zs[np.where(zs<-90)]+180
        
        #make a dictionary of strikes with keys as period
        mdictinv=dict([(ff,jj) for ff,jj in zip(z1.period,zs)])
        invlst.append(mdictinv)
    
        #------------get strike from phase tensor strike angle---------------
        pt=z1.getPhaseTensor(thetar=thetar)
        az=pt.azimuth
        azerr=pt.azimuthvar
        #put an error max on the estimation of strike angle
        if pterr:
            az[np.where(azerr>pterr)]=0.0
        #don't need to add 90 because pt assumes 90 is north.
        az[np.where(az>90)]=az[np.where(az>90)]-180
        az[np.where(az<-90)]=az[np.where(az<-90)]+180
        #make a dictionary of strikes with keys as period
        mdictpt=dict([(ff,jj) for ff,jj in zip(z1.period,az)])
        ptlst.append(mdictpt)
        
        #-----------get tipper strike------------------------------------
        if tipper=='y':
            tip=z1.getTipper(thetar=thetar)
            tipr=-tip.anglereal
            
            tipr[np.where(tipr>90)]=tipr[np.where(tipr>90)]-180
            tipr[np.where(tipr<-90)]=tipr[np.where(tipr<-90)]+180
            
            tiprdict=dict([(ff,jj) for ff,jj in zip(z1.period,tipr)])
            tiprlst.append(tiprdict)
            
        
        
    
    #get min and max period
    maxper=np.max([np.max(mm.keys()) for mm in invlst])
    minper=np.min([np.min(mm.keys()) for mm in ptlst])
    
    #make empty arrays to put data into for easy manipulation
    medinv=np.zeros((nt,nc))
    medpt=np.zeros((nt,nc))
    if tipper=='y':
        medtipr=np.zeros((nt,nc))
    
    #make a list of periods from the longest period list
    plst=np.logspace(np.log10(minper),np.log10(maxper),num=nt,base=10)
    pdict=dict([(ii,jj) for jj,ii in enumerate(plst)])
    
    #put data into arrays
    for ii,mm in enumerate(invlst):
        mperiod=mm.keys()
        for jj,mp in enumerate(mperiod):
            for kk in pdict.keys():
                if mp>kk*(1-ptol) and mp<kk*(1+ptol):
                    ll=pdict[kk]
                    medinv[ll,ii]=invlst[ii][mp]
                    medpt[ll,ii]=ptlst[ii][mp]
                    if tipper=='y':
                        medtipr[ll,ii]=tiprlst[ii][mp]
                else:
                    pass
        
    #-----Plot Histograms of the strike angles-----------------------------
    if prange=='data':
        brange=np.arange(np.floor(np.log10(minper)),
                         np.ceil(np.log10(maxper)),1)
    else:
        brange=np.arange(np.floor(prange[0]),np.ceil(prange[1]),1)
    
    
    #font dictionary
    fd={'size':fs,'weight':'normal'}
    
    #------------------plot indivdual decades---------------------------------
    if plottype==1:
        #plot specs
        plt.rcParams['figure.subplot.hspace']=.3
        plt.rcParams['figure.subplot.wspace']=.3
        
        fig3=plt.figure(fignum,dpi=dpi)
        plt.clf()
        nb=len(brange)
        for jj,bb in enumerate(brange,1):
            #make subplots for invariants and phase tensor azimuths
            if tipper=='n':
                axhinv=fig3.add_subplot(2,nb,jj,polar=True)
                axhpt=fig3.add_subplot(2,nb,jj+nb,polar=True)
                axlst=[axhinv,axhpt]
            if tipper=='y':
                axhinv=fig3.add_subplot(3,nb,jj,polar=True)
                axhpt=fig3.add_subplot(3,nb,jj+nb,polar=True)
                axhtip=fig3.add_subplot(3,nb,jj+2*nb,polar=True)
                axlst=[axhinv,axhpt,axhtip]
            
            #make a list of indicies for each decades    
            binlst=[]
            for ii,ff in enumerate(plst):
                if ff>10**bb and ff<10**(bb+1):
                    binlst.append(ii)
            
            #extract just the subset for each decade
            hh=medinv[binlst,:]
            gg=medpt[binlst,:]
            if tipper=='y':
                tr=medtipr[binlst,:]
                trhist=np.histogram(tr[np.nonzero(tr)].flatten(),bins=72,
                                       range=(-180,180))
                bartr=axhtip.bar((trhist[1][:-1])*np.pi/180,trhist[0],
                                 width=5*np.pi/180)
                for cc,bar in enumerate(bartr):
                    fc=float(trhist[0][cc])/trhist[0].max()*.9
                    bar.set_facecolor((0,1-fc/2,fc))
                        
            
            #estimate the histogram for the decade for invariants and pt
            invhist=np.histogram(hh[np.nonzero(hh)].flatten(),bins=72,
                                    range=(-180,180))
            pthist=np.histogram(gg[np.nonzero(gg)].flatten(),bins=72,
                                   range=(-180,180))
            
            #plot the histograms    
            barinv=axhinv.bar((invhist[1][:-1])*np.pi/180,invhist[0],width=5*np.pi/180)
            barpt=axhpt.bar((pthist[1][:-1])*np.pi/180,pthist[0],width=5*np.pi/180)
            
            for cc,bar in enumerate(barinv):
                fc=float(invhist[0][cc])/invhist[0].max()*.8
                bar.set_facecolor((fc,0,1-fc))
            for cc,bar in enumerate(barpt):
                fc=float(pthist[0][cc])/pthist[0].max()*.8
                bar.set_facecolor((fc,1-fc,0))
                
            #make axis look correct with N to the top at 90.
            for aa,axh in enumerate(axlst):
                axh.xaxis.set_major_locator(MultipleLocator(30*np.pi/180))
                axh.xaxis.set_ticklabels(['E','','',
                                          'N','','',
                                          'W','','',
                                          'S','',''])
                axh.grid(alpha=galpha)
                if aa==0:
                    axh.set_xlim(-90*np.pi/180,270*np.pi/180)
                    axh.text(np.pi,axh.get_ylim()[1]*tpad,
                             '{0:.1f}$^o$'.format(90-np.median(hh[np.nonzero(hh)])),
                              horizontalalignment='center',
                              verticalalignment='baseline',
                              fontdict={'size':fs-nb},
                              bbox={'facecolor':(.9,0,.1),'alpha':.25})
                    print '-----Period Range {0:.3g} to {1:.3g} (s)-----'.format(10**bb,
                          10**(bb+1))
                         
                    print '   *Z-Invariants:  median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                            90-np.median(hh[np.nonzero(hh)]),
                            90-invhist[1][np.where(invhist[0]==invhist[0].max())[0][0]],
                            90-np.mean(hh[np.nonzero(hh)])) 
                    if bb==-5:
                        axh.set_title('10$^{-5}$-10$^{-4}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==-4:
                        axh.set_title('10$^{-4}$-10$^{-3}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==-3:
                        axh.set_title('10$^{-3}$-10$^{-2}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==-2:
                        axh.set_title('10$^{-2}$-10$^{-1}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==-1:
                        axh.set_title('10$^{-1}$-10$^{0}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==0:
                        axh.set_title('10$^{0}$-10$^{1}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==1:
                        axh.set_title('10$^{1}$-10$^{2}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==2:
                        axh.set_title('10$^{2}$-10$^{3}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==3:
                        axh.set_title('10$^{3}$-10$^{4}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==4:
                        axh.set_title('10$^{4}$-10$^{5}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    elif bb==5:
                        axh.set_title('10$^{5}$-10$^{6}$s',fontdict=fd,
                                      bbox={'facecolor':'white','alpha':galpha})
                    axh.titleOffsetTrans._t=(0,.1)
        
                elif aa==1:
                    axh.set_xlim(-180*np.pi/180,180*np.pi/180)
                    axh.text(np.pi,axh.get_ylim()[1]*tpad,
                             '{0:.1f}$^o$'.format(90-np.median(gg[np.nonzero(gg)])),
                              horizontalalignment='center',
                              verticalalignment='baseline',
                              fontdict={'size':fs-nb},
                              bbox={'facecolor':(.9,.9,0),'alpha':galpha})
                    print '   *PT Strike:     median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                            90-np.median(gg[np.nonzero(gg)]),
                            90-pthist[1][np.where(pthist[0]==pthist[0].max())[0][0]],
                            90-np.mean(gg[np.nonzero(gg)]))
                    if tipper!='y':
                        print '\n'
                    if nb>5: 
                        if bb==-5:
                            axh.set_title('10$^{-5}$-10$^{-4}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-4:
                            axh.set_title('10$^{-4}$-10$^{-3}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-3:
                            axh.set_title('10$^{-3}$-10$^{-2}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-2:
                            axh.set_title('10$^{-2}$-10$^{-1}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-1:
                            axh.set_title('10$^{-1}$-10$^{0}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==0:
                            axh.set_title('10$^{0}$-10$^{1}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==1:
                            axh.set_title('10$^{1}$-10$^{2}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==2:
                            axh.set_title('10$^{2}$-10$^{3}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==3:
                            axh.set_title('10$^{3}$-10$^{4}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==4:
                            axh.set_title('10$^{4}$-10$^{5}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==5:
                            axh.set_title('10$^{5}$-10$^{6}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                elif aa==2:
                    axh.set_xlim(-180*np.pi/180,180*np.pi/180)
                    axh.text(np.pi,axh.get_ylim()[1]*tpad,
                             '{0:.1f}$^o$'.format(90-np.median(tr[np.nonzero(tr)])),
                              horizontalalignment='center',
                              verticalalignment='baseline',
                              fontdict={'size':fs-nb},
                              bbox={'facecolor':(0,.1,.9),'alpha':galpha})
                    print '   *Tipper Strike: median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                            90-np.median(tr[np.nonzero(tr)]),
                            90-trhist[1][np.where(trhist[0]==trhist[0].max())[0][0]],
                            90-np.mean(tr[np.nonzero(tr)])) 
                    print '\n'
                    if nb>5: 
                        if bb==-5:
                            axh.set_title('10$^{-5}$-10$^{-4}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-4:
                            axh.set_title('10$^{-4}$-10$^{-3}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-3:
                            axh.set_title('10$^{-3}$-10$^{-2}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-2:
                            axh.set_title('10$^{-2}$-10$^{-1}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==-1:
                            axh.set_title('10$^{-1}$-10$^{0}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==0:
                            axh.set_title('10$^{0}$-10$^{1}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==1:
                            axh.set_title('10$^{1}$-10$^{2}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==2:
                            axh.set_title('10$^{2}$-10$^{3}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==3:
                            axh.set_title('10$^{3}$-10$^{4}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==4:
                            axh.set_title('10$^{4}$-10$^{5}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                        elif bb==5:
                            axh.set_title('10$^{5}$-10$^{6}$s',fontdict=fd,
                                          bbox={'facecolor':'white',
                                                'alpha':galpha})
                                  
                if jj==1:
                    if aa==0:
                        axh.set_ylabel('Strike (Z)',fontdict=fd,labelpad=5000./dpi,
                                       bbox={'facecolor':(.9,0,.1),'alpha':galpha})
                    elif aa==1:
                        axh.set_ylabel('PT Azimuth',fontdict=fd,labelpad=5000./dpi,
                                       bbox={'facecolor':(.9,.9,0),'alpha':galpha})
                    elif aa==2:
                        axh.set_ylabel('Tipper Strike',fd,labelpad=5000./dpi,
                                       bbox={'facecolor':(0,.1,.9),'alpha':galpha})
                
                plt.setp(axh.yaxis.get_ticklabels(),visible=False)
                
        print 'Note: North is assumed to be 0 and the strike angle is measured'+\
              'clockwise positive.'
        
        plt.show()
    
    #------------------Plot strike angles for all period ranges--------------------
    elif plottype==2:
        #plot specs
        plt.rcParams['figure.subplot.left']=.07
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.100
        plt.rcParams['figure.subplot.top']=.88
        plt.rcParams['figure.subplot.hspace']=.3
        plt.rcParams['figure.subplot.wspace']=.2
        
        fig3=plt.figure(fignum,dpi=dpi)
        plt.clf()
        #make subplots for invariants and phase tensor azimuths
        if tipper=='n':
            axhinv=fig3.add_subplot(1,2,1,polar=True)
            axhpt=fig3.add_subplot(1,2,2,polar=True)
            axlst=[axhinv,axhpt]
        else:
            axhinv=fig3.add_subplot(1,3,1,polar=True)
            axhpt=fig3.add_subplot(1,3,2,polar=True)
            axhtip=fig3.add_subplot(1,3,3,polar=True)
            axlst=[axhinv,axhpt,axhtip]
        
        #make a list of indicies for each decades    
        binlst=[pdict[ff] for ff in plst 
                if ff>10**brange.min() and ff<10**brange.max()]
        
        #extract just the subset for each decade
        hh=medinv[binlst,:]
        gg=medpt[binlst,:]
        
        #estimate the histogram for the decade for invariants and pt
        invhist=np.histogram(hh[np.nonzero(hh)].flatten(),bins=72,range=(-180,180))
        pthist=np.histogram(gg[np.nonzero(gg)].flatten(),bins=72,range=(-180,180))
        
        #plot the histograms    
        barinv=axhinv.bar((invhist[1][:-1])*np.pi/180,invhist[0],width=5*np.pi/180)
        barpt=axhpt.bar((pthist[1][:-1])*np.pi/180,pthist[0],width=5*np.pi/180)
        
        for cc,bar in enumerate(barinv):
            fc=float(invhist[0][cc])/invhist[0].max()*.8
            bar.set_facecolor((fc,0,1-fc))
        for cc,bar in enumerate(barpt):
            fc=float(pthist[0][cc])/pthist[0].max()*.8
            bar.set_facecolor((fc,1-fc,0))
            
        if tipper=='y':
                tr=medtipr[binlst,:]
                trhist=np.histogram(tr[np.nonzero(tr)].flatten(),bins=72,
                                       range=(-180,180))
                bartr=axhtip.bar((trhist[1][:-1])*np.pi/180,trhist[0],
                                 width=5*np.pi/180)
                for cc,bar in enumerate(bartr):
                    fc=float(trhist[0][cc])/trhist[0].max()*.9
                    bar.set_facecolor((0,1-fc/2,fc))
                    
        #make axis look correct with N to the top at 90.
        for aa,axh in enumerate(axlst):
            axh.xaxis.set_major_locator(MultipleLocator(30*np.pi/180))
            axh.grid(alpha=galpha)
            plt.setp(axh.yaxis.get_ticklabels(),visible=False)
            axh.xaxis.set_ticklabels(['E','','',
                                      'N','','',
                                      'W','','',
                                      'S','',''])
                                      
            #put median angle in a box at bottom of polar diagram
            if aa==0:
                axh.set_ylim(0,invhist[0].max())
                axh.text(170*np.pi/180,axh.get_ylim()[1]*.65,
                         '{0:.1f}$^o$'.format(90-np.median(hh[np.nonzero(hh)])),
                          horizontalalignment='center',
                          verticalalignment='baseline',
                          fontdict={'size':fs-2},
                          bbox={'facecolor':(.9,0,.1),'alpha':.25})
                print '-----Period Range {0:.3g} to {1:.3g} (s)-----'.format(min(plst),
                          max(plst))
                         
                print '   *Z-Invariants:  median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                        90-np.median(hh[np.nonzero(hh)]),
                        90-invhist[1][np.where(invhist[0]==invhist[0].max())[0][0]],
                        90-np.mean(hh[np.nonzero(hh)])) 
    
            elif aa==1:
                axh.set_ylim(0,pthist[0].max())
                axh.text(170*np.pi/180,axh.get_ylim()[1]*.65,
                         '{0:.1f}$^o$'.format(90-np.median(gg[np.nonzero(gg)])),
                          horizontalalignment='center',
                          verticalalignment='baseline',
                          fontdict={'size':fs-2},
                          bbox={'facecolor':(.9,.9,0),'alpha':galpha})
                print '   *PT Strike:     median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                        90-np.median(gg[np.nonzero(gg)]),
                        90-pthist[1][np.where(pthist[0]==pthist[0].max())[0][0]],
                        90-np.mean(gg[np.nonzero(gg)]))
                if tipper!='y':
                    print '\n'
            elif aa==2:
                axh.set_ylim(0,trhist[0].max())
                axh.text(170*np.pi/180,axh.get_ylim()[1]*.65,
                         '{0:.1f}$^o$'.format(90-np.median(tr[np.nonzero(tr)])),
                          horizontalalignment='center',
                          verticalalignment='baseline',
                          fontdict={'size':fs-2},
                          bbox={'facecolor':(0,.1,.9),'alpha':galpha})
                print '   *Tipper Stike:  median={0:.1f} mode={1:.1f} mean={2:.1f}\n'.format(
                        90-np.median(tr[np.nonzero(tr)]),
                        90-trhist[1][np.where(trhist[0]==trhist[0].max())[0][0]],
                        90-np.mean(tr[np.nonzero(tr)]))
            #set the title of the diagrams
            if aa==0:
                axh.set_title('Strike (Z)',fontdict=fd,
                               bbox={'facecolor':(.9,0,.1),'alpha':galpha})
            elif aa==1:
                axh.set_title('PT Azimuth',fontdict=fd,
                               bbox={'facecolor':(.9,.9,0),'alpha':galpha})
            elif aa==2:
                axh.set_title('Tipper Strike',fontdict=fd,
                               bbox={'facecolor':(0,.1,.9),'alpha':galpha})
            axh.titleOffsetTrans._t=(0,.15)
            
                
                
        print 'Note: North is assumed to be 0 and the strike angle is measured '+\
              'clockwise positive.'
        
        plt.show()
