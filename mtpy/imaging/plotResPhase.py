'''
This module contains functions to plot different MT things.

J Peacock
'''

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
import mtpy1.core.z as Z
import mtpy1.utils.latlongutmconversion as utm2ll
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec

#==============================================================================
# Make some color maps for plotting
#==============================================================================
#yellow to red
ptcmapdict = {'red':((0.0,1.0,1.0), (1.0,1.0,1.0)),
              'green':((0.0,0.0,1.0), (1.0,0.0,1.0)),
              'blue':((0.0,0.0,0.0), (1.0,0.0,0.0))}
mt_yl2rd=colors.LinearSegmentedColormap('mt_yl2rd',ptcmapdict,256)

#blue to yellow to red
skcmapdict = {'red':((0.0,0.0,0.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,1.0,1.0)),
              'green':((0.0,1.0,0.0),(.5,1.0,0.0),(.5,0.0,1.0),(1.0,0.0,1.0)),
              'blue':((0.0,0.0,1.0),(.5,0.0,1.0),(0.5,0.1,0.1),(1.0,0.1,0.1))}
mt_bl2yl2rd=colors.LinearSegmentedColormap('mt_bl2yl2rd',skcmapdict,256)

#blue to white to red
skcmapdict2 = {'red':((0.0,0.0,0.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,1.0,1.0)),
               'green':((0.0,1.0,0.0),(.5,1.0,0.0),(.5,0.0,1.0),(1.0,0.0,1.0)),
               'blue':((0.0,0.0,1.0),(.5,1.0,1.0),(0.5,0.0,1.0),(1.0,0.0,0.0))}
mt_bl2wh2rd=colors.LinearSegmentedColormap('mt_bl2wh2rd',skcmapdict2,256)

            
#blue to white to red in segmented colors
mt_seg_bl2wh2rd = colors.ListedColormap(((0,0,1),(.5,.5,1),(.75,.75,1),
                                         (.9,.9,1),(1,1,1),(1.0,.9,.9),
                                         (1,.75,.75),(1,.5,.5),
                                         (1,0,0)))

#white to blue
ptcmapdict3 = {'red':((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'green':((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'blue':((0.0,1.0,1.0),(1.0,1.0,1.0))}
mt_wh2bl = colors.LinearSegmentedColormap('mt_wh2bl',ptcmapdict3,256)

#red to blue
rtcmapdict = {'red':((0.0,0.0,1.0),(1.0,0.0,1.0)),
              'green':((0.0,0.0,0.0),(1.0,0.0,0.0)),
              'blue':((0.0,1.0,0.0),(1.0,1.0,0.0))}
mt_rd2bl = colors.LinearSegmentedColormap('mt_rd2bl',rtcmapdict,256)

#define text formating for plotting
ckdict = {'phiminang':'$\Phi_{min}$ (deg)','phimin':'$\Phi_{min}$ (deg)',
          'phimaxang':'$\Phi_{max}$ (deg)','phimax':'$\Phi_{max}$ (deg)',
          'phidet':'Det{$\Phi$} (deg)','beta':'Skew (deg)',
          'ellipticity':'Ellipticity','beta_seg':'Skew (deg)'}
          
cmapdict = {'mt_yl2rd':mt_yl2rd,
            'mt_bl2yl2rd':mt_bl2yl2rd,
            'mt_wh2bl':mt_wh2bl,
            'mt_rd2bl':mt_rd2bl,
            'mt_bl2wh2rd':mt_bl2wh2rd,
            'mt_seg_bl2wh2rd':mt_seg_bl2wh2rd}

#define a dictionary for different zones of utm to put stations in correct 
#location
zonedict = dict([(a,ii) for ii,a in enumerate(['a','b','c','d','e','f','g','h',
               'i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x',
               'y','z'])])

class PlotResPhase(object):
    """
    Plots Resistivity and phase for the different modes of the MT response.  At
    the moment is supports the input of an .edi file. Other formats that will
    be supported are the impedance tensor and errors with an array of periods
    and .j format.
    
    
    Arguments:
        ----------
            **filename** : string
                           filename containing impedance (.edi) is the only 
                           format supported at the moment
                          
            **fignum** : int
                         figure number
                         *default* is 1
                         
            **ffactor** : float
                          scaling factor for computing resistivity from 
                          impedances.
                          *Default* is 1
            
            **rotz** : float
                       rotation angle of impedance tensor (deg or radians), 
                       *Note* : rotaion is clockwise positive
                       *default* is 0
            
            **plotnum** : [ 1 | 2 | 3 ]
                            * 1 for just Ex/By and Ey/Bx *default*
                            * 2 for all 4 components
                            * 3 for off diagonal plus the determinant
        
            **title** : string
                        title of plot
                        *default* is station name
                        
            **marker_dict** : dictionary
                              Dictionary for the marker properties of the plot
                              Keys include:
                                  * 
                        
        :Example: ::
            
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> edifile = r"/home/MT01/MT01.edi"
            >>> p1 = mtplot.PlotResPhase(edifile, plotnum=2)
            >>> # plots all 4 components
    
    """   
    
    def __init__(self, filename, fignum=1, plotnum=1, title=None, dpi=300, 
                 rotz=0, plot_yn='y'):
        
        #set some of the properties as attributes much to Lars' discontent
        self.fn = filename
        self.fignum = fignum
        self.plotnum = plotnum
        self.title = title
        self.dpi = dpi
        self.rotz = rotz
        
        #-->line properties
        #line style between points
        self.xy_ls='None'        
        self.yx_ls='None'        
        self.det_ls='None'        
        
        #outline color
        self.xy_color = 'b'
        self.yx_color = 'r'
        self.det_color = 'g'
        
        #face color
        self.xy_mfc='None'
        self.yx_mfc='None'
        self.det_mfc='None'
        
        #maker
        self.xy_marker = 's'
        self.yx_marker = 'o'
        self.det_marker = 'd'
        
        #size
        self.marker_size = 2
        
        #set plot limits
        self.xlimits = None
        self.res_limits = None
        self.phase_limits = None
        
        #set font parameters
        self.font_size = 7
        
        #set scaling factor of the apparent resistivity
        self.ffactor = 1

        #plot on initializing
        if plot_yn=='y':
            self.plot()
        
        
    def plot(self):
        """
        plotResPhase(filename,fignum) will plot the apparent resistivity and 
        phase for TE and TM modes 
        from a .dat file produced by writedat.  If you want to save the plot 
        use the save button on top left.
        
        
            
        """
        
        #set figure size according to what the plot will be.
        if self.plotnum==1 or self.plotnum==3:
            self.figsize = [4,6]
        elif self.plotnum==2:
            self.figsize = [6,6]
            
            
        #check to see if the file given is an edi file
        if self.fn.find('edi',-4,len(self.fn))>=0:
            print 'Reading '+self.fn
            #create a Z object from the file            
            self.fn_z = Z.Z(self.fn)
            
            #get the reistivity and phase object
            self.rp = self.fn_z.getResPhase(thetar=self.rotz)
            
            #set x-axis limits from short period to long period
            if self.xlimits==None:
                self.xlimits = (10**(np.floor(np.log10(self.fn_z.period[0]))),
                                10**(np.ceil(np.log10((self.fn_z.period[-1])))))
            if self.phase_limits==None:
                self.phase_limits = (0,89.9)
                
            if self.res_limits == None:
                self.res_limits = (10**(np.floor(np.log10(min[self.rp.resxy,
                                                              self.rp.resyx]))),
                                  10**(np.ceil(np.log10(max[self.rp.resxy,
                                                            self.rp.resyx]))))
                     
        # ==> will add in other file types and the ability to put in just
        #     the impedance tensor and its error
            
            
        else:
            raise ValueError('Could not read file: '+self.fn)
        
        
        #set some parameters of the figure and subplot spacing
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .93
        
        #set the font properties for the axis labels
        fontdict={'size':self.font_size+2, 'weight':'bold'}
        
        #create a grid to place the figures into, set to have 2 rows and 2 
        #columns to put any of the 4 components.  Make the phase plot
        #slightly shorter than the apparent resistivity plot and have the two
        #close to eachother vertically.
        gs = gridspec.GridSpec(2, 2, height_ratios=[2,1.5], hspace=.01)
        
        #--> make figure for xy,yx components
        if self.plotnum==1 or self.plotnum==3:
            #make figure instance
            self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
            
            #space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.1)
            
            #--> create the axes instances
            #apparent resistivity axis
            self.axr = plt.subplot(gs[0,:])
            
            #phase axis which shares the x axis (period) with the resistivity
            self.axp = plt.subplot(gs[1,:], sharex=self.axr)
            
            #place the y-coordinate labels in the same location for each axis
            self.axr.yaxis.set_label_coords(-.0565, 0.5)
            self.axp.yaxis.set_label_coords(-.0565, 0.5)
            
        #--> make figure for all 4 components
        elif self.plotnum==2:
            #make figure instance
            self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
            
            #space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.07)
            
            #--> create axes instance
            #apparent resistivity
            self.axr = plt.subplot(gs[0,0])
            
            #phase axis which shares the x axis (period) with the resistivity
            self.axp = plt.subplot(gs[1,0],sharex=self.axr)
            
            #place the y-coordinate labels in the same location for each axis
            self.axr.yaxis.set_label_coords(-.075, 0.5)
            self.axp.yaxis.set_label_coords(-.075, 0.5)
        
        #---------plot the apparent resistivity--------------------------------
        #--> plot as error bars and just as points xy-blue, yx-red
        #res_xy
        self.ebxyr = self.axr.errorbar(self.fn_z.period, 
                                       self.rp.resxy, 
                                       marker=self.xy_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.xy_mfc, 
                                       mec=self.xy_color, 
                                       mew=100./self.dpi, 
                                       ls=self.xy_ls, 
                                       yerr=self.rp.resxyerr, 
                                       ecolor=self.xy_color,
                                       capsize=self.marker_size)
        
        #res_yx                              
        self.ebyxr = self.axr.errorbar(self.fn_z.period, 
                                       self.rp.resyx, 
                                       marker=self.yx_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.yx_mfc,
                                       mec=self.yx_color, 
                                       mew=100./self.dpi,
                                       ls=self.yx_ls, 
                                       yerr=self.rp.resyxerr, 
                                       ecolor=self.yx_color,
                                       capsize=self.marker_size)
                                      
        #--> set axes properties
        plt.setp(self.axr.get_xticklabels(), visible=False)
        self.axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                            fontdict=fontdict)
        self.axr.set_yscale('log')
        self.axr.set_xscale('log')
        self.axr.set_xlim(self.xlimits)
        self.axr.set_ylim(self.res_limits)
        self.axr.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                      lw=.25)
        self.axr.legend((self.ebxyr[0],self.ebyxr[0]), 
                        ('$E_x/B_y$','$E_y/B_x$'),
                        loc=3, 
                        markerscale=1, 
                        borderaxespad=.01,
                        labelspacing=.07, 
                        handletextpad=.2, 
                        borderpad=.02)
        
        #-----Plot the phase---------------------------------------------------
        #phase_xy
        self.ebxyp = self.axp.errorbar(self.fn_z.period, 
                                       self.rp.phasexy, 
                                       marker=self.xy_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.xy_mfc,
                                       mec=self.xy_color, 
                                       mew=100./self.dpi,
                                       ls=self.xy_ls,
                                       yerr=self.rp.phasexyerr, 
                                       ecolor=self.xy_color,
                                       capsize=self.marker_size)
                                      
        #phase_yx: Note add 180 to place it in same quadrant as phase_xy
        self.ebyxp = self.axp.errorbar(self.fn_z.period, 
                                       self.rp.phaseyx+180, 
                                       marker=self.yx_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.yx_mfc, 
                                       mec=self.yx_color, 
                                       mew=100./self.dpi,
                                       ls=self.yx_ls, 
                                       yerr=self.rp.phaseyxerr, 
                                       ecolor=self.yx_color,
                                       capsize=self.marker_size)

        #check the phase to see if any point are outside of [0:90]
        if self.phase_limits==None:
            if min(self.rp.phasexy)<0 or min(self.rp.phaseyx+180)<0:
                 pymin = min([min(self.rp.phasexy), min(self.rp.phaseyx+180)])
                 if pymin>0:
                    pymin = 0
            else:
                pymin = 0
            
            if max(self.rp.phasexy)>90 or max(self.rp.phaseyx+180)>90:
                pymax = min([max(self.rp.phasexy), max(self.rp.phaseyx+180)])
                if pymax<91:
                    pymax = 89.9
            else:
                pymax = 89.9
                
            self.phase_limits = (pymin, pymax)
        
        #--> set axes properties
        self.axp.set_xlabel('Period (s)', fontdict)
        self.axp.set_ylabel('Phase (deg)', fontdict)
        self.axp.set_xscale('log')
        self.axp.set_ylim(self.phase_limits)        
        self.axp.yaxis.set_major_locator(MultipleLocator(15))
        self.axp.yaxis.set_minor_locator(MultipleLocator(5))
        self.axp.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                      lw=.25)
        
        #===Plot the xx, yy components if desired==============================
        if self.plotnum==2:
            #---------plot the apparent resistivity----------------------------
            self.axr2 = plt.subplot(gs[0,1])
            self.axr2.yaxis.set_label_coords(-.1, 0.5)
            
            #res_xx
            self.ebxxr = self.axr2.errorbar(self.fn_z.period, 
                                            self.rp.resxx, 
                                            marker=self.xy_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.xy_mfc,
                                            mec=self.xy_color, 
                                            mew=100./self.dpi,
                                            ls=self.xy_ls, 
                                            yerr=self.rp.resxxerr, 
                                            ecolor=self.xy_color,
                                            capsize=self.marker_size)
            
            #res_yy                              
            self.ebyyr = self.axr2.errorbar(self.fn_z.period, 
                                            self.rp.resyy, 
                                            marker=self.yx_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.yx_mfc, 
                                            mec=self.yx_color, 
                                            mew=100./self.dpi, 
                                            ls=self.yx_ls, 
                                            yerr=self.rp.resyyerr, 
                                            ecolor=self.yx_color,
                                            capsize=self.marker_size)

            #--> set axes properties
            plt.setp(self.axr2.get_xticklabels(), visible=False)
            self.axr2.set_yscale('log')
            self.axr2.set_xscale('log')
            self.axr2.set_xlim(self.xlimits)
            self.axr2.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                           lw=.25)
            self.axr2.legend((self.ebxxr[0], self.ebyyr[0]), 
                            ('$E_x/B_x$','$E_y/B_y$'),
                            loc=3, markerscale=1, 
                            borderaxespad=.01,
                            labelspacing=.07, 
                            handletextpad=.2, 
                            borderpad=.02)
            
            #-----Plot the phase-----------------------------------------------
            self.axp2=plt.subplot(gs[1,1],sharex=self.axr2)
            
            self.axp2.yaxis.set_label_coords(-.1, 0.5)
            
            #phase_xx
            self.ebxxp = self.axp2.errorbar(self.fn_z.period, 
                                            self.rp.phasexx, 
                                            marker=self.xy_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.xy_mfc, 
                                            mec=self.xy_color, 
                                            mew=100./self.dpi, 
                                            ls=self.xy_ls, 
                                            yerr=self.rp.phasexxerr,
                                            ecolor=self.xy_color,
                                            capsize=self.marker_size)
                                            
            #phase_yy
            self.ebyyp = self.axp2.errorbar(self.fn_z.period,
                                            self.rp.phaseyy, 
                                            marker=self.yx_marker,
                                            ms=self.marker_size, 
                                            mfc=self.yx_mfc,
                                            mec=self.yx_color, 
                                            mew=100./self.dpi, 
                                            ls=self.yx_ls, 
                                            yerr=self.rp.phaseyyerr, 
                                            ecolor=self.yx_color,
                                            capsize=self.marker_size)
            
            #--> set axes properties
            self.axp2.set_xlabel('Period (s)', fontdict)
            self.axp2.set_xscale('log')
            self.axp2.set_ylim(ymin=-179.9, ymax=179.9)        
            self.axp2.yaxis.set_major_locator(MultipleLocator(30))
            self.axp2.yaxis.set_minor_locator(MultipleLocator(5))
            self.axp2.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                           lw=.25)       
        
        #===Plot the Determinant if desired====================================                          
        if self.plotnum==3:
                
            #res_det
            self.ebdetr = self.axr.errorbar(self.fn_z.period, 
                                            self.rp.resdet, 
                                            marker=self.det_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.det_mfc,
                                            mec=self.det_color, 
                                            mew=100./self.dpi, 
                                            ls=self.det_ls, 
                                            yerr=self.rp.resdeterr, 
                                            ecolor=self.det_color,
                                            capsize=self.marker_size)
        
            #phase_det
            self.ebdetp = self.axp.errorbar(self.fn_z.period, 
                                            self.rp.phasedet, 
                                            marker=self.det_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.det_mfc, 
                                            mec=self.det_color, 
                                            mew=100./self.dpi, 
                                            ls=self.det_ls, 
                                            yerr=self.rp.phasedeterr, 
                                            ecolor=self.det_color,
                                            capsize=self.marker_size)
                                            
            self.axr.legend((self.ebxyr[0],self.ebyxr[0],self.ebdetr[0]),
                            ('$E_x/B_y$','$E_y/B_x$','$\det(\mathbf{\hat{Z}})$'),
                            loc=3,
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.02)
        
        
        #make title and show
        if self.title!=None:
            self.fig.suptitle(self.title, fontsize=self.font_size+2, 
                              fontweight='bold')
        else:
            self.fig.suptitle(self.fn_z.station, fontsize=self.font_size+2, 
                         fontweight='bold')
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

        if fig_dpi==None:
            fig_dpi = self.dpi
            
        if os.path.isdir(save_fn)==False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            plt.clf()
            plt.close(self.fig)
            
        else:
            save_fn = os.path.join(save_fn,self.fn_z.station+'_ResPhase.'+
                                    file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation)
        
        if close_plot=='y':
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
        
        return "Plots Resistivity and phase for the different modes of the MT \n" +\
              "response.  At the moment it supports the input of an .edi \n"+\
              "file. Other formats that will be supported are the impedance\n"+\
              "tensor and errors with an array of periods and .j format.\n"
              
class PlotMultipleFiles(object):
    """
    PlotMultipleFiles has methods to plot different representations of your 
    data includeing:
        
        **resistivityPhaseModes** : plots the modes in one plot or multiple 
                                    plots to compare the responses
                                    
        **resistivityPhasePseudoSections** : plots the different componets in
                                             a pseudo section format
                                             
        **resistivityPhaseMaps** : plots the different components in map view
                                   for different periods
        
        **phaseTensorPseudoSections** : plots the phase tensors in a pseudo
                                        section format, with the option of
                                        coloring and tipper arrows
                                        
        **phaseTensorMaps** : plots the phase tensors in map view for different
                              periods with the options for coloring and tipper
                              arrows
                              
        **plotStrikeRoseDiagrams** : plots the estimated strike direction in as
                                     Rose diagrams for the strike estimated 
                                     from the invariants of Weaver et al., 2003
                                     and the phase tensor strike of Caldwell et
                                     al., 2004.  In the works is the strike
                                     angle estimated from the Groom Bailey 
                                     decomposition of McNiece & Jones, 2001.
    """

    def resistivityPhaseModes(self):
        pass
    
    def resistivityPhasePseudoSections(self):
        pass
    
    def resistivityPhaseMaps(self):
        pass
    
    def phaseTensorPseudoSections(self):
        pass
    
    def phaseTensorMaps(self):
        pass
    
    def plotStrikeRoseDiagrams(self):
        pass

class MultipleResistivityPhasePlots(object):
    """
    plots multiple MT responses simultaneously either in single plots or in 
    one plot of subfigures or in a single plot with subfigures for each 
    component.
    
    
    """

    def __init__(self, edilst, plotnum=1, plotstyle='1', fignum=1, dpi=300, 
                 rotz=0, plot_yn='y'):
        """
        Initialize parameters
        """
        
        #set some of the properties as attributes much to Lars' discontent
        self.edilst = edilst
        self.plotstyle = plotstyle
        self.fignum = fignum
        self.plotnum = plotnum
        self.dpi = dpi
        self.rotz = rotz
        
        #-->line properties
        #line style between points
        self.xy_ls='None'        
        self.yx_ls='None'        
        self.det_ls='None'        
        
        #outline color
        self.xy_color = 'b'
        self.yx_color = 'r'
        self.det_color = 'g'
        
        #face color
        self.xy_mfc='None'
        self.yx_mfc='None'
        self.det_mfc='None'
        
        #maker
        self.xy_marker = 's'
        self.yx_marker = 'o'
        self.det_marker = 'd'
        
        #size
        self.marker_size = 2
        
        #set plot limits
        self.xlimits = None
        self.res_limits = (10.**0,10.**3)
        self.phase_limits = None
        
        #set font parameters
        self.font_size = 5
        
        #set scaling factor of the apparent resistivity
        self.ffactor = 1
        
        self.label_dict = {-6:'$10^{-6}$',
                           -5:'$10^{-5}$',
                           -4:'$10^{-4}$',
                           -3:'$10^{-3}$',
                           -2:'$10^{-2}$', 
                           -1:'$10^{-1}$', 
                            0:'$10^{0}$',
                            1:'$10^{1}$',
                            2:'$10^{2}$',
                            3:'$10^{3}$',
                            4:'$10^{4}$',
                            5:'$10^{5}$',
                            6:'$10^{6}$',
                            7:'$10^{7}$',
                            8:'$10^{8}$'}

        #plot on initializing
        if plot_yn=='y':
            self.plot()
        
        
    def plot(self):
        """
        plot the apparent resistivity and phase
        """
        
        if self.plotstyle=='1':
            self.plotlst=[]
            for ii,edi in enumerate(self.edilst,1):
                p1 = PlotResPhase(edi, fignum=ii, plotnum=self.plotnum, 
                                  dpi=self.dpi, rotz=self.rotz, plot_yn='n')
                
                #make sure all the properties are set to match the users
                #line style between points
                p1.xy_ls = self.xy_ls        
                p1.yx_ls = self.yx_ls        
                p1.det_ls = self.det_ls        
                
                #outline color
                p1.xy_color = self.xy_color 
                p1.yx_color = self.yx_color 
                p1.det_color = self.det_color
                
                #face color
                p1.xy_mfc = self.xy_mfc
                p1.yx_mfc = self.yx_mfc
                p1.det_mfc = self.det_mfc
                
                #maker
                p1.xy_marker = self.xy_marker
                p1.yx_marker = self.yx_marker
                p1.det_marker = self.det_marker 
                
                #size
                p1.marker_size = 2
                
                #set plot limits
                p1.xlimits = self.xlimits
                p1.res_limits = self.res_limits
                p1.phase_limits = self.phase_limits

                #set font parameters
                p1.font_size = self.font_size
                
                #set scaling factor of the apparent resistivity
                p1.ffactor = self.ffactor
                
                #--> plot the apparent resistivity and phase
                self.plotlst.append(p1)
                
                p1.plot()
        
        #-----Plot All in one figure with each plot as a subfigure------------        
        if self.plotstyle=='all':

            ns = len(self.edilst)

            #set some parameters of the figure and subplot spacing
            plt.rcParams['font.size'] = self.font_size
            plt.rcParams['figure.subplot.right'] = .98
            plt.rcParams['figure.subplot.bottom'] = .1
            plt.rcParams['figure.subplot.top'] = .93
            
            #set the font properties for the axis labels
            fontdict={'size':self.font_size+2, 'weight':'bold'}
            
            #set figure size according to what the plot will be.
            if self.plotnum==1 or self.plotnum==3:
                self.figsize = [ns*4,6]
                
            elif self.plotnum==2:
                self.figsize = [ns*8,6]
                
            #make a figure instance
            self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
                
            #make subplots as columns for all stations that need to be plotted
            gs0 = gridspec.GridSpec(1, ns)
                
            #space out the subplots
            gs0.update(hspace=.05, wspace=.05, left=.085)

#            gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs0[0])   
            ii = -1    
            for edi in self.edilst:
                #check to see if the file given is an edi file
                if edi.find('edi',-4,len(edi))>=0:
                    print 'Reading ' + edi
                    #create a Z object from the file            
                    zz = Z.Z(edi)
                    
                    #get the reistivity and phase object
                    rp = zz.getResPhase(thetar=self.rotz)
                    
                    #set x-axis limits from short period to long period
                    if self.xlimits==None:
                        self.xlimits = (10**(np.floor(np.log10(zz.period[0]))),
                                        10**(np.ceil(np.log10((zz.period[-1])))))
                    if self.phase_limits==None:
                        self.phase_limits = (0,89.9)
                    
                    ii+=1

                else:
                    raise ValueError('Could not read file: ' + edi)

                
                #create a grid to place the figures into, set to have 2 rows and 2 
                #columns to put any of the 4 components.  Make the phase plot
                #slightly shorter than the apparent resistivity plot and have the two
                #close to eachother vertically.
                
                #--> make figure for xy,yx components
                if self.plotnum==1 or self.plotnum==3:
                    #make subplots for each subplot in the figure
                    ax1 = gridspec.GridSpecFromSubplotSpec(2, 1, 
                                                           subplot_spec=gs0[ii],
                                                           height_ratios=[2,1.5], 
                                                           hspace=.01)
                    axr = self.fig.add_subplot(ax1[0,0])
                    axp = self.fig.add_subplot(ax1[1,0], sharex=axr)
                    
                    #place the y-coordinate labels in the same location for each axis
                    if ii==0:
                        axr.yaxis.set_label_coords(-.195, 0.5)
                        axp.yaxis.set_label_coords(-.195, 0.5)
                    
                #--> make figure for all 4 components
                elif self.plotnum==2:
                    #make subplots for each subplot in the figure
                    ax1 = gridspec.GridSpecFromSubplotSpec(2, 2, 
                                                           subplot_spec=gs0[ii],
                                                           height_ratios=[2,1.5], 
                                                           hspace=.01)
                    axr = self.fig.add_subplot(ax1[0,0])
                    axp = self.fig.add_subplot(ax1[1,0], sharex=axr)
                    axr2 = self.fig.add_subplot(ax1[0,1], sharex=axr)
                    axp2 = self.fig.add_subplot(ax1[1,1], sharex=axr)
                    
                    #place the y-coordinate labels in the same location for each axis
                    if ii==0:
                        axr.yaxis.set_label_coords(-.195, 0.5)
                        axp.yaxis.set_label_coords(-.195, 0.5)
                
                #=========Plot Z_xy and Z_yx=================================== 
                if self.plotnum==1 or self.plotnum==2:
                    #---------plot the apparent resistivity--------------------
                    #--> plot as error bars and just as points xy-blue, yx-red
                    #res_xy
                    ebxyr = axr.errorbar(zz.period, 
                                         rp.resxy, 
                                         marker=self.xy_marker, 
                                         ms=self.marker_size, 
                                         mfc=self.xy_mfc, 
                                         mec=self.xy_color, 
                                         mew=100./self.dpi, 
                                         ls=self.xy_ls, 
                                         yerr=rp.resxyerr, 
                                         ecolor=self.xy_color,
                                         capsize=self.marker_size)
                    
                    #res_yx                              
                    ebyxr = axr.errorbar(zz.period, 
                                         rp.resyx, 
                                         marker=self.yx_marker, 
                                         ms=self.marker_size, 
                                         mfc=self.yx_mfc,
                                         mec=self.yx_color, 
                                         mew=100./self.dpi,
                                         ls=self.yx_ls, 
                                         yerr=rp.resyxerr, 
                                         ecolor=self.yx_color,
                                         capsize=self.marker_size)
                                                  
                    #--> set axes properties
                    plt.setp(axr.get_xticklabels(), visible=False)
                    if ii==0:
                        axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                            fontdict=fontdict)
                        axr.legend((ebxyr[0],ebyxr[0]), 
                                    ('$E_x/B_y$','$E_y/B_x$'),
                                    loc=3, 
                                    markerscale=1, 
                                    borderaxespad=.01,
                                    labelspacing=.07, 
                                    handletextpad=.2, 
                                    borderpad=.02)
                    else:
                        plt.setp(axr.get_yticklabels(), visible=False)
                                    
                    axr.set_yscale('log')
                    axr.set_xscale('log')
                    axr.set_ylim(self.res_limits)
                    axr.set_xlim(self.xlimits)
                    axr.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                                  lw=.25)
                    
                    
                    #-----Plot the phase---------------------------------------
                    #phase_xy
                    ebxyp = axp.errorbar(zz.period, 
                                         rp.phasexy, 
                                         marker=self.xy_marker, 
                                         ms=self.marker_size, 
                                         mfc=self.xy_mfc,
                                         mec=self.xy_color, 
                                         mew=100./self.dpi,
                                         ls=self.xy_ls,
                                         yerr=rp.phasexyerr, 
                                         ecolor=self.xy_color,
                                         capsize=self.marker_size)
                                                  
                    #phase_yx: Note add 180 to place it in same quadrant as
                    #phase_xy
                    ebyxp = axp.errorbar(zz.period, 
                                         rp.phaseyx+180, 
                                         marker=self.yx_marker, 
                                         ms=self.marker_size, 
                                         mfc=self.yx_mfc, 
                                         mec=self.yx_color, 
                                         mew=100./self.dpi,
                                         ls=self.yx_ls, 
                                         yerr=rp.phaseyxerr, 
                                         ecolor=self.yx_color,
                                         capsize=self.marker_size)
            
                    #check the phase to see if any point are outside of [0:90]
                    if self.phase_limits==None:
                        if min(self.rp.phasexy)<0 or min(self.rp.phaseyx+180)<0:
                             pymin = min([min(self.rp.phasexy), 
                                          min(self.rp.phaseyx+180)])
                             if pymin>0:
                                pymin = 0
                        else:
                            pymin = 0
                        
                        if max(self.rp.phasexy)>90 or max(self.rp.phaseyx+180)>90:
                            pymax = min([max(self.rp.phasexy), 
                                         max(self.rp.phaseyx+180)])
                            if pymax<91:
                                pymax = 89.9
                        else:
                            pymax = 89.9
                            
                        self.phase_limits = (pymin, pymax)
                    
                    #--> set axes properties
                    axp.set_xlabel('Period (s)', fontdict)
                   
                    if ii==0:
                        axp.set_ylabel('Phase (deg)', fontdict)
                    
                    else:
                        plt.setp(axp.get_yticklabels(), visible=False)
                        
                    axp.set_xscale('log')
                    axp.set_ylim(self.phase_limits)        
                    axp.yaxis.set_major_locator(MultipleLocator(15))
                    axp.yaxis.set_minor_locator(MultipleLocator(5))
                    axp.grid(True, alpha=.25, 
                                  which='both', 
                                  color=(.25,.25,.25),
                                  lw=.25)
                                  
                    tklabels=[self.label_dict[tt] 
                              for tt in np.arange(np.log10(self.xlimits[0]),
                                                  np.log10(self.xlimits[1])+1)]
                    tklabels[0]=''
                    tklabels[-1]=''
                    
                    axp.set_xticklabels(tklabels,
                                        fontdict={'size':self.font_size})
                    
                    #====Plot the Z_xx, Z_yy components if desired=============
                    if self.plotnum==2:
                        #---------plot the apparent resistivity----------------
                        axr2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #res_xx
                        ebxxr = axr2.errorbar(zz.period, 
                                              rp.resxx, 
                                              marker=self.xy_marker, 
                                              ms=self.marker_size, 
                                              mfc=self.xy_mfc,
                                              mec=self.xy_color, 
                                              mew=100./self.dpi,
                                              ls=self.xy_ls, 
                                              yerr=rp.resxxerr, 
                                              ecolor=self.xy_color,
                                              capsize=self.marker_size)
                        
                        #res_yy                              
                        ebyyr = axr2.errorbar(zz.period, 
                                              rp.resyy, 
                                              marker=self.yx_marker, 
                                              ms=self.marker_size, 
                                              mfc=self.yx_mfc, 
                                              mec=self.yx_color, 
                                              mew=100./self.dpi, 
                                              ls=self.yx_ls, 
                                              yerr=rp.resyyerr, 
                                              ecolor=self.yx_color,
                                              capsize=self.marker_size)
            
                        #--> set axes properties
                        plt.setp(axr2.get_xticklabels(), visible=False)
                        
                        axr2.set_yscale('log')
                        axr2.set_xscale('log')
                        axr2.set_xlim(self.xlimits)
                        axr2.grid(True, 
                                  alpha=.25, 
                                  which='both', 
                                  color=(.25,.25,.25),
                                  lw=.25)
                        if ii==0:
                            axr2.legend((ebxxr[0], ebyyr[0]), 
                                        ('$E_x/B_x$','$E_y/B_y$'),
                                        loc=3, markerscale=1, 
                                        borderaxespad=.01,
                                        labelspacing=.07, 
                                        handletextpad=.2, 
                                        borderpad=.02)
                                            
    #                    else:
    #                        plt.setp(axr2.get_yticklabels(), visible=False)
                        
                        #-----Plot the phase-----------------------------------
                        
                        axp2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #phase_xx
                        ebxxp = axp2.errorbar(zz.period, 
                                              rp.phasexx, 
                                              marker=self.xy_marker, 
                                              ms=self.marker_size, 
                                              mfc=self.xy_mfc, 
                                              mec=self.xy_color, 
                                              mew=100./self.dpi, 
                                              ls=self.xy_ls, 
                                              yerr=rp.phasexxerr,
                                              ecolor=self.xy_color,
                                              capsize=self.marker_size)
                                                        
                        #phase_yy
                        ebyyp = axp2.errorbar(zz.period,
                                              rp.phaseyy, 
                                              marker=self.yx_marker,
                                              ms=self.marker_size, 
                                              mfc=self.yx_mfc,
                                              mec=self.yx_color, 
                                              mew=100./self.dpi, 
                                              ls=self.yx_ls, 
                                              yerr=rp.phaseyyerr, 
                                              ecolor=self.yx_color,
                                              capsize=self.marker_size)
                        
                        #--> set axes properties
                        axp2.set_xlabel('Period (s)', fontdict)
                        axp2.set_xscale('log')
                        axp2.set_ylim(ymin=-179.9, ymax=179.9)        
                        axp2.yaxis.set_major_locator(MultipleLocator(30))
                        axp2.yaxis.set_minor_locator(MultipleLocator(5))
                        axp2.set_xticklabels(tklabels,
                                            fontdict={'size':self.font_size})
                        axp2.grid(True, 
                                  alpha=.25, 
                                  which='both', 
                                  color=(.25,.25,.25),
                                  lw=.25) 
                                   
                
                #===Plot the Determinant if desired============================                          
                if self.plotnum==3:
                        
                    #res_det
                    ebdetr = axr.errorbar(zz.period, 
                                          rp.resdet, 
                                          marker=self.det_marker, 
                                          ms=self.marker_size, 
                                          mfc=self.det_mfc,
                                          mec=self.det_color, 
                                          mew=100./self.dpi, 
                                          ls=self.det_ls, 
                                          yerr=rp.resdeterr, 
                                          ecolor=self.det_color,
                                          capsize=self.marker_size)
                
                    #phase_det
                    ebdetp = axp.errorbar(zz.period, 
                                          rp.phasedet, 
                                          marker=self.det_marker, 
                                          ms=self.marker_size, 
                                          mfc=self.det_mfc, 
                                          mec=self.det_color, 
                                          mew=100./self.dpi, 
                                          ls=self.det_ls, 
                                          yerr=rp.phasedeterr, 
                                          ecolor=self.det_color,
                                          capsize=self.marker_size)
                    
#                    if ii==0:                                
#                        axr.legend((ebdetr[0]),
#                                    ('$\det(\mathbf{\hat{Z}})$'),
#                                    loc=3,
#                                    markerscale=1,
#                                    borderaxespad=.01,
#                                    labelspacing=.07,
#                                    handletextpad=.2,
#                                    borderpad=.02)
                                    
                    #--> set axes properties
                    plt.setp(axr.get_xticklabels(), visible=False)
                    if ii==0:
                        axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                            fontdict=fontdict)
                    else:
                        plt.setp(axr.get_yticklabels(), visible=False)
                                    
                    axr.set_yscale('log')
                    axr.set_xscale('log')
                    axr.set_ylim(self.res_limits)
                    axr.set_xlim(self.xlimits)
                    axr.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                                  lw=.25)
                                  
                    #--> set axes properties
                    axp.set_xlabel('Period (s)', fontdict)
                   
                    if ii==0:
                        axp.set_ylabel('Phase (deg)', fontdict)
                    
                    else:
                        plt.setp(axp.get_yticklabels(), visible=False)
                        
                    axp.set_xscale('log')
                    axp.set_ylim(self.phase_limits)        
                    axp.yaxis.set_major_locator(MultipleLocator(15))
                    axp.yaxis.set_minor_locator(MultipleLocator(5))
                    tklabels=[self.label_dict[tt] 
                              for tt in np.arange(np.log10(self.xlimits[0]),
                                                  np.log10(self.xlimits[1])+1)]
                    tklabels[0]=''
                    tklabels[-1]=''
                    
                    axp.set_xticklabels(tklabels,
                                        fontdict={'size':self.font_size})
                    axp.grid(True, alpha=.25, 
                                  which='both', 
                                  color=(.25,.25,.25),
                                  lw=.25)
                
                
                #make title and show
                axr.set_title(zz.station, fontsize=self.font_size, 
                              fontweight='bold')
                plt.show()
        
        #=====Plot all responses into one plot to compare changes=============
        if self.plotstyle=='compare':
            ns = len(self.edilst)
            
            cxy=[(0,0+float(cc)/ns,1-float(cc)/ns) for cc in range(ns)]
            cyx=[(1,float(cc)/ns,0) for cc in range(ns)]
            cdet=[(0,1-float(cc)/ns,0) for cc in range(ns)]
            mxy=['s','D','x','+','*','1','3','4']*5
            myx=['o','h','8','p','H',7,4,6]*5
            
            legendlst=[]
            stationlst=[]

            #set some parameters of the figure and subplot spacing
            plt.rcParams['font.size'] = self.font_size
            plt.rcParams['figure.subplot.right'] = .98
            plt.rcParams['figure.subplot.bottom'] = .1
            plt.rcParams['figure.subplot.top'] = .93
            
            #set the font properties for the axis labels
            fontdict={'size':self.font_size+2, 'weight':'bold'}
            
            #set figure size according to what the plot will be.
            if self.plotnum==1 or self.plotnum==3:
                self.figsize = [4,6]
                
            elif self.plotnum==2:
                self.figsize = [6,6]
                
            #make a figure instance
            self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
            
            #make subplots
            gs = gridspec.GridSpec(2, 2, height_ratios=[2,1.5], hspace=.01) 
            
             #--> make figure for xy,yx components
            if self.plotnum==1 or self.plotnum==3:
                
                #space out the subplots
                gs.update(hspace=.05, wspace=.15, left=.1)
                
                #--> create the axes instances
                #apparent resistivity axis
                axr = plt.subplot(gs[0,:])
                
                #phase axis which shares the x axis (period) with the resistivity
                axp = plt.subplot(gs[1,:], sharex=axr)
                
                #place the y-coordinate labels in the same location for each axis
                axr.yaxis.set_label_coords(-.0565, 0.5)
                axp.yaxis.set_label_coords(-.0565, 0.5)
                
            #--> make figure for all 4 components
            elif self.plotnum==2:
                #make figure instance
                self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
                
                #space out the subplots
                gs.update(hspace=.05, wspace=.15, left=.07)
                
                #--> create axes instance
                #apparent resistivity
                axr = plt.subplot(gs[0,0])
                axr2 = plt.subplot(gs[0,1])
                
                #phase axis which shares the x axis (period) with the resistivity
                axp = plt.subplot(gs[1,0],sharex=axr)
                axp2 = plt.subplot(gs[1,1],sharex=axr)
                
                #place the y-coordinate labels in the same location for each axis
                axr.yaxis.set_label_coords(-.075, 0.5)
                axp.yaxis.set_label_coords(-.075, 0.5)
    
            for ii,edi in enumerate(self.edilst):
                #check to see if the file given is an edi file
                if edi.find('edi',-4,len(edi))>=0:
                    print 'Reading ' + edi
                    #create a Z object from the file            
                    zz = Z.Z(edi)
                    
                    #get the reistivity and phase object
                    rp = zz.getResPhase(thetar=self.rotz)
                    
                    #set x-axis limits from short period to long period
                    if self.xlimits==None:
                        self.xlimits = (10**(np.floor(np.log10(zz.period[0]))),
                                        10**(np.ceil(np.log10((zz.period[-1])))))
                    if self.phase_limits==None:
                        self.phase_limits = (0,89.9)
                    
                    stationlst.append(zz.station)

                else:
                    raise ValueError('Could not read file: ' + edi)
                
                #=========Plot Z_xy and Z_yx=================================== 
                if self.plotnum==1 or self.plotnum==2:
                    #---------plot the apparent resistivity--------------------
                    #--> plot as error bars and just as points xy-blue, yx-red
                    #res_xy
                    ebxyr = axr.errorbar(zz.period, 
                                         rp.resxy, 
                                         marker=mxy[ii], 
                                         ms=self.marker_size, 
                                         mfc='None', 
                                         mec=cxy[ii], 
                                         mew=100./self.dpi, 
                                         ls=self.xy_ls, 
                                         yerr=rp.resxyerr, 
                                         ecolor=cxy[ii],
                                         capsize=self.marker_size)
                    
                    #res_yx                              
                    ebyxr = axr.errorbar(zz.period, 
                                         rp.resyx, 
                                         marker=myx[ii], 
                                         ms=self.marker_size, 
                                         mfc='None',
                                         mec=cyx[ii], 
                                         mew=100./self.dpi,
                                         ls=self.yx_ls, 
                                         yerr=rp.resyxerr, 
                                         ecolor=cyx[ii],
                                         capsize=self.marker_size)
                                                  

                    #-----Plot the phase---------------------------------------
                    #phase_xy
                    ebxyp = axp.errorbar(zz.period, 
                                         rp.phasexy, 
                                         marker=mxy[ii], 
                                         ms=self.marker_size, 
                                         mfc='None',
                                         mec=cxy[ii], 
                                         mew=100./self.dpi,
                                         ls=self.xy_ls,
                                         yerr=rp.phasexyerr, 
                                         ecolor=cxy[ii],
                                         capsize=self.marker_size)
                                                  
                    #phase_yx: Note add 180 to place it in same quadrant as
                    #phase_xy
                    ebyxp = axp.errorbar(zz.period, 
                                         rp.phaseyx+180, 
                                         marker=myx[ii], 
                                         ms=self.marker_size, 
                                         mfc='None',
                                         mec=cyx[ii], 
                                         mew=100./self.dpi,
                                         ls=self.yx_ls, 
                                         yerr=rp.phaseyxerr, 
                                         ecolor=cyx[ii],
                                         capsize=self.marker_size)
            

                    legendlst.append([ebxyr,ebyxr])             
                    
                    #====Plot the Z_xx, Z_yy components if desired=============
                    if self.plotnum==2:
                        #---------plot the apparent resistivity----------------
                        axr2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #res_xx
                        ebxxr = axr2.errorbar(zz.period, 
                                              rp.resxx, 
                                              marker=mxy[ii], 
                                              ms=self.marker_size, 
                                              mfc='None',
                                              mec=cxy[ii], 
                                              mew=100./self.dpi,
                                              ls=self.xy_ls, 
                                              yerr=rp.resxxerr, 
                                              ecolor=cxy[ii],
                                              capsize=self.marker_size)
                        
                        #res_yy                              
                        ebyyr = axr2.errorbar(zz.period, 
                                              rp.resyy, 
                                              marker=myx[ii], 
                                              ms=self.marker_size, 
                                              mfc='None', 
                                              mec=cyx[ii], 
                                              mew=100./self.dpi, 
                                              ls=self.yx_ls, 
                                              yerr=rp.resyyerr, 
                                              ecolor=cyx[ii],
                                              capsize=self.marker_size)
            
                                            
                        #-----Plot the phase-----------------------------------
                        
                        axp2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #phase_xx
                        ebxxp = axp2.errorbar(zz.period, 
                                              rp.phasexx, 
                                              marker=mxy[ii], 
                                              ms=self.marker_size, 
                                              mfc='None', 
                                              mec=cxy[ii], 
                                              mew=100./self.dpi, 
                                              ls=self.xy_ls, 
                                              yerr=rp.phasexxerr,
                                              ecolor=cxy[ii],
                                              capsize=self.marker_size)
                                                        
                        #phase_yy
                        ebyyp = axp2.errorbar(zz.period,
                                              rp.phaseyy, 
                                              marker=myx[ii],
                                              ms=self.marker_size, 
                                              mfc='None',
                                              mec=cyx[ii], 
                                              mew=100./self.dpi, 
                                              ls=self.yx_ls, 
                                              yerr=rp.phaseyyerr, 
                                              ecolor=cyx[ii],
                                              capsize=self.marker_size)
                        

                                   
                
                #===Plot the Determinant if desired============================                          
                if self.plotnum==3:
                        
                    #res_det
                    ebdetr = axr.errorbar(zz.period, 
                                          rp.resdet, 
                                          marker=mxy[ii], 
                                          ms=self.marker_size, 
                                          mfc='None',
                                          mec=cdet[ii], 
                                          mew=100./self.dpi, 
                                          ls=self.det_ls, 
                                          yerr=rp.resdeterr, 
                                          ecolor=cdet[ii],
                                          capsize=self.marker_size)
                
                    #phase_det
                    ebdetp = axp.errorbar(zz.period, 
                                          rp.phasedet, 
                                          marker=mxy[ii], 
                                          ms=self.marker_size, 
                                          mfc='None', 
                                          mec=cdet[ii], 
                                          mew=100./self.dpi, 
                                          ls=self.det_ls, 
                                          yerr=rp.phasedeterr, 
                                          ecolor=cdet[ii],
                                          capsize=self.marker_size)
                    
                    legendlst.append(ebdetr)
                                    
            #-------set axis properties----------------------------------------
            axr.set_yscale('log')
            axr.set_xscale('log')
            axr.set_ylim(self.res_limits)
            axr.set_xlim(self.xlimits)
            axr.grid(True, alpha=.25, which='both', color=(.25,.25,.25),
                          lw=.25)
            plt.setp(axr.get_xticklabels(), visible=False)
            
            axr.set_ylabel('App. Resistivity ($\Omega \cdot$m)',
                           fontdict=fontdict)
                          
            #check the phase to see if any point are outside of [0:90]
            if self.phase_limits==None:
                self.phase_limits = (0, 89.99)
            
            #--> set axes properties
            axp.set_xlabel('Period (s)', fontdict=fontdict)
            axp.set_ylabel('Phase (deg)',fontdict=fontdict)
            axp.set_xscale('log')
            axp.set_ylim(self.phase_limits)        
            axp.yaxis.set_major_locator(MultipleLocator(15))
            axp.yaxis.set_minor_locator(MultipleLocator(5))
            axp.grid(True, alpha=.25, 
                          which='both', 
                          color=(.25,.25,.25),
                          lw=.25)
            
            #make legend
            if self.plotnum==1:
                llst=[ll[0] for ll in legendlst]+[ll[1] for ll in legendlst]
                slst=[ss+'_xy' for ss in stationlst]+\
                     [ss+'_yx' for ss in stationlst]
                     
                axr.legend(llst, 
                        slst,
                        loc=3,
                        ncol=2,
                        markerscale=1, 
                        borderaxespad=.01,
                        labelspacing=.07, 
                        handletextpad=.2, 
                        borderpad=.02)
            elif self.plotnum==3:
                llst=[ll[0] for ll in legendlst]
                slst=[ss+'_det' for ss in stationlst]
                     
                axr.legend(llst, 
                        slst,
                        loc=3,
                        markerscale=1, 
                        borderaxespad=.01,
                        labelspacing=.07, 
                        handletextpad=.2, 
                        borderpad=.02)

            
            if self.plotnum==2:
                llst=[ll[0] for ll in legendlst]+[ll[1] for ll in legendlst]
                slst=[ss+'_xx' for ss in stationlst]+\
                     [ss+'_yy' for ss in stationlst]
                     
                axr2.legend(llst, 
                        slst,
                        loc=3,
                        ncol=2,
                        markerscale=1, 
                        borderaxespad=.01,
                        labelspacing=.07, 
                        handletextpad=.2, 
                        borderpad=.02)
                
                axr2.set_yscale('log')
                axr2.set_xscale('log')
                axr2.set_xlim(self.xlimits)
                axr2.grid(True, 
                          alpha=.25, 
                          which='both', 
                          color=(.25,.25,.25),
                          lw=.25)
                plt.setp(axr2.get_xticklabels(), visible=False)
                                
                #--> set axes properties
                axp2.set_xlabel('Period (s)', fontdict)
                axp2.set_xscale('log')
                axp2.set_ylim(ymin=-179.9, ymax=179.9)        
                axp2.yaxis.set_major_locator(MultipleLocator(30))
                axp2.yaxis.set_minor_locator(MultipleLocator(5))
                axp2.grid(True, 
                          alpha=.25, 
                          which='both', 
                          color=(.25,.25,.25),
                          lw=.25) 
        
        
            plt.show()
            
            

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
        
        return "Plots resistivity and phase for the different modes of the MT \n" +\
              "response for multiple sites. At the moment it supports the \n"+\
              "input of an .edi file. Other formats that will be supported\n"+\
              "are the impedance tensor and errors with an array of periods\n"+\
              "and .j format.\n"        

    
class PlotPhaseTensor(object):
    """
    PlotPhaseTensor will plot the components of the phase tensor as ellipses
    as well as the azimuth, phi_min and phi_max, ellipticity and skew.
    """
    pass

class PlotPhaseTensorMaps(object):
    """
    PlotPhaseTensorMaps will plot the phase tensor in map view.
    
    *The options should be to plot the induction arrows
    *Set a background image
    
    
    """
    
    pass

class PlotPhaseTensorPseudoSection(object):
    """
    PlotPhaseTensorPseudoSection will plot the phase tensor ellipses in a 
    pseudo section format 
    
    *The options should be to plot the induction arrows
    
    
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
    
    
    def __init__(self, filenamelst, ellipse_dict={}, offsetscaling=500, 
                 stationid=(0,4), title=None, cb_dict={}, linedir='ns', 
                 fignum=1, rotz=0, figsize=[6,6], dpi=300, indarrows='n', 
                 arrow_dict={}, tscale='period', font_size=7, plot_yn='y'):

        #----set attributes for the class-------------------------
        self.fn_list = filenamelst
        
        #--> set the ellipse properties
        
        #set default size to 2
        try:
            self.ellipse_size = ellipse_dict['size']
        except KeyError:
            self.ellipse_size = 2
        
        #set default colorby to phimin
        try:
            self.ellipse_colorby = ellipse_dict['colorby']
        except KeyError:
            self.ellipse_colorby = 'phimin'
        
        #set color range to 0-90
        try:
            self.ellipse_range = ellipse_dict['range']
        except KeyError:
            if self.ellipse_colorby=='beta' or self.ellipse_colorby=='beta_seg':
                self.ellipse_range = (-9,9,3)
            elif self.ellipse_colorby=='ellipticity':
                self.ellipse_range = (0,1,.1)
            else:
                self.ellipse_range = (0,90,5)
            
        #set colormap to yellow to red
        try:
            self.ellipse_cmap = ellipse_dict['cmap']
        except KeyError:
            if self.ellipse_colorby=='beta':
                self.ellipse_cmap = 'mt_bl2wh2rd'
            else:
                self.ellipse_cmap = 'mt_yl2rd'
            
        #--> set colorbar properties
        #set orientation to horizontal
        try:
            self.cb_orientation = cb_dict['orientation']
        except KeyError:
            self.cb_orientation = 'vertical'
        
        #set the position to middle outside the plot            
        try:
            self.cb_poistion = cb_dict['position']
        except KeyError:
            self.cb_position = (.92, .375, .025, .3)
            
        #--> set plot properties
        self.offsetscaling = offsetscaling
        self.dpi = dpi
        self.font_size = font_size
        self.tscale = tscale
        self.rotz = rotz
        self.figsize = figsize
        self.fignum = fignum
        self.linedir = linedir
        self.stationid = stationid
        self.title = title
        
        #--> set induction arrow properties 
        self.indarrows = indarrows
        
        if self.indarrows.find('y')==0:
            #set arrow length
            try:
                self.arrow_size = arrow_dict['size']
            except KeyError:
                self.arrow_size = 5
                
            #set head length
            try:
                self.arrow_head_length = arrow_dict['head_length']
            except KeyError:
                self.arrow_head_length = .5*self.ellipse_size
                
            #set head width
            try:
                self.arrow_head_width = arrow_dict['head_width']
            except KeyError:
                self.arrow_head_width = .2*self.ellipse_size
                
            #set line width
            try:
                self.arrow_lw = arrow_dict['lw']
            except KeyError:
                self.arrow_lw = .5*self.ellipse_size
                
            #set real color to black
            try:
                self.arrow_color_real = arrow_dict['color'][0]
            except KeyError:
                self.arrow_color_real = 'k'
                
            #set imaginary color to black
            try:
                self.arrow_color_imag = arrow_dict['color'][1]
            except KeyError:
                self.arrow_color_imag = 'b'
                
            #set threshold of induction arrows to plot
            try:
                self.arrow_threshold = arrow_dict['threshold']
            except KeyError:
                self.arrow_threshold = 1
            
        #--> plot if desired
        self.plot_yn = plot_yn
        if self.plot_yn=='y':
            self.plot()
            
        
    def plot(self):
        """
        plots the phase tensor pseudo section
        """
            
        plt.rcParams['font.size']=self.font_size
        plt.rcParams['figure.subplot.left']=.08
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.06
        plt.rcParams['figure.subplot.top']=.96
        plt.rcParams['figure.subplot.wspace']=.55
        plt.rcParams['figure.subplot.hspace']=.70
        
        #create a plot instance
        self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
        self.ax = self.fig.add_subplot(1,1,1,aspect='equal')
        
        #create empty lists to put things into
        self.stationlst = []
        self.offsetlst = []
        minlst = []
        maxlst = []
        
        #plot phase tensor ellipses
        for ii,fn in enumerate(self.fn_list):
            imp = Z.Z(fn)
            self.stationlst.append(imp.station[self.stationid[0]:self.stationid[1]])
            
            #set the an arbitrary origin to compare distance to all other 
            #stations.
            if ii==0:
                east0 = imp.lon
                north0 = imp.lat
                offset = 0.0
            else:
                east = imp.lon
                north = imp.lat
                if self.linedir=='ew': 
                    if east0<east:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif east0>east:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset = 0
                elif self.linedir=='ns':
                    if north0<north:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif north0>north:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset=0
                        
            self.offsetlst.append(offset)
            
            #get phase tensor elements and flip so the top is small periods/high 
            #frequencies
            pt = imp.getPhaseTensor(thetar=self.rotz)
            periodlst = imp.period[::-1]
            phimax = pt.phimax[::-1]
            phimin = pt.phimin[::-1]
            azimuth = pt.azimuth[::-1]
            
            #if there are induction arrows, flip them as pt
            if self.indarrows.find('y')==0:
                tip = imp.getTipper(thetar=self.rotz)
                tmr = tip.magreal[::-1]
                tmi = tip.magimag[::-1]
                tar = tip.anglereal[::-1]
                tai = tip.angleimag[::-1]
                
                aheight = self.arrow_head_length 
                awidth = self.arrow_head_width
                alw = self.arrow_lw
                
            #get the properties to color the ellipses by
            if self.ellipse_colorby=='phiminang' or \
               self.ellipse_colorby=='phimin':
                colorarray = pt.phiminang[::-1]
                
            elif self.ellipse_colorby=='phidet':
                 colorarray = np.sqrt(abs(pt.phidet[::-1]))*(180/np.pi)
                
            elif self.ellipse_colorby=='beta' or\
                 self.ellipse_colorby=='beta_seg':
                colorarray = pt.beta[::-1]
                
            elif self.ellipse_colorby=='ellipticity':
                colorarray = pt.ellipticity[::-1]
                
            else:
                raise NameError(self.ellipse_colorby+' is not supported')
            
            #get the number of periods
            n = len(periodlst)
            
            #get min and max of the color array for scaling later
            minlst.append(min(colorarray))
            maxlst.append(max(colorarray))
            
            #set local parameters with shorter names
            es = self.ellipse_size
            ck = self.ellipse_colorby
            cmap = self.ellipse_cmap
            ckmin = self.ellipse_range[0]
            ckmax = self.ellipse_range[1]
            ckseg = self.ellipse_range[2]
            nseg = (ckmax-ckmin)/(2*ckseg)
            
    
            for jj in range(n):
                
                #make sure the ellipses will be visable
                eheight = phimin[jj]/phimax[jj]*es
                ewidth = phimax[jj]/phimax[jj]*es
            
                #create an ellipse scaled by phimin and phimax and oriented along
                #the azimuth    
                ellipd=patches.Ellipse((offset*self.offsetscaling, 3*jj),
                                       width=ewidth,
                                       height=eheight,
                                       angle=azimuth[jj])
                
                #get face color info
                if ck=='phiminang' or  ck=='phimin':
                    cvar = (colorarray[jj]-ckmin)/(ckmax-ckmin)
                    
                elif ck=='phidet':
                    cvar = (colorarray[jj]-ckmin)/(ckmax-ckmin)
                    
                elif ck=='beta':
                    cvar = 2*colorarray[jj]/(ckmax-ckmin)
                    
                elif ck=='beta_seg':
                    cvar = 2*nseg*(np.round(colorarray[jj]/ckseg)/(ckmax-ckmin))
                    if jj>70:
                        print imp.station,imp.frequency[-jj],cvar,colorarray[jj]

                elif ck=='ellipticity':
                    cvar = (colorarray[jj]-ckmin)/(ckmax-ckmin)
                    
                else:
                    raise NameError('color key '+ck+' not supported')
                    
                
                #set facecolor depending on the colormap
                #yellow to red
                if cmap=='mt_yl2rd':
                    if abs(cvar)>1:
                        ellipd.set_facecolor((1,0,0))
                    elif cvar<0:
                        ellipd.set_facecolor((0,0,0))
                    else:
                        ellipd.set_facecolor((1,1-abs(cvar),.1))
                        
                #white to blue
                elif cmap=='mt_wh2bl':
                    if abs(cvar)>1:
                        ellipd.set_facecolor((0,0,0))
                    else:
                        ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                        
                #blue to white to red
                elif cmap=='mt_bl2wh2rd':
                    if cvar<0 and cvar>-1:
                        ellipd.set_facecolor((1+cvar,1+cvar,1))
                    elif cvar<-1:
                        ellipd.set_facecolor((0,0,1))
                    elif cvar>=0 and cvar<1:
                        ellipd.set_facecolor((1,1-cvar,1-cvar))
                    elif cvar>1:
                        ellipd.set_facecolor((1,0,0))
                        
                #blue to white to red
                elif cmap=='mt_seg_bl2wh2rd':
                    if cvar<0 and cvar>-1:
                        ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                    elif cvar<=-1:
                        ellipd.set_facecolor((0,0,1))
                    elif cvar>0 and cvar<1:
                        ellipd.set_facecolor((1,1-cvar,1-cvar))
                    elif cvar>=1:
                        ellipd.set_facecolor((1,0,0))
                        
                #blue to yellow to red
                elif cmap=='mt_bl2yl2rd':
                    if cvar<0 and cvar>-1:
                        ellipd.set_facecolor((abs(cvar),abs(cvar),1-abs(cvar)))
                    elif cvar<-1:
                        ellipd.set_facecolor((0,0,1))
                    elif cvar>0 and cvar<1:
                        ellipd.set_facecolor((1,1-abs(cvar),.01))
                    elif cvar>1:
                        ellipd.set_facecolor((1,0,0))
               
                else:
                    raise NameError('Colormap '+cmap+' is not supported')
                    
                #===add the ellipse to the plot==========
                self.ax.add_artist(ellipd)
                
                
                #--------- Add induction arrows if desired --------------------
                if self.indarrows.find('y')==0:
                    
                    
                    #--> plot real tipper
                    if self.indarrows=='yri' or self.indarrows=='yr':
                        txr = tmr[jj]*np.cos(tar[jj]*np.pi/180)*self.arrow_size
                        tyr = tmr[jj]*np.sin(tar[jj]*np.pi/180)*self.arrow_size
                        
                        maxlength = np.sqrt((txr/self.arrow_size)**2+\
                                            (tyr/self.arrow_size)**2)
                        if maxlength>self.arrow_threshold:
                            pass
                        else:
                            self.ax.arrow(offset*self.offsetscaling, 
                                          3*jj, 
                                          txr,
                                          tyr,
                                          lw=alw,
                                          facecolor=self.arrow_color_real,
                                          edgecolor=self.arrow_color_real,
                                          length_includes_head=False,
                                          head_width=awidth,
                                          head_length=aheight)
                                      
                    #--> plot imaginary tipper
                    if self.indarrows=='yri' or self.indarrows=='yi':
                        txi = tmi[jj]*np.cos(tai[jj]*np.pi/180)*self.arrow_size
                        tyi = tmi[jj]*np.sin(tai[jj]*np.pi/180)*self.arrow_size
                        
                        maxlength = np.sqrt((txi/self.arrow_size)**2+\
                                            (tyi/self.arrow_size)**2)
                        if maxlength>self.arrow_threshold:
                            pass
                        else:
                            self.ax.arrow(offset*self.offsetscaling,
                                          3*jj,
                                          txi,
                                          tyi,
                                          lw=alw,
                                          facecolor=self.arrow_color_imag,
                                          edgecolor=self.arrow_color_imag,
                                          length_includes_head=False,
                                          head_width=awidth,
                                          head_length=aheight)
        
        #--> Set plot parameters                
        self.offsetlst = np.array(self.offsetlst)
        
        #set x-limits
        self.ax.set_xlim(self.offsetlst.min()*self.offsetscaling-es*2,
                         self.offsetlst.max()*self.offsetscaling+es*2)
        #set y-limits    
        self.ax.set_ylim(-5,n*3+5)
        
        #set y-axis labels for either period or frequency
        if self.tscale=='period':
            yticklabels=['{0:.5g}'.format(periodlst[ii]) 
                        for ii in np.arange(start=0,stop=n,step=2)]
            self.ax.set_ylabel('Period (s)',
                               fontsize=self.font_size,
                               fontweight='bold')
                               
        elif self.tscale=='frequency':
            yticklabels=['{0:>5}'.format('{0:.5g}'.format(1./periodlst[ii])) 
                         for ii in np.arange(start=0,stop=n,step=2)]
            self.ax.set_ylabel('Frequency (Hz)',
                               fontsize=self.font_size,
                               fontweight='bold')
        #set x-axis label                       
        self.ax.set_xlabel('Station',
                           fontsize=self.font_size+2,
                           fontweight='bold')
        
        #set title of the plot
        if self.title==None:
            pass
        else:
            self.ax.set_title(self.title, fontsize=self.font_size+2)
         
        #set tick locations
        plt.yticks(np.arange(start=0,stop=3*n,step=6), yticklabels)
        plt.xticks(self.offsetlst*self.offsetscaling, self.stationlst)
        
        #make a legend for the induction arrows
        if self.indarrows.find('y')==0:
            if self.indarrows=='yri':
                treal=self.ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,
                                   self.arrow_dict['color'])
                timag=self.ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,
                                   'b')
                self.ax.legend([treal[0],timag[0]],
                               ['Tipper_real','Tipper_imag'],
                               loc='lower right',
                               prop={'size':self.font_size-1,'weight':'bold'},
                               ncol=2,
                               markerscale=.5,
                               borderaxespad=.005,
                               borderpad=.25)
                          
            elif self.indarrows=='yr':
                treal = self.ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,
                                     self.arrow_color_real)
                self.ax.legend([treal[0]],
                               ['Tipper_real'],
                               loc='lower right',
                               prop={'size':self.font_size-1,'weight':'bold'},
                               ncol=2,
                               markerscale=.5,
                               borderaxespad=.005,
                               borderpad=.25)
                          
            elif self.indarrows=='yi':
                timag = self.ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,
                                     self.arrow_color_imag)
                self.ax.legend([timag[0]],
                               ['Tipper_imag'],
                               loc='lower right',
                               prop={'size':self.font_size-1,'weight':'bold'},
                               ncol=2,
                               markerscale=.5,
                               borderaxespad=.005,
                               borderpad=.25)
        
        #put a grid on the plot
        self.ax.grid(alpha=.25,which='both',color=(.25,.25,.25))
        
        print '-'*25
        print ck+' min = {0:.2f}'.format(min(minlst))
        print ck+' max = {0:.2f}'.format(max(maxlst))
        print '-'*25
        
        #make a colorbar with appropriate colors             
        self.ax2 = self.fig.add_axes(self.cb_position)
        
        if cmap=='mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax+ckseg, ckseg)
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)
            self.cb=mcb.ColorbarBase(self.ax2,
                                     cmap=mt_seg_bl2wh2rd,
                                     norm=norms,
                                     orientation=self.cb_orientation,
                                     ticks=bounds)
        else:
            self.cb=mcb.ColorbarBase(self.ax2,
                                     cmap=cmapdict[cmap],
                                     norm=colors.Normalize(vmin=ckmin,
                                                           vmax=ckmax),
                                     orientation=self.cb_orientation)

        #label the color bar accordingly
        self.cb.set_label(ckdict[ck],
                          fontdict={'size':self.font_size,'weight':'bold'})
            
        #place the label in the correct location                   
        if self.cb_orientation=='horizontal':
            self.cb.ax.xaxis.set_label_position('top')
            self.cb.ax.xaxis.set_label_coords(.5,1.3)
            
            
        elif self.cb_orientation=='vertical':
            self.cb.ax.yaxis.set_label_position('right')
            self.cb.ax.yaxis.set_label_coords(1.25,.5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis='y',direction='in')
        
        plt.show()
        
    def writeTextFiles(self):
        """
        This will write text files for all the phase tensor parameters
        """
        pass
    
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
        
        return "Plots pseudo section of phase tensor ellipses" 

class PlotStrike(object):
    """
    PlotStrike will plot the strike estimated from the invariants, phase tensor
    and the tipper in either a rose diagram of xy plot
    
    """    

    pass

class PlotTimeSeries(object):
    """
    plot the time series
    """
  

#    
#def plotPTpseudoSection(filenamelst,colorkey='phimin',esize=2,
#                        offsetscaling=.005,colorkeymm=(0,90),stationid=[0,4],
#                        title=None,cbshrink=.8,linedir='ns',fignum=1,rotz=0,
#                        pxy=[8,8],dpi=300,indarrows='n',ascale=5,
#                        cmap='ptcmap',tscale='period'):
#    
#    """
#    plotPTpseudoSection(filenamelst,colorkey='beta',esize=2,offsetscaling=
#    .005) will plot a pseudo section of phase tensor ellipses given a list of 
#    full path filenames. 
#    
#    Arguments:
#    ----------
#    
#        **filenamelst** : list of strings
#                          full paths to .edi files to plot
#        
#        **colorkey** : [ 'phimin' | 'beta' | 'ellipticity' | 'phidet' ]
#                       fill color of the ellipses and can be:
#                       * 'phimin'      -> minimum phase
#                       * 'beta'        -> phase tensor skew angle
#                       * 'ellipticity' -> phase tensor ellipticity
#                       * 'phidet'      -> determinant of the phase tensor
#                       * *Default* is 'phimin'
#                       
#        **colorkeymm** : tuple (min,max)
#                        min and max of colorkey to which colorbar is normalized
#                        to.  In degrees except for ellipticity which is already
#                        normalized.
#                        
#        **esize** : float
#                    size of ellipse in map units 
#        
#        **offsetscaling** : float
#                            is a factor that scales the distance from one 
#                            station to the next to make the plot readable.
#                            *Default* is 0.005
#        **linedir** : [ 'ns' | 'ew' ]
#                      predominant direction of profile line
#                      * 'ns' -> North-South Line
#                      * 'ew' -> East-West line
#                      * *Default* is 'ns'
#        
#        **stationid** : tuple or list 
#                        start and stop of station name indicies.  
#                        ex: for MT01dr stationid=(0,4) will be MT01
#        
#        **rotz** : float
#                   angle in degrees to rotate the data clockwise positive.
#                   *Default* is 0
#        
#        **title** : string
#                    figure title added to Phase Tensor + title
#                    
#        **cbshrink** : float
#                       percent to shrink the color bar to fit the image better
#                       *Default* is 0.8 -> 80 percent
#                       
#        **fignum** : int
#                     figure number.  *Default* is 1
#        
#        **indarrow** : [ 'yri' | 'yr' | 'yi' | 'n' ]
#                        * 'yri' to plot induction both real and imaginary 
#                           induction arrows 
#                        * 'yr' to plot just the real induction arrows
#                        * 'yi' to plot the imaginary induction arrows
#                        * 'n' to not plot them
#                        * *Default* is 'n'                        
#                        **Note: convention is to point towards a conductor **
#                         
#        **ascale** : float
#                     scaling factor to make induction arrows bigger
#        
#        **cmap** :[ 'ptcmap' | 'ptcmap3' | 'skcmap' | 'skcmap2' | 'rtcmap' ]
#                  color map of ellipse facecolor.  So far the colormaps are:
#                      * 'ptcmap'  -> yellow (low phase) to red (high phase)
#                      * 'ptcmap3' -> white (low numbers) to blue (high numbers)
#                      * 'skcmap'  -> blue to yellow to red
#                      * 'skcmap2' -> blue to white to red
#                      * 'rtcmap'  -> blue to purple to red
#                      
#        **tscale** : [ 'period' | 'frequency' ]
#                     * 'period'    -> plot vertical scale in period
#                     * 'frequency' -> plot vertical scale in frequency
#                     
#    :Example: ::
#        
#        >>> import mtpy.imaging.mtplottools as mtplot
#        >>> import os
#        >>> edipath = r"/home/EDIfiles"
#        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
#        >>> ...       if edi.find('.edi')]
#        >>> # color by phimin with a range of 20-70 deg
#        >>> mtplot.plotPTpseudoSection(edilst,colorkeymm=(20,70))
#        >>> 
#        >>> # add induction arrows to skew angle with range (-5,5)
#        >>> mtplot.plotPTpseudoSection(edilst,colorkey='beta',
#        >>> ...                        colorkeymm=(-5,5),indarrows='yri')
#
#    """
#    
#    fs=int(dpi/30.)
#    plt.rcParams['font.size']=fs
#    plt.rcParams['figure.subplot.left']=.08
#    plt.rcParams['figure.subplot.right']=.98
#    plt.rcParams['figure.subplot.bottom']=.06
#    plt.rcParams['figure.subplot.top']=.96
#    plt.rcParams['figure.subplot.wspace']=.55
#    plt.rcParams['figure.subplot.hspace']=.70
#    #plt.rcParams['font.family']='helvetica'
#    
#    ckmin=colorkeymm[0]
#    ckmax=colorkeymm[1]
#    
#    #create a plot instance
#    fig=plt.figure(fignum,pxy,dpi=300)
#    ax=fig.add_subplot(1,1,1,aspect='equal')
#    stationlst=[]
#    offsetlst=[]
#    minlst=[]
#    maxlst=[]
#    #plot phase tensor ellipses
#    for ii,fn in enumerate(filenamelst):
#        imp=Z.Z(fn)
#        stationlst.append(imp.station[stationid[0]:stationid[1]])
#        zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
#        
#        if ii==0:
#            east0=east
#            north0=north
#            offset=0.0
#        else:
#            if linedir=='ew': 
#                if east0<east:
#                    offset=np.sqrt((east0-east)**2+(north0-north)**2)
#                elif east0>east:
#                    offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
#                else:
#                    offset=0
#            elif linedir=='ns':
#                if north0<north:
#                    offset=np.sqrt((east0-east)**2+(north0-north)**2)
#                elif north0>north:
#                    offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
#                else:
#                    offset=0
#        offsetlst.append(offset)
#        #get phase tensor elements and flip so the top is small periods/high 
#        #frequencies
#        pt=imp.getPhaseTensor(thetar=rotz)
#        periodlst=imp.period[::-1]
#        phimax=pt.phimax[::-1]
#        phimin=pt.phimin[::-1]
#        azimuth=pt.azimuth[::-1]
#        if colorkey=='phiminang' or colorkey=='phimin':
#            colorarray=pt.phiminang[::-1]
#        if colorkey=='phidet':
#            colorarray=np.sqrt(abs(pt.phidet[::-1]))*(180/np.pi)
#        if colorkey=='beta':
#            colorarray=pt.beta[::-1]
#        if colorkey=='ellipticity':
#            colorarray=pt.ellipticity[::-1]
#            
#        n=len(periodlst)
#        minlst.append(min(colorarray))
#        maxlst.append(max(colorarray))
#
#        for jj in range(n):
#            
#            #make sure the ellipses will be visable
#            eheight=phimin[jj]/phimax[jj]*esize
#            ewidth=phimax[jj]/phimax[jj]*esize
#        
#                
#            #create an ellipse scaled by phimin and phimax and oriented along
#            #the azimuth    
#            ellipd=Ellipse((offset*offsetscaling,3*jj),width=ewidth,
#                          height=eheight,
#                          angle=azimuth[jj])
#            ax.add_artist(ellipd)
#            
#            #get face color info
#            if colorkey=='phiminang' or  colorkey=='phimin':
#                cvar=(pt.phiminang[jj]-ckmin)/(ckmax-ckmin)
#            elif colorkey=='phidet':
#                cvar=(pt.phidet[jj]-ckmin)/(ckmax-ckmin)
#            elif colorkey=='beta':
#                cvar=(pt.beta[jj]-abs(ckmin))/(ckmax-ckmin)
#            elif colorkey=='ellipticity':
#                cvar=(pt.ellipticity[jj]-ckmin)/(ckmax-ckmin)
#            else:
#                raise NameError('color key '+colorkey+' not supported')
#                
#            
#            #set facecolor depending on the colormap
#            #yellow to red
#            if cmap=='ptcmap':
#                if abs(cvar)>1:
#                    ellipd.set_facecolor((1,0,0))
#                elif cvar<0:
#                    ellipd.set_facecolor((1-abs(cvar),1,abs(cvar)))
#                else:
#                    ellipd.set_facecolor((1,1-abs(cvar),.1))
#            #white to blue
#            elif cmap=='ptcmap3':
#                if abs(cvar)>1:
#                    ellipd.set_facecolor((0,0,0))
#                else:
#                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
#            #blue to yellow to red
#            elif cmap=='skcmap2':
#                if cvar<0 and cvar>-1:
#                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
#                elif cvar<-1:
#                    ellipd.set_facecolor((0,0,1))
#                elif cvar>0 and cvar<1:
#                    ellipd.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
#                elif cvar>1:
#                    ellipd.set_facecolor((1,0,0))
#            #blue to white to red
#            elif cmap=='skcmap':
#                if cvar<0 and cvar>-1:
#                    ellipd.set_facecolor((abs(cvar),abs(cvar),1-abs(cvar)))
#                elif cvar<-1:
#                    ellipd.set_facecolor((0,0,1))
#                elif cvar>0 and cvar<1:
#                    ellipd.set_facecolor((1,1-abs(cvar),.01))
#                elif cvar>1:
#                    ellipd.set_facecolor((1,0,0))
#                    
#            if indarrows.find('y')==0:
#                tip=imp.getTipper(thetar=rotz)
#                aheight=.5*esize
#                awidth=.2*esize
#                #plot real tipper
#                if indarrows=='yri' or indarrows=='yr':
#                    txr=tip.magreal[jj]*np.cos(tip.anglereal[jj]*np.pi/180)*\
#                        esize*5
#                    tyr=tip.magreal[jj]*np.sin(tip.anglereal[jj]*np.pi/180)*\
#                        esize*5
#
#                    ax.arrow(offset*offsetscaling,3*jj,txr,tyr,lw=.75*awidth,
#                         facecolor='k',edgecolor='k',
#                         length_includes_head=False,
#                         head_width=awidth,head_length=aheight)
#                #plot imaginary tipper
#                if indarrows=='yri' or indarrows=='yi':
#                    txi=tip.magimag[jj]*np.cos(tip.angleimag[jj]*np.pi/180)*\
#                        esize*5
#                    tyi=tip.magimag[jj]*np.sin(tip.angleimag[jj]*np.pi/180)*\
#                        esize*5
#
#                    ax.arrow(offset*offsetscaling,3*jj,txi,tyi,lw=.75*awidth,
#                             facecolor='b',edgecolor='b',
#                             length_includes_head=False,
#                             head_width=awidth,head_length=aheight)
#                    
#    offsetlst=np.array(offsetlst)
#    ax.set_xlim(min(offsetlst)*offsetscaling-4,max(offsetlst)*offsetscaling+4)
#    ax.set_ylim(-5,n*3+5)
#    if tscale=='period':
#        yticklabels=['{0:.3g}'.format(periodlst[ii]) for ii in np.arange(start=0,stop=n,
#                     step=2)]
#        ax.set_ylabel('Period (s)',fontsize=fs+5,fontweight='bold')
#    elif tscale=='frequency':
#        yticklabels=['{0:.4g}'.format(1./periodlst[ii]) for ii in np.arange(start=0,stop=n,
#                     step=2)]
#        ax.set_ylabel('Frequency (Hz)',fontsize=fs+5,fontweight='bold')
#    ax.set_xlabel('Station',fontsize=fs+5,fontweight='bold')
#    
#    if title==None:
#        pass
#    else:
#        ax.set_title(title,fontsize=fs+4)
#    plt.yticks(np.arange(start=0,stop=3*n,step=6),yticklabels)
#    plt.xticks(np.array(offsetlst)*offsetscaling,stationlst)
#    
#    if indarrows.find('y')==0:
#        if indarrows=='yri':
#            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
#            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
#            ax.legend([treal[0],timag[0]],['Tipper_real','Tipper_imag'],
#                      loc='lower right',
#                      prop={'size':10,'weight':'bold'},
#                      ncol=2,markerscale=.5,borderaxespad=.005,
#                      borderpad=.25)
#        elif indarrows=='yr':
#            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
#            ax.legend([treal[0]],['Tipper_real'],
#                      loc='lower right',
#                      prop={'size':10,'weight':'bold'},
#                      ncol=2,markerscale=.5,borderaxespad=.005,
#                      borderpad=.25)
#        elif indarrows=='yi':
#            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
#            ax.legend([timag[0]],['Tipper_imag'],
#                      loc='lower right',
#                      prop={'size':10,'weight':'bold'},
#                      ncol=2,markerscale=.5,borderaxespad=.005,
#                      borderpad=.25)
#    
#    ax.grid(alpha=.25,which='both')
#    
#    print 'Colorkey min = ',min(minlst)
#    print 'Colorkey max = ',max(maxlst)
#    
#    #make a colorbar with appropriate colors             
#    ax2=make_axes(ax,shrink=.8)
#    if cmap=='ptcmap': 
#        cb=ColorbarBase(ax2[0],cmap=ptcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='ptcmap3': 
#        cb=ColorbarBase(ax2[0],cmap=ptcmap3,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='skcmap': 
#        cb=ColorbarBase(ax2[0],cmap=skcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='skcmap2': 
#        cb=ColorbarBase(ax2[0],cmap=skcmap2,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='rtcmap': 
#        cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#
#    #label the color bar accordingly
#    if colorkey=='phimin' or colorkey=='phiminang':
#        cb.set_label('$\Phi_{min}$ (deg)',fontdict={'size':10,'weight':'bold'})
#    elif colorkey=='beta':
#        cb.set_label('Skew (deg)',fontdict={'size':10,'weight':'bold'})
#    elif colorkey=='phidet':
#        cb.set_label('Det{$\Phi$} (deg)',fontdict={'size':10,'weight':'bold'})
#    elif colorkey=='ellipticity':
#        cb.set_label('Ellipticity',fontdict={'size':10,'weight':'bold'})
#    
#    plt.show()
#
#def plotRTpseudoSection(filenamelst,colorkey='rhodet',esize=2,
#                        offsetscaling=.005,colorkeymm=[0,90],stationid=[0,4],
#                        title=None,cbshrink=.8,linedir='ns',fignum=1,rotz=0,
#                        yscale='period',pxy=[8,8],dpi=300):
#    
#    """
#    plotRTpseudoSection(filenamelst,colorkey='beta',esize=2,offsetscaling=
#    .005) will plot a pseudo section of resistivity tensor ellipses given a list of 
#    full path filenames. (Weckmann et al. 2002)
#    
#    **filenamelst** : list of strings
#                          full paths to .edi files to plot
#        
#        **colorkey** : [ 'rhodet' ]
#                       fill color of the ellipses and can be:
#                       * 'rhodet'      -> minimum phase
#                       * more to come
#                       * *Default* is 'rhodet'
#                       
#        **colorkeymm** : tuple (min,max)
#                        min and max of colorkey to which colorbar is normalized
#                        to.  In log10 resistivity
#                        
#        **esize** : float
#                    size of ellipse in map units 
#        
#        **offsetscaling** : float
#                            is a factor that scales the distance from one 
#                            station to the next to make the plot readable.
#                            *Default* is 0.005
#        **linedir** : [ 'ns' | 'ew' ]
#                      predominant direction of profile line
#                      * 'ns' -> North-South Line
#                      * 'ew' -> East-West line
#                      * *Default* is 'ns'
#        
#        **stationid** : tuple or list 
#                        start and stop of station name indicies.  
#                        ex: for MT01dr stationid=(0,4) will be MT01
#        
#        **rotz** : float
#                   angle in degrees to rotate the data clockwise positive.
#                   *Default* is 0
#        
#        **title** : string
#                    figure title added to Phase Tensor + title
#                    
#        **cbshrink** : float
#                       percent to shrink the color bar to fit the image better
#                       *Default* is 0.8 -> 80 percent
#                       
#        **fignum** : int
#                     figure number.  *Default* is 1
#                      
#        **yscale** : [ 'period' | 'frequency' ]
#                     * 'period'    -> plot vertical scale in period
#                     * 'frequency' -> plot vertical scale in frequency
#                     
#    :Example: ::
#        
#        >>> import mtpy.imaging.mtplottools as mtplot
#        >>> import os
#        >>> edipath = r"/home/EDIfiles"
#        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
#        >>> ...       if edi.find('.edi')]
#        >>> # color by rhodet with a range of 0-4 log10 resistivity
#        >>> mtplot.plotRTpseudoSection(edilst,colorkeymm=(0,4))
#    
#    """
#    
#    fs=int(dpi/30.)
#    plt.rcParams['font.size']=fs
#    plt.rcParams['figure.subplot.left']=.08
#    plt.rcParams['figure.subplot.right']=.98
#    plt.rcParams['figure.subplot.bottom']=.06
#    plt.rcParams['figure.subplot.top']=.96
#    plt.rcParams['figure.subplot.wspace']=.55
#    plt.rcParams['figure.subplot.hspace']=.70
#    #plt.rcParams['font.family']='helvetica'
#    
#    
#    #create a plot instance
#    fig=plt.figure(fignum,pxy,dpi=dpi)
#    ax=fig.add_subplot(1,1,1,aspect='equal')
#    stationlst=[]
#    offsetlst=[]
#    minlst=[]
#    maxlst=[]
#    #plot phase tensor ellipses
#    for ii,fn in enumerate(filenamelst):
#        imp=Z.Z(fn)
#        stationlst.append(imp.station[stationid[0]:stationid[1]])
#        zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
#        
#        if ii==0:
#            east0=east
#            north0=north
#            offset=0.0
#        else:
#            if linedir=='ew': 
#                if east0<east:
#                    offset=np.sqrt((east0-east)**2+(north0-north)**2)
#                elif east0>east:
#                    offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
#                else:
#                    offset=0
#            elif linedir=='ns':
#                if north0<north:
#                    offset=np.sqrt((east0-east)**2+(north0-north)**2)
#                elif north0>north:
#                    offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
#                else:
#                    offset=0
#        offsetlst.append(offset)
#        #get phase tensor elements and flip so the top is small periods/high 
#        #frequencies
#        rt=imp.getResTensor(thetar=rotz)
#        periodlst=imp.period[::-1]
#        phimax=rt.rhomax[::-1]
#        phimin=rt.rhomin[::-1]
#        azimuth=rt.rhoazimuth[::-1]
#        if colorkey=='rhomin':
#            colorarray=phimin
#        elif colorkey=='rhomax':
#            colorarray=phimax
#        elif colorkey=='rhobeta':
#            colorarray=rt.rhobeta[::-1]
#        elif colorkey=='rhodet':
#            colorarray=rt.rhodet[::-1]
#        n=len(periodlst)
#        minlst.append(min(colorarray))
#        maxlst.append(max(colorarray))
#
#        for jj in range(n):
#            #make sure the ellipses will be visable
#            eheight=phimax[jj]/phimax[jj]*esize
#            ewidth=phimin[jj]/phimax[jj]*esize
#        
#            #create an ellipse scaled by phimin and phimax and oriented along
#            #the azimuth 
#            
#            emax=10*esize
#            if eheight>emax or ewidth>emax:
#                pass
#            else:
#                ellip=Ellipse((offset*offsetscaling,3*jj),width=ewidth,
#                              height=eheight,
#                              angle=azimuth[jj])
#                #color the ellipse as red high conductivity, blue low
#                cvars=abs(np.log10(colorarray[jj]))/colorkeymm[1]
#                if cvars>1:
#                    ellip.set_facecolor((0,0,1))
#                elif cvars<float(colorkeymm[0])/colorkeymm[1] or cvars<0:
#                    ellip.set_facecolor((1,0,0))
#                else:
#                    ellip.set_facecolor((1-cvars,0,cvars))
#                ax.add_artist(ellip)
#            
#    offsetlst=np.array(offsetlst)
#    ax.set_xlim(min(offsetlst)*offsetscaling-4,max(offsetlst)*offsetscaling+4)
#    ax.set_ylim(-5,n*3+5)
#    if yscale=='period':
#        yticklabels=['{0:.3g}'.format(periodlst[ii]) for ii in np.arange(start=0,stop=n,
#                     step=2)]
#        ax.set_ylabel('Period (s)',fontsize=fs+5,fontweight='bold')
#    elif yscale=='frequency':
#        yticklabels=['{0:.4g}'.format(1./periodlst[ii]) for ii in np.arange(start=0,stop=n,
#                     step=2)]
#        ax.set_ylabel('Frequency (Hz)',fontsize=fs+5,fontweight='bold')
#    ax.set_xlabel('Station',fontsize=fs+5,fontweight='bold')
#    
#    if title==None:
#        pass
##        plt.title('Phase Tensor Pseudo Section '+title,fontsize=16)
#    else:
#        ax.set_title(title,fontsize=fs+4)
#    plt.yticks(np.arange(start=0,stop=3*n,step=6),yticklabels)
#    plt.xticks(np.array(offsetlst)*offsetscaling,stationlst)
#    
#    ax.grid(alpha=.25,which='both')
#    
#    print 'Colorkey min = ',min(minlst)
#    print 'Colorkey max = ',max(maxlst)
#    
#    #make colorbar
#    ax2=make_axes(ax,shrink=cbshrink)
#    cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=colorkeymm[0],
#                        vmax=colorkeymm[1]))
#    cb.set_label(colorkey,fontdict={'size':14,'weight':'bold'})
#    
#    plt.show()
#
#def plotPTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='phimin',xpad=.2,
#               ypad=.2,tickstrfmt='%2.2f',cborientation='vertical',
#               colorkeymm=[0,90.],figsave='y',fmt=['png'],rotz=0,pxy=[10,12],
#               galpha=.25,stationid=None,stationpad=.0005,
#               sfdict={'size':12,'weight':'bold'},indarrows='n',
#               cmap='ptcmap',tscale='period',mapscale='latlon',fignum=1,
#               imagefile=None,image_extent=None,refpoint=(0,0),cbshrink=.8,
#               arrowprop={'headheight':0.25,'headwidth':0.25,'linewidth':0.5,
#                          'arrowscale':1},
#               arrowlegend={'placement':'lower right','xborderpad':.2,
#                            'yborderpad':.2,'fontpad':.05,
#                            'fontdict':{'size':10,'weight':'bold'}}):
#    """  
#    Plots phase tensor ellipses in map view from a list of edifiles with full 
#    path.
#    
#    Arguments:
#    ----------
#        **edilst** : list of strings
#                          full paths to .edi files to plot
#                          
#        **freqspot** : position in frequency list for plotting, an integer
#        
#        **esize** : float
#                    size of ellipse in map units 
#        
#        **colorkey** : [ 'phimin' | 'beta' | 'ellipticity' | 'phidet' ]
#                       fill color of the ellipses and can be:
#                       * 'phimin'      -> minimum phase
#                       * 'beta'        -> phase tensor skew angle
#                       * 'ellipticity' -> phase tensor ellipticity
#                       * 'phidet'      -> determinant of the phase tensor
#                       * *Default* is 'phimin'
#                       
#        **colorkeymm** : tuple (min,max)
#                        min and max of colorkey to which colorbar is normalized
#                        to.  In degrees except for ellipticity which is already
#                        normalized.
#                    
#        **xpad** : float
#                   padding in x-direction in map units from furthest ellipses
#        
#        **ypad** : float
#                   padding in y-direction in map units from furthest ellipses
#        
#        **tickstrfmt** : string
#                        format of tick strings needs to be a string format.
#                        ex: '%.2f' for 2 decimal place in floating number
#        
#        **cborientation** : [ 'horizontal' | 'vertical' ]
#                            colorbar orientation horizontal or vertical
#        
#        **figsave** :  [ 'y' | 'n' ]
#                        * 'y' will save figure to edifilelist path in a folder 
#                           called PTfigures.  The figure will also close once
#                           saved.
#                        * 'n' will not save the figure
#                
#        **fmt** : list of formats 
#                 can be pdf,eps or any other formats supported by matplotlib. 
#                 Can be a list of multiple formats.
#                 **Note: that pdf and eps sometimes do not work properly.
#            
#        **rotz** : float
#                   angle in degrees to rotate the data clockwise positive.
#                   *Default* is 0
#        
#        **pxy** : tuple (width,height)
#                  dimensions of the figure in inches.  *Default* is (10,12)
#        
#        **galpha** : float [0:1]
#                     opacity of the grid.  0 for transparent, 1 for opaque
#                     *Default* is 0.25
#        
#        **stationid** : tuple or list 
#                        start and stop of station name indicies.  
#                        ex: for MT01dr stationid=(0,4) will be MT01.
#                        *Default* is None
#                    
#        **stationpad** : float
#                         padding for station name in the y-direction in map
#                         units
#        
#        **sfdict** : dictionary
#                     font dictionary for station name. Keys can be
#                     matplotlib.text properties, common ones are:
#                         * 'size'   -> for font size
#                         * 'weight' -> for font weight
#                         * 'color'  -> for color of font
#                                     
#        **indarrow** : [ 'yri' | 'yr' | 'yi' | 'n' ]
#                        * 'yri' to plot induction both real and imaginary 
#                           induction arrows 
#                        * 'yr' to plot just the real induction arrows
#                        * 'yi' to plot the imaginary induction arrows
#                        * 'n' to not plot them
#                        * *Default* is 'n'                        
#                        **Note: convention is to point towards a conductor **
#       
#       **arrowprop** : dictionary of arrow properties with keys:
#                        * 'linewidth'  -> width of the arrow line
#                        * 'headheight' -> height of arrow head
#                        * 'headwidth'  -> width of the arrow head
#                        * 'arrowscale' -> scale size of the arrow
#        
#        
#        **arrowlegend** : dictionary of properties for legend with keys:
#                        * 'placement -> placement of arrow legend can be:
#                            - 'upper right'
#                            - 'lower right'
#                            - 'upper left'
#                            - 'lower left'
#                        * 'xborderpad' -> padding from x axis
#                        * yborderpad'  -> padding from y axis
#                        * fontpad'     -> padding between arrow and legend text
#                        * fontdict'    -> dictionary of font properties
#        
#        **cmap** :[ 'ptcmap' | 'ptcmap3' | 'skcmap' | 'skcmap2' | 'rtcmap' ]
#                  color map of ellipse facecolor.  So far the colormaps are:
#                      * 'ptcmap'  -> yellow (low phase) to red (high phase)
#                      * 'ptcmap3' -> white (low numbers) to blue (high numbers)
#                      * 'skcmap'  -> blue to yellow to red
#                      * 'skcmap2' -> blue to white to red
#                      * 'rtcmap'  -> blue to purple to red
#            
#        **tscale** : [ 'period' | 'frequency' ]
#                     * 'period'    -> plot vertical scale in period
#                     * 'frequency' -> plot vertical scale in frequency
#        
#        **mapscale** : [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
#                      scale of map
#                          * 'latlon'     -> lats and longs
#                          * 'eastnorth'  -> for easting and northing, this is 
#                             recomended if you want to plot tipper data for 
#                             small surveys.
#                          * 'eastnorthkm' -> for kilometer scaling
#                          ** Note: if the ellipses are plotting in the wrong 
#                          spot in east-north scale, check to see if your survey
#                          crosses UTM grids.  Currently there is no way of 
#                          dealing with this.
#                   
#        **imagefile** : string
#                        path to an image file jpg or png or svg
#        
#        **image_extent** : tuple (xmin,xmax,ymin,ymax)
#                           coordinates according to mapscale. Must be input if
#                           image file is not None.
#        
#        **refpoint** : tuple (x,y)
#                       reference point estimate relative distance to.  This 
#                       point will be (0,0) on the map and everything else is 
#                       referenced to this point
#         
#    :Example: ::
#        
#        >>> import mtpy.imaging.mtplottools as mtplot
#        >>> import os
#        >>> edipath = r"/home/EDIfiles"
#        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
#        >>> ...       if edi.find('.edi')]
#        >>> # color by phimin with a range of 20-70 deg in meter scale
#        >>> mtplot.plotPTMaps(edilst,colorkeymm=(20,70),mapscale='eastnorth')
#        >>> 
#        >>> # add induction arrows to skew angle with range (-5,5)
#        >>> mtplot.plotPTMaps(edilst,colorkey='beta',
#        >>> ...               colorkeymm=(-5,5),indarrows='yri')
#        >>>
#        >>> # put an underlying image as a basemap in km scale
#        >>> mtplot.plotPTMaps(edilst,imagefile=r"/home/Maps/Basemap.jpg",
#        >>> ...               image_extent=(0,20,0,20),mapscale='eastnorthkm')
#
#
#    """
#    jj=freqspot
#    fig=plt.figure(fignum,pxy,dpi=200)
#    plt.clf()
#    ax=fig.add_subplot(1,1,1,aspect='equal')
#    
#    if imagefile!=None:
#        if image_extent==None:
#            raise ValueError('Need to put in image extent')
#        im=plt.imread(imagefile)
#        ax.imshow(im,origin='lower',extent=image_extent,aspect='auto')
#        
#    elliplst=[]
#    latlst=[]
#    lonlst=[]
#    if not esize is float:
#        esize=float(esize)
#    
#    if mapscale=='latlon':
#        tickstrfmt='%.3f'
#    elif mapscale=='eastnorth' or mapscale=='eastnorthkm':
#        tickstrfmt='%.0f'
#    
#    ckmin=colorkeymm[0]
#    ckmax=colorkeymm[1]
#    
#    plt.rcParams['font.size']=8
#    plt.rcParams['figure.subplot.left']=.1
#    plt.rcParams['figure.subplot.right']=.98
#    plt.rcParams['figure.subplot.bottom']=.1
#    plt.rcParams['figure.subplot.top']=.93
#    plt.rcParams['figure.subplot.wspace']=.55
#    plt.rcParams['figure.subplot.hspace']=.70
#    #plt.rcParams['font.family']='helvetica'
#    
#    for ii,filename in enumerate(edifilelst):
#        #get phase tensor info
#        imp=Z.Z(filename)
#        #get phase tensor
#        pt=imp.getPhaseTensor(thetar=rotz)
#        pt.phimax=np.nan_to_num(pt.phimax)
#        pt.phimin=np.nan_to_num(pt.phimin)
#        #check to see if the period is there
#        try:
#            freq=1./imp.period[jj]
#            if mapscale=='latlon':
#                latlst.append(imp.lat)
#                lonlst.append(imp.lon)
#                plotx=imp.lon-refpoint[0]
#                ploty=imp.lat-refpoint[1]
#            elif mapscale=='eastnorth':
#                zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
#                if ii==0:
#                    zone1=zone
#                    plotx=east-refpoint[0]
#                    ploty=north-refpoint[1]
#                else:
#                    if zone1!=zone:
#                        if zone1[0:2]==zone[0:2]:
#                            pass
#                        elif int(zone1[0:2])<int(zone[0:2]):
#                            east=east+500000
#                        else:
#                            east=east-500000
#                        latlst.append(north-refpoint[1])
#                        lonlst.append(east-refpoint[0])
#                        plotx=east-refpoint[0]
#                        ploty=north-refpoint[1]
#                    else:
#                        latlst.append(north-refpoint[1])
#                        lonlst.append(east-refpoint[0])
#                        plotx=east-refpoint[0]
#                        ploty=north-refpoint[1]
#            elif mapscale=='eastnorthkm':
#                zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
#                if ii==0:
#                    zone1=zone
#                    plotx=(east-refpoint[0])/1000.
#                    ploty=(north-refpoint[1])/1000.
#                else:
#                    if zone1!=zone:
#                        if zone1[0:2]==zone[0:2]:
#                            pass
#                        elif int(zone1[0:2])<int(zone[0:2]):
#                            east=east+500000
#                        else:
#                            east=east-500000
#                        latlst.append((north-refpoint[1])/1000.)
#                        lonlst.append((east-refpoint[0])/1000.)
#                        plotx=(east-refpoint[0])/1000.
#                        ploty=(north-refpoint[1])/1000.
#                    else:
#                        latlst.append((north-refpoint[1])/1000.)
#                        lonlst.append((east-refpoint[0])/1000.)
#                        plotx=(east-refpoint[0])/1000.
#                        ploty=(north-refpoint[1])/1000.
#            else:
#                raise NameError('mapscale not recognized')
#                    
#            
#            phimin=pt.phimin[jj]
#            phimax=pt.phimax[jj]
#            eangle=pt.azimuth[jj]
#            #create an ellipse object
#            if phimax==0 or phimax>100 or phimin==0 or pt.phiminang[jj]<0 or \
#                pt.phiminang[jj]>100:
#                eheight=.0000001*esize
#                ewidth=.0000001*esize
#            else:
#                scaling=esize/phimax
#                eheight=phimin*scaling
#                ewidth=phimax*scaling
#            ellipd=Ellipse((plotx,ploty),width=ewidth,height=eheight,
#                          angle=eangle)
#            #print imp.lon,imp.lat,scaling,ewidth,eheight,phimin,phimax
#            elliplst.append(ellipd)
#            ax.add_artist(ellipd)
#            
#            #get face color info
#            if colorkey=='phiminang' or  colorkey=='phimin':
#                cvar=(pt.phiminang[jj]-ckmin)/(ckmax-ckmin)
#            elif colorkey=='phidet':
#                cvar=(pt.phidet[jj]-ckmin)/(ckmax-ckmin)
#            elif colorkey=='beta':
#                cvar=(pt.beta[jj]-abs(ckmin))/(ckmax-ckmin)
#            elif colorkey=='ellipticity':
#                cvar=(pt.ellipticity[jj]-ckmin)/(ckmax-ckmin)
#            else:
#                raise NameError('color key '+colorkey+' not supported')
#            
#            #set facecolor depending on the colormap
#            #yellow to red
#            if cmap=='ptcmap':
#                if abs(cvar)>1:
#                    ellipd.set_facecolor((1,0,0))
#                elif cvar<0:
#                    ellipd.set_facecolor((1-abs(cvar),1,abs(cvar)))
#                else:
#                    ellipd.set_facecolor((1,1-abs(cvar),.1))
#            #white to blue
#            elif cmap=='ptcmap3':
#                if abs(cvar)>1:
#                    ellipd.set_facecolor((0,0,0))
#                else:
#                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
#            #blue to yellow to red
#            elif cmap=='skcmap2':
#                if cvar<0 and cvar>-1:
#                    ellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
#                elif cvar<-1:
#                    ellipd.set_facecolor((0,0,1))
#                elif cvar>0 and cvar<1:
#                    ellipd.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
#                elif cvar>1:
#                    ellipd.set_facecolor((1,0,0))
#            #blue to white to red
#            elif cmap=='skcmap':
#                if cvar<0 and cvar>-1:
#                    ellipd.set_facecolor((abs(cvar),abs(cvar),1-abs(cvar)))
#                elif cvar<-1:
#                    ellipd.set_facecolor((0,0,1))
#                elif cvar>0 and cvar<1:
#                    ellipd.set_facecolor((1,1-abs(cvar),.01))
#                elif cvar>1:
#                    ellipd.set_facecolor((1,0,0))
#                    
#            #-----------Plot Induction Arrows---------------------------
#            if indarrows.find('y')==0:
##                if mapscale=='latlon':
##                    print 'Might try mapscale=latlon for better scale of arrows'
#                    
#                tip=imp.getTipper(thetar=rotz)
#                aheight=arrowprop['headheight']
#                awidth=arrowprop['headwidth']
#                ascale=arrowprop['arrowscale']
#                #plot real tipper
#                if indarrows=='yri' or indarrows=='yr':
#                    if tip.magreal[jj]<=1.0:
#                        txr=tip.magreal[jj]*ascale*\
#                            np.sin((tip.anglereal[jj])*np.pi/180)
#                        tyr=tip.magreal[jj]*ascale*\
#                            np.cos((tip.anglereal[jj])*np.pi/180)
#    
#                        ax.arrow(plotx,ploty,txr,tyr,lw=arrowprop['linewidth'],
#                             facecolor='k',edgecolor='k',
#                             length_includes_head=False,
#                             head_width=awidth,head_length=aheight)
#                    else:
#                        pass
#                #plot imaginary tipper
#                if indarrows=='yri' or indarrows=='yi':
#                    if tip.magimag[jj]<=1.0:
#                        txi=tip.magimag[jj]*\
#                            np.sin((tip.angleimag[jj])*np.pi/180)*scaling
#                        tyi=tip.magimag[jj]*\
#                            np.cos((tip.angleimag[jj])*np.pi/180)*scaling
#    
#                        ax.arrow(plotx,ploty,txi,tyi,lw=arrowprop['linewidth'],
#                                 facecolor='b',edgecolor='b',
#                                 length_includes_head=False,
#                                 head_width=awidth,head_length=aheight)
#                         
#            
#            #------------Plot station name------------------------------
#            if stationid!=None:
#                ax.text(plotx,ploty+stationpad,
#                        imp.station[stationid[0]:stationid[1]],
#                        horizontalalignment='center',
#                        verticalalignment='baseline',
#                        fontdict=sfdict)
#                
#        #if the period is not there 
#        except IndexError:
#            print 'Did not find index for station'.format(jj)+imp.station
#    
#    if mapscale=='latlon':    
#        ax.set_xlabel('longitude',fontsize=10,fontweight='bold')
#        ax.set_ylabel('latitude',fontsize=10,fontweight='bold')
#        ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
#        ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#        ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
#        ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#    elif mapscale=='eastnorth':
#        ax.set_xlabel('Easting (m)',fontsize=10,fontweight='bold')
#        ax.set_ylabel('Northing (m)',fontsize=10,fontweight='bold')
#        ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
#        ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#        ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
#        ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#    elif mapscale=='eastnorthkm':
#        ax.set_xlabel('Easting (km)',fontsize=10,fontweight='bold')
#        ax.set_ylabel('Northing (km)',fontsize=10,fontweight='bold')
#        ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
#        ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#        ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
#        ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#    if tscale=='period':
#        titlefreq='{0:.5g} (s)'.format(1./freq)
#    else:
#        titlefreq='{0:.5g} (Hz)'.format(freq)
#    ax.set_title('Phase Tensor Map for '+titlefreq,
#                 fontsize=10,fontweight='bold')
#    #plot induction arrow scale bar
#    if indarrows.find('y')==0:
#        parrx=ax.get_xlim()
#        parry=ax.get_ylim()
#        try:
#            axpad=arrowlegend['xborderpad']
#        except KeyError:
#            axpad=xpad+arrowprop['arrowscale']
#        try:
#            aypad=arrowlegend['yborderpad']
#        except KeyError:
#            aypad=ypad
#        try:
#            txtpad=arrowlegend['fontpad']
#        except KeyError:
#            txtpad=.25*esize
#            
#        
#        if arrowlegend['placement']=='lower right':
#            pax=parrx[1]-axpad
#            pay=parry[0]+aypad
#            ptx=arrowprop['arrowscale']
#            pty=0
#            txa=parrx[1]-axpad+arrowprop['arrowscale']/2.
#            txy=pay+txtpad
#        elif arrowlegend['placement']=='upper right':
#            pax=parrx[1]-axpad
#            pay=parry[1]-aypad
#            ptx=arrowprop['arrowscale']
#            pty=0
#            txa=parrx[1]-axpad+arrowprop['arrowscale']/2.
#            txy=pay+txtpad
#        elif arrowlegend['placement']=='lower left':
#            pax=parrx[0]+axpad
#            pay=parry[0]+aypad
#            ptx=arrowprop['arrowscale']
#            pty=0
#            txa=parrx[0]+axpad+arrowprop['arrowscale']/2.
#            txy=pay+txtpad
#        elif arrowlegend['placement']=='upper left':
#            pax=parrx[0]+axpad
#            pay=parry[1]-aypad
#            ptx=arrowprop['arrowscale']
#            pty=0
#            txa=parrx[0]+axpad+arrowprop['arrowscale']/2.
#            txy=pay+txtpad
#        else:
#            raise NameError('arrowlegend not supported.')
#            
#        ax.arrow(pax,pay,ptx,pty,lw=arrowprop['linewidth'],
#                             facecolor='k',edgecolor='k',
#                             length_includes_head=False,
#                             head_width=arrowprop['headwidth'],
#                             head_length=arrowprop['headheight'])
#        
#        ax.text(txa,txy,'|T|=1',
#                horizontalalignment='center',
#                verticalalignment='baseline',
#                fontdict={'size':10,'weight':'bold'})
#                
#    ax.grid(alpha=galpha)
#    
##    if indarrows.find('y')==0:
##        if indarrows=='yri':
##            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
##            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
##            ax.legend([treal[0],timag[0]],['Tipper_real','Tipper_imag'],
##                      loc='upper center',
##                      prop={'size':10,'weight':'bold'},
##                      ncol=2,markerscale=.5,borderaxespad=.005,
##                      borderpad=.25)
##        elif indarrows=='yr':
##            treal=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'k')
##            ax.legend([treal[0]],['Tipper_real'],
##                      loc='upper center',
##                      prop={'size':10,'weight':'bold'},
##                      ncol=2,markerscale=.5,borderaxespad=.005,
##                      borderpad=.25)
##        elif indarrows=='yi':
##            timag=ax.plot(np.arange(10)*.000005,np.arange(10)*.00005,'b')
##            ax.legend([timag[0]],['Tipper_imag'],
##                      loc='upper center',
##                      prop={'size':10,'weight':'bold'},
##                      ncol=2,markerscale=.5,borderaxespad=.005,
##                      borderpad=.25)
#    
#    
#    
#    #make a colorbar with appropriate colors             
#    ax2=make_axes(ax,shrink=cbshrink)
#    if cmap=='ptcmap': 
#        cb=ColorbarBase(ax2[0],cmap=ptcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='ptcmap3': 
#        cb=ColorbarBase(ax2[0],cmap=ptcmap3,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='skcmap': 
#        cb=ColorbarBase(ax2[0],cmap=skcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='skcmap2': 
#        cb=ColorbarBase(ax2[0],cmap=skcmap2,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    elif cmap=='rtcmap': 
#        cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#
#    #label the color bar accordingly
#    if colorkey=='phimin' or colorkey=='phiminang':
#        cb.set_label('$\Phi_{min}$ (deg)',fontdict={'size':10,'weight':'bold'})
#    elif colorkey=='beta':
#        cb.set_label('Skew (deg)',fontdict={'size':10,'weight':'bold'})
#    elif colorkey=='phidet':
#        cb.set_label('Det{$\Phi$} (deg)',fontdict={'size':10,'weight':'bold'})
#    elif colorkey=='ellipticity':
#        cb.set_label('Ellipticity',fontdict={'size':10,'weight':'bold'})
#        
#    plt.show()
#    
#    #save the figure if desired
#    if figsave=='y':
#        sf='_{0:.5g}'.format(freq)
#        savepath=os.path.join(os.path.dirname(edifilelst[0]),'PTfigures')
#        if not os.path.exists(savepath):
#            os.mkdir(savepath)
#            print 'Made directory: '+savepath
#        else:
#            pass
#        for f in fmt:
#            fig.savefig(os.path.join(savepath,
#                                 'PTmap_'+colorkey+sf+'Hz.'+f),
#                                 format=f)
#            print 'Saved file figures to: '+ os.path.join(savepath, \
#                                     'PTmap_'+colorkey+sf+'Hz.'+f)
#        plt.close()
#        
#
#def plotResPhasePseudoSection(edifilelst,stationid=[0,4],ffactor=1,
#                              maxperiod=60,aspect=4,cmap='jet_r',xtickspace=2,
#                              linedir='ns',rotz=0,dpi=300):
#    """
#    plotResPhasePseudoSection(edifilelst,stationid=4,ffactor=10E3,df=100.,
#                              maxperiod=24,aspect=4,cmap='jet_r') plots a 
#    pseudo section from a list of edifiles with full path with descending 
#    frequency or ascending period on the y axis and relative position on the x. 
#    keyword arguments can be:
#        stationid = start and finish of station string for plotting
#        ffactor = fudge factor to make sure resistivities plot correctly
#        df = sampling frequency (Hz)
#        maxperiod = maximum period to plot
#        aspect = horizontal to vertical ratio of plot
#        cmap = colormap of image
#    """
#    
#    fs=int(dpi/40)
#    plt.rcParams['font.size']=fs
#    plt.rcParams['figure.subplot.left']=.07
#    plt.rcParams['figure.subplot.right']=.98
#    plt.rcParams['figure.subplot.bottom']=.06
#    plt.rcParams['figure.subplot.top']=.94
#    plt.rcParams['figure.subplot.wspace']=.01
#    plt.rcParams['figure.subplot.hspace']=.20
#    #plt.rcParams['font.family']='helvetica'
#    
#    
#    #create a plot instance
##    ax1=fig.add_subplot(2,2,1)
##    ax2=fig.add_subplot(2,2,2,sharex=ax1)
##    ax3=fig.add_subplot(2,2,3,sharex=ax1)
##    ax4=fig.add_subplot(2,2,4,sharex=ax1)
#    
#    #create empty lists to put things into
#    stationlst=[]
#    offsetlst=[]
#    minperiodlst=[]
#    maxperiodlst=[]
#    periodlst=[]
#    nlst=[]
#    
#    #create empty arrays to put data into
#    n=len(edifilelst)
#    resxx=np.zeros((maxperiod,n))
#    resxy=np.zeros((maxperiod,n))
#    resyx=np.zeros((maxperiod,n))
#    resyy=np.zeros((maxperiod,n))
#    
#    phasexx=np.zeros((maxperiod,n))
#    phasexy=np.zeros((maxperiod,n))
#    phaseyx=np.zeros((maxperiod,n))
#    phaseyy=np.zeros((maxperiod,n))
#    
#    #get data into arrays
#    for ii,fn in enumerate(edifilelst):
#        #get offsets between stations
#        imp=Z.Z(fn)
#        stationlst.append(imp.station[stationid[0]:stationid[1]])
#        zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
#        
#        if ii==0:
#            east0=east
#            north0=north
#            offset=0.0
#        else:
#            if linedir=='ew': 
#                if east0<east:
#                    offset=np.sqrt((east0-east)**2+(north0-north)**2)
#                elif east0>east:
#                    offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
#                else:
#                    offset=0
#            elif linedir=='ns':
#                if north0<north:
#                    offset=np.sqrt((east0-east)**2+(north0-north)**2)
#                elif north0>north:
#                    offset=-1*np.sqrt((east0-east)**2+(north0-north)**2)
#                else:
#                    offset=0
#        offsetlst.append(offset)
#        #get resisitivity and phase in a dictionary and append to a list
#        rp=imp.getResPhase(ffactor=ffactor,thetar=rotz)
#        m=len(imp.period)
#        resxx[0:m,ii]=np.log10(rp.resxx)#[::-1]
#        resxy[0:m,ii]=np.log10(rp.resxy)#[::-1]
#        resyx[0:m,ii]=np.log10(rp.resyx)#[::-1]
#        resyy[0:m,ii]=np.log10(rp.resyy)#[::-1]
#        
#        phasexx[0:m,ii]=rp.phasexx#[::-1]
#        phasexy[0:m,ii]=rp.phasexy#[::-1]
#        phaseyx[0:m,ii]=rp.phaseyx+180#[::-1]
#        phaseyy[0:m,ii]=rp.phaseyy#[::-1]
#        
#        #get min and max of the period, just in case some edi files have 
#        #different period information
#        minperiodlst.append(min(imp.period))
#        maxperiodlst.append(max(imp.period))
#        nlst.append(len(imp.period))
#        periodlst.append(imp.period)
#        
#    nperiod=max(nlst)
#    for pp in range(len(periodlst)):
#        if nlst[pp]==nperiod:
#            period=np.log10(periodlst[pp])#[::-1]
#
#    #plot data
#    fig=plt.figure(1,[12,4],dpi=150)
#    extent=(0,n,period[-1],period[0])
#    aspect=aspect
#    cmap=cmap
#    pad=.2
#    kwargs={'aspect':aspect,'cmap':cmap,'extent':extent}
#    cbarkwargs={'shrink':.7,'orientation':'horizontal','pad':pad}
#
#    
#    ax2=plt.subplot(2,4,2)
#    plt.imshow(resxy[0:nperiod,:],**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
##    plt.ylabel('Log$_{10}$ Period',fontsize=10,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('Log$_{10}\mathbf{(1/\sigma_{xy})}$',fontsize=fs+4,fontweight='bold')
#    ax2.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
#    ax3=plt.subplot(2,4,3)
#    plt.imshow(resyx[0:nperiod,:],**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
##    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('Log$_{10}\mathbf{(1/\sigma_{yx})}$',fontsize=fs+4,fontweight='bold')
#    ax3.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
#    ax6=plt.subplot(2,4,6)
#    plt.imshow(phasexy[0:nperiod,:],vmin=0,vmax=90,**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
##    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('$\mathbf{\phi_{xy}}$',fontsize=fs+4)
#    ax6.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
#    ax7=plt.subplot(2,4,7)
#    plt.imshow(phaseyx[0:nperiod,:],vmin=0,vmax=90,**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
##    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('$\mathbf{\phi_{yx}}$',fontsize=fs+4)
#    ax7.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
##    fig2=plt.figure(2,dpi=150)    
#
#    ax1=plt.subplot(2,4,1)
#    plt.imshow(resxx[0:nperiod,:],**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
#    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('Log$_{10}\mathbf{(1/\sigma_{xx})}$',fontsize=fs+4,fontweight='bold')
#    ax1.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
#    ax4=plt.subplot(2,4,4)
#    plt.imshow(resyy[0:nperiod,:],**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
##    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('Log$_{10}\mathbf{(1/\sigma_{yy})}$',fontsize=fs+4,fontweight='bold')
#    ax4.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
#    ax5=plt.subplot(2,4,5)
#    plt.imshow(phasexx[0:nperiod,:],**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
#    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('$\mathbf{\phi_{xx}}$',fontsize=fs+4)
#    ax5.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    
#    ax8=plt.subplot(2,4,8)
#    plt.imshow(phaseyy[0:nperiod,:],**kwargs)
#    plt.colorbar(**cbarkwargs)
#    plt.xticks(np.arange(0,n,xtickspace),
#               [stationlst[st] for st in range(0,n,xtickspace)])
##    plt.ylabel('Log$_{10}$ Period',fontsize=fs+4,fontweight='bold')
#    plt.xlabel('Station',fontsize=fs+4,fontweight='bold')
#    plt.title('$\mathbf{\phi_{yy}}$',fontsize=fs+4)
#    ax8.yaxis.set_minor_locator(MultipleLocator(.2))
#    plt.grid()
#    plt.show()
#    
#def plotRTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='rhodet',xpad=.2,
#               ypad=.2,tickstrfmt='%2.4f',cborientation='vertical',
#               colorkeymm=[0,4],figsave='y',fmt=['png'],rotz=0,pxy=[10,12],
#                galpha=.25):
#    """ 
#    plotPTMaps(edifilelst,freqspot=10,esize=2.0,colorkey='phimin',xpad=.2,
#               ypad=.2,tickstrfmt='%2.4f',cborientation='vertical',
#               colorkeymax=90.,figsave='y',fmt=['png']) 
#    plots phase tensor ellipses in map view from a list of edifiles with full 
#    path.  Parameters are:
#        
#        freqspot = position in frequency list for plotting, an integer
#        esize = size of ellipse, float
#        colorkey =  the fill color of the ellipses and can be any of the 
#        dictionary keys returned by Z.getPhaseTensor(), note skew is beta:
#        'phimin','phi', 'phiminvar', 'azimuthvar', 'azimuth', 'betavar', 
#        'phivar', 'alphavar', 'beta', 'ellipticityvar', 'phiminangvar', 
#        'ellipticity', 'phimaxangvar', 'alpha', 'phiminang', 'phimaxvar', 
#        'phimaxang', 'phimax'
#        colorkeymm = [min,max] min and max of colorkey to which colorbar is
#                    normalized to.
#        xpad = pad from xmin and xmax on plot
#        ypad = pad from ymin and ymax on plot
#        tickstrfmt = format of tick strings needs to be a string format
#        cborientation = colorbar orientation horizontal or vertical
#        figsave = y or n, if yes figure will be saved to edifilelist path in 
#                a folder called PTfigures.
#        fmt = ['png'] is format of save figure can be pdf,eps or any other 
#            formats supported by matplotlib. Can be a list of multiple formats.
#            Note that pdf and eps do not properly yet.
#        
#    """
#    jj=freqspot
#    fig=plt.figure(1,pxy,dpi=150)
#    ax=fig.add_subplot(1,1,1,aspect='equal')
#    elliplst=[]
#    latlst=[]
#    lonlst=[]
#    if not esize is float:
#        esize=float(esize)
#    
#    ckmin=colorkeymm[0]
#    ckmax=colorkeymm[1]
#    
#    plt.rcParams['font.size']=8
#    plt.rcParams['figure.subplot.left']=.1
#    plt.rcParams['figure.subplot.right']=.98
#    plt.rcParams['figure.subplot.bottom']=.1
#    plt.rcParams['figure.subplot.top']=.93
#    plt.rcParams['figure.subplot.wspace']=.55
#    plt.rcParams['figure.subplot.hspace']=.70
#    #plt.rcParams['font.family']='helvetica'
#    
#    for ii,filename in enumerate(edifilelst):
#        #get phase tensor info
#        imp=Z.Z(filename)
#        try:
#            freq=1./imp.period[jj]
#            latlst.append(imp.lat)
#            lonlst.append(imp.lon)
#            rt=imp.getResTensor(thetar=rotz)
#            phimax=rt.rhomax[jj]
#            phimin=rt.rhomin[jj]
#            eangle=rt.rhoazimuth[jj]
#            if colorkey=='rhomin':
#                colorarray=phimin
#            elif colorkey=='rhomax':
#                colorarray=phimax
#            elif colorkey=='rhobeta':
#                colorarray=rt.rhobeta[jj]
#            elif colorkey=='rhodet':
#                colorarray=rt.rhodet[jj]
#            eangle=rt.rhoazimuth[jj]
#            #create an ellipse object
#            scaling=esize/phimax
#            eheight=phimin*scaling
#            ewidth=phimax*scaling
#            
##            ellipd=Ellipse((imp.lon,imp.lat),width=ewidth,height=eheight,
##                          angle=eangle)
##            #print imp.lon,imp.lat,scaling,ewidth,eheight,phimin,phimax
##            elliplst.append(ellipd)
##            ax.add_artist(ellipd)
##            #get face color info
#        
#            #create an ellipse scaled by phimin and phimax and oriented along
#            #the azimuth 
#            emax=10*esize
##            print eheight,ewidth,emax
#            if eheight>emax or ewidth>emax:
#                pass
#            else:
#                ellip=Ellipse((imp.lon,imp.lat),width=ewidth,
#                              height=eheight,
#                              angle=eangle)
#                #color the ellipse as red high conductivity, blue low
#                cvars=(np.log10(abs(colorarray))-ckmin)/(ckmax-ckmin)
#                #print cvars
#                if cvars>1:
#                    ellip.set_facecolor((0,0,1))
#                elif cvars<float(ckmin)/ckmax or cvars<0:
#                    ellip.set_facecolor((1,0,0))
#                else:
#                    ellip.set_facecolor((1-cvars,0,cvars))
#                ax.add_artist(ellip)            
#            
#        except IndexError:
#            pass
#        
#    ax.set_xlabel('longitude',fontsize=10,fontweight='bold')
#    ax.set_ylabel('latitude',fontsize=10,fontweight='bold')
#    ax.set_xlim(min(lonlst)-xpad,max(lonlst)+xpad)
#    ax.xaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#    ax.set_ylim(min(latlst)-xpad,max(latlst)+xpad)
#    ax.yaxis.set_major_formatter(FormatStrFormatter(tickstrfmt))
#    titlefreq='%2.3f' % (1/freq)
#    ax.set_title('Resistivity Tensor for '+titlefreq+'(s)',
#                 fontsize=10,fontweight='bold')
#    ax.grid(alpha=galpha)
#                 
#    ax2=make_axes(ax,shrink=.8)
#    if colorkey=='rhodet' or colorkey=='phidet':
#        cb=ColorbarBase(ax2[0],cmap=rtcmap,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label('Log$_{10}$ '+colorkey+' ($\Omega \cdot$m)')
#    elif colorkey=='beta' or colorkey=='ellipticity':
#        cb=ColorbarBase(ax2[0],cmap=ptcmap3,norm=Normalize(vmin=ckmin,vmax=ckmax))
#        cb.set_label(colorkey+' (deg)')
#    plt.show()
#    
#    if figsave=='y':
#        sf='%2.3g' % freq
#        savepath=os.path.join(os.path.dirname(edifilelst[0]),'RTfigures')
#        if not os.path.exists(savepath):
#            os.mkdir(savepath)
#            print 'Made directory: '+savepath
#        else:
#            pass
#        for f in fmt:
#            fig.savefig(os.path.join(savepath,
#                                 'RTmap'+sf+'Hz.'+f),
#                                 format=f)
#            print 'Saved file figures to: '+ os.path.join(savepath, \
#                                     'RTmap'+sf+'Hz.'+f)
#        plt.close()
#    
#def plotRoseStrikeAngles(edilst,fignum=1,fs=10,dpi=300,thetar=0,ptol=.05,
#                         tpad=.60,galpha=.25,prange='data',plottype=1,
#                         tipper='n',pterr=None):
#    """
#    plots the strike angle as determined by phase tensor azimuth (Caldwell et 
#    al. [2004]) and invariants of the impedance tensor (Weaver et al. [2003]).
#    
#    The data is split into decades where the histogram for each is plotted in 
#    the form of a rose diagram with a range of 0 to 180 degrees.
#    Where 0 is North and 90 is East.   The median angle of the period band is 
#    set in polar diagram.  The top row is the strike estimated from
#    the invariants of the impedance tensor.  The bottom row is the azimuth
#    estimated from the phase tensor.  If tipper is 'y' then the 3rd row is the
#    strike determined from the tipper, which is orthogonal to the induction
#    arrow direction.  
#    
#    Arguments:
#    ----------
#        **edilst** : list of full paths to edifiles to be used
#        
#        **fignum** : int
#                     figure number to be plotted. *Default* is 1
#        
#        **fs** : float
#                 font size for labels of plotting. *Default* is 10
#                 
#        **dpi** : int
#                  dots-per-inch resolution of figure, 300 is needed for 
#                  publications. *Default* is 300
#                  
#        **thetar** : float
#                     angle of rotation clockwise positive. *Default* is 0
#                     
#        **ptol** : float
#                   Tolerance level to match periods from different edi files.
#                   *Default* is 0.05
#                   
#        **tpad** : float
#                   padding of the angle label at the bottom of each polar 
#                   diagram.  *Default* is 1.65
#                   
#        **galpha** : float (0:1)
#                     transparency of the font boxes and grid lines.  0 is fully
#                     tranparent and 1 is opaque.  *Default* is 0.25
#                     
#        **prange** : [ 'data' | (period_min,period_max) ]
#                    period range to estimate the strike angle. Options are:
#                        * *'data'* for estimating the strike for all periods
#                            in the data.
#                        * (pmin,pmax) for period min and period max, input as
#                          (log10(pmin),log10(pmax))
#        
#        **plottype** : [ 1 | 2 ]
#                        * *1* to plot individual decades in one plot
#                        * *2* to plot all period ranges into one polar diagram
#                              for each strike angle estimation
#        
#        **tipper** : [ 'y' | 'n' ]
#                      * *'y'* to plot the tipper strike
#                      * *'n'* to not plot tipper strike
#                      
#        **pterr** : float   
#                    Maximum error in degrees that is allowed to estimate strike.
#                    *Default* is None allowing all estimates to be used.
#                
#
#    :Example: ::
#        
#        >>> import os
#        >>> import mtpy.imaging.mtplottools as mtplot
#        >>> edipath = r"/home/EDIFiles"
#        >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
#        >>> ...       if edi.find('.edi')>0]
#        >>> # plot rose plots in decades with tipper and an error floor on pt
#        >>> mtplot.plotRoseStrikeAngles(edilst,plottype=1,pterr=5)
#        >>> # plot all decades into one rose plot for each estimation
#        >>> mtplot.plotRoseStrikeAngles(edilst,plottype=2,pterr=5)
#    """
#    plt.rcParams['font.size']=fs-2
#    plt.rcParams['figure.subplot.left']=.07
#    plt.rcParams['figure.subplot.right']=.98
#    plt.rcParams['figure.subplot.bottom']=.09
#    plt.rcParams['figure.subplot.top']=.90
#    plt.rcParams['figure.subplot.wspace']=.2
#    plt.rcParams['figure.subplot.hspace']=.4   
#    
#    invlst=[]
#    ptlst=[]
#    
#    if tipper=='y':
#        tiprlst=[]
#    
#    nc=len(edilst)
#    nt=0
#    kk=0
#    
#    for dd,edi in enumerate(edilst):
#        z1=Z.Z(edi)
#        mm=np.remainder(dd,4)
#        period=z1.period
#    
#        #get maximum length of periods and which index that is
#        if len(period)>nt:
#            nt=len(period)
#            kk=dd
#        #-----------get strike angle from invariants--------------------------
#        zinv=z1.getInvariants(thetar=thetar)
#        
#        #add 90 degrees because invariants assume 0 is north, but plotting assumes
#        #that 90 is north and measures clockwise, thus the negative
#        zs=90-zinv.strike
#        
#        #for plotting put the NW angles into the SE quadrant 
#        zs[np.where(zs>90)]=zs[np.where(zs>90)]-180
#        zs[np.where(zs<-90)]=zs[np.where(zs<-90)]+180
#        
#        #make a dictionary of strikes with keys as period
#        mdictinv=dict([(ff,jj) for ff,jj in zip(z1.period,zs)])
#        invlst.append(mdictinv)
#    
#        #------------get strike from phase tensor strike angle---------------
#        pt=z1.getPhaseTensor(thetar=thetar)
#        az=pt.azimuth
#        azerr=pt.azimuthvar
#        #put an error max on the estimation of strike angle
#        if pterr:
#            az[np.where(azerr>pterr)]=0.0
#        #don't need to add 90 because pt assumes 90 is north.
#        az[np.where(az>90)]=az[np.where(az>90)]-180
#        az[np.where(az<-90)]=az[np.where(az<-90)]+180
#        #make a dictionary of strikes with keys as period
#        mdictpt=dict([(ff,jj) for ff,jj in zip(z1.period,az)])
#        ptlst.append(mdictpt)
#        
#        #-----------get tipper strike------------------------------------
#        if tipper=='y':
#            tip=z1.getTipper(thetar=thetar)
#            tipr=-tip.anglereal
#            
#            tipr[np.where(tipr>90)]=tipr[np.where(tipr>90)]-180
#            tipr[np.where(tipr<-90)]=tipr[np.where(tipr<-90)]+180
#            
#            tiprdict=dict([(ff,jj) for ff,jj in zip(z1.period,tipr)])
#            tiprlst.append(tiprdict)
#            
#        
#        
#    
#    #get min and max period
#    maxper=np.max([np.max(mm.keys()) for mm in invlst])
#    minper=np.min([np.min(mm.keys()) for mm in ptlst])
#    
#    #make empty arrays to put data into for easy manipulation
#    medinv=np.zeros((nt,nc))
#    medpt=np.zeros((nt,nc))
#    if tipper=='y':
#        medtipr=np.zeros((nt,nc))
#    
#    #make a list of periods from the longest period list
#    plst=np.logspace(np.log10(minper),np.log10(maxper),num=nt,base=10)
#    pdict=dict([(ii,jj) for jj,ii in enumerate(plst)])
#    
#    #put data into arrays
#    for ii,mm in enumerate(invlst):
#        mperiod=mm.keys()
#        for jj,mp in enumerate(mperiod):
#            for kk in pdict.keys():
#                if mp>kk*(1-ptol) and mp<kk*(1+ptol):
#                    ll=pdict[kk]
#                    medinv[ll,ii]=invlst[ii][mp]
#                    medpt[ll,ii]=ptlst[ii][mp]
#                    if tipper=='y':
#                        medtipr[ll,ii]=tiprlst[ii][mp]
#                else:
#                    pass
#        
#    #-----Plot Histograms of the strike angles-----------------------------
#    if prange=='data':
#        brange=np.arange(np.floor(np.log10(minper)),
#                         np.ceil(np.log10(maxper)),1)
#    else:
#        brange=np.arange(np.floor(prange[0]),np.ceil(prange[1]),1)
#    
#    
#    #font dictionary
#    fd={'size':fs,'weight':'normal'}
#    
#    #------------------plot indivdual decades---------------------------------
#    if plottype==1:
#        #plot specs
#        plt.rcParams['figure.subplot.hspace']=.3
#        plt.rcParams['figure.subplot.wspace']=.3
#        
#        fig3=plt.figure(fignum,dpi=dpi)
#        plt.clf()
#        nb=len(brange)
#        for jj,bb in enumerate(brange,1):
#            #make subplots for invariants and phase tensor azimuths
#            if tipper=='n':
#                axhinv=fig3.add_subplot(2,nb,jj,polar=True)
#                axhpt=fig3.add_subplot(2,nb,jj+nb,polar=True)
#                axlst=[axhinv,axhpt]
#            if tipper=='y':
#                axhinv=fig3.add_subplot(3,nb,jj,polar=True)
#                axhpt=fig3.add_subplot(3,nb,jj+nb,polar=True)
#                axhtip=fig3.add_subplot(3,nb,jj+2*nb,polar=True)
#                axlst=[axhinv,axhpt,axhtip]
#            
#            #make a list of indicies for each decades    
#            binlst=[]
#            for ii,ff in enumerate(plst):
#                if ff>10**bb and ff<10**(bb+1):
#                    binlst.append(ii)
#            
#            #extract just the subset for each decade
#            hh=medinv[binlst,:]
#            gg=medpt[binlst,:]
#            if tipper=='y':
#                tr=medtipr[binlst,:]
#                trhist=np.histogram(tr[np.nonzero(tr)].flatten(),bins=72,
#                                       range=(-180,180))
#                bartr=axhtip.bar((trhist[1][:-1])*np.pi/180,trhist[0],
#                                 width=5*np.pi/180)
#                for cc,bar in enumerate(bartr):
#                    fc=float(trhist[0][cc])/trhist[0].max()*.9
#                    bar.set_facecolor((0,1-fc/2,fc))
#                        
#            
#            #estimate the histogram for the decade for invariants and pt
#            invhist=np.histogram(hh[np.nonzero(hh)].flatten(),bins=72,
#                                    range=(-180,180))
#            pthist=np.histogram(gg[np.nonzero(gg)].flatten(),bins=72,
#                                   range=(-180,180))
#            
#            #plot the histograms    
#            barinv=axhinv.bar((invhist[1][:-1])*np.pi/180,invhist[0],width=5*np.pi/180)
#            barpt=axhpt.bar((pthist[1][:-1])*np.pi/180,pthist[0],width=5*np.pi/180)
#            
#            for cc,bar in enumerate(barinv):
#                fc=float(invhist[0][cc])/invhist[0].max()*.8
#                bar.set_facecolor((fc,0,1-fc))
#            for cc,bar in enumerate(barpt):
#                fc=float(pthist[0][cc])/pthist[0].max()*.8
#                bar.set_facecolor((fc,1-fc,0))
#                
#            #make axis look correct with N to the top at 90.
#            for aa,axh in enumerate(axlst):
#                axh.xaxis.set_major_locator(MultipleLocator(30*np.pi/180))
#                axh.xaxis.set_ticklabels(['E','','',
#                                          'N','','',
#                                          'W','','',
#                                          'S','',''])
#                axh.grid(alpha=galpha)
#                if aa==0:
#                    axh.set_xlim(-90*np.pi/180,270*np.pi/180)
#                    axh.text(np.pi,axh.get_ylim()[1]*tpad,
#                             '{0:.1f}$^o$'.format(90-np.median(hh[np.nonzero(hh)])),
#                              horizontalalignment='center',
#                              verticalalignment='baseline',
#                              fontdict={'size':fs-nb},
#                              bbox={'facecolor':(.9,0,.1),'alpha':.25})
#                    print '-----Period Range {0:.3g} to {1:.3g} (s)-----'.format(10**bb,
#                          10**(bb+1))
#                         
#                    print '   *Z-Invariants:  median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
#                            90-np.median(hh[np.nonzero(hh)]),
#                            90-invhist[1][np.where(invhist[0]==invhist[0].max())[0][0]],
#                            90-np.mean(hh[np.nonzero(hh)])) 
#                    if bb==-5:
#                        axh.set_title('10$^{-5}$-10$^{-4}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==-4:
#                        axh.set_title('10$^{-4}$-10$^{-3}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==-3:
#                        axh.set_title('10$^{-3}$-10$^{-2}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==-2:
#                        axh.set_title('10$^{-2}$-10$^{-1}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==-1:
#                        axh.set_title('10$^{-1}$-10$^{0}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==0:
#                        axh.set_title('10$^{0}$-10$^{1}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==1:
#                        axh.set_title('10$^{1}$-10$^{2}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==2:
#                        axh.set_title('10$^{2}$-10$^{3}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==3:
#                        axh.set_title('10$^{3}$-10$^{4}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==4:
#                        axh.set_title('10$^{4}$-10$^{5}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    elif bb==5:
#                        axh.set_title('10$^{5}$-10$^{6}$s',fontdict=fd,
#                                      bbox={'facecolor':'white','alpha':galpha})
#                    axh.titleOffsetTrans._t=(0,.1)
#        
#                elif aa==1:
#                    axh.set_xlim(-180*np.pi/180,180*np.pi/180)
#                    axh.text(np.pi,axh.get_ylim()[1]*tpad,
#                             '{0:.1f}$^o$'.format(90-np.median(gg[np.nonzero(gg)])),
#                              horizontalalignment='center',
#                              verticalalignment='baseline',
#                              fontdict={'size':fs-nb},
#                              bbox={'facecolor':(.9,.9,0),'alpha':galpha})
#                    print '   *PT Strike:     median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
#                            90-np.median(gg[np.nonzero(gg)]),
#                            90-pthist[1][np.where(pthist[0]==pthist[0].max())[0][0]],
#                            90-np.mean(gg[np.nonzero(gg)]))
#                    if tipper!='y':
#                        print '\n'
#                    if nb>5: 
#                        if bb==-5:
#                            axh.set_title('10$^{-5}$-10$^{-4}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-4:
#                            axh.set_title('10$^{-4}$-10$^{-3}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-3:
#                            axh.set_title('10$^{-3}$-10$^{-2}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-2:
#                            axh.set_title('10$^{-2}$-10$^{-1}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-1:
#                            axh.set_title('10$^{-1}$-10$^{0}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==0:
#                            axh.set_title('10$^{0}$-10$^{1}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==1:
#                            axh.set_title('10$^{1}$-10$^{2}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==2:
#                            axh.set_title('10$^{2}$-10$^{3}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==3:
#                            axh.set_title('10$^{3}$-10$^{4}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==4:
#                            axh.set_title('10$^{4}$-10$^{5}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==5:
#                            axh.set_title('10$^{5}$-10$^{6}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                elif aa==2:
#                    axh.set_xlim(-180*np.pi/180,180*np.pi/180)
#                    axh.text(np.pi,axh.get_ylim()[1]*tpad,
#                             '{0:.1f}$^o$'.format(90-np.median(tr[np.nonzero(tr)])),
#                              horizontalalignment='center',
#                              verticalalignment='baseline',
#                              fontdict={'size':fs-nb},
#                              bbox={'facecolor':(0,.1,.9),'alpha':galpha})
#                    print '   *Tipper Strike: median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
#                            90-np.median(tr[np.nonzero(tr)]),
#                            90-trhist[1][np.where(trhist[0]==trhist[0].max())[0][0]],
#                            90-np.mean(tr[np.nonzero(tr)])) 
#                    print '\n'
#                    if nb>5: 
#                        if bb==-5:
#                            axh.set_title('10$^{-5}$-10$^{-4}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-4:
#                            axh.set_title('10$^{-4}$-10$^{-3}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-3:
#                            axh.set_title('10$^{-3}$-10$^{-2}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-2:
#                            axh.set_title('10$^{-2}$-10$^{-1}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==-1:
#                            axh.set_title('10$^{-1}$-10$^{0}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==0:
#                            axh.set_title('10$^{0}$-10$^{1}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==1:
#                            axh.set_title('10$^{1}$-10$^{2}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==2:
#                            axh.set_title('10$^{2}$-10$^{3}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==3:
#                            axh.set_title('10$^{3}$-10$^{4}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==4:
#                            axh.set_title('10$^{4}$-10$^{5}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                        elif bb==5:
#                            axh.set_title('10$^{5}$-10$^{6}$s',fontdict=fd,
#                                          bbox={'facecolor':'white',
#                                                'alpha':galpha})
#                                  
#                if jj==1:
#                    if aa==0:
#                        axh.set_ylabel('Strike (Z)',fontdict=fd,labelpad=5000./dpi,
#                                       bbox={'facecolor':(.9,0,.1),'alpha':galpha})
#                    elif aa==1:
#                        axh.set_ylabel('PT Azimuth',fontdict=fd,labelpad=5000./dpi,
#                                       bbox={'facecolor':(.9,.9,0),'alpha':galpha})
#                    elif aa==2:
#                        axh.set_ylabel('Tipper Strike',fd,labelpad=5000./dpi,
#                                       bbox={'facecolor':(0,.1,.9),'alpha':galpha})
#                
#                plt.setp(axh.yaxis.get_ticklabels(),visible=False)
#                
#        print 'Note: North is assumed to be 0 and the strike angle is measured'+\
#              'clockwise positive.'
#        
#        plt.show()
#    
#    #------------------Plot strike angles for all period ranges--------------------
#    elif plottype==2:
#        #plot specs
#        plt.rcParams['figure.subplot.left']=.07
#        plt.rcParams['figure.subplot.right']=.98
#        plt.rcParams['figure.subplot.bottom']=.100
#        plt.rcParams['figure.subplot.top']=.88
#        plt.rcParams['figure.subplot.hspace']=.3
#        plt.rcParams['figure.subplot.wspace']=.2
#        
#        fig3=plt.figure(fignum,dpi=dpi)
#        plt.clf()
#        #make subplots for invariants and phase tensor azimuths
#        if tipper=='n':
#            axhinv=fig3.add_subplot(1,2,1,polar=True)
#            axhpt=fig3.add_subplot(1,2,2,polar=True)
#            axlst=[axhinv,axhpt]
#        else:
#            axhinv=fig3.add_subplot(1,3,1,polar=True)
#            axhpt=fig3.add_subplot(1,3,2,polar=True)
#            axhtip=fig3.add_subplot(1,3,3,polar=True)
#            axlst=[axhinv,axhpt,axhtip]
#        
#        #make a list of indicies for each decades    
#        binlst=[pdict[ff] for ff in plst 
#                if ff>10**brange.min() and ff<10**brange.max()]
#        
#        #extract just the subset for each decade
#        hh=medinv[binlst,:]
#        gg=medpt[binlst,:]
#        
#        #estimate the histogram for the decade for invariants and pt
#        invhist=np.histogram(hh[np.nonzero(hh)].flatten(),bins=72,range=(-180,180))
#        pthist=np.histogram(gg[np.nonzero(gg)].flatten(),bins=72,range=(-180,180))
#        
#        #plot the histograms    
#        barinv=axhinv.bar((invhist[1][:-1])*np.pi/180,invhist[0],width=5*np.pi/180)
#        barpt=axhpt.bar((pthist[1][:-1])*np.pi/180,pthist[0],width=5*np.pi/180)
#        
#        for cc,bar in enumerate(barinv):
#            fc=float(invhist[0][cc])/invhist[0].max()*.8
#            bar.set_facecolor((fc,0,1-fc))
#        for cc,bar in enumerate(barpt):
#            fc=float(pthist[0][cc])/pthist[0].max()*.8
#            bar.set_facecolor((fc,1-fc,0))
#            
#        if tipper=='y':
#                tr=medtipr[binlst,:]
#                trhist=np.histogram(tr[np.nonzero(tr)].flatten(),bins=72,
#                                       range=(-180,180))
#                bartr=axhtip.bar((trhist[1][:-1])*np.pi/180,trhist[0],
#                                 width=5*np.pi/180)
#                for cc,bar in enumerate(bartr):
#                    fc=float(trhist[0][cc])/trhist[0].max()*.9
#                    bar.set_facecolor((0,1-fc/2,fc))
#                    
#        #make axis look correct with N to the top at 90.
#        for aa,axh in enumerate(axlst):
#            axh.xaxis.set_major_locator(MultipleLocator(30*np.pi/180))
#            axh.grid(alpha=galpha)
#            plt.setp(axh.yaxis.get_ticklabels(),visible=False)
#            axh.xaxis.set_ticklabels(['E','','',
#                                      'N','','',
#                                      'W','','',
#                                      'S','',''])
#                                      
#            #put median angle in a box at bottom of polar diagram
#            if aa==0:
#                axh.set_ylim(0,invhist[0].max())
#                axh.text(170*np.pi/180,axh.get_ylim()[1]*.65,
#                         '{0:.1f}$^o$'.format(90-np.median(hh[np.nonzero(hh)])),
#                          horizontalalignment='center',
#                          verticalalignment='baseline',
#                          fontdict={'size':fs-2},
#                          bbox={'facecolor':(.9,0,.1),'alpha':.25})
#                print '-----Period Range {0:.3g} to {1:.3g} (s)-----'.format(min(plst),
#                          max(plst))
#                         
#                print '   *Z-Invariants:  median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
#                        90-np.median(hh[np.nonzero(hh)]),
#                        90-invhist[1][np.where(invhist[0]==invhist[0].max())[0][0]],
#                        90-np.mean(hh[np.nonzero(hh)])) 
#    
#            elif aa==1:
#                axh.set_ylim(0,pthist[0].max())
#                axh.text(170*np.pi/180,axh.get_ylim()[1]*.65,
#                         '{0:.1f}$^o$'.format(90-np.median(gg[np.nonzero(gg)])),
#                          horizontalalignment='center',
#                          verticalalignment='baseline',
#                          fontdict={'size':fs-2},
#                          bbox={'facecolor':(.9,.9,0),'alpha':galpha})
#                print '   *PT Strike:     median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
#                        90-np.median(gg[np.nonzero(gg)]),
#                        90-pthist[1][np.where(pthist[0]==pthist[0].max())[0][0]],
#                        90-np.mean(gg[np.nonzero(gg)]))
#                if tipper!='y':
#                    print '\n'
#            elif aa==2:
#                axh.set_ylim(0,trhist[0].max())
#                axh.text(170*np.pi/180,axh.get_ylim()[1]*.65,
#                         '{0:.1f}$^o$'.format(90-np.median(tr[np.nonzero(tr)])),
#                          horizontalalignment='center',
#                          verticalalignment='baseline',
#                          fontdict={'size':fs-2},
#                          bbox={'facecolor':(0,.1,.9),'alpha':galpha})
#                print '   *Tipper Stike:  median={0:.1f} mode={1:.1f} mean={2:.1f}\n'.format(
#                        90-np.median(tr[np.nonzero(tr)]),
#                        90-trhist[1][np.where(trhist[0]==trhist[0].max())[0][0]],
#                        90-np.mean(tr[np.nonzero(tr)]))
#            #set the title of the diagrams
#            if aa==0:
#                axh.set_title('Strike (Z)',fontdict=fd,
#                               bbox={'facecolor':(.9,0,.1),'alpha':galpha})
#            elif aa==1:
#                axh.set_title('PT Azimuth',fontdict=fd,
#                               bbox={'facecolor':(.9,.9,0),'alpha':galpha})
#            elif aa==2:
#                axh.set_title('Tipper Strike',fontdict=fd,
#                               bbox={'facecolor':(0,.1,.9),'alpha':galpha})
#            axh.titleOffsetTrans._t=(0,.15)
#            
#                
#                
#        print 'Note: North is assumed to be 0 and the strike angle is measured '+\
#              'clockwise positive.'
#        
#        plt.show()
