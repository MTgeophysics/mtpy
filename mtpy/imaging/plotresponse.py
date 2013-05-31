# -*- coding: utf-8 -*-
"""
============
plotresponse
============

Plots the resistivity and phase for different modes and components

Created on Thu May 30 16:54:08 2013

@author: jpeacock-pr
"""
#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import mtpy.utils.exceptions as mtex
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl

#==============================================================================


#==============================================================================
#  Plot apparent resistivity and phase
#==============================================================================
class PlotResponse(object):
    """
    Plots Resistivity and phase for the different modes of the MT response.  At
    the moment is supports the input of an .edi file. Other formats that will
    be supported are the impedance tensor and errors with an array of periods
    and .j format.
    
    The normal use is to input an .edi file, however it would seem that not
    everyone uses this format, so you can input the data and put it into 
    arrays or objects like class mtpy.core.z.Z.  Or if the data is in 
    resistivity and phase format they can be input as arrays or a class
    mtpy.imaging.mtplot.ResPhase.  Or you can put it into a class
    mtpy.imaging.mtplot.MTplot.
    
    The plot places the apparent resistivity in log scale in the top panel(s), 
    depending on the plotnum.  The phase is below this, note that 180 degrees
    has been added to the yx phase so the xy and yx phases plot in the same
    quadrant.  Both the resistivity and phase share the same x-axis which is in
    log period, short periods on the left to long periods on the right.  So
    if you zoom in on the plot both plots will zoom in to the same 
    x-coordinates.  If there is tipper information, you can plot the tipper
    as a third panel at the bottom, and also shares the x-axis.  The arrows are
    in the convention of pointing towards a conductor.  The xx and yy 
    components can be plotted as well, this adds two panels on the right.  
    Here the phase is left unwrapped.
    
    To manipulate the plot you can change any of the attributes listed below
    and call redraw_plot().  If you know more aout matplotlib and want to 
    change axes parameters, that can be done by changing the parameters in the
    axes attributes and then call update_plot(), note the plot must be open.
    
    
    Arguments:
    ----------
        **filename** : string
                       filename containing impedance (.edi) is the only 
                       format supported at the moment
                       
        **z_array** : np.array((nf, 2, 2), dtype='complex')
                impedance tensor with length of nf -> the number of freq
                *default* is None
                
        **z_err_array** : np.array((nf, 2, 2), dtype='real')
                    impedance tensor error estimates, same shape as z.
                    *default* is None
                    
        **res_array** : np.array((nf, 2, 2))
                        array of resistivity values in linear scale.
                        *default* is None
                        
        **res_err_array** : np.array((nf, 2, 2))
                            array of resistivity error estimates, same shape 
                            as res_array. *default* is None
                            
        **phase_array** : np.array((nf, 2, 2))
                          array of phase values in degrees, same shape as 
                          res_array. *default* is None
                          
        **phase_err_array** : np.array((nf, 2, 2))
                              array of phase error estimates, same shape as 
                              phase_array. *default* is None
                              
        **tipper_array** : np.array((nf, 1, 2), dtype='complex')
                           array of tipper values for tx, ty. *default* is None
                           
        **tipper_err_array** : np.array((nf, 1, 2))
                               array of tipper error estimates, same shape as
                               tipper_array. *default* is None
                               
        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None
                      
        **tipper_object** : class mtpy.core.z.Tipper
                            object of mtpy.core.z. If this is input be sure the
                            attribute z.freq is filled.  
                            *default* is None 
                            
        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None
                      
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
                    
        **plot_tipper** : [ 'yri' | 'yr' | 'yi' | 'n' ]
                          Plots the tipper in a bottom pannel
                          * 'yri'  --> plots the real and imaginar parts
                          * 'yr'   --> plots just the real part
                          * 'yi'   --> plots just the imaginary part
                          
                          **Note:** the convention is to point towards a 
                          conductor.  Can change this by setting the
                          parameter arrow_direction = 1.
                          
        **plot_strike** : [ 'y' | 1 | 2 | 3 | 'n' ]
                          Plots the strike angle from different parameters:
                              * 'y'  --> plots strike angle determined from 
                                         the invariants of Weaver et al. [2000]
                                         and the phase tensor of
                                         Caldwell et al. [2004], if Tipper is 
                                         plotted the strike of the tipper is
                                         also plotted.
                                         
                               * 1  --> plots strike angle determined from 
                                        the invariants of Weaver et al. [2000]
                               * 2  --> plots strike angle determined from 
                                        the phase tensor of 
                                        Caldwell et al. [2004]
                               * 3  --> plots strike angle determined from 
                                        the tipper
                               * 'n' --> doesn't plot the strike, *default*
                               
        **plot_skew** : [ 'y' | 'n' ]
                        plots the skew angle calculated from the phase tensor
                        in the same subfigure as the strike angle.  The axes
                        are labelled on the right hand side.
                            * 'y' --> plots skew angle
                            * 'n' --> does not plot skew angle *default*
                          
        **dpi** : int
                 dots-per-inch resolution, *default* is 300
                    
                        
        :Example: ::
            
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> edifile = r"/home/MT01/MT01.edi"
            >>> rp1 = mtplot.PlotResPhase(edifile, plotnum=2)
            >>> # plots all 4 components
            >>> rp1 = mtplot.PlotResPhase(edifile, plot_tipper='yr')
            >>> # plots the real part of the tipper
            
    Attributes:
    -----------
        -fn             filename to be plotted (only supports .edi so far) 
        -fignum         figure number for plotting
        -plotnum        plot type, see arguments for details 
        -title          title of the plot, *default* is station name
        -dpi            Dots-per-inch resolution of plot, *default* is 300
        -rotz           Rotate impedance tensor by this angle (deg) assuming
                        that North is 0 and angle is positive clockwise
                        
        -plot_tipper    string to tell the program to plot tipper arrows or 
                        not, see accepted values above in arguments
        
        -plot_strike    string or integer telling the program to plot the 
                        strike angle, see values above in arguments
                        
        -plot_skew      string to tell the program to plot skew angle.
                        The skew is plotted in the same subplot as the strike
                        angle at the moment
                        
                
        -period          period array cooresponding to the impedance tensor
        -font_size       size of font for the axis ticklabels, note that the 
                         axis labels will be font_size+2
        
        -axr             matplotlib.axes object for the xy,yx resistivity plot.  
        -axp             matplotlib.axes object for the xy,yx phase plot
        -axt             matplotlib.axes object for the tipper plot
        -ax2r            matplotlib.axes object for the xx,yy resistivity plot
        -ax2p            matplotlib.axes object for the xx,yy phase plot
        -axs             matplotlib.axes object for the strike plot
        -axs2            matplotlib.axes object for the skew plot       
        ..
        
             **Note:** that from these axes object you have control of the
             plot.  You can do this by changing any parameter in the 
             axes object and then calling update_plot()
        
        -erxyr          class matplotlib.container.ErrorbarContainer for 
                        xy apparent resistivity.
        -erxyp          class matplotlib.container.ErrorbarContainer for 
                        xy.
        -eryxr          class matplotlib.container.ErrorbarContainer for 
                        yx apparent resistivity.
        -eryxp          class matplotlib.container.ErrorbarContainer for 
                        yx phase.
        
        ..
        
            **Note:** that from these line objects you can manipulate the 
            error bar properties and then call update_plot()
                     
        -xy_ls           line style for xy and xx components, *default* is None       
        -yx_ls           line style for yx and yy components, *default* is None        
        -det_ls          line style for determinant, *default* is None
        
        -xy_marker       marker for xy and xx, *default* is squares
        -yx_marker       marker for yx and yy, *default* is circles
        -det_marker      marker for determinant, *default* is diamonds
        
        -xy_color        marker color for xy and xx, *default* is blue
        -yx_color        marker color for yx and yy, *default* is red
        -det_color       marker color for determinant, *default* is green
        
        -xy_mfc          marker face color for xy and xx, *default* is None 
        -yx_mfc          marker face color for yx and yy, *default* is None
        -det_mfc         marker face color for determinant, *default* is None
        
        -skew_marker     marker for skew angle, *default* is 'd'
        -skew_color      color for skew angle, *default* is 'orange'
        
        -strike_inv_marker  marker for strike angle determined by invariants
                            *default* is '^'
        -strike_inv_color   color for strike angle determined by invaraiants 
                            *default* is (.2, .2, .7)
        -strike_pt_marker  marker for strike angle determined by pt, 
                           *default* is'v'
        -strike_pt_color   color for strike angle determined by pt
                           *default* is (.7, .2, .2)
        
        -strike_tip_marker  marker for strike angle determined by tipper
                            *default* is '>'
        -strike_tip_color   color for strike angle determined by tipper
                            *default* is (.2, .7, .2)
        
        -marker_size     size of marker in relative dimenstions, *default* is 2
        -marker_lw       line width of marker, *default* is 100./dpi 
        ..
        
         *For more on line and marker styles see matplotlib.lines.Line2D*
        
        -arrow_lw          line width of the arrow, *default* is 0.75
        -arrow_head_width  head width of the arrow, *default* is 0 for no arrow
                           head.  Haven't found a good way to scale the arrow
                           heads in a log scale.
                         
        -arrow_head_height  head width of the arrow, *default* is 0 for no arrow
                            head.  Haven't found a good way to scale the arrow
                            heads in a log scale.
                          
        -arrow_color_real  color of the real arrows, *default* is black
        -arrow_color_imag  color of the imaginary arrows, *default* is blue
        
        -arrow_direction   0 for pointing towards a conductor and -1 for 
                           pointing away from a conductor.
                          

        -xlimits        limits on the x-limits (period), *default* is None
                        which will estimate the min and max from the data, 
                        setting the min as the floor(min(period)) and the max
                        as ceil(max(period)).  Input in linear scale if you
                        want to change the period limits, ie. (.1,1000)
                        
        -res_limits     limits on the resistivity, *default* is None, which 
                        will estimate the min and max from the data, rounding
                        to the lowest and highest increments to the power of 10
                        Input in linear scale if you want to change them, 
                        ie. (1,10000). Note this only sets the xy and yx 
                        components, not the xx and yy.
                        
        -phase_limits   limits on the phase, *default* is (0,90) but will 
                        adapt to the data if there is phase above 90 or below
                        0.  Input in degrees.  Note this only changes the xy
                        and yx components.
                        
        -tipper_limits  limits of the y-axis, *default* is (-1,1)
        
        -skew_limits    limits for skew angle, *default* is (-9,9)
        
        -strike_limits  limits for strike angle, *default* is (-90,90)   
        


    """   
    
    def __init__(self, filename=None, z_array=None, z_err_array=None, 
                 period=None, fignum=1, plotnum=1, title=None, dpi=300, 
                 rotz=0, plot_yn='y', plot_tipper='n', plot_strike='n',
                 plot_skew='n', tipper_array=None, tipper_err_array=None, 
                 tipper_object=None, res_array=None, res_err_array=None,
                 phase_array=None, phase_err_array=None,
                 res_phase_object=None, z_object=None, mt_object=None):
        
        #--> initialize an MTplot object
        if mt_object is None:
            self._mt = mtpl.MTplot(filename=filename, 
                                   z=z_array, 
                                   z_err=z_err_array,
                                   z_object=z_object,
                                   phase_array=phase_array, 
                                   res_array=res_array, 
                                   phase_err_array=phase_err_array,
                                   res_err_array=res_err_array,
                                   tipper=tipper_array, 
                                   tipper_err=tipper_err_array,
                                   tipper_object=tipper_object,
                                   period=period)
                              
        else:
            self._mt = mt_object

        
        #if you have a resistivity and phase object 
        if res_phase_object is not None:
            if period is None:
                raise mtex.MTpyError_Z('Need to input period array for '+\
                                           'plotting')
            self.rp = res_phase_object
            
        
        #set some of the properties as attributes much to Lars' discontent
        self.fn = filename
        self.fignum = fignum
        self.plotnum = plotnum
        self.title = title
        self.dpi = dpi
        self.rotz = rotz
        
        #-->line properties
        #line style between points
        self.xy_ls = 'None'        
        self.yx_ls = 'None'        
        self.det_ls = 'None'        
        
        #outline color
        self.xy_color = 'b'
        self.yx_color = 'r'
        self.det_color = 'g'
        
        #face color
        self.xy_mfc = 'None'
        self.yx_mfc = 'None'
        self.det_mfc = 'None'
        
        #maker
        self.xy_marker = 's'
        self.yx_marker = 'o'
        self.det_marker = 'd'
        
        #size
        self.marker_size = 2
        
        #marker line width
        self.marker_lw = 100./self.dpi
        
        #set plot limits
        self.xlimits = None
        self.res_limits = None
        self.phase_limits = None
        self.tipper_limits = None
        self.strike_limits = None
        self.skew_limits = None
        
        #set font parameters
        self.font_size = 7
        
        #set period
        self.period = self._mt.period
        
        #set plot tipper or not
        self.plot_tipper = plot_tipper
        
        #plot strike angle or not
        self.plot_strike = plot_strike
        
        #plot skew angle
        self.plot_skew = plot_skew
        
        #set arrow properties
        self.arrow_lw = .75
        self.arrow_head_width = 0.0
        self.arrow_head_height = 0.0
        self.arrow_color_real = 'k'
        self.arrow_color_imag = 'b'
        self.arrow_direction = 0
        
        #skew properties
        self.skew_color = (.85, .35, 0)
        self.skew_marker = 'd'
        
        #strike properties
        self.strike_inv_marker = '^'
        self.strike_inv_color = (.2, .2, .7)
        
        self.strike_pt_marker = 'v'
        self.strike_pt_color = (.7, .2, .2)
        
        self.strike_tip_marker = '>'
        self.strike_tip_color = (.2, .7, .2)


        #plot on initializing
        if plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plotResPhase(filename,fignum) will plot the apparent resistivity and 
        phase for a single station. 
 
        """
        
        #set figure size according to what the plot will be.
        if self.plotnum == 1 or self.plotnum == 3:
            self.figsize = [5, 7]
        
        elif self.plotnum == 2:
            self.figsize = [7, 7]
            
        #--> rotate the impedance tensor if desired
        if self.rotz != 0:
            self._mt.rot_z = self.rotz
        
        #get the reistivity and phase object
        try:
            self.rp = self._mt.get_ResPhase()
        except AttributeError:
            pass
        
        #set x-axis limits from short period to long period
        if self.xlimits == None:
            self.xlimits = (10**(np.floor(np.log10(self.period[0]))),
                            10**(np.ceil(np.log10((self.period[-1])))))
        if self.phase_limits == None:
            pass
            
        if self.res_limits == None:
            self.res_limits = (10**(np.floor(
                                    np.log10(min([self.rp.resxy.min(),
                                                  self.rp.resyx.min()])))),
                              10**(np.ceil(
                                    np.log10(max([self.rp.resxy.max(),
                                                  self.rp.resyx.max()])))))

        #set some parameters of the figure and subplot spacing
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .93
        plt.rcParams['figure.subplot.left'] = .80
        if self.plot_skew == 'y':
            plt.rcParams['figure.subplot.right'] = .90
        else:
            plt.rcParams['figure.subplot.right'] = .98
        
        #set the font properties for the axis labels
        fontdict = {'size':self.font_size+2, 'weight':'bold'}
        
        #create a grid to place the figures into, set to have 2 rows and 2 
        #columns to put any of the 4 components.  Make the phase plot
        #slightly shorter than the apparent resistivity plot and have the two
        #close to eachother vertically.  If there is tipper add a 3rd row and
        #if there is strike add another row
        if self.plot_tipper.find('y') == 0:
            if self.plot_strike != 'n' or self.plot_skew == 'y':
                gs = gridspec.GridSpec(4, 2, height_ratios=[2, 1.5, 1, 1], 
                                       hspace=.05)
            else:
                gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1.5, 1], 
                                       hspace=.05)
        else:
            if self.plot_strike != 'n' or self.plot_skew == 'y':
                gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1.5, 1],
                                       hspace=.05)
            else:
                gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1.5], 
                                       hspace=.01)
        
        #make figure instance
        self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
        
        #--> make figure for xy,yx components
        if self.plotnum == 1 or self.plotnum == 3:
            #space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.1)
            
            #--> create the axes instances
            if self.plot_tipper.find('y') == 0:
                #apparent resistivity axis
                self.axr = self.fig.add_subplot(gs[0, :])
                
                #phase axis that shares period axis with resistivity
                self.axp = self.fig.add_subplot(gs[1, :], sharex=self.axr)
                
                #tipper axis  
                self.axt = self.fig.add_subplot(gs[2, :], sharex=self.axr)
                
                #place y coordinate labels in the same location                
                self.axr.yaxis.set_label_coords(-.075, 0.5)
                self.axp.yaxis.set_label_coords(-.075, 0.5)
                self.axt.yaxis.set_label_coords(-.075, 0.5)
                if self.plot_strike != 'n' or self.plot_skew != 'n':
                    self.axs = self.fig.add_subplot(gs[3, :], sharex=self.axr)
                    self.axs.yaxis.set_label_coords(-.075, 0.5)
            else:
                #apparent resistivity axis
                self.axr = plt.subplot(gs[0, :])
                
                #phase axis which shares the x axis with the resistivity
                self.axp = plt.subplot(gs[1, :], sharex=self.axr)
                
                #place the y-coordinate labels in the same location
                self.axr.yaxis.set_label_coords(-.075, 0.5)
                self.axp.yaxis.set_label_coords(-.075, 0.5)
                
                #add stike axis if desired
                if self.plot_strike != 'n' or self.plot_skew != 'n':
                    self.axs = self.fig.add_subplot(gs[2, :], sharex=self.axr)
                    self.axs.yaxis.set_label_coords(-.075, 0.5)
            
        #--> make figure for all 4 components
        elif self.plotnum == 2:
            
            #space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.07)
            
            #--> create the axes instances
            if self.plot_tipper.find('y') == 0:
                #apparent resistivity axis
                self.axr = self.fig.add_subplot(gs[0, 0])
                
                #phase axis that shares period axis with resistivity
                self.axp = self.fig.add_subplot(gs[1, 0], sharex=self.axr)
                
                #tipper axis that shares period with resistivity 
                self.axt = self.fig.add_subplot(gs[2, :], sharex=self.axr)
                
                #place y coordinate labels in the same location                
                self.axr.yaxis.set_label_coords(-.095, 0.5)
                self.axp.yaxis.set_label_coords(-.095, 0.5)
                self.axt.yaxis.set_label_coords(-.095, 0.5)
                
                #add strike axis if desired
                if self.plot_strike != 'n' or self.plot_skew != 'n':
                    self.axs = self.fig.add_subplot(gs[3, :], sharex=self.axr)
                    self.axs.yaxis.set_label_coords(-.095, 0.5)
            else:
                #apparent resistivity axis
                self.axr = plt.subplot(gs[0, 0])
                
                #phase axis which shares the x axis with the resistivity
                self.axp = plt.subplot(gs[1, 0], sharex=self.axr)
                
                #place the y-coordinate labels in the same location
                self.axr.yaxis.set_label_coords(-.095, 0.5)
                self.axp.yaxis.set_label_coords(-.095, 0.5)
                
                #add strike axis if desired
                if self.plot_strike != 'n' or self.plot_skew != 'n':
                    self.axs = self.fig.add_subplot(gs[2, :], sharex=self.axr)
                    self.axs.yaxis.set_label_coords(-.095, 0.5)
        
        #---------plot the apparent resistivity--------------------------------
        #--> plot as error bars and just as points xy-blue, yx-red
        #res_xy
        self.ebxyr = self.axr.errorbar(self.period, 
                                       self.rp.resxy, 
                                       marker=self.xy_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.xy_mfc, 
                                       mec=self.xy_color, 
                                       mew=self.marker_lw, 
                                       ls=self.xy_ls, 
                                       yerr=self.rp.resxy_err, 
                                       ecolor=self.xy_color,
                                       capsize=self.marker_size,
                                       elinewidth=self.marker_lw)
        
        #res_yx                              
        self.ebyxr = self.axr.errorbar(self.period, 
                                       self.rp.resyx, 
                                       marker=self.yx_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.yx_mfc,
                                       mec=self.yx_color, 
                                       mew=self.marker_lw,
                                       ls=self.yx_ls, 
                                       yerr=self.rp.resyx_err, 
                                       ecolor=self.yx_color,
                                       capsize=self.marker_size,
                                       elinewidth=self.marker_lw)
                                      
        #--> set axes properties
        plt.setp(self.axr.get_xticklabels(), visible=False)
        self.axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                            fontdict=fontdict)
        self.axr.set_yscale('log')
        self.axr.set_xscale('log')
        self.axr.set_xlim(self.xlimits)
        self.axr.set_ylim(self.res_limits)
        self.axr.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)
        self.axr.legend((self.ebxyr[0], self.ebyxr[0]), 
                        ('$Z_{xy}$', '$Z_{yx}$'),
                        loc=3, 
                        markerscale=1, 
                        borderaxespad=.01,
                        labelspacing=.07, 
                        handletextpad=.2, 
                        borderpad=.02)
        
        #-----Plot the phase---------------------------------------------------
        #phase_xy
        self.ebxyp = self.axp.errorbar(self.period, 
                                       self.rp.phasexy, 
                                       marker=self.xy_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.xy_mfc,
                                       mec=self.xy_color, 
                                       mew=self.marker_lw,
                                       ls=self.xy_ls,
                                       yerr=self.rp.phasexy_err, 
                                       ecolor=self.xy_color,
                                       capsize=self.marker_size,
                                       elinewidth=self.marker_lw)
                                      
        #phase_yx: Note add 180 to place it in same quadrant as phase_xy
        self.ebyxp = self.axp.errorbar(self.period, 
                                       self.rp.phaseyx, 
                                       marker=self.yx_marker, 
                                       ms=self.marker_size, 
                                       mfc=self.yx_mfc, 
                                       mec=self.yx_color, 
                                       mew=self.marker_lw,
                                       ls=self.yx_ls, 
                                       yerr=self.rp.phaseyx_err, 
                                       ecolor=self.yx_color,
                                       capsize=self.marker_size,
                                       elinewidth=self.marker_lw)

        #check the phase to see if any point are outside of [0:90]
        if self.phase_limits == None:
            if min(self.rp.phasexy)<0 or min(self.rp.phaseyx)<0:
                pymin = min([min(self.rp.phasexy), min(self.rp.phaseyx)])
                if pymin > 0:
                   pymin = 0
            else:
                pymin = 0
            
            if max(self.rp.phasexy)>90 or max(self.rp.phaseyx)>90:
                pymax = min([max(self.rp.phasexy), max(self.rp.phaseyx)])
                if pymax < 91:
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
        self.axp.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)
        
        #-----plot tipper----------------------------------------------------              
        if self.plot_tipper.find('y') == 0:
            plt.setp(self.axp.xaxis.get_ticklabels(), visible=False)
            
            tp = self._mt.get_Tipper()
            
            txr = tp.mag_real*np.cos(tp.ang_real*np.pi/180+\
                                     np.pi*self.arrow_direction)
            tyr = tp.mag_real*np.sin(tp.ang_real*np.pi/180+\
                                     np.pi*self.arrow_direction)
    
            txi = tp.mag_imag*np.cos(tp.ang_imag*np.pi/180+\
                                     np.pi*self.arrow_direction)
            tyi = tp.mag_imag*np.sin(tp.ang_imag*np.pi/180+\
                                     np.pi*self.arrow_direction)
            
            nt = len(txr)
            
            tiplst = []
            tiplabel = []
            
            for aa in range(nt):
                xlenr = txr[aa]*self.period[aa]
                xleni = txi[aa]*self.period[aa]
                
                #scale the arrow head height and width to fit in a log scale
                if np.log10(self.period[aa]) < 0:
                    hwidth = self.arrow_head_width*\
                                10**(np.floor(np.log10(self.period[aa])))
                    hheight = self.arrow_head_height*\
                                 10**(np.floor(np.log10(self.period[aa])))
                else:
                    hwidth = self.arrow_head_width/\
                                 10**(np.floor(np.log10(self.period[aa])))
                    hheight = self.arrow_head_height/\
                                 10**(np.floor(np.log10(self.period[aa]))) 
                if np.log10(self.period[aa])<0:
                    alw = self.arrow_lw*self.period[aa]
                else:
                    alw = self.arrow_lw
                #--> plot real arrows
                if self.plot_tipper.find('r') > 0:
                    self.axt.arrow(self.period[aa],
                                   0,
                                   xlenr,
                                   tyr[aa],
                                   lw=alw,
                                   facecolor=self.arrow_color_real,
                                   edgecolor=self.arrow_color_real,
                                   head_width=hwidth,
                                   head_length=hheight,
                                   length_includes_head=False)
                    
                    if aa == 0:
                        line1 = self.axt.plot(0, 0, self.arrow_color_real)
                        tiplst.append(line1[0])
                        tiplabel.append('real')
                                   
                #--> plot imaginary arrows
                if self.plot_tipper.find('i') > 0:               
                    self.axt.arrow(self.period[aa],
                                   0,
                                   xleni,
                                   tyi[aa],
                                   lw=alw,
                                   facecolor=self.arrow_color_imag,
                                   edgecolor=self.arrow_color_imag,
                                   length_includes_head=False)
                    if aa == 0:              
                        line2 = self.axt.plot(0, 0, self.arrow_color_imag)
                        tiplst.append(line2[0])
                        tiplabel.append('imag')
                
            #make a line at 0 for reference
            self.axt.plot(self.period, [0]*nt, 'k', lw=.5)
        
          
            self.axt.legend(tiplst, tiplabel,
                            loc='upper left',
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.1,
                            prop={'size':self.font_size})

            #set axis properties            
            self.axt.yaxis.set_major_locator(MultipleLocator(.2))               
            self.axt.yaxis.set_minor_locator(MultipleLocator(.1))               
            self.axt.set_xlabel('Period (s)', fontdict=fontdict)
            self.axt.set_ylabel('Tipper', fontdict=fontdict)    
            
            self.axt.set_xscale('log')
            if self.tipper_limits is None:
                tmax = max([np.sqrt(txr.max()**2+tyr.max()**2),
                            np.sqrt(txi.max()**2+tyi.max()**2)])
                if tmax > 1:
                    tmax = .899
                            
                tmin = -min([np.sqrt(txr.min()**2+tyr.min()**2),
                            np.sqrt(txi.min()**2+tyi.min()**2)])
                if tmin < -1:
                    tmin = -.899
                            
                self.tipper_limits = (tmin-.1, tmax+.1)
            
            self.axt.set_ylim(self.tipper_limits)
            self.axt.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
                          
        #------plot strike angles----------------------------------------------
        if self.plot_strike != 'n' or self.plot_skew == 'y':
            try:
                plt.setp(self.axp.xaxis.get_ticklabels(), visible=False)
                plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
            except AttributeError:
                pass
            
            stlst = []
            stlabel = []
            st_maxlst = []
            st_minlst = []
            
            if self.plot_strike == 'y' or self.plot_strike == 1:
                #strike from invariants
                zinv = self._mt.get_Zinvariants()
                s1 = zinv.strike
                
                #fold angles so go from -90 to 90
                s1[np.where(s1>90)] = s1[np.where(s1>90)]-180
                s1[np.where(s1<-90)] = s1[np.where(s1<-90)]+180
                
                #plot strike with error bars
                ps1 = self.axs.errorbar(self.period, 
                                        s1, 
                                        marker=self.strike_inv_marker, 
                                        ms=self.marker_size, 
                                        mfc=self.strike_inv_color, 
                                        mec=self.strike_inv_color, 
                                        mew=self.marker_lw,
                                        ls='none', 
                                        yerr=zinv.strike_err, 
                                        ecolor=self.strike_inv_color,
                                        capsize=self.marker_size,
                                        elinewidth=self.marker_lw)
                                        
                stlst.append(ps1[0])
                stlabel.append('Z_inv')
                st_maxlst.append(s1.max())
                st_minlst.append(s1.min())
                                        
            if self.plot_strike == 'y' or self.plot_strike == 2:
                
                #strike from phase tensor
                pt = self._mt.get_PhaseTensor()
                s2, s2_err = pt.azimuth
                
                #fold angles to go from -90 to 90
                s2[np.where(s2>90)] = s2[np.where(s2>90)]-180
                s2[np.where(s2<-90)] = s2[np.where(s2<-90)]+180
                
                #plot strike with error bars
                ps2 = self.axs.errorbar(self.period, 
                                        s2, 
                                        marker=self.strike_pt_marker, 
                                        ms=self.marker_size, 
                                        mfc=self.strike_pt_color, 
                                        mec=self.strike_pt_color, 
                                        mew=self.marker_lw,
                                        ls='none', 
                                        yerr=s2_err, 
                                        ecolor=self.strike_pt_color,
                                        capsize=self.marker_size,
                                        elinewidth=self.marker_lw)
                                        
                stlst.append(ps2[0])
                stlabel.append('PT')
                st_maxlst.append(s2.max())
                st_minlst.append(s2.min())
            
            if self.plot_strike == 'y' or self.plot_strike == 3:
                #strike from tipper
                tp = self._mt.get_Tipper()
                s3 = tp.ang_real+90
                
                #fold to go from -90 to 90
                s3[np.where(s3 > 90)] -= 180
                s3[np.where(s3 < -90)] += 180
                
                #plot strike with error bars
                ps3 = self.axs.errorbar(self.period, 
                                        s3, 
                                        marker=self.strike_tip_marker, 
                                        ms=self.marker_size, 
                                        mfc=self.strike_tip_color, 
                                        mec=self.strike_tip_color, 
                                        mew=self.marker_lw,
                                        ls='none', 
                                        yerr=np.zeros_like(s3), 
                                        ecolor=self.strike_tip_color,
                                        capsize=self.marker_size,
                                        elinewidth=self.marker_lw)
                                        
                stlst.append(ps3[0])
                stlabel.append('Tip')
                st_maxlst.append(s3.max())
                st_minlst.append(s3.min())
            
            #------plot skew angle---------------------------------------------
            if self.plot_skew == 'y':
                #strike from phase tensor
                pt = self._mt.get_PhaseTensor()
                sk, sk_err = pt.beta
                
                
                self.axs2 = self.axs.twinx()
                ps4 = self.axs2.errorbar(self.period, 
                                        sk, 
                                        marker=self.skew_marker, 
                                        ms=self.marker_size, 
                                        mfc=self.skew_color, 
                                        mec=self.skew_color, 
                                        mew=self.marker_lw,
                                        ls='none', 
                                        yerr=sk_err, 
                                        ecolor=self.skew_color,
                                        capsize=self.marker_size,
                                        elinewidth=self.marker_lw)
                stlst.append(ps4[0])
                stlabel.append('Skew')
                if self.skew_limits is None:
                    self.skew_limits = (-9, 9)
                
                self.axs2.set_ylim(self.skew_limits)
                self.axs2.yaxis.set_major_locator(MultipleLocator(3))
                self.axs2.yaxis.set_minor_locator(MultipleLocator(1))
                self.axs2.set_ylabel('Skew', color=self.skew_color)
                self.axs2.set_xscale('log')
                for tl in self.axs2.get_yticklabels():
                    tl.set_color(self.skew_color)
                    
                st_minlst.append(0.0)
                st_maxlst.append(0.0)
                
            #--> set axes properties
            if self.strike_limits is None:
                stmin = min(st_minlst)
                if stmin-3 < -90:
                    stmin -= 3
                else:
                    stmin = -89.99
                    
                stmax = max(st_maxlst)
                if stmin+3 < 90:
                    stmin += 3
                else:
                    stmin = 89.99
                self.strike_limits = (-max([abs(stmin), abs(stmax)]),
                                       max([abs(stmin), abs(stmax)]))
                                        
                
            self.axs.plot(self.axr.get_xlim(), [0, 0], color='k', lw=.5)
            
            self.axs.set_ylabel('Strike',
                                fontdict=fontdict)
            self.axs.set_xlabel('Period (s)',
                                fontdict=fontdict)
            self.axs.set_ylim(self.strike_limits)
            self.axs.yaxis.set_major_locator(MultipleLocator(30))
            self.axs.yaxis.set_minor_locator(MultipleLocator(5))
            self.axs.set_xscale('log')
            self.axs.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
            try:
                self.axs.legend(stlst, 
                                stlabel,
                                loc=3, 
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02,
                                prop={'size':self.font_size-1})
            except:
                pass

        #===Plot the xx, yy components if desired==============================
        if self.plotnum == 2:
            #---------plot the apparent resistivity----------------------------
            self.axr2 = self.fig.add_subplot(gs[0, 1], sharex=self.axr)
            self.axr2.yaxis.set_label_coords(-.1, 0.5)
            
            #res_xx
            self.ebxxr = self.axr2.errorbar(self.period, 
                                            self.rp.resxx, 
                                            marker=self.xy_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.xy_mfc,
                                            mec=self.xy_color, 
                                            mew=self.marker_lw,
                                            ls=self.xy_ls, 
                                            yerr=self.rp.resxx_err, 
                                            ecolor=self.xy_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
            
            #res_yy                              
            self.ebyyr = self.axr2.errorbar(self.period, 
                                            self.rp.resyy, 
                                            marker=self.yx_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.yx_mfc, 
                                            mec=self.yx_color, 
                                            mew=self.marker_lw, 
                                            ls=self.yx_ls, 
                                            yerr=self.rp.resyy_err, 
                                            ecolor=self.yx_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)

            #--> set axes properties
            plt.setp(self.axr2.get_xticklabels(), visible=False)
            self.axr2.set_yscale('log')
            self.axr2.set_xscale('log')
            self.axr2.set_xlim(self.xlimits)
            self.axr2.grid(True, alpha=.25, 
                           which='both', 
                           color=(.25, .25, .25),
                           lw=.25)
                           
            self.axr2.legend((self.ebxxr[0], self.ebyyr[0]), 
                            ('$Z_{xx}$','$Z_{yy}$'),
                            loc=3, markerscale=1, 
                            borderaxespad=.01,
                            labelspacing=.07, 
                            handletextpad=.2, 
                            borderpad=.02)
            
            #-----Plot the phase-----------------------------------------------
            self.axp2 = self.fig.add_subplot(gs[1, 1], sharex=self.axr)
            
            self.axp2.yaxis.set_label_coords(-.1, 0.5)
            
            #phase_xx
            self.ebxxp = self.axp2.errorbar(self.period, 
                                            self.rp.phasexx, 
                                            marker=self.xy_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.xy_mfc, 
                                            mec=self.xy_color, 
                                            mew=self.marker_lw, 
                                            ls=self.xy_ls, 
                                            yerr=self.rp.phasexx_err,
                                            ecolor=self.xy_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
                                            
            #phase_yy
            self.ebyyp = self.axp2.errorbar(self.period,
                                            self.rp.phaseyy, 
                                            marker=self.yx_marker,
                                            ms=self.marker_size, 
                                            mfc=self.yx_mfc,
                                            mec=self.yx_color, 
                                            mew=self.marker_lw, 
                                            ls=self.yx_ls, 
                                            yerr=self.rp.phaseyy_err, 
                                            ecolor=self.yx_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
            
            #--> set axes properties
            self.axp2.set_xlabel('Period (s)', fontdict)
            self.axp2.set_xscale('log')
            self.axp2.set_ylim(ymin=-179.9, ymax=179.9)        
            self.axp2.yaxis.set_major_locator(MultipleLocator(30))
            self.axp2.yaxis.set_minor_locator(MultipleLocator(5))
            self.axp2.grid(True, alpha=.25, 
                           which='both', 
                           color=(.25, .25, .25),
                           lw=.25) 
                           
            if self.plot_tipper.find('y')==0 or self.plot_skew=='y':
                plt.setp(self.axp2.xaxis.get_ticklabels(), visible=False)
                plt.setp(self.axp2.xaxis.get_label(), visible=False)
        
        #===Plot the Determinant if desired==================================                          
        if self.plotnum == 3:
                
            #res_det
            self.ebdetr = self.axr.errorbar(self.period, 
                                            self.rp.resdet, 
                                            marker=self.det_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.det_mfc,
                                            mec=self.det_color, 
                                            mew=self.marker_lw, 
                                            ls=self.det_ls, 
                                            yerr=self.rp.resdet_err, 
                                            ecolor=self.det_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
        
            #phase_det
            self.ebdetp = self.axp.errorbar(self.period, 
                                            self.rp.phasedet, 
                                            marker=self.det_marker, 
                                            ms=self.marker_size, 
                                            mfc=self.det_mfc, 
                                            mec=self.det_color, 
                                            mew=self.marker_lw, 
                                            ls=self.det_ls, 
                                            yerr=self.rp.phasedet_err, 
                                            ecolor=self.det_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
                                            
            self.axr.legend((self.ebxyr[0], self.ebyxr[0], self.ebdetr[0]),
                            ('$Z_{xy}$','$Z_{yx}$','$\det(\mathbf{\hat{Z}})$'),
                            loc=3,
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.02)
        
        
        #make title and show
        if self.title != None:
            self.fig.suptitle(self.title, fontdict=fontdict)
        else:
            if self._mt.station != None:
                self.fig.suptitle(self._mt.station, fontdict=fontdict)
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
            save_fn = os.path.join(save_fn, self._mt.station+'_ResPhase.'+
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
        
        return "Plots Resistivity and phase for the different modes of the" +\
              "MT response."
