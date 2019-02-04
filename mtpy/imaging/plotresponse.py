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
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl

#==============================================================================


#==============================================================================
#  Plot apparent resistivity and phase
#==============================================================================
class PlotResponse(mtpl.PlotSettings):
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
    depending on the plot_num.  The phase is below this, note that 180 degrees
    has been added to the yx phase so the xy and yx phases plot in the same
    quadrant.  Both the resistivity and phase share the same x-axis which is in
    log period, short periods on the left to long periods on the right.  So
    if you zoom in on the plot both plots will zoom in to the same 
    x-coordinates.  If there is tipper information, you can plot the tipper
    as a third panel at the bottom, and also shares the x-axis.  The arrows are
    in the convention of pointing towards a conductor.  The xx and yy 
    components can be plotted as well, this adds two panels on the right.  
    Here the phase is left unwrapped.  Other parameters can be added as 
    subplots such as strike, skew and phase tensor ellipses.
    
    To manipulate the plot you can change any of the attributes listed below
    and call redraw_plot().  If you know more aout matplotlib and want to 
    change axes parameters, that can be done by changing the parameters in the
    axes attributes and then call update_plot(), note the plot must be open.

    
    Arguments:
    ----------
        **fn**: string
                       filename containing impedance (.edi) is the only 
                       format supported at the moment
                       
        **z_array**: np.ndarray((nf, 2, 2), dtype='complex')
                impedance tensor with length of nf -> the number of freq
                *default* is None
                
        **z_err_array**: np.ndarray((nf, 2, 2), dtype='real')
                    impedance tensor error estimates, same shape as z.
                    *default* is None
                    
        **res_array**: np.ndarray((nf, 2, 2))
                        array of resistivity values in linear scale.
                        *default* is None
                        
        **res_err_array**: np.ndarray((nf, 2, 2))
                            array of resistivity error estimates, same shape 
                            as res_array. *default* is None
                            
        **phase_array**: np.ndarray((nf, 2, 2))
                          array of phase values in degrees, same shape as 
                          res_array. *default* is None
                          
        **phase_err_array**: np.ndarray((nf, 2, 2))
                              array of phase error estimates, same shape as 
                              phase_array. *default* is None
                              
        **tipper_array**: np.ndarray((nf, 1, 2), dtype='complex')
                           array of tipper values for tx, ty. *default* is None
                           
        **tipper_err_array**: np.ndarray((nf, 1, 2))
                               array of tipper error estimates, same shape as
                               tipper_array. *default* is None
                               
        **z_object**: class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None
                      
        **tipper_object**: class mtpy.core.z.Tipper
                            object of mtpy.core.z. If this is input be sure the
                            attribute z.freq is filled.  
                            *default* is None 
                            
        **mt_object** : class mtpy.imaging.mtplottools.MTplot
                        object of mtpy.imaging.mtplottools.MTplot
                        *default* is None
                        
    Optional Key Words:
    --------------------
                  
        *fig_num*: int
                     figure number
                     *default* is 1
        
        *fig_size*: [width, height] in inches of actual figure size
                     
        *ffactor*: float
                      scaling factor for computing resistivity from 
                      impedances.
                      *Default* is 1
        
        *rotation_angle*: float
                   rotation angle of impedance tensor (deg or radians), 
                   *Note* : rotaion is clockwise positive
                   *default* is 0
        
        *plot_num*: [ 1 | 2 | 3 ]
                        * 1 for just Ex/By and Ey/Bx *default*
                        * 2 for all 4 components
                        * 3 for off diagonal plus the determinant
    
        *plot_title*: string
                         plot_title of plot
                         *default* is station name
                    
        *plot_tipper*: [ 'yri' | 'yr' | 'yi' | 'n' ]
                          Plots the tipper in a bottom pannel
                          * 'yri'  --> plots the real and imaginar parts
                          * 'yr'   --> plots just the real part
                          * 'yi'   --> plots just the imaginary part
                          
                          *Note:* the convention is to point towards a 
                          conductor.  Can change this by setting the
                          parameter arrow_direction = 1.
                          
        *plot_strike*: [ 'y' | 1 | 2 | 3 | 'n' ]
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
                               
        *plot_skew*: [ 'y' | 'n' ]
                        plots the skew angle calculated from the phase tensor
                            * 'y' --> plots skew angle
                            * 'n' --> does not plot skew angle *default*
                            
        *plot_pt*: [ 'y' | 'n' ]
                    plots the phase tensor ellipses which have the properties
                    of ellipse_
                        * 'y' --> plots phase tensor ellipses
                        * 'n' --> does not plot ellipses *default*
                          
        *fig_dpi*: int
                 dots-per-inch resolution, *default* is 300
                    
                        
        :Example: ::
            
            >>> import mtpy.imaging.mtplot as mtplot
            >>> edifile = r"/home/MT01/MT01.edi"
            >>> rp1 = mtplot.PlotResPhase(fn=edifile, plot_num=2)
            >>> # plots all 4 components
            >>> rp1 = mtplot.PlotResPhase(fn=edifile, plot_tipper='yr')
            >>> # plots the real part of the tipper
            
    Attributes:
    -----------
        -fn              filename to be plotted (only supports .edi so far) 
        -fig_num         figure number for plotting
        -fig_size        size of figure in inches [width, height]
        -plot_num        plot type, see arguments for details 
        -plot_title      title of the plot, *default* is station name
        -fig_dpi         Dots-per-inch resolution of plot, *default* is 300
        -rotation_angle           Rotate impedance tensor by this angle (deg) assuming
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
        -marker_lw       line width of marker, *default* is 100./fig_dpi
        -lw              line width of line and errorbar lines
        ..
        
         *For more on line and marker styles see matplotlib.lines.Line2D*
        
        -arrow_lw          line width of the arrow, *default* is 0.75
        -arrow_head_width  head width of the arrow, *default* is 0 for no arrow
                           head.  Haven't found a good way to scale the arrow
                           heads in a log scale.
                         
        -arrow_head_length  head width of the arrow, *default* is 0 for no arrow
                            head.  Haven't found a good way to scale the arrow
                            heads in a log scale.
                          
        -arrow_color_real  color of the real arrows, *default* is black
        -arrow_color_imag  color of the imaginary arrows, *default* is blue
        
        -arrow_direction   0 for pointing towards a conductor and -1 for 
                           pointing away from a conductor.
                          

        -x_limits        limits on the x-limits (period), *default* is None
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
                        
        -phase_quadrant [ 1 | 3 ] 
                        * 1 for both phases to be in 0 to 90, 
                        * 3 for xy to be in 0-90 and yx to be in -180 to 270  
                        
        -tipper_limits  limits of the y-axis, *default* is (-1,1)
        
        -skew_limits    limits for skew angle, *default* is (-9,9)
        
        -strike_limits  limits for strike angle, *default* is (-90,90) 
        
        
    Methods:
    --------
        * *plot*: plots the pseudosection according to keywords
        * *redraw_plot*: redraws the plot, use if you change some of the 
                         attributes.
        * *update_plot*: updates the plot, use if you change some of the 
                         axes attributes, figure needs to be open to update.
        * *save_plot*: saves the plot to given filepath.
        

    """   
    
    def __init__(self, **kwargs):
        super(PlotResponse, self).__init__()

        fn = kwargs.pop('fn', None)
        z_object = kwargs.pop('z_object', None)
        tipper_object = kwargs.pop('tipper_object', None)
        mt_object = kwargs.pop('mt_object', None)

        #--> initialize an MTplot object
        if mt_object is None:
            self.mt = mtpl.MTplot(fn=fn, 
                                   z_object=z_object,
                                   tipper_object=tipper_object)
                              
        else:
            self.mt = mt_object 
            
        self.phase_quadrant =  1
        
        #mtpl.PlotSettings.__init__(self)
        
        self.plot_num = kwargs.pop('plot_num', 1)
        self.rotation_angle = kwargs.pop('rotation_angle', 0)
        
        #set plot tipper or not
        self._plot_tipper = kwargs.pop('plot_tipper', 'n')
        
        #plot strike angle or not
        self._plot_strike = kwargs.pop('plot_strike', 'n')
        
        #plot skew angle
        self._plot_skew = kwargs.pop('plot_skew', 'n')
        
        #plot phase tensor ellipses
        self._plot_pt = kwargs.pop('plot_pt', 'n')
        
        #order of plots
        self.plot_order = kwargs.pop('plot_order', 
                                     ['tip' , 'pt', 'strike', 'skew'])
                                     
        self.plot_dict = dict([(kk, vv) for kk, vv in zip(['tip' , 'pt', 
                                                           'strike', 'skew'],
                                                           [self._plot_tipper, 
                                                            self._plot_pt,
                                                            self._plot_strike, 
                                                            self._plot_skew])])
        
        #set arrow properties
        self.arrow_head_length = 0.03
        self.arrow_head_width = 0.03
        self.arrow_lw = .5
        
        #ellipse_properties
        self.ellipse_size = 0.25
        self.ellipse_spacing = kwargs.pop('ellipse_spacing', 1)
        if self.ellipse_size == 2 and self.ellipse_spacing == 1:
            self.ellipse_size = 0.25

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        #plot on initializing
        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plotResPhase(filename,fig_num) will plot the apparent resistivity and 
        phase for a single station. 
 
        """
        
        #set figure size according to what the plot will be.
        if self.fig_size is None:
            if self.plot_num == 1 or self.plot_num == 3:
                self.fig_size = [5, 7]
            
            elif self.plot_num == 2:
                self.fig_size = [7, 7]
            
        #--> rotate the impedance tensor if desired
        if self.rotation_angle != 0:
            self.mt.rotation_angle = self.rotation_angle
        
        if self._plot_tipper.find('y') == 0:
            if np.all(self.mt.Tipper.tipper == 0+0j):
                print('No Tipper data for station {0}'.format(self.mt.station))
                self.plot_tipper = 'n'
        
        #set x-axis limits from short period to long period
        if self.x_limits == None:
            self.x_limits = (10**(np.floor(np.log10(self.mt.period[0]))),
                            10**(np.ceil(np.log10((self.mt.period[-1])))))
        if self.phase_limits == None:
            pass
            
        if self.res_limits == None:
            self.res_limits = (10**(np.floor(
                                    np.log10(min([self.mt.Z.res_xy.min(),
                                                  self.mt.Z.res_yx.min()])))),
                              10**(np.ceil(
                                    np.log10(max([self.mt.Z.res_xy.max(),
                                                  self.mt.Z.res_yx.max()])))))

        #set some parameters of the figure and subplot spacing
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .93
        plt.rcParams['figure.subplot.left'] = .80
        plt.rcParams['figure.subplot.right'] = .98
        
        #set the font properties for the axis labels
        fontdict = {'size':self.font_size+2, 'weight':'bold'}
        
        #create a dictionary for the number of subplots needed
        pdict = {'res' : 0, 
                 'phase' : 1}
        #start the index at 2 because resistivity and phase is permanent for 
        #now 
        index = 2
        for key in self.plot_order:
            if self.plot_dict[key].find('y')==0:
                pdict[key] = index
                index += 1
        
        #get number of rows needed
        nrows = index
        
        #set height ratios of the subplots
        hr = [2, 1.5]+[1]*(len(list(pdict.keys()))-2)
        
        #create a grid to place the figures into, set to have 2 rows and 2 
        #columns to put any of the 4 components.  Make the phase plot
        #slightly shorter than the apparent resistivity plot and have the two
        #close to eachother vertically.  If there is tipper add a 3rd row and
        #if there is strike add another row
        gs = gridspec.GridSpec(nrows, 2, height_ratios=hr,hspace=.05)

        #make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        
        #--> make figure for xy,yx components
        if self.plot_num == 1 or self.plot_num == 3:
            #set label coordinates
            labelcoords = (-0.075, 0.5)
            
            #space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.1)
            
            #--> create the axes instances
            #apparent resistivity axis
            self.axr = self.fig.add_subplot(gs[0, :])
                
            #phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs[1, :], sharex=self.axr)
            

        
        #--> make figure for all 4 components
        elif self.plot_num == 2:
            #set label coordinates
            labelcoords = (-0.095, 0.5)            
            
            #space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.07)
            
            #--> create the axes instances
            #apparent resistivity axis
            self.axr = self.fig.add_subplot(gs[0, 0])
            
            #phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs[1, 0], sharex=self.axr)
        
        #place y coordinate labels in the same location                
        self.axr.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        self.axp.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        
        #--> plot tipper
        try:
            self.axt = self.fig.add_subplot(gs[pdict['tip'], :], )
            self.axt.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        except KeyError:
            pass
        
        #--> plot phase tensors
        try:
            #can't share axis because not on the same scale
            self.axpt = self.fig.add_subplot(gs[pdict['pt'], :], 
                                             aspect='equal')
            self.axpt.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        except KeyError:
            pass
        
        #--> plot strike
        try:
            self.axst = self.fig.add_subplot(gs[pdict['strike'], :], 
                                             sharex=self.axr)
            self.axst.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        except KeyError:
            pass
        
        #--> plot skew
        try:
            self.axsk = self.fig.add_subplot(gs[pdict['skew'], :], 
                                             sharex=self.axr)
            self.axsk.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        except KeyError:
            pass
    
        
        #---------plot the apparent resistivity--------------------------------
        #--> plot as error bars and just as points xy, yx
        #res_xy
        self.ebxyr = mtpl.plot_errorbar(self.axr, 
                                       self.mt.period, 
                                       self.mt.Z.res_xy, 
                                       marker=self.xy_marker, 
                                       ms=self.marker_size,
                                       color = self.xy_color, 
                                       ls=self.xy_ls,
                                       lw=self.lw,
                                       y_error=self.mt.Z.res_err_xy, 
                                       e_capsize=self.marker_size)
        
        #res_yx                              
        self.ebyxr = mtpl.plot_errorbar(self.axr, 
                                       self.mt.period, 
                                       self.mt.Z.res_yx, 
                                       marker=self.yx_marker, 
                                       ms=self.marker_size, 
                                       color=self.yx_color, 
                                       ls=self.yx_ls,
                                       lw=self.lw, 
                                       y_error=self.mt.Z.res_err_yx, 
                                       e_capsize=self.marker_size)
                                      
        #--> set axes properties
        plt.setp(self.axr.get_xticklabels(), visible=False)
        self.axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                            fontdict=fontdict)
        self.axr.set_yscale('log')
        self.axr.set_xscale('log')
        self.axr.set_xlim(self.x_limits)
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
        self.ebxyp = mtpl.plot_errorbar(self.axp, 
                                       self.mt.period, 
                                       self.mt.Z.phase_xy, 
                                       marker=self.xy_marker, 
                                       ms=self.marker_size, 
                                       color=self.xy_color, 
                                       ls=self.xy_ls,
                                       lw=self.lw,
                                       y_error=self.mt.Z.phase_err_xy, 
                                       e_capsize=self.marker_size)
                                      
        #phase_yx: Note add 180 to place it in same quadrant as phase_xy
        self.ebyxp = mtpl.plot_errorbar(self.axp, 
                                       self.mt.period, 
                                       self.mt.Z.phase_yx+180*self.phase_quadrant, 
                                       marker=self.yx_marker, 
                                       ms=self.marker_size,  
                                       color=self.yx_color, 
                                       ls=self.yx_ls, 
                                       lw=self.lw,
                                       y_error=self.mt.Z.phase_err_yx, 
                                       e_capsize=self.marker_size)

        #check the phase to see if any point are outside of [0:90]
        if self.phase_limits == None:
            if min(self.mt.Z.phase_xy)<0 or min(self.mt.Z.phase_yx)<0:
                pymin = min([min(self.mt.Z.phase_xy), min(self.mt.Z.phase_yx)])
                if pymin > 0:
                   pymin = 0
            else:
                pymin = 0
            
            if max(self.mt.Z.phase_xy)>90 or max(self.mt.Z.phase_yx)>90:
                pymax = min([max(self.mt.Z.phase_xy), max(self.mt.Z.phase_yx)])
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
        #set th xaxis tick labels to invisible
        if pdict['phase'] != nrows-1:
            plt.setp(self.axp.xaxis.get_ticklabels(), visible=False)
            self.axp.set_xlabel('')
            
        #-----plot tipper----------------------------------------------------              
        if self._plot_tipper.find('y') == 0:
            
            txr = self.mt.Tipper.mag_real*np.sin(self.mt.Tipper.angle_real*np.pi/180+\
                                     np.pi*self.arrow_direction)
            tyr = self.mt.Tipper.mag_real*np.cos(self.mt.Tipper.angle_real*np.pi/180+\
                                     np.pi*self.arrow_direction)
    
            txi = self.mt.Tipper.mag_imag*np.sin(self.mt.Tipper.angle_imag*np.pi/180+\
                                     np.pi*self.arrow_direction)
            tyi = self.mt.Tipper.mag_imag*np.cos(self.mt.Tipper.angle_imag*np.pi/180+\
                                     np.pi*self.arrow_direction)
            
            nt = len(txr)
            
            tiplist = []
            tiplabel = []
            
            for aa in range(nt):
                if type(txr) == np.ma.core.MaskedArray:
                    if txr.mask[aa]:
                        continue
                xlenr = txr[aa]*np.log10(self.mt.period[aa])
                xleni = txi[aa]*np.log10(self.mt.period[aa])
                
                #--> plot real arrows
                if self._plot_tipper.find('r') > 0:
                    self.axt.arrow(np.log10(self.mt.period[aa]),
                                   0,
                                   xlenr,
                                   tyr[aa],
                                   lw=self.arrow_lw,
                                   facecolor=self.arrow_color_real,
                                   edgecolor=self.arrow_color_real,
                                   head_width=self.arrow_head_width,
                                   head_length=self.arrow_head_length,
                                   length_includes_head=False)
                    
                    if aa == 0:
                        line1 = self.axt.plot(0, 0, self.arrow_color_real)
                        tiplist.append(line1[0])
                        tiplabel.append('real')
                                   
                #--> plot imaginary arrows
                if self._plot_tipper.find('i') > 0:               
                    self.axt.arrow(np.log10(self.mt.period[aa]),
                                   0,
                                   xleni,
                                   tyi[aa],
                                   lw=self.arrow_lw,
                                   facecolor=self.arrow_color_imag,
                                   edgecolor=self.arrow_color_imag,
                                   head_width=self.arrow_head_width,
                                   head_length=self.arrow_head_length,
                                   length_includes_head=False)
                    if aa == 0:              
                        line2 = self.axt.plot(0, 0, self.arrow_color_imag)
                        tiplist.append(line2[0])
                        tiplabel.append('imag')
                
            #make a line at 0 for reference
            self.axt.plot(np.log10(self.mt.period), [0]*nt, 'k', lw=.5)
        
            
            self.axt.legend(tiplist, tiplabel,
                            loc='upper left',
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.1,
                            prop={'size':self.font_size})

            #set axis properties  
            
            self.axt.set_xlim(np.log10(self.x_limits[0]),
                               np.log10(self.x_limits[1]))
                               
            tklabels = []
            xticks = []

            for tk in self.axt.get_xticks():
                try:
                    tklabels.append(mtpl.labeldict[tk])
                    xticks.append(tk)
                except KeyError:
                    pass
            self.axt.set_xticks(xticks)
            self.axt.set_xticklabels(tklabels, 
                                      fontdict={'size':self.font_size})
            self.axt.set_xlabel('Period (s)', fontdict=fontdict)
            #need to reset the x_limits caouse they get reset when calling
            #set_ticks for some reason
            self.axt.set_xlim(np.log10(self.x_limits[0]),
                               np.log10(self.x_limits[1]))
                               
            self.axt.yaxis.set_major_locator(MultipleLocator(.2))               
            self.axt.yaxis.set_minor_locator(MultipleLocator(.1))               
            self.axt.set_xlabel('Period (s)', fontdict=fontdict)
            self.axt.set_ylabel('Tipper', fontdict=fontdict)    
            
            #self.axt.set_xscale('log')
            if self.tipper_limits is None:
                tmax = max([tyr.max(), tyi.max()])
                if tmax > 1:
                    tmax = .899
                            
                tmin = min([tyr.min(), tyi.min()])
                if tmin < -1:
                    tmin = -.899
                            
                self.tipper_limits = (tmin-.1, tmax+.1)
            
            self.axt.set_ylim(self.tipper_limits)
            self.axt.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
            
            #set th xaxis tick labels to invisible
            if pdict['tip'] != nrows-1:
                plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
        
        #------plot strike angles----------------------------------------------
        if self._plot_strike.find('y') == 0:
            
            stlist = []
            stlabel = []
            st_maxlist = []
            st_minlist = []
            
#            if self._plot_strike.find('i') > 0:
#                #strike from invariants
#                zinv = self.mt.get_Zinvariants()
#                s1 = zinv.strike
#                
#                #fold angles so go from -90 to 90
#                s1[np.where(s1>90)] -= -180
#                s1[np.where(s1<-90)] += 180
#                
#                #plot strike with error bars
#                ps1 = mtpl.plot_errorbar(self.axst,
#                                         self.mt.period, 
#                                        s1, 
#                                        marker=self.strike_inv_marker, 
#                                        ms=self.marker_size,  
#                                        color=self.strike_inv_color, 
#                                        ls='none',
#                                        lw=self.lw, 
#                                        y_error=zinv.strike_err, 
#                                        e_capsize=self.marker_size)
#                                        
#                stlist.append(ps1[0])
#                stlabel.append('Z_inv')
#                st_maxlist.append(s1.max())
#                st_minlist.append(s1.min())
                                        
            if self._plot_strike.find('p') > 0:
                
                #strike from phase tensor
                s2 = self.mt.pt.azimuth
                s2_err = self.mt.pt.azimuth_err
                #fold angles to go from -90 to 90
                s2[np.where(s2>90)] -= 180
                s2[np.where(s2<-90)] += 180
                
                #plot strike with error bars
                ps2 = mtpl.plot_errorbar(self.axst, 
                                         self.mt.period, 
                                         s2, 
                                         marker=self.strike_pt_marker, 
                                         ms=self.marker_size, 
                                         color=self.strike_pt_color, 
                                         ls='none',
                                         lw=self.lw, 
                                         y_error=s2_err, 
                                         e_capsize=self.marker_size)
                                        
                stlist.append(ps2[0])
                stlabel.append('PT')
                st_maxlist.append(s2.max())
                st_minlist.append(s2.min())
            
            if self._plot_strike.find('t') > 0:
                #strike from tipper
                s3 = self.mt.Tipper.angle_real+90
                
                #fold to go from -90 to 90
                s3[np.where(s3 > 90)] -= 180
                s3[np.where(s3 < -90)] += 180
                
                #plot strike with error bars
                ps3 = mtpl.plot_errorbar(self.axst,
                                         self.mt.period, 
                                        s3, 
                                        marker=self.strike_tip_marker, 
                                        ms=self.marker_size,  
                                        color=self.strike_tip_color, 
                                        ls='none', 
                                        lw=self.lw,
                                        y_error=None, 
                                        e_capsize=self.marker_size)
                                        
                stlist.append(ps3[0])
                stlabel.append('Tip')
                st_maxlist.append(s3.max())
                st_minlist.append(s3.min())
                
            #--> set axes properties
            if self.strike_limits is None:
                try:
                    stmin = min(st_minlist)
                except ValueError:
                    stmin = -89.99
                if stmin-3 < -90:
                    stmin -= 3
                else:
                    stmin = -89.99
                
                try:
                    stmax = min(st_maxlist)
                except ValueError:
                    stmax = 89.99    
                if stmax+3 < 90:
                    stmax += 3
                else:
                    stmax = 89.99
                self.strike_limits = (-max([abs(stmin), abs(stmax)]),
                                       max([abs(stmin), abs(stmax)]))
                                        
                
            self.axst.plot(self.axr.get_xlim(), [0, 0], color='k', lw=.5)
            
            self.axst.set_ylabel('Strike',
                                fontdict=fontdict)
            self.axst.set_xlabel('Period (s)',
                                fontdict=fontdict)
            self.axst.set_ylim(self.strike_limits)
            self.axst.yaxis.set_major_locator(MultipleLocator(30))
            self.axst.yaxis.set_minor_locator(MultipleLocator(5))
            self.axst.set_xscale('log')
            self.axst.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
            try:
                self.axst.legend(stlist, 
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
            
            #set th xaxis tick labels to invisible
            if pdict['strike'] != nrows-1:
                plt.setp(self.axst.xaxis.get_ticklabels(), visible=False)
            
        #------plot skew angle---------------------------------------------
        if self._plot_skew == 'y':
            #strike from phase tensor
            sk = self.mt.pt.beta
            sk_err = self.mt.pt.beta_err
            
            ps4 = mtpl.plot_errorbar(self.axsk, 
                                     self.mt.period, 
                                    sk, 
                                    marker=self.skew_marker, 
                                    ms=self.marker_size,  
                                    color=self.skew_color, 
                                    ls='none',
                                    lw=self.lw, 
                                    y_error=sk_err, 
                                    e_capsize=self.marker_size)
                                    
            self.axsk.plot(self.mt.period, [3]*len(self.mt.period), '-k',
                           lw = self.lw)
            self.axsk.plot(self.mt.period, [-3]*len(self.mt.period), '-k',
                           lw = self.lw)
                                    
            if self.skew_limits is None:
                self.skew_limits = (-9, 9)
            
            self.axsk.set_ylim(self.skew_limits)
            self.axsk.yaxis.set_major_locator(MultipleLocator(3))
            self.axsk.yaxis.set_minor_locator(MultipleLocator(1))
            self.axsk.grid(True, alpha=.25, color=(.25, .25, .25))
            self.axsk.set_ylabel('Skew', fontdict)
            self.axsk.set_xlabel('Period (s)', fontdict)
            self.axsk.set_xscale('log')
           
            #set th xaxis tick labels to invisible
            if pdict['skew'] != nrows-1:
                plt.setp(self.axst.xaxis.get_ticklabels(), visible=False)
        
        #----plot phase tensor ellipse---------------------------------------    
        if self._plot_pt == 'y':        
            
            cmap = self.ellipse_cmap
            ckmin = self.ellipse_range[0]
            ckmax = self.ellipse_range[1]
            try:
                ckstep = float(self.ellipse_range[2])
            except IndexError:
                ckstep = 3
                
            if cmap == 'mt_seg_bl2wh2rd':
                bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
                nseg = float((ckmax-ckmin)/(2*ckstep))
    
            #get the properties to color the ellipses by
            if self.ellipse_colorby == 'phiminang' or \
               self.ellipse_colorby == 'phimin':
                colorarray = self.mt.pt.phimin
                
            elif self.ellipse_colorby == 'phimaxang' or \
               self.ellipse_colorby == 'phimax':
                colorarray = self.mt.pt.phimax
        
                                               
            elif self.ellipse_colorby == 'phidet':
                colorarray = np.sqrt(abs(self.mt.pt.det))*(180/np.pi)
                 
                
            elif self.ellipse_colorby == 'skew' or\
                 self.ellipse_colorby == 'skew_seg':
                colorarray = self.mt.pt.beta
                
            elif self.ellipse_colorby == 'ellipticity':
                colorarray = self.mt.pt.ellipticity
                
            else:
                raise NameError(self.ellipse_colorby+' is not supported')
         
            #-------------plot ellipses-----------------------------------
            for ii, ff in enumerate(self.mt.period):
                #make sure the ellipses will be visable
                eheight = self.mt.pt.phimin[ii]/self.mt.pt.phimax[ii]*\
                                                            self.ellipse_size
                ewidth = self.mt.pt.phimax[ii]/self.mt.pt.phimax[ii]*\
                                                            self.ellipse_size
            
                #create an ellipse scaled by phimin and phimax and oriented 
                #along the azimuth which is calculated as clockwise but needs 
                #to be plotted counter-clockwise hence the negative sign.
                ellipd = patches.Ellipse((np.log10(ff)*self.ellipse_spacing,
                                          0),
                                         width=ewidth,
                                         height=eheight,
                                         angle=90-self.mt.pt.azimuth[ii])
                                         
                self.axpt.add_patch(ellipd)
                
            
                #get ellipse color
                if cmap.find('seg') > 0:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[ii],
                                                         self.ellipse_colorby,
                                                         cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
                else:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[ii],
                                                         self.ellipse_colorby,
                                                         cmap,
                                                         ckmin,
                                                         ckmax))
                
        
            #----set axes properties-----------------------------------------------
            #--> set tick labels and limits
            self.axpt.set_xlim(np.log10(self.x_limits[0]),
                               np.log10(self.x_limits[1]))
            
            tklabels = []
            xticks = []
            for tk in self.axpt.get_xticks():
                try:
                    tklabels.append(mtpl.labeldict[tk])
                    xticks.append(tk)
                except KeyError:
                    pass
            self.axpt.set_xticks(xticks)
            self.axpt.set_xticklabels(tklabels, 
                                      fontdict={'size':self.font_size})
            self.axpt.set_xlabel('Period (s)', fontdict=fontdict)
            self.axpt.set_ylim(ymin=-1.5*self.ellipse_size, 
                               ymax=1.5*self.ellipse_size)
            #need to reset the x_limits caouse they get reset when calling
            #set_ticks for some reason
            self.axpt.set_xlim(np.log10(self.x_limits[0]),
                               np.log10(self.x_limits[1]))
            self.axpt.grid(True, 
                         alpha=.25, 
                         which='major', 
                         color=(.25,.25,.25),
                         lw=.25)
            
            plt.setp(self.axpt.get_yticklabels(), visible=False)
            if pdict['pt'] != nrows-1:
                plt.setp(self.axpt.get_xticklabels(), visible=False)
                
            #add colorbar for PT
            axpos = self.axpt.get_position()
            cb_position = (axpos.bounds[0]-.0575,
                           axpos.bounds[1]+.02,
                           .01,
                           axpos.bounds[3]*.75)
            self.cbax = self.fig.add_axes(cb_position)
            if cmap == 'mt_seg_bl2wh2rd':
                #make a color list
                clist = [(cc, cc, 1) 
                        for cc in np.arange(0,1+1./(nseg),1./(nseg))]+\
                       [(1, cc, cc) 
                        for cc in np.arange(1,-1./(nseg),-1./(nseg))]
                
                #make segmented colormap
                mt_seg_bl2wh2rd = colors.ListedColormap(clist)
    
                #make bounds so that the middle is white
                bounds = np.arange(ckmin-ckstep, ckmax+2*ckstep, ckstep)
                
                #normalize the colors
                norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)
                
                #make the colorbar
                self.cbpt = mcb.ColorbarBase(self.cbax,
                                           cmap=mt_seg_bl2wh2rd,
                                           norm=norms,
                                           orientation='vertical',
                                           ticks=bounds[1:-1])
            else:
                self.cbpt = mcb.ColorbarBase(self.cbax,
                                           cmap=mtcl.cmapdict[cmap],
                                           norm=colors.Normalize(vmin=ckmin,
                                                                 vmax=ckmax),
                                            orientation='vertical')
            self.cbpt.set_ticks([ckmin, (ckmax-ckmin)/2, ckmax])
            self.cbpt.set_ticklabels(['{0:.0f}'.format(ckmin),
                                      '{0:.0f}'.format((ckmax-ckmin)/2),
                                      '{0:.0f}'.format(ckmax)])
            self.cbpt.ax.yaxis.set_label_position('left')
            self.cbpt.ax.yaxis.set_label_coords(-1.05, .5)
            self.cbpt.ax.yaxis.tick_right()
            self.cbpt.ax.tick_params(axis='y', direction='in')
            self.cbpt.set_label(mtpl.ckdict[self.ellipse_colorby], 
                                fontdict={'size':self.font_size})
                        
            
        #===Plot the xx, yy components if desired==============================
        if self.plot_num == 2:
            #---------plot the apparent resistivity----------------------------
            self.axr2 = self.fig.add_subplot(gs[0, 1], sharex=self.axr)
            self.axr2.yaxis.set_label_coords(-.1, 0.5)
            
            #res_xx
            self.ebxxr = mtpl.plot_errorbar(self.axr2, 
                                            self.mt.period, 
                                            self.mt.Z.res_xx, 
                                            marker=self.xy_marker, 
                                            ms=self.marker_size, 
                                            color=self.xy_color, 
                                            ls=self.xy_ls,
                                            lw=self.lw, 
                                            y_error=self.mt.Z.res_err_xx, 
                                            e_capsize=self.marker_size)
            
            #res_yy                              
            self.ebyyr = mtpl.plot_errorbar(self.axr2,
                                            self.mt.period, 
                                            self.mt.Z.res_yy, 
                                            marker=self.yx_marker, 
                                            ms=self.marker_size,  
                                            color=self.yx_color,  
                                            ls=self.yx_ls,
                                            lw=self.lw, 
                                            y_error=self.mt.Z.res_err_yy, 
                                            e_capsize=self.marker_size)

            #--> set axes properties
            plt.setp(self.axr2.get_xticklabels(), visible=False)
            self.axr2.set_yscale('log')
            self.axr2.set_xscale('log')
            self.axr2.set_xlim(self.x_limits)
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
            self.ebxxp = mtpl.plot_errorbar(self.axp2, 
                                            self.mt.period, 
                                            self.mt.Z.phase_xx, 
                                            marker=self.xy_marker, 
                                            ms=self.marker_size,
                                            color=self.xy_color,  
                                            ls=self.xy_ls,
                                            lw=self.lw, 
                                            y_error=self.mt.Z.phase_err_xx,
                                            e_capsize=self.marker_size)
                                            
            #phase_yy
            self.ebyyp = mtpl.plot_errorbar(self.axp2, 
                                            self.mt.period,
                                            self.mt.Z.phase_yy, 
                                            marker=self.yx_marker,
                                            ms=self.marker_size, 
                                            color=self.yx_color, 
                                            ls=self.yx_ls,
                                            lw=self.lw, 
                                            y_error=self.mt.Z.phase_err_yy, 
                                            e_capsize=self.marker_size)
            
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
                           
            if len(list(pdict.keys()))>2:
                plt.setp(self.axp2.xaxis.get_ticklabels(), visible=False)
                plt.setp(self.axp2.xaxis.get_label(), visible=False)
        
        #===Plot the Determinant if desired==================================                          
        if self.plot_num == 3:
                
            #res_det
            self.ebdetr = mtpl.plot_errorbar(self.axr, 
                                            self.mt.period, 
                                            self.mt.Z.res_det, 
                                            marker=self.det_marker, 
                                            ms=self.marker_size, 
                                            color=self.det_color,  
                                            ls=self.det_ls, 
                                            lw=self.lw,
                                            y_error=self.mt.Z.res_det_err, 
                                            e_capsize=self.marker_size)
        
            #phase_det
            self.ebdetp = mtpl.plot_errorbar(self.axp, 
                                            self.mt.period, 
                                            self.mt.Z.phase_det, 
                                            marker=self.det_marker, 
                                            ms=self.marker_size,  
                                            color=self.det_color,  
                                            ls=self.det_ls,
                                            lw=self.lw, 
                                            y_error=self.mt.Z.phase_det_err, 
                                            e_capsize=self.marker_size)
                                            
            self.axr.legend((self.ebxyr[0], self.ebyxr[0], self.ebdetr[0]),
                            ('$Z_{xy}$','$Z_{yx}$','$\det(\mathbf{\hat{Z}})$'),
                            loc=3,
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.02)
        
        
        #make plot_title and show
        if self.plot_title is None:
            self.plot_title = self.mt.station
        
        self.fig.suptitle(self.plot_title, fontdict=fontdict)
        
        # be sure to show
        plt.show()
        
    def _set_plot_tipper(self, plot_tipper):
        """
        If plotting tipper make arrow attributes

        """
        
        self._plot_tipper = plot_tipper
        
        self.plot_dict['tip'] = self._plot_tipper
        
    def _get_plot_tipper(self):
        self._plot_tipper
        
    plot_tipper = property(fget=_get_plot_tipper, fset=_set_plot_tipper, 
                           doc="""string to plot tipper""")
                           
    def _set_plot_pt(self, plot_pt):
        """
        If plotting tipper make arrow attributes

        """
        
        self._plot_pt = plot_pt
        
        self.plot_dict['pt'] = self._plot_pt
        
    def _get_plot_pt(self):
        self._plot_pt
        
    plot_pt = property(fget=_get_plot_pt, fset=_set_plot_pt, 
                       doc="""string to plot phase tensor ellipses""")
                       
    def _set_plot_strike(self, plot_strike):
        """
        change plot_dict when changing plot_strike

        """
        
        self._plot_strike = plot_strike
        
        self.plot_dict['strike'] = self._plot_strike
    
    def _get_plot_strike(self):    
        return self._plot_strike
        
    plot_strike = property(fget=_get_plot_strike, fset=_set_plot_strike, 
                           doc="""string to plot strike""")
    def _set_plot_skew(self, plot_skew):
        """
        change plot_dict when changing plot_strike

        """
        
        self._plot_skew = plot_skew
        
        self.plot_dict['skew'] = self._plot_skew
    
    def _get_plot_skew(self):
        return self._plot_skew
        
    plot_skew = property(fget=_get_plot_skew, fset=_set_plot_skew, 
                           doc="""string to plot skew""")
        

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
                              
            **fig_fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the fig_dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than fig_dpi of the figure.
                          
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
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, fig_dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            plt.clf()
            plt.close(self.fig)
            
        else:
            save_fn = os.path.join(save_fn, self.mt.station+'_ResPhase.'+
                                    file_format)
            self.fig.savefig(save_fn, fig_dpi=fig_dpi, format=file_format,
                        orientation=orientation)
        
        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)
        
        else:
            pass
        
        self.fig_fn = save_fn
        print('Saved figure to: '+self.fig_fn)

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
