# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:02:39 2013

@author: jpeacock-pr
"""

#============================================================================

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import mtpy.imaging.mtplottools as mtpl
from mtpy.imaging.plotresponse import PlotResponse as plotresponse

#============================================================================

class PlotMultipleResponses(object):
    """
    plots multiple MT responses simultaneously either in single plots or in 
    one plot of subfigures or in a single plot with subfigures for each 
    component.
    
    expecting only one type of input --> can be:
        **fn_lst** : list of filenames to plot
           
         **z_object_lst** : list of mtpy.core.z.Z objects
         
         **res_object_lst** : list of mtpy.imaging.mtplot.ResPhase objects
         
         **tipper_object_lst** : list of mtpy.imaging.mtplot.Tipper objects
         
         **mt_object_lst** : list of mtpy.imaging.mtplot.MTplot objects
         
        
    Arguments:
    ----------
        **fn_lst** : list of filenames to plot
                     ie. [fn_1, fn_2, ...], *default* is None
           
         **z_object_lst** : list of mtpy.core.z.Z objects
                            *default* is None
         
         **res_object_lst** : list of mtpy.imaging.mtplot.ResPhase objects
                              *default* is None
         
         **tipper_object_lst** : list of mtpy.imaging.mtplot.Tipper objects
                                 *default* is None
         
         **mt_object_lst** : list of mtpy.imaging.mtplot.MTplot objects
                             *default* is None
                      
        **fignum** : int
                     figure number
                     *default* is 1
        
        **rot_z** : float or np.ndarray
                   rotation angle of impedance tensor (deg or radians), 
                   *Note* : rotaion is clockwise positive
                   *default* is 0
                   Can input so each station is rotated at a constant angle or
                   each period is rotated differently, or both.
        
        **plotnum** : [ 1 | 2 | 3 ]
                        * 1 for just Ex/By and Ey/Bx *default*
                        * 2 for all 4 components
                        * 3 for off diagonal plus the determinant
                        
        **plot_style** : [ '1' | 'all' | 'compare' ]
                        determines the plotting style:
                            * '1' for plotting each station in a different 
                                  figure. *default*
                                  
                            * 'all' for plotting each station in a subplot
                                    all in the same figure
                            
                            * 'compare' for comparing the responses all in 
                                        one plot.  Here the responses are 
                                        colored from dark to light.  This 
                                        plot can get messy if too many stations
                                        are plotted.  
                                    
    
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
                       string for plotting skew angle.  This is plotted in 
                       the same plot as strike angle at the moment.  
                           * 'y' for plotting the skew
                           * 'n' for not plotting skew *default*
                          
        **dpi** : int
                 dots-per-inch resolution, *default* is 300
                    
                        
        :Example: ::
            
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> import os
            >>> edipath = r"/home/Edifiles"
            >>> edilst = [os.path.join(edipath,edi) 
            >>> ...       for edi in os.listdir(edipath)
            >>> ...       if edi.find('.edi')>0]
            >>> plot each station in a subplot all in one figure with tipper
            >>> rp1 = mtplot.PlotMultipleResPhase(edilst, plotnum=1, 
            >>> ...                                plot_tipper='yr,
            >>> ...                                plot_style='all')

            
    Attributes:
    -----------
        -mt_lst         list of mtplot.MTplot objects made from inputs 
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
    
    """

    def __init__(self, fn_lst=None, res_object_lst=None,
                 z_object_lst=None, tipper_object_lst=None, mt_object_lst=None,
                 fignum=1, dpi=300, rot_z=0, plot_num=1, plot_style='1', 
                 plot_yn='y', plot_tipper='n',plot_strike='n', plot_skew='n',
                 plotnum=1,title=None):
        """
        Initialize parameters
        """
        
        #--> get the inputs into a list of mt objects
        self.mt_lst = mtpl.get_mtlst(fn_lst=fn_lst, 
                                     res_object_lst=res_object_lst,
                                     z_object_lst=z_object_lst, 
                                     tipper_object_lst=tipper_object_lst, 
                                     mt_object_lst=mt_object_lst)
        
        #set some of the properties as attributes much to Lars' discontent
        self.fignum = fignum
        self.plotnum = plotnum
        self.plot_style = plot_style
        self.title = title
        self.dpi = dpi
        
        #if rotation angle is an int or float make an array the length of 
        #mt_lst for plotting purposes
        if type(rot_z) is float or type(rot_z) is int:
            self.rot_z = np.array([rot_z]*len(self.mt_lst))
        
        #if the rotation angle is an array for rotation of different 
        #freq than repeat that rotation array to the len(mt_lst)
        elif type(rot_z) is np.ndarray:
            if rot_z.shape[0]  !=  len(self.mt_lst):
                self.rot_z = np.repeat(rot_z, len(self.mt_lst))
                
        else:
            self.rot_z = rot_z

        
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
        self.font_size = 6
        
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
        if plot_yn == 'y':
            self.plot()
            
    #---rotate data on setting rot_z
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
            if rot_z.shape[0]  !=  len(self.mt_lst):
                rot_z = np.repeat(rot_z, len(self.mt_lst))
                
        else:
            pass
            
        for ii, mt in enumerate(self.mt_lst):
            mt.rot_z = rot_z[ii]
            
    rot_z = property(fset=_set_rot_z, doc="rotation angle(s)")
                            
    #---plot the resistivity and phase
    def plot(self):
        """
        plot the apparent resistivity and phase
        """
        
        if self.plot_style == '1':
            self.plotlst = []
            
            #--> plot from edi's if given
            for ii, mt in enumerate(self.mt_lst, 1):
                p1 = plotresponse(mt_object=mt, 
                                  fignum=ii, 
                                  plotnum=self.plotnum, 
                                  dpi=self.dpi, 
                                  rotz=self.rot_z[ii-1], 
                                  plot_yn='n',
                                  plot_tipper=self.plot_tipper,
                                  plot_strike=self.plot_strike,
                                  plot_skew=self.plot_skew)
                
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
                
                #set arrow properties
                p1.arrow_lw = self.arrow_lw
                p1.arrow_head_width = self.arrow_head_width 
                p1.arrow_head_height = self.arrow_head_height
                p1.arrow_color_real = self.arrow_color_real 
                p1.arrow_color_imag = self.arrow_color_imag 
                p1.arrow_direction = self.arrow_direction
                p1.tipper_limits = self.tipper_limits 
                
                #skew properties
                p1.skew_color = self.skew_color
                p1.skew_marker = self.skew_marker
                
                #strike properties
                p1.strike_inv_marker = self.strike_inv_marker
                p1.strike_inv_color = self.strike_inv_color
                
                p1.strike_pt_marker = self.strike_pt_marker
                p1.strike_pt_color = self.strike_pt_color
                
                p1.strike_tip_marker = self.strike_tip_marker
                p1.strike_tip_color = self.strike_tip_color
                
                #--> plot the apparent resistivity and phase
                self.plotlst.append(p1)
                
                p1.plot()
                    
        
        #-----Plot All in one figure with each plot as a subfigure------------        
        if self.plot_style == 'all':

            ns = len(self.mt_lst)

            #set some parameters of the figure and subplot spacing
            plt.rcParams['font.size'] = self.font_size
            if self.plot_skew == 'y':
                plt.rcParams['figure.subplot.right'] = .94
            else:
                plt.rcParams['figure.subplot.right'] = .98
            plt.rcParams['figure.subplot.bottom'] = .1
            plt.rcParams['figure.subplot.top'] = .93
            
            #set the font properties for the axis labels
            fontdict = {'size':self.font_size+2, 'weight':'bold'}
            
            #set figure size according to what the plot will be.
            if self.plotnum == 1 or self.plotnum == 3:
                self.figsize = [ns*4, 6]
                
            elif self.plotnum == 2:
                self.figsize = [ns*8, 6]
                
            #make a figure instance
            self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
                
            #make subplots as columns for all stations that need to be plotted
            gs0 = gridspec.GridSpec(1, ns)
                
            #space out the subplots
            gs0.update(hspace=.025, wspace=.025, left=.085)
     
            for ii, mt in enumerate(self.mt_lst):
                #get the reistivity and phase object
                rp = mt.get_ResPhase()
                
                #set x-axis limits from short period to long period
                if self.xlimits == None:
                    self.xlimits = (10**(np.floor(np.log10(mt.period[0]))),
                                    10**(np.ceil(np.log10((mt.period[-1])))))
                if self.phase_limits == None:
                    pass
                    
                if self.res_limits == None:
                    self.res_limits = (10**(np.floor(
                                            np.log10(min([rp.resxy.min(),
                                                          rp.resyx.min()])))),
                                      10**(np.ceil(
                                          np.log10(max([rp.resxy.max(),
                                                        rp.resyx.max()])))))

                # create a grid to place the figures into, set to have 2 rows 
                # and 2 columns to put any of the 4 components.  Make the phase
                # plot slightly shorter than the apparent resistivity plot and 
                # have the two close to eachother vertically.
                
                #--> make figure for xy,yx components
                if self.plotnum == 1 or self.plotnum == 3:
                    #--> plot tipper
                    if self.plot_tipper.find('y') == 0:
                        #--> plot strike and or skew
                        if self.plot_strike == 'y' or self.plot_skew == 'y':
                            #make subplots for each subplot in the figure
                            ax1 = gridspec.GridSpecFromSubplotSpec(4, 1, 
                                                       subplot_spec=gs0[ii],
                                                       height_ratios=[2,
                                                                      1.5,
                                                                      1,
                                                                      1], 
                                                       hspace=.00)
                            axr = self.fig.add_subplot(ax1[0, 0])
                            axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                            axt = self.fig.add_subplot(ax1[2, 0], sharex=axr)
                            axs = self.fig.add_subplot(ax1[3, 0], sharex=axr)
                       
                       #--> don't plot strike or skew    
                        else:
                            #make subplots for each subplot in the figure
                            ax1 = gridspec.GridSpecFromSubplotSpec(3, 1, 
                                                       subplot_spec=gs0[ii],
                                                       height_ratios=[2,
                                                                      1.5,
                                                                      1], 
                                                       hspace=.00)
                            axr = self.fig.add_subplot(ax1[0, 0])
                            axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                            axt = self.fig.add_subplot(ax1[2, 0], sharex=axr)
                    
                    #--> don't plot tipper but plot skew and or strike
                    elif self.plot_strike == 'y' or self.plot_skew == 'y':
                        #make subplots for each subplot in the figure
                        ax1 = gridspec.GridSpecFromSubplotSpec(3, 1, 
                                                         subplot_spec=gs0[ii],
                                                         height_ratios=[2,
                                                                        1.5,
                                                                        1], 
                                                         hspace=.00)
                        axr = self.fig.add_subplot(ax1[0, 0])
                        axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                        axs = self.fig.add_subplot(ax1[2, 0], sharex=axr)
                   
                    #--> just plot resistivity and phase 
                    else:
                        #make subplots for each subplot in the figure
                        ax1 = gridspec.GridSpecFromSubplotSpec(2, 1, 
                                                         subplot_spec=gs0[ii],
                                                         height_ratios=[2,
                                                                        1.5], 
                                                         hspace=.00)
                        axr = self.fig.add_subplot(ax1[0, 0])
                        axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                        
                    #place the y-coordinate labels in the same location for
                    #each axis
                    if ii == 0:
                        axr.yaxis.set_label_coords(-.145, 0.5)
                        axp.yaxis.set_label_coords(-.145, 0.5)
                        if self.plot_tipper.find('y') == 0:
                            axt.yaxis.set_label_coords(-.145, 0.5)
                        if self.plot_skew == 'y' or self.plot_strike == 'y':
                            axs.yaxis.set_label_coords(-.145, 0.5)
                    
                #--> make figure for all 4 components
                elif self.plotnum == 2:
                    #--> plot tipper
                    if self.plot_tipper.find('y') == 0:
                        #--> plot strike and or skew
                        if self.plot_strike == 'y' or self.plot_skew == 'y':
                            #make subplots for each subplot in the figure
                            ax1 = gridspec.GridSpecFromSubplotSpec(4, 2, 
                                                       subplot_spec=gs0[ii],
                                                       height_ratios=[2,
                                                                      1.5,
                                                                      1,
                                                                      1], 
                                                       hspace=.00)
                            axr = self.fig.add_subplot(ax1[0, 0])
                            axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                            axr2 = self.fig.add_subplot(ax1[0, 1], sharex=axr)
                            axp2 = self.fig.add_subplot(ax1[1, 1], sharex=axr)
                            axt = self.fig.add_subplot(ax1[2, :], sharex=axr)
                            axs = self.fig.add_subplot(ax1[3, :], sharex=axr)
                       
                       #--> don't plot strike or skew    
                        else:
                            #make subplots for each subplot in the figure
                            ax1 = gridspec.GridSpecFromSubplotSpec(3, 2, 
                                                       subplot_spec=gs0[ii],
                                                       height_ratios=[2,
                                                                      1.5,
                                                                      1], 
                                                       hspace=.00)
                            axr = self.fig.add_subplot(ax1[0, 0])
                            axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                            axr2 = self.fig.add_subplot(ax1[0, 1], sharex=axr)
                            axp2 = self.fig.add_subplot(ax1[1, 1], sharex=axr)
                            axt = self.fig.add_subplot(ax1[2, :], sharex=axr)
                    
                    #--> don't plot tipper but plot skew and or strike
                    elif self.plot_strike == 'y' or self.plot_skew == 'y':
                        #make subplots for each subplot in the figure
                        ax1 = gridspec.GridSpecFromSubplotSpec(3, 2, 
                                                         subplot_spec=gs0[ii],
                                                         height_ratios=[2,
                                                                        1.5,
                                                                        1], 
                                                         hspace=.00)
                        axr = self.fig.add_subplot(ax1[0, 0])
                        axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                        axr2 = self.fig.add_subplot(ax1[0, 1], sharex=axr)
                        axp2 = self.fig.add_subplot(ax1[1, 1], sharex=axr)
                        axs = self.fig.add_subplot(ax1[2, :], sharex=axr)
                   
                    #--> just plot resistivity and phase 
                    else:
                        #make subplots for each subplot in the figure
                        ax1 = gridspec.GridSpecFromSubplotSpec(2, 2, 
                                                         subplot_spec=gs0[ii],
                                                         height_ratios=[2,
                                                                        1.5], 
                                                         hspace=.00)
                    
                        axr = self.fig.add_subplot(ax1[0, 0])
                        axp = self.fig.add_subplot(ax1[1, 0], sharex=axr)
                        axr2 = self.fig.add_subplot(ax1[0, 1], sharex=axr)
                        axp2 = self.fig.add_subplot(ax1[1, 1], sharex=axr)
                    #place the y-coordinate labels in the same location for
                    #each axis
                    if ii == 0:
                        axr.yaxis.set_label_coords(-.145, 0.5)
                        axp.yaxis.set_label_coords(-.145, 0.5)
                        if self.plot_tipper.find('y') == 0:
                            axt.yaxis.set_label_coords(-.145, 0.5)
                        if self.plot_skew == 'y' or self.plot_strike == 'y':
                            axs.yaxis.set_label_coords(-.145, 0.5)
                 
                #---------plot the apparent resistivity----------------------
                #--> plot as error bars and just as points xy-blue, yx-red
                #res_xy
                ebxyr = axr.errorbar(mt.period, 
                                     rp.resxy, 
                                     marker=self.xy_marker, 
                                     ms=self.marker_size, 
                                     mfc=self.xy_mfc, 
                                     mec=self.xy_color, 
                                     mew=self.marker_lw, 
                                     ls=self.xy_ls, 
                                     yerr=rp.resxy_err, 
                                     ecolor=self.xy_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)
                
                #res_yx                              
                ebyxr = axr.errorbar(mt.period, 
                                     rp.resyx, 
                                     marker=self.yx_marker, 
                                     ms=self.marker_size, 
                                     mfc=self.yx_mfc,
                                     mec=self.yx_color, 
                                     mew=self.marker_lw,
                                     ls=self.yx_ls, 
                                     yerr=rp.resyx_err, 
                                     ecolor=self.yx_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)
                                              
                #--> set axes properties
                plt.setp(axr.get_xticklabels(), visible=False)
                axr.set_yscale('log')
                axr.set_xscale('log')
                axr.set_xlim(self.xlimits)
                axr.set_ylim(self.res_limits)
                axr.grid(True, alpha=.25, which='both', 
                              color=(.25,.25,.25),
                              lw=.25)
                if ii == 0:
                    axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                    fontdict=fontdict)
                    axr.legend((ebxyr[0], ebyxr[0]), 
                                ('$Z_{xy}$', '$Z_{yx}$'),
                                loc=3, 
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)
                else:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    
                    
                #-----Plot the phase----------------------------------------
                #phase_xy
                ebxyp = axp.errorbar(mt.period, 
                                     rp.phasexy, 
                                     marker=self.xy_marker, 
                                     ms=self.marker_size, 
                                     mfc=self.xy_mfc,
                                     mec=self.xy_color, 
                                     mew=self.marker_lw,
                                     ls=self.xy_ls,
                                     yerr=rp.phasexy_err, 
                                     ecolor=self.xy_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)
                                              
                #phase_yx:
                ebyxp = axp.errorbar(mt.period, 
                                     rp.phaseyx, 
                                     marker=self.yx_marker, 
                                     ms=self.marker_size, 
                                     mfc=self.yx_mfc, 
                                     mec=self.yx_color, 
                                     mew=self.marker_lw,
                                     ls=self.yx_ls, 
                                     yerr=rp.phaseyx_err, 
                                     ecolor=self.yx_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)
        
                #check the phase to see if any point are outside of [0:90]
                if self.phase_limits == None:
                    if min(rp.phasexy)<0 or min(rp.phaseyx)<0:
                        pymin = min([min(rp.phasexy), 
                                     min(rp.phaseyx)])
                        if pymin > 0:
                            pymin = 0
                    else:
                        pymin = 0
                    
                    if max(rp.phasexy) > 90 or max(rp.phaseyx)  >90:
                        pymax = min([max(rp.phasexy), 
                                     max(rp.phaseyx)])
                        if pymax < 91:
                            pymax = 89.9
                    else:
                        pymax = 89.9
                        
                    self.phase_limits = (pymin, pymax)
                
                #--> set axes properties
                if ii == 0:
                    axp.set_ylabel('Phase (deg)', fontdict)
                else:
                    plt.setp(axp.get_yticklabels(), visible=False)
                    
                if self.plot_tipper == 'n' and self.plot_skew == 'n' and \
                        self.plot_strike == 'n':
                    axp.set_xlabel('Period (s)', fontdict)
                    
                axp.set_xscale('log')
                axp.set_ylim(self.phase_limits)        
                axp.yaxis.set_major_locator(MultipleLocator(15))
                axp.yaxis.set_minor_locator(MultipleLocator(5))
                axp.grid(True, alpha=.25, which='both', 
                              color=(.25,.25,.25),
                              lw=.25)
                                  
                tklabels = [self.label_dict[tt] 
                            for tt in np.arange(np.log10(self.xlimits[0]),
                                              np.log10(self.xlimits[1])+1)]
                tklabels[0] = ''
                tklabels[-1] = ''
                
                axp.set_xticklabels(tklabels,
                                    fontdict={'size':self.font_size})
                                    
                #-----plot tipper--------------------------------------------              
                if self.plot_tipper.find('y') == 0:
                    plt.setp(axp.xaxis.get_ticklabels(), visible=False)
                    
                    tp = mt.get_Tipper()
                    
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
                        xlenr = txr[aa]*mt.period[aa]
                        xleni = txi[aa]*mt.period[aa]
                        
                        #scale the arrow head height and width to fit in a
                        #log scale
                        if np.log10(mt.period[aa])<0:
                            hwidth = self.arrow_head_width*\
                                     10**(np.floor(np.log10(mt.period[aa])))
                            hheight = self.arrow_head_height*\
                                     10**(np.floor(np.log10(mt.period[aa])))
                        else:
                            hwidth = self.arrow_head_width/\
                                    10**(np.floor(np.log10(mt.period[aa])))
                            hheight = self.arrow_head_height/\
                                    10**(np.floor(np.log10(mt.period[aa]))) 
                        if np.log10(mt.period[aa])<0:
                            alw = self.arrow_lw*mt.period[aa]
                        else:
                            alw = self.arrow_lw
                        #--> plot real arrows
                        if self.plot_tipper.find('r')>0:
                            axt.arrow(mt.period[aa],
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
                                line1 = axt.plot(0, 0, self.arrow_color_real)
                                tiplst.append(line1[0])
                                tiplabel.append('real')
                                           
                        #--> plot imaginary arrows
                        if self.plot_tipper.find('i')>0:               
                            axt.arrow(mt.period[aa],
                                      0,
                                      xleni,
                                      tyi[aa],
                                      lw=alw,
                                      facecolor=self.arrow_color_imag,
                                      edgecolor=self.arrow_color_imag,
                                      length_includes_head=False)
                            if aa == 0:              
                                line2 = axt.plot(0, 0, self.arrow_color_imag)
                                tiplst.append(line2[0])
                                tiplabel.append('imag')
                        
                    #make a line at 0 for reference
                    axt.plot(mt.period, [0]*nt, 'k', lw=.5)
                
                  
                    if ii == 0:
                        axt.legend(tiplst, tiplabel,
                                    loc='upper left',
                                    markerscale=1,
                                    borderaxespad=.01,
                                    labelspacing=.07,
                                    handletextpad=.2,
                                    borderpad=.1,
                                    prop={'size':self.font_size})
                        
                        axt.set_ylabel('Tipper', fontdict=fontdict) 
                    else:
                        plt.setp(axt.get_yticklabels(), visible=False)
        
                    #set axis properties            
                    axt.yaxis.set_major_locator(MultipleLocator(.2))               
                    axt.yaxis.set_minor_locator(MultipleLocator(.1))
                    axt.set_xlabel('Period (s)', fontdict=fontdict)
   
                    
                    axt.set_xscale('log')
                    if self.tipper_limits is None:
                        tmax = max([np.sqrt(txr.max()**2+tyr.max()**2),
                                    np.sqrt(txi.max()**2+tyi.max()**2)])
                        if tmax > 1:
                            tmax = .999
                                    
                        tmin = -min([np.sqrt(txr.min()**2+tyr.min()**2),
                                    np.sqrt(txi.min()**2+tyi.min()**2)])
                        if tmin < -1:
                            tmin = -.999
                                    
                        self.tipper_limits = (tmin-.1, tmax+.1)
                    
                    axt.set_ylim(self.tipper_limits)
                    axt.grid(True, alpha=.25, which='both', 
                             color=(.25,.25,.25),
                             lw=.25)
                             
                    tklabels = [self.label_dict[tt] 
                                for tt in np.arange(np.log10(self.xlimits[0]),
                                              np.log10(self.xlimits[1])+1)]
                    tklabels[0] = ''
                    tklabels[-1] = ''
                
                    axt.set_xticklabels(tklabels,
                                        fontdict={'size':self.font_size})
                    
                #------plot strike angles---------------------------------
                if self.plot_strike != 'n' or self.plot_skew == 'y':
                    try:
                        plt.setp(axp.xaxis.get_ticklabels(), visible=False)
                    except UnboundLocalError:
                        pass
                    try:
                        plt.setp(axt.xaxis.get_ticklabels(), visible=False)
                    except UnboundLocalError:
                        pass
                    
                    stlst = []
                    stlabel = []
                    st_maxlst = []
                    st_minlst = []
                    
                    if self.plot_strike == 'y' or self.plot_strike == 1:
                        #strike from invariants
                        zinv = mt.get_Zinvariants()
                        s1 = zinv.strike
                        
                        #fold angles so go from -90 to 90
                        s1[np.where(s1>90)] -= 180
                        s1[np.where(s1<-90)] += 180
                        
                        #plot strike with error bars
                        ps1 = axs.errorbar(mt.period, 
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
                        pt = mt.get_PhaseTensor()
                        s2, s2_err = pt.azimuth
                        
                        #fold angles to go from -90 to 90
                        s2[np.where(s2>90)] -= 180
                        s2[np.where(s2<-90)] += 180
                        
                        #plot strike with error bars
                        ps2 = axs.errorbar(mt.period, 
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
                        tp = mt.get_Tipper()
                        s3 = tp.ang_real+90
                        
                        #fold to go from -90 to 90
                        s3[np.where(s3>90)] -= -180
                        s3[np.where(s3<-90)] += 180
                        
                        #plot strike with error bars
                        ps3 = axs.errorbar(mt.period, 
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
                        
                    #------plot skew angle-------------------------------------
                    if self.plot_skew == 'y':
                        #strike from phase tensor
                        pt = mt.get_PhaseTensor()
                        sk, sk_err = pt.beta
                        
                        
                        axs2 = axs.twinx()
                        ps4 = axs2.errorbar(mt.period, 
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
                        
                        
                        axs2.set_ylim(self.skew_limits)
                        axs2.yaxis.set_major_locator(MultipleLocator(3))
                        axs2.yaxis.set_minor_locator(MultipleLocator(1))
                        if ii == len(self.mt_lst)-1:
                            axs2.set_ylabel('Skew', color=self.skew_color)
                        else:
                            plt.setp(axs2.get_yticklabels(), visible=False)
                            
                        axs2.set_xscale('log')
                        for tl in axs2.get_yticklabels():
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
                                                
                        
                    axs.plot(axr.get_xlim(), [0, 0], color='k', lw=.5)
                    
                    if ii == 0:
                        axs.set_ylabel('Strike', fontdict=fontdict)
                        try:
                            axs.legend(stlst, 
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
                    else:
                        plt.setp(axs.get_yticklabels(), visible=False)
                    
                    axs.set_xlabel('Period (s)', fontdict=fontdict)
                    axs.set_ylim(self.strike_limits)
                    axs.yaxis.set_major_locator(MultipleLocator(30))
                    axs.yaxis.set_minor_locator(MultipleLocator(5))
                    axs.set_xscale('log')
                    axs.grid(True, 
                             alpha=.25, 
                             which='both', 
                             color=(.25, .25, .25),
                             lw=.25)
                             
                    tklabels = [self.label_dict[tt] 
                                for tt in np.arange(np.log10(self.xlimits[0]),
                                              np.log10(self.xlimits[1])+1)]
                    tklabels[0] = ''
                    tklabels[-1] = ''
                
                    axs.set_xticklabels(tklabels,
                                    fontdict={'size':self.font_size})
                    
                    # ==  == Plot the Z_xx, Z_yy components if desired ==  
                    if self.plotnum == 2:
                        #---------plot the apparent resistivity----------------
                        axr2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #res_xx
                        ebxxr = axr2.errorbar(mt.period, 
                                              rp.resxx, 
                                              marker=self.xy_marker, 
                                              ms=self.marker_size, 
                                              mfc=self.xy_mfc,
                                              mec=self.xy_color, 
                                              mew=self.marker_lw,
                                              ls=self.xy_ls, 
                                              yerr=rp.resxx_err, 
                                              ecolor=self.xy_color,
                                              capsize=self.marker_size,
                                              elinewidth=self.marker_lw)
                        
                        #res_yy                              
                        ebyyr = axr2.errorbar(mt.period, 
                                              rp.resyy, 
                                              marker=self.yx_marker, 
                                              ms=self.marker_size, 
                                              mfc=self.yx_mfc, 
                                              mec=self.yx_color, 
                                              mew=self.marker_lw, 
                                              ls=self.yx_ls, 
                                              yerr=rp.resyy_err, 
                                              ecolor=self.yx_color,
                                              capsize=self.marker_size,
                                              elinewidth=self.marker_lw)
            
                        #--> set axes properties
                        plt.setp(axr2.get_xticklabels(), visible=False)
                        
                        axr2.set_yscale('log')
                        axr2.set_xscale('log')
                        axr2.set_xlim(self.xlimits)
                        axr2.grid(True, 
                                  alpha=.25, 
                                  which='both', 
                                  color=(.25, .25, .25),
                                  lw=.25)
                        if ii == 0:
                            axr2.legend((ebxxr[0], ebyyr[0]), 
                                        ('$Z_{xx}$', '$Z_{yy}$'),
                                        loc=3, 
                                        markerscale=1, 
                                        borderaxespad=.01,
                                        labelspacing=.07, 
                                        handletextpad=.2, 
                                        borderpad=.02)
                                            
                        
                        #-----Plot the phase-----------------------------------
                        
                        axp2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #phase_xx
                        ebxxp = axp2.errorbar(mt.period, 
                                              rp.phasexx, 
                                              marker=self.xy_marker, 
                                              ms=self.marker_size, 
                                              mfc=self.xy_mfc, 
                                              mec=self.xy_color, 
                                              mew=self.marker_lw, 
                                              ls=self.xy_ls, 
                                              yerr=rp.phasexx_err,
                                              ecolor=self.xy_color,
                                              capsize=self.marker_size,
                                              elinewidth=self.marker_lw)
                                                        
                        #phase_yy
                        ebyyp = axp2.errorbar(mt.period,
                                              rp.phaseyy, 
                                              marker=self.yx_marker,
                                              ms=self.marker_size, 
                                              mfc=self.yx_mfc,
                                              mec=self.yx_color, 
                                              mew=self.marker_lw, 
                                              ls=self.yx_ls, 
                                              yerr=rp.phaseyy_err, 
                                              ecolor=self.yx_color,
                                              capsize=self.marker_size,
                                              elinewidth=self.marker_lw)
                        
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
                                  color=(.25, .25, .25),
                                  lw=.25) 
                                   
                
                # == =Plot the Determinant if desired ==  ==  ==  == 
                if self.plotnum == 3:
                        
                    #res_det
                    ebdetr = axr.errorbar(mt.period, 
                                          rp.resdet, 
                                          marker=self.det_marker, 
                                          ms=self.marker_size, 
                                          mfc=self.det_mfc,
                                          mec=self.det_color, 
                                          mew=self.marker_lw, 
                                          ls=self.det_ls, 
                                          yerr=rp.resdet_err, 
                                          ecolor=self.det_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)
                
                    #phase_det
                    ebdetp = axp.errorbar(mt.period, 
                                          rp.phasedet, 
                                          marker=self.det_marker, 
                                          ms=self.marker_size, 
                                          mfc=self.det_mfc, 
                                          mec=self.det_color, 
                                          mew=self.marker_lw, 
                                          ls=self.det_ls, 
                                          yerr=rp.phasedet_err, 
                                          ecolor=self.det_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)
                    
                    #--> set axes properties
                    plt.setp(axr.get_xticklabels(), visible=False)
                    if ii == 0:
                        axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                            fontdict=fontdict)
                    else:
                        plt.setp(axr.get_yticklabels(), visible=False)
                                    
                    axr.set_yscale('log')
                    axr.set_xscale('log')
                    axr.set_ylim(self.res_limits)
                    axr.set_xlim(self.xlimits)
                    axr.grid(True, alpha=.25,
                             which='both', 
                             color=(.25, .25, .25),
                             lw=.25)
                                  
                    #--> set axes properties
                    axp.set_xlabel('Period (s)', fontdict)
                   
                    if ii == 0:
                        axp.set_ylabel('Phase (deg)', fontdict)
                    
                    else:
                        plt.setp(axp.get_yticklabels(), visible=False)
                        
                    axp.set_xscale('log')
                    axp.set_ylim(self.phase_limits)        
                    axp.yaxis.set_major_locator(MultipleLocator(15))
                    axp.yaxis.set_minor_locator(MultipleLocator(5))
                    tklabels = [self.label_dict[tt] 
                                for tt in np.arange(np.log10(self.xlimits[0]),
                                                  np.log10(self.xlimits[1])+1)]
                    tklabels[0] = ''
                    tklabels[-1] = ''
                    
                    axp.set_xticklabels(tklabels,
                                        fontdict={'size':self.font_size})
                    axp.grid(True, alpha=.25, 
                                  which='both', 
                                  color=(.25, .25, .25),
                                  lw=.25)
                
                
                #make title and show
                axr.set_title(mt.station, fontsize=self.font_size, 
                              fontweight='bold')
                plt.show()
        
        #===Plot all responses into one plot to compare changes ==
        if self.plot_style == 'compare':
            ns = len(self.mt_lst)
            
            #make color lists for the plots going light to dark
            cxy = [(0, 0+float(cc)/ns, 1-float(cc)/ns) for cc in range(ns)]
            cyx = [(1, float(cc)/ns, 0) for cc in range(ns)]
            cdet = [(0, 1-float(cc)/ns, 0) for cc in range(ns)]
            ctipr = [(.75*cc/ns, .75*cc/ns, .75*cc/ns) for cc in range(ns)]
            ctipi = [(float(cc)/ns, 1-float(cc)/ns, .25) for cc in range(ns)]
            cst = [(.5*cc/ns, 0, .5*cc/ns) for cc in range(ns)]
            csk = [(.75*cc/ns, .25*cc/ns, 0) for cc in range(ns)]
            
            #make marker lists for the different components
            mxy = ['s', 'D', 'x', '+', '*', '1', '3', '4']*5
            myx = ['o', 'h', '8', 'p', 'H', 7, 4, 6]*5
            
            legendlst = []
            stationlst = []
            tiplst = []
            stlst = []
            sklst = []

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
            
            #set figure size according to what the plot will be.
            if self.plotnum == 1 or self.plotnum == 3:
                self.figsize = [5, 7]
                
            elif self.plotnum == 2:
                self.figsize = [6, 6]
                
            #make a figure instance
            self.fig = plt.figure(self.fignum, self.figsize, dpi=self.dpi)
            
            #make subplots
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
                        self.axs = self.fig.add_subplot(gs[2, :], 
                                                        sharex=self.axr)
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
                        self.axs = self.fig.add_subplot(gs[3, :], 
                                                        sharex=self.axr)
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
                        self.axs = self.fig.add_subplot(gs[2, :], 
                                                        sharex=self.axr)
                        self.axs.yaxis.set_label_coords(-.095, 0.5)
                        
            if self.plot_skew == 'y':
                self.axs2 = self.axs.twinx()
    
            for ii,mt in enumerate(self.mt_lst):
                    
                #get the reistivity and phase object
                rp = mt.get_ResPhase()
                rp.rotate(self.rot_z[ii])
                
                #set x-axis limits from short period to long period
                if self.xlimits == None:
                    self.xlimits = (10**(np.floor(np.log10(mt.period[0]))),
                                    10**(np.ceil(np.log10((mt.period[-1])))))
                if self.phase_limits == None:
                    self.phase_limits = (0, 89.9)
                
                stationlst.append(mt.station)
                
                # ==  ==  ==  == =Plot Z_xy and Z_yx ==
                if self.plotnum == 1 or self.plotnum == 2:
                    #---------plot the apparent resistivity--------------------
                    #--> plot as error bars and just as points xy-blue, yx-red
                    #res_xy
                    ebxyr = self.axr.errorbar(mt.period, 
                                             rp.resxy, 
                                             marker=mxy[ii], 
                                             ms=self.marker_size, 
                                             mfc='None', 
                                             mec=cxy[ii], 
                                             mew=self.marker_lw, 
                                             ls=self.xy_ls, 
                                             yerr=rp.resxy_err, 
                                             ecolor=cxy[ii],
                                             capsize=self.marker_size,
                                             elinewidth=self.marker_lw)
                    
                    #res_yx                              
                    ebyxr = self.axr.errorbar(mt.period, 
                                             rp.resyx, 
                                             marker=myx[ii], 
                                             ms=self.marker_size, 
                                             mfc='None',
                                             mec=cyx[ii], 
                                             mew=self.marker_lw,
                                             ls=self.yx_ls, 
                                             yerr=rp.resyx_err, 
                                             ecolor=cyx[ii],
                                             capsize=self.marker_size,
                                             elinewidth=self.marker_lw)
                                                  

                    #-----Plot the phase---------------------------------------
                    #phase_xy
                    ebxyp = self.axp.errorbar(mt.period, 
                                             rp.phasexy, 
                                             marker=mxy[ii], 
                                             ms=self.marker_size, 
                                             mfc='None',
                                             mec=cxy[ii], 
                                             mew=self.marker_lw,
                                             ls=self.xy_ls,
                                             yerr=rp.phasexy_err, 
                                             ecolor=cxy[ii],
                                             capsize=self.marker_size,
                                             elinewidth=self.marker_lw)
                                                  
                    #phase_yx: Note add 180 to place it in same quadrant as
                    #phase_xy
                    ebyxp = self.axp.errorbar(mt.period, 
                                             rp.phaseyx, 
                                             marker=myx[ii], 
                                             ms=self.marker_size, 
                                             mfc='None',
                                             mec=cyx[ii], 
                                             mew=self.marker_lw,
                                             ls=self.yx_ls, 
                                             yerr=rp.phaseyx_err, 
                                             ecolor=cyx[ii],
                                             capsize=self.marker_size,
                                             elinewidth=self.marker_lw)
            

                    legendlst.append([ebxyr, ebyxr])             
                    
                    # ==== Plot the Z_xx, Z_yy components if desired ==               
                    if self.plotnum == 2:
                        #---------plot the apparent resistivity----------------
                        self.axr2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #res_xx
                        ebxxr = self.axr2.errorbar(mt.period, 
                                                  rp.resxx, 
                                                  marker=mxy[ii], 
                                                  ms=self.marker_size, 
                                                  mfc='None',
                                                  mec=cxy[ii], 
                                                  mew=self.marker_lw,
                                                  ls=self.xy_ls, 
                                                  yerr=rp.resxx_err, 
                                                  ecolor=cxy[ii],
                                                  capsize=self.marker_size,
                                                  elinewidth=self.marker_lw)
                        
                        #res_yy                              
                        ebyyr = self.axr2.errorbar(mt.period, 
                                                  rp.resyy, 
                                                  marker=myx[ii], 
                                                  ms=self.marker_size, 
                                                  mfc='None', 
                                                  mec=cyx[ii], 
                                                  mew=self.marker_lw, 
                                                  ls=self.yx_ls, 
                                                  yerr=rp.resyy_err, 
                                                  ecolor=cyx[ii],
                                                  capsize=self.marker_size,
                                                  elinewidth=self.marker_lw)
            
                                            
                        #-----Plot the phase-----------------------------------
                        
                        self.axp2.yaxis.set_label_coords(-.1, 0.5)
                        
                        #phase_xx
                        ebxxp = self.axp2.errorbar(mt.period, 
                                                  rp.phasexx, 
                                                  marker=mxy[ii], 
                                                  ms=self.marker_size, 
                                                  mfc='None', 
                                                  mec=cxy[ii], 
                                                  mew=self.marker_lw, 
                                                  ls=self.xy_ls, 
                                                  yerr=rp.phasexx_err,
                                                  ecolor=cxy[ii],
                                                  capsize=self.marker_size,
                                                  elinewidth=self.marker_lw)
                                                        
                        #phase_yy
                        ebyyp = self.axp2.errorbar(mt.period,
                                                  rp.phaseyy, 
                                                  marker=myx[ii],
                                                  ms=self.marker_size, 
                                                  mfc='None',
                                                  mec=cyx[ii], 
                                                  mew=self.marker_lw, 
                                                  ls=self.yx_ls, 
                                                  yerr=rp.phaseyy_err, 
                                                  ecolor=cyx[ii],
                                                  capsize=self.marker_size,
                                                  elinewidth=self.marker_lw)
                        

                                   
                
                #===Plot the Determinant if desired ==                           
                if self.plotnum == 3:
                        
                    #res_det
                    ebdetr = self.axr.errorbar(mt.period, 
                                              rp.resdet, 
                                              marker=mxy[ii], 
                                              ms=self.marker_size, 
                                              mfc='None',
                                              mec=cdet[ii], 
                                              mew=self.marker_lw, 
                                              ls=self.det_ls, 
                                              yerr=rp.resdet_err, 
                                              ecolor=cdet[ii],
                                              capsize=self.marker_size,
                                              elinewidth=self.marker_lw)
                
                    #phase_det
                    ebdetp = self.axp.errorbar(mt.period, 
                                              rp.phasedet, 
                                              marker=mxy[ii], 
                                              ms=self.marker_size, 
                                              mfc='None', 
                                              mec=cdet[ii], 
                                              mew=self.marker_lw, 
                                              ls=self.det_ls, 
                                              yerr=rp.phasedet_err, 
                                              ecolor=cdet[ii],
                                              capsize=self.marker_size,
                                              elinewidth=self.marker_lw)
                    
                    legendlst.append(ebdetr)
                    
                #-----plot tipper----------------------------------------------              
                if self.plot_tipper.find('y') == 0:
                    
                    tp = mt.get_Tipper()
                    
                    txr = tp.mag_real*np.cos(tp.ang_real*np.pi/180+\
                                             np.pi*self.arrow_direction)
                    tyr = tp.mag_real*np.sin(tp.ang_real*np.pi/180+\
                                             np.pi*self.arrow_direction)
            
                    txi = tp.mag_imag*np.cos(tp.ang_imag*np.pi/180+\
                                             np.pi*self.arrow_direction)
                    tyi = tp.mag_imag*np.sin(tp.ang_imag*np.pi/180+\
                                             np.pi*self.arrow_direction)
                    
                    nt = len(txr)
                    
                    for aa in range(nt):
                        xlenr = txr[aa]*mt.period[aa]
                        xleni = txi[aa]*mt.period[aa]
                        
                        #scale the arrow head height and width to fit in a log
                        #scale
                        if np.log10(mt.period[aa])<0:
                            hwidth = self.arrow_head_width*\
                                      10**(np.floor(np.log10(mt.period[aa])))
                            hheight = self.arrow_head_height*\
                                      10**(np.floor(np.log10(mt.period[aa])))
                        else:
                            hwidth = self.arrow_head_width/\
                                      10**(np.floor(np.log10(mt.period[aa])))
                            hheight = self.arrow_head_height/\
                                      10**(np.floor(np.log10(mt.period[aa]))) 
                            
                        #--> plot real arrows
                        if self.plot_tipper.find('r') > 0:
                            self.axt.arrow(mt.period[aa],
                                           0,
                                           xlenr,
                                           tyr[aa],
                                           lw=self.arrow_lw,
                                           facecolor=ctipr[ii],
                                           edgecolor=ctipr[ii],
                                           head_width=hwidth,
                                           head_length=hheight,
                                           length_includes_head=False)
                            
                                           
                        #--> plot imaginary arrows
                        if self.plot_tipper.find('i')>0:               
                            self.axt.arrow(mt.period[aa],
                                           0,
                                           xleni,
                                           tyi[aa],
                                           lw=alw,
                                           facecolor=ctipi[ii],
                                           edgecolor=ctipi[ii],
                                           length_includes_head=False)
                        
                    lt = self.axt.plot(0, 0, lw=1, color=ctipr[ii])
                    tiplst.append(lt[0])
                        
                    
            
                #---plot strike------------------------------------------------
                if self.plot_strike == 'y' or self.plot_strike == 2:
                    
                    #strike from phase tensor
                    pt = mt.get_PhaseTensor()
                    s2, s2_err = pt.azimuth
                    
                    #fold angles to go from -90 to 90
                    s2[np.where(s2>90)] = s2[np.where(s2>90)]-180
                    s2[np.where(s2<-90)] = s2[np.where(s2<-90)]+180
                    
                    #plot strike with error bars
                    ps2 = self.axs.errorbar(mt.period, 
                                            s2, 
                                            marker=mxy[ii], 
                                            ms=self.marker_size, 
                                            mfc=cst[ii], 
                                            mec=cst[ii], 
                                            mew=self.marker_lw,
                                            ls='none', 
                                            yerr=s2_err, 
                                            ecolor=cst[ii],
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
                                            
                    stlst.append(ps2[0])
                                        
                #------plot skew angle-----------------------------------------
                if self.plot_skew == 'y':
                    #strike from phase tensor
                    pt = mt.get_PhaseTensor()
                    sk, sk_err = pt.beta
                    
                    ps4 = self.axs2.errorbar(mt.period, 
                                            sk, 
                                            marker=myx[ii], 
                                            ms=self.marker_size, 
                                            mfc=csk[ii], 
                                            mec=csk[ii], 
                                            mew=self.marker_lw,
                                            ls='none', 
                                            yerr=sk_err, 
                                            ecolor=csk[ii],
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)
    
                    if self.skew_limits is None:
                        self.skew_limits = (-9, 9)
                    sklst.append(ps4[0])
                                    
            #-------set axis properties----------------------------------------
            self.axr.set_yscale('log')
            self.axr.set_xscale('log')
            self.axr.set_ylim(self.res_limits)
            self.axr.set_xlim(self.xlimits)
            self.axr.grid(True, alpha=.25, 
                          which='both', 
                          color=(.25, .25, .25),
                          lw=.25)
                          
            plt.setp(self.axr.get_xticklabels(), visible=False)
            
            self.axr.set_ylabel('App. Resistivity ($\Omega \cdot$m)',
                                fontdict=fontdict)
                          
            #check the phase to see if any point are outside of [0:90]
            if self.phase_limits == None:
                self.phase_limits = (0, 89.99)
            
            #--> set axes properties
            self.axp.set_xlabel('Period (s)', fontdict=fontdict)
            self.axp.set_ylabel('Phase (deg)',fontdict=fontdict)
            self.axp.set_xscale('log')
            self.axp.set_ylim(self.phase_limits)        
            self.axp.yaxis.set_major_locator(MultipleLocator(15))
            self.axp.yaxis.set_minor_locator(MultipleLocator(5))
            self.axp.grid(True, alpha=.25, 
                          which='both', 
                          color=(.25, .25, .25),
                          lw=.25)
            
            #make legend
            if self.plotnum == 1:
                llst = [ll[0] for ll in legendlst]+[ll[1] for ll in legendlst]
                slst = [ss+'_xy' for ss in stationlst]+\
                       [ss+'_yx' for ss in stationlst]
                     
                self.axr.legend(llst, 
                                slst,
                                loc=3,
                                ncol=2,
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)
            elif self.plotnum == 3:
                llst = [ll[0] for ll in legendlst]
                slst = [ss+'_det' for ss in stationlst]
                     
                self.axr.legend(llst, 
                                slst,
                                loc=3,
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)

            
            if self.plotnum == 2:
                llst = [ll[0] for ll in legendlst]+[ll[1] for ll in legendlst]
                slst = [ss+'_xx' for ss in stationlst]+\
                       [ss+'_yy' for ss in stationlst]
                     
                self.axr2.legend(llst, 
                                 slst,
                                 loc=3,
                                 ncol=2,
                                 markerscale=1, 
                                 borderaxespad=.01,
                                 labelspacing=.07, 
                                 handletextpad=.2, 
                                 borderpad=.02)
                
                self.axr2.set_yscale('log')
                self.axr2.set_xscale('log')
                self.axr2.set_xlim(self.xlimits)
                self.axr2.grid(True, 
                               alpha=.25, 
                               which='both', 
                               color=(.25,.25,.25),
                              lw=.25)
                plt.setp(self.axr2.get_xticklabels(), visible=False)
                                
                #--> set axes properties Phaseyx
                self.axp2.set_xlabel('Period (s)', fontdict)
                self.axp2.set_xscale('log')
                self.axp2.set_ylim(ymin=-179.9, ymax=179.9)        
                self.axp2.yaxis.set_major_locator(MultipleLocator(30))
                self.axp2.yaxis.set_minor_locator(MultipleLocator(5))
                self.axp2.grid(True, 
                               alpha=.25, 
                               which='both', 
                               color=(.25, .25, .25),
                              lw=.25) 
                          
            if self.plot_tipper.find('y') == 0:
                self.axt.plot(self.axr.get_xlim(), [0, 0], color='k', lw=.5)
                #--> set axis properties Tipper 
                plt.setp(self.axp.get_xticklabels(), visible=False)
                self.axp.set_xlabel('')
                if self.plotnum == 2:
                    plt.setp(self.axp2.get_xticklabels(), visible=False)
                    self.axp2.set_xlabel('')
                    
                self.axt.yaxis.set_major_locator(MultipleLocator(.2))               
                self.axt.yaxis.set_minor_locator(MultipleLocator(.1))               
                self.axt.set_xlabel('Period (s)', fontdict=fontdict)
                self.axt.set_ylabel('Tipper', fontdict=fontdict)    
                
                self.axt.set_xscale('log')
                if self.tipper_limits is None:
                    tmax = max([np.sqrt(txr.max()**2+tyr.max()**2),
                                np.sqrt(txi.max()**2+tyi.max()**2)])
                    if tmax > 1:
                        tmax = .99
                                
                    tmin = -min([np.sqrt(txr.min()**2+tyr.min()**2),
                                np.sqrt(txi.min()**2+tyi.min()**2)])
                    if tmin < -1:
                        tmin = -.99
                                
                    self.tipper_limits = (tmin-.1, tmax+.1)
                
                self.axt.set_ylim(self.tipper_limits)
                self.axt.grid(True, alpha=.25, 
                              which='both', 
                              color=(.25,.25,.25),
                              lw=.25)
                              
                self.axt.legend(tiplst, 
                                stationlst,
                                loc=3,
                                ncol=2,
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)
                              
            #--> set axes properties for strike and skew
            if self.plot_strike == 'y':
                plt.setp(self.axt.get_xticklabels(), visible=False)
                self.axt.set_xlabel('')
                                        
                if self.strike_limits is None:
                    self.strike_limits = (-89.99, 89.99)
                self.axs.plot(self.axr.get_xlim(), [0, 0], color='k', lw=.5)
                
                self.axs.set_ylabel('Strike',
                                    fontdict=fontdict)
                self.axs.set_xlabel('Period (s)',
                                    fontdict=fontdict)
                self.axs.set_ylim(self.strike_limits)
                self.axs.yaxis.set_major_locator(MultipleLocator(30))
                self.axs.yaxis.set_minor_locator(MultipleLocator(5))
                self.axs.set_xscale('log')
                self.axs.grid(True, alpha=.25, 
                              which='both', 
                              color=(.25, .25, .25),
                              lw=.25)
                self.axs.legend(stlst, 
                                stationlst,
                                loc=3,
                                ncol=2,
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)
           
           #--> set axes properties for skew
            if self.plot_skew == 'y':
                self.axs2.set_ylim(self.skew_limits)
                self.axs2.yaxis.set_major_locator(MultipleLocator(3))
                self.axs2.yaxis.set_minor_locator(MultipleLocator(1))
                self.axs2.set_ylabel('Skew', color=self.skew_color)
                self.axs2.set_xscale('log')
                for tl in self.axs2.get_yticklabels():
                    tl.set_color(self.skew_color)
                    
                self.axs2.legend(sklst, 
                                 stationlst,
                                 loc=4,
                                 ncol=2,
                                 markerscale=1, 
                                 borderaxespad=.01,
                                 labelspacing=.07, 
                                 handletextpad=.2, 
                                 borderpad=.02)
                    
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
