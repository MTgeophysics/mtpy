# -*- coding: utf-8 -*-
"""
PlotResidualPhaseTensorPseudoSection
=======================

    *plots the residual phase tensor for two different sets of measurments
    
    
Created on Wed Oct 16 14:56:04 2013

@author: jpeacock-pr
"""
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter
import mtpy.utils.conversions as utm2ll
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl
import mtpy.analysis.pt as mtpt
import scipy.signal as sps

#==============================================================================

class PlotResidualPTps(mtpl.MTEllipse):
    """
    plot residual phase tensor pseudo section. 
    
    
    Arguments:
    ----------
    
        **fn_list1** : list of strings
                        full paths to .edi files for survey 1

        **fn_list2** : list of strings
                        full paths to .edi files for survey 2
                        
        **Note** it is assumed that the edi file lists have the same number
                 of edi files and have the same station names.  If you get
                 an error check this first.
                 
        
        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor 
                          ellipses with keys:
                              * 'size' -> size of ellipse in points 
                                         *default* is 2
                              
                              * 'colorby' : [ 'phimin' | 'phimax' | 'skew' | 
                                              'skew_seg' | 'phidet' | 
                                              'ellipticity' ]
                                        
                                        - 'phimin' -> colors by minimum phase
                                        - 'phimax' -> colors by maximum phase
                                        - 'skew' -> colors by beta (skew)
                                        - 'skew_seg' -> colors by beta in 
                                                       discrete segments 
                                                       defined by the range
                                        - 'phidet' -> colors by determinant of
                                                     the phase tensor
                                        - 'ellipticity' -> colors by ellipticity
                                        *default* is 'phimin'
                                
                               * 'range' : tuple (min, max, step)
                                     Need to input at least the min and max
                                     and if using 'skew_seg' to plot
                                     discrete values input step as well
                                     *default* depends on 'colorby'
                                     
                          * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 
                                      'mt_wh2bl' | 'mt_rd2bl' | 
                                      'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                                      'mt_rd2gr2bl' ]
                                      
                                   - 'mt_yl2rd' -> yellow to red
                                   - 'mt_bl2yl2rd' -> blue to yellow to red
                                   - 'mt_wh2bl' -> white to blue
                                   - 'mt_rd2bl' -> red to blue
                                   - 'mt_bl2wh2rd' -> blue to white to red
                                   - 'mt_bl2gr2rd' -> blue to green to red
                                   - 'mt_rd2gr2bl' -> red to green to blue
                                   - 'mt_seg_bl2wh2rd' -> discrete blue to 
                                                         white to red
        
                                         
        
        **stretch** : float or tuple (xstretch, ystretch)
                        is a factor that scales the distance from one 
                        station to the next to make the plot readable.
                        *Default* is 200
                        
        **linedir** : [ 'ns' | 'ew' ]
                      predominant direction of profile line
                      * 'ns' -> North-South Line
                      * 'ew' -> East-West line
                      *Default* is 'ns'
        
        **station_id** : tuple or list 
                        start and stop of station name indicies.  
                        ex: for MT01dr station_id=(0,4) will be MT01
        
        **rotz** : float or np.ndarray
                   angle in degrees to rotate the data clockwise positive.
                   Can be an array of angle to individually rotate stations or
                   periods or both. 
                       - If rotating each station by a constant
                         angle the array needs to have a shape of 
                         (# of stations)
                        - If rotating by period needs to have shape 
                           # of periods
                        - If rotating both individually shape=(ns, nf)
                  *Default* is 0
        
        **title** : string
                    figure title
                    
        **dpi** : int 
                  dots per inch of the resolution. *default* is 300
                    
                       
        **fignum** : int
                     figure number.  *Default* is 1
        
                         
        **tscale** : [ 'period' | 'frequency' ]
        
                     * 'period'    -> plot vertical scale in period
                     
                     * 'frequency' -> plot vertical scale in frequency
                     
        **cb_dict** : dictionary to control the color bar
        
                      * 'orientation' : [ 'vertical' | 'horizontal' ]
                                       orientation of the color bar 
                                       *default* is vertical
                                       
                      * 'position' : tuple (x,y,dx,dy)
                                    - x -> lateral position of left hand corner 
                                          of the color bar in figure between 
                                          [0,1], 0 is left side
                                          
                                    - y -> vertical position of the bottom of 
                                          the color bar in figure between 
                                          [0,1], 0 is bottom side.
                                          
                                    - dx -> width of the color bar [0,1]
                                    
                                    - dy -> height of the color bar [0,1]
        **font_size** : float
                        size of the font that labels the plot, 2 will be added
                        to this number for the axis labels.
                        
        **plot_yn** : [ 'y' | 'n' ]
                      * 'y' to plot on creating an instance
                      
                      * 'n' to not plot on creating an instance
                      
        **xlim** : tuple(xmin, xmax)
                   min and max along the x-axis in relative distance of degrees
                   and multiplied by xstretch
                   
        **ylim** : tuple(ymin, ymax)
                   min and max period to plot, note that the scaling will be
                   done in the code.  So if you want to plot from (.1s, 100s)
                   input ylim=(.1,100)
    
    To get a list of .edi files that you want to plot -->
    :Example: ::
        
        >>> import mtpy.imaging.mtplot as mtplot
        >>> import os
        >>> edipath1 = r"/home/EDIfiles1"
        >>> edilist1 = [os.path.join(edipath1,edi) for edi in os.listdir(edipath1)
        >>> ...       if edi.find('.edi')>0]
        >>> edipath2 = r"/home/EDIfiles2"
        >>> edilist2 = [os.path.join(edipath2,edi) for edi in os.listdir(edipath2)
        >>> ...       if edi.find('.edi')>0]
        >>> # color by phimin with a range of 0-5 deg
    
    * If you want to plot minimum phase colored from blue to red in a range of
     20 to 70 degrees you can do it one of two ways--> 
    
    1)          
    :Example: ::
        
        >>> edict = {'range':(20,70), 'cmap':'mt_bl2gr2rd','colorby':'phimin'}
        >>> pt1 = mtplot.residual_pt_ps(edilist1, edilst2, ellipse_dict=edict)
     
    2)
    :Example: ::
        
        >>> pt1 = mtplot.residual_pt_ps(edilist1, edilst2, ellipse_dict=edict,\
                                        plot_yn='n')
        >>> pt1.ellipse_colorby = 'phimin'
        >>> pt1.ellipse_cmap = 'mt_bl2gr2rd'
        >>> pt1.ellipse_range = (20,70)
        >>> pt1.plot()
        
    * If you want to add real induction arrows that are scaled by 10 and point
     away from a conductor --> 
    :Example: ::
        
        >>> pt1.plot_tipper = 'yr'
        >>> pt1.arrow_size = 10
        >>> pt1.arrow_direction = -1
        >>> pt1.redraw_plot()
    
    * If you want to save the plot as a pdf with a generic name -->
    :Example: ::
        >>> pt1.save_figure(r"/home/PTFigures", file_format='pdf', dpi=300)
        File saved to '/home/PTFigures/PTPseudoSection.pdf'
        
    Attributes:
    -----------
        -arrow_color_imag     color of imaginary induction arrow
        -arrow_color_real     color of real induction arrow
        -arrow_direction      convention of arrows pointing to or away from 
                              conductors, see above.
        -arrow_head_length    length of arrow head in relative points
        -arrow_head_width     width of arrow head in relative points
        -arrow_lw             line width of arrows
        -arrow_size           scaling factor to multiple arrows by to be visible
        -arrow_threshold      threshold for plotting arrows, anything above 
                              this number will not be plotted.
        
        -ax                   matplotlib.axes instance for the main plot
        -ax2                  matplotlib.axes instance for the color bar
        -cb                   matplotlib.colors.ColorBar instance for color bar
        -cb_orientation       color bar orientation ('vertical' | 'horizontal')
        -cb_position          color bar position (x, y, dx, dy)
        
        -dpi                  dots-per-inch resolution
        
        -ellipse_cmap         ellipse color map, see above for options
        -ellipse_colorby      parameter to color ellipse by
        -ellipse_range        (min, max, step) values to color ellipses
        -ellipse_size         scaling factor to make ellipses visible
        
        -fig                  matplotlib.figure instance for the figure 
        -fignum               number of figure being plotted
        -figsize              size of figure in inches
        -font_size            font size of axes tick label, axes labels will be
                              font_size + 2
        
        -linedir              prominent direction of profile being plotted 
             
        -mt_list               list of mtplot.MTplot instances containing all
                              the important information for each station
        -offsetlist            array of relative offsets of each station
        
        -plot_tipper          string to inform program to plot induction arrows
        -plot_yn              plot the pseudo section on instance creation
        
        -rot_z                rotates the data by this angle assuming North is
                              0 and angle measures clockwise
                              
        -station_id            index [min, max] to reaad station name
        -stationlist           list of stations plotted
        -title                title of figure
        -tscale               temporal scale of y-axis ('frequency' | 'period')
        
        -xlimits              limits on x-axis (xmin, xmax)
        -xstretch             scaling factor to stretch x offsets
        
        -ylimits              limits on y-axis (ymin, ymax)
        -ystep                step to set major ticks on y-axis
        -ystretch             scaling factor to strech axes in y direction
        
    Methods:
    --------

        -plot                 plots the pseudo section
        -redraw_plot          on call redraws the plot from scratch
        -save_figure          saves figure to a file of given format
        -update_plot          updates the plot while still active
        -writeTextFiles       writes parameters of the phase tensor and tipper
                              to text files.

    
    """
    
    def __init__(self, fn_list1, fn_list2, **kwargs):
        assert len(fn_list1) == len(fn_list2)
        self.fn_list1 = fn_list1
        self.fn_list2 = fn_list2
        
        self.mt_list1 = mtpl.get_mtlist(fn_list=fn_list1)
        self.mt_list2 = mtpl.get_mtlist(fn_list=fn_list2)
        
        self.residual_pt_list = []
        self.med_filt_kernel = kwargs.pop('med_filt_kernel', None)
        
        #--> set the ellipse properties
        self._ellipse_dict = kwargs.pop('ellipse_dict', {'cmap':'mt_wh2or',
                                                         'range':(0, 2)})
        self._read_ellipse_dict()
        
        #--> set colorbar properties---------------------------------
        #set orientation to horizontal
        cb_dict = kwargs.pop('cb_dict', {})
        try:
            self.cb_orientation = cb_dict['orientation']
        except KeyError:
            self.cb_orientation = 'vertical'
        
        #set the position to middle outside the plot            
        try:
            self.cb_position = cb_dict['position']
        except KeyError:
            self.cb_position = None
            
        #--> set plot properties ------------------------------
        self.fig_num = kwargs.pop('fig_num', 'residual_PT')
        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.tscale = kwargs.pop('tscale', 'period')
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.linedir = kwargs.pop('linedir', 'ew')
        self.font_size = kwargs.pop('font_size', 7)
        self.station_id = kwargs.pop('station_id', [0,4])
        self.ystep = kwargs.pop('ystep', 4)
        self.xstep = kwargs.pop('xstep', 1)
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)
        
        #--> set the freq to plot
        self.ftol = kwargs.pop('ftol', .1)
        
        #--> set spacing of plot
        self.subplot_wspace = .1
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1
        
        #set the stretching in each direction
        stretch = kwargs.pop('stretch', (50, 25))
        if type(stretch) == float or type(stretch) == int:
            self.xstretch = stretch
            self.ystretch = stretch
        else:
            self.xstretch = stretch[0]
            self.ystretch = stretch[1]
        
        #if rotation angle is an int or float make an array the length of 
        #mt_list for plotting purposes
        
        self._rot_z = kwargs.pop('rot_z', 0)
        if type(self._rot_z) is float or type(self._rot_z) is int:
            self._rot_z = np.array([self._rot_z]*len(self.mt_list1))
        
        #if the rotation angle is an array for rotation of different 
        #freq than repeat that rotation array to the len(mt_list)
        elif type(self._rot_z) is np.ndarray:
            if self._rot_z.shape[0]  !=  len(self.mt_list1):
                self._rot_z = np.repeat(self._rot_z, len(self.mt_list1))
                
        else:
            pass
        
            
        #--> set station name properties
        station_dict = kwargs.pop('station_dict', {})
        self.station_id = station_dict.pop('id', (0, 3))
        self.station_pad = station_dict.pop('pad', .0005)
        self.station_font_dict = station_dict.pop('font_dict', 
                                                  {'size':self.font_size,
                                                   'weight':'bold'})     
        
        #--> plot if desired ------------------------
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()
        
    #---need to rotate data on setting rotz
    def _set_rot_z(self, rot_z):
        """
        need to rotate data when setting z
        """
        
        #if rotation angle is an int or float make an array the length of 
        #mt_list for plotting purposes
        if type(rot_z) is float or type(rot_z) is int:
            self._rot_z = np.array([rot_z]*len(self.mt_list))
        
        #if the rotation angle is an array for rotation of different 
        #freq than repeat that rotation array to the len(mt_list)
        elif type(rot_z) is np.ndarray:
            if rot_z.shape[0]!=len(self.mt_list):
                self._rot_z = np.repeat(rot_z, len(self.mt_list))
                
        else:
            pass
            
        for ii,rpt in enumerate(self.residual_pt_list):
            rpt.rot_z = self._rot_z[ii]
    def _get_rot_z(self):
        return self._rot_z
        
    rot_z = property(fget=_get_rot_z, fset=_set_rot_z, 
                     doc="""rotation angle(s)""")
                     
    def _compute_residual_pt(self):
        """
        compute residual phase tensor so the result is something useful to 
        plot
        """
        self.residual_pt_list = []
        for mt1 in self.mt_list1:
            station_find = False
            for mt2 in self.mt_list2:
                if mt2.station == mt1.station:
                    pt1 = mt1.get_PhaseTensor()
                    pt2 = mt2.get_PhaseTensor()
                    rpt = mtpt.ResidualPhaseTensor(pt1, pt2)
                    rpt.compute_residual_pt(pt1, pt2)
                    rpt.station = mt1.station
                    rpt.lat = mt1.lat
                    rpt.lon = mt1.lon
                    rpt.freq = mt1.freq
                    self.residual_pt_list.append(rpt)
                    station_find = True
                    break
                else:
                    pass
            if station_find == False:
                print 'Did not find {0} from list 1 in list 2'.format(mt1.station) 
                
    def _apply_median_filter(self, kernel=(3, 3)):
        """
        apply a median filter to the data to remove extreme outliers
        
        kernel is (station, frequency)
        
        """
        
        if self.ellipse_colorby == 'phimin':
            color_array = np.array([rpt.residual_pt.phimin[0] for 
                                    rpt in self.residual_pt_list])
            filt_array = sps.medfilt2d(color_array, kernel_size=kernel)
            for ss in range(filt_array.shape[0]):
                self.residual_pt_list[ss].residual_pt.phimin[0] = filt_array[ss]
        
        elif self.ellipse_colorby == 'phimax':
            color_array = np.array([rpt.residual_pt.phimax[0] for 
                                    rpt in self.residual_pt_list])
            filt_array = sps.medfilt2d(color_array, kernel_size=kernel)
            for ss in range(filt_array.shape[0]):
                self.residual_pt_list[ss].residual_pt.phimax[0] = filt_array[ss]
                
        elif self.ellipse_colorby == 'skew':
            color_array = np.array([rpt.residual_pt.beta[0] for 
                                    rpt in self.residual_pt_list])
            filt_array = sps.medfilt2d(color_array, kernel_size=kernel)
            for ss in range(filt_array.shape[0]):
                self.residual_pt_list[ss].residual_pt.beta[0] = filt_array[ss]
        
        #--> need to do azimuth for all
        color_array = np.array([rpt.residual_pt.azimuth[0] for 
                                rpt in self.residual_pt_list])
        filt_array = sps.medfilt2d(color_array, kernel_size=kernel)
        for ss in range(filt_array.shape[0]):
            self.residual_pt_list[ss].residual_pt.azimuth[0] = filt_array[ss] 
            
    def _get_ellipse_size_max(self):
        """
        get the maximum size of ellipses for the given period to normalize by
        suggesting size of change
        
        """
        color_array = np.array([rpt.residual_pt.phimax[0] for 
                                    rpt in self.residual_pt_list])
                             
        return color_array.max()
      
    def plot(self):
        """
        plot residual phase tensor
        """                            
        #get residual phase tensor for plotting        
        self._compute_residual_pt()
        
        #filter data if desired
        if self.med_filt_kernel is not None:
            self._apply_median_filter(kernel=self.med_filt_kernel)
        
        
        #set position properties for the plot
        plt.rcParams['font.size']=self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace
        
        #make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        
        self.ax = self.fig.add_subplot(1, 1, 1, aspect='equal')
        
        #create empty lists to put things into
        self.stationlist = []
        self.offsetlist = []
        minlist = []
        maxlist = []
        plot_periodlist = None
        
        #set local parameters with shorter names
        es = self.ellipse_size
        ck = self.ellipse_colorby
        cmap = self.ellipse_cmap
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            ckstep = 3
                
        nseg = float((ckmax-ckmin)/(2*ckstep))

        if cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
            
        #get largest ellipse
        emax = self._get_ellipse_size_max()
            
        #plot phase tensor ellipses
        for ii, rpt in enumerate(self.residual_pt_list):
            self.stationlist.append(
                              rpt.station[self.station_id[0]:self.station_id[1]])
            
            #set the an arbitrary origin to compare distance to all other 
            #stations.
            if ii == 0:
                east0 = rpt.lon
                north0 = rpt.lat
                offset = 0.0
            else:
                east = rpt.lon
                north = rpt.lat
                if self.linedir == 'ew': 
                    if east0 < east:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif east0 > east:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset = 0
                elif self.linedir == 'ns':
                    if north0 < north:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif north0 > north:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset = 0
                        
            self.offsetlist.append(offset)
            
            periodlist = 1./rpt.freq[::-1]
            phimax = rpt.residual_pt.phimax[0][::-1]
            phimin = rpt.residual_pt.phimin[0][::-1]
            azimuth = rpt.residual_pt.azimuth[0][::-1]
                
            #get the properties to color the ellipses by
            if self.ellipse_colorby == 'phimin':
                colorarray = rpt.residual_pt.phimin[0][::-1]
                
            elif self.ellipse_colorby == 'phimax':
                colorarray = rpt.residual_pt.phimin[0][::-1]
                
            elif self.ellipse_colorby == 'phidet':
                colorarray = np.sqrt(abs(rpt.residual_pt.det[::-1]))*(180/np.pi)
                
            elif self.ellipse_colorby == 'skew' or\
                 self.ellipse_colorby == 'skew_seg':
                colorarray = rpt.residual_pt.beta[0][::-1]
                
            elif self.ellipse_colorby == 'ellipticity':
                colorarray = rpt.residual_pt.ellipticity[::-1]
                
            else:
                raise NameError(self.ellipse_colorby+' is not supported')
            
            #get the number of periods
            n = len(periodlist)
            
            if ii == 0:
                plot_periodlist = periodlist
            
            else:
                if n > len(plot_periodlist):
                    plot_periodlist = periodlist
            
            #get min and max of the color array for scaling later
            minlist.append(min(colorarray))
            maxlist.append(max(colorarray))

            for jj, ff in enumerate(periodlist):
                
                #make sure the ellipses will be visable
                eheight = phimin[jj]/emax*es
                ewidth = phimax[jj]/emax*es
            
                #create an ellipse scaled by phimin and phimax and orient
                #the ellipse so that north is up and east is right
                #need to add 90 to do so instead of subtracting
                ellipd = patches.Ellipse((offset*self.xstretch,
                                          np.log10(ff)*self.ystretch),
                                            width=ewidth,
                                            height=eheight,
                                            angle=azimuth[jj]-90)
                                            
                #get ellipse color
                if cmap.find('seg')>0:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[jj],
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                else:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[jj],
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax))
                    
                # == =add the ellipse to the plot == ========
                self.ax.add_artist(ellipd)
                
                
        #--> Set plot parameters 
        self._plot_periodlist = plot_periodlist
        n = len(plot_periodlist)
        
        
        #calculate minimum period and maximum period with a stretch factor
        pmin = np.log10(plot_periodlist.min())*self.ystretch
        pmax = np.log10(plot_periodlist.max())*self.ystretch
        
        #need to sort the offsets and station labels so they plot correctly
        sdtype = [('offset', np.float), ('station','|S10')]
        slist = np.array([(oo, ss) for oo, ss in zip(self.offsetlist, 
                         self.stationlist)], dtype=sdtype)
        offset_sort = np.sort(slist, order='offset')
     
        self.offsetlist = offset_sort['offset']
        self.stationlist = offset_sort['station']
        
        #set y-ticklabels
        if self.tscale == 'period':
            yticklabels = ['{0:>4}'.format('{0: .1e}'.format(plot_periodlist[ll])) 
                            for ll in np.arange(0, n, self.ystep)]+\
                        ['{0:>4}'.format('{0: .1e}'.format(plot_periodlist[-1]))]
            
            self.ax.set_ylabel('Period (s)',
                               fontsize=self.font_size+2,
                               fontweight='bold')
                               
        elif self.tscale == 'frequency':
            yticklabels = ['{0:>4}'.format('{0: .1e}'.format(1./plot_periodlist[ll])) 
                            for ll in np.arange(0, n, self.ystep)]+\
                            ['{0:>4}'.format('{0: .1e}'.format(1./plot_periodlist[-1]))]
            
            self.ax.set_ylabel('Frequency (Hz)',
                               fontsize=self.font_size+2,
                               fontweight='bold')
        #set x-axis label                       
        self.ax.set_xlabel('Station',
                           fontsize=self.font_size+2,
                           fontweight='bold')
         
        #--> set tick locations and labels
        #set y-axis major ticks
        self.ax.yaxis.set_ticks([np.log10(plot_periodlist[ll])*self.ystretch 
                             for ll in np.arange(0, n, self.ystep)])
        
        #set y-axis minor ticks                     
        self.ax.yaxis.set_ticks([np.log10(plot_periodlist[ll])*self.ystretch 
                             for ll in np.arange(0, n, 1)],minor=True)
        #set y-axis tick labels
        self.ax.set_yticklabels(yticklabels)
        
        #set x-axis ticks
        self.ax.set_xticks(self.offsetlist*self.xstretch)
        
        #set x-axis tick labels as station names
        xticklabels = self.stationlist
        if self.xstep != 1:
            xticklabels = np.zeros(len(self.stationlist), 
                                   dtype=self.stationlist.dtype)
            for xx in range(0,len(self.stationlist),self.xstep):
                xticklabels[xx] = self.stationlist[xx]
        self.ax.set_xticklabels(xticklabels)
        
        #--> set x-limits
        if self.xlimits == None:
            self.ax.set_xlim(self.offsetlist.min()*self.xstretch-es*2,
                             self.offsetlist.max()*self.xstretch+es*2)
        else:
            self.ax.set_xlim(self.xlimits)
            
        #--> set y-limits
        if self.ylimits == None:
            self.ax.set_ylim(pmax+es*2, pmin-es*2)
        else:
            pmin = np.log10(self.ylimits[0])*self.ystretch
            pmax = np.log10(self.ylimits[1])*self.ystretch
            self.ax.set_ylim(pmax+es*2, pmin-es*2)
            
        #--> set title of the plot
        if self.plot_title == None:
            pass
        else:
            self.ax.set_title(self.plot_title, fontsize=self.font_size+2)
        
        #put a grid on the plot
        self.ax.grid(alpha=.25, which='both', color=(.25, .25, .25))
        
        #print out the min an max of the parameter plotted
        print '-'*25
        print ck+' min = {0:.2f}'.format(min(minlist))
        print ck+' max = {0:.2f}'.format(max(maxlist))
        print '-'*25

        #==> make a colorbar with appropriate colors
        if self.cb_position == None:
            self.ax2, kw = mcb.make_axes(self.ax,
                                         orientation=self.cb_orientation,
                                         shrink=.35)
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)
        
        if cmap == 'mt_seg_bl2wh2rd':
            #make a color list
            self.clist = [(cc, cc, 1) 
                         for cc in np.arange(0, 1+1./(nseg), 1./(nseg))]+\
                        [(1, cc, cc) 
                         for cc in np.arange(1, -1./(nseg), -1./(nseg))]
            
            #make segmented colormap
            mt_seg_bl2wh2rd = colors.ListedColormap(self.clist)

            #make bounds so that the middle is white
            bounds = np.arange(ckmin-ckstep, ckmax+2*ckstep, ckstep)
            
            #normalize the colors
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)
            
            #make the colorbar
            self.cb = mcb.ColorbarBase(self.ax2,
                                       cmap=mt_seg_bl2wh2rd,
                                       norm=norms,
                                       orientation=self.cb_orientation,
                                       ticks=bounds[1:-1])
        else:
            self.cb = mcb.ColorbarBase(self.ax2,
                                       cmap=mtcl.cmapdict[cmap],
                                       norm=colors.Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                       orientation=self.cb_orientation)

        #label the color bar accordingly
        self.cb.set_label(mtpl.ckdict[ck],
                          fontdict={'size':self.font_size,'weight':'bold'})
            
        #place the label in the correct location                   
        if self.cb_orientation == 'horizontal':
            self.cb.ax.xaxis.set_label_position('top')
            self.cb.ax.xaxis.set_label_coords(.5, 1.3)
            
            
        elif self.cb_orientation == 'vertical':
            self.cb.ax.yaxis.set_label_position('right')
            self.cb.ax.yaxis.set_label_coords(1.5, .5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis='y', direction='in')
        
        plt.show()

    def save_figure(self, save_fn, file_format='pdf', 
                  orientation='portrait', fig_dpi=None, close_plot='y'):
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
            >>> p1 = mtplot.PlotPhaseTensorMaps(edilist,freqspot=10)
            >>> p1.save_plot(r'/home/MT', file_format='jpg')
            'Figure saved to /home/MT/PTMaps/PTmap_phimin_10Hz.jpg'
            
        """

        sf='_{0:.6g}'.format(self.plot_freq)
        
        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation,  bbox_inches='tight')
            plt.clf()
            plt.close(self.fig)
            
        else:
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)
            if not os.path.exists(os.path.join(save_fn, 'PTMaps')):
                os.mkdir(os.path.join(save_fn, 'PTMaps'))
                save_fn = os.path.join(save_fn, 'PTMaps')
                
            save_fn = os.path.join(save_fn, 'PTmap_'+self.ellipse_colorby+sf+
                                    'Hz.'+file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')
                        
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