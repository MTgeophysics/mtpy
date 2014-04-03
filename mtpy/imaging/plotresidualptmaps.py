# -*- coding: utf-8 -*-
"""
PlotResidualPhaseTensor
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

class PlotResidualPTMaps(mtpl.MTEllipse):
    """
    plot residual phase tensor maps
    
     Arguments:
    ----------
    
        **fn_list1** : list of strings
                        full paths to .edi files for survey 1

        **fn_list2** : list of strings
                        full paths to .edi files for survey 2
                        
        **Note** it is assumed that the edi file lists have the same number
                 of edi files and have the same station names.  If you get
                 an error check this first.
                          
        **plot_freq** : float
                             freq to plot in Hz
                             *default* is 1
                             
        **ftol** : float
                   tolerance in freq range to look for in each file.
                   *default* is 0.1 (10 percent)
                             
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
        
        **xpad** : float
                   padding in the east-west direction of plot boundaries.  Note
                   this is by default set to lat and long units, so if you use
                   easting/northing put it the respective units. 
                   *default* is 0.2
                   
        **pad** : float
                   padding in the north-south direction of plot boundaries.  
                   Note this is by default set to lat and long units, so if you
                   use easting/northing put it the respective units.
                   *default* is 0.2
        
        **rotz** : float
                   angle in degrees to rotate the data clockwise positive.
                   *Default* is 0

        **figsize** : tuple or list (x, y) in inches
                      dimensions of the figure box in inches, this is a default
                      unit of matplotlib.  You can use this so make the plot
                      fit the figure box to minimize spaces from the plot axes
                      to the figure box.  *default* is [8, 8]
                      
        **station_dict** : dictionary
                            * 'id' --> for station id index.  Ex: If you want 
                                      'S01' from 'S01dr' input as (0,3).
                                      
                            * 'pad' --> pad from the center of the ellipse to 
                                       the station label
                                       
                            * 'font_dict'--> dictionary of font properties
                                           font dictionary for station name. 
                                           Keys can be matplotlib.text 
                                           properties, common ones are:
                                           * 'size'   -> for font size
                                           * 'weight' -> for font weight
                                           * 'color'  -> for color of font
                                           * 'angle'  -> for angle of text
                               
        **tscale** : [ 'period' | 'freq' ]
        
                     * 'period'    -> plot vertical scale in period
                     
                     * 'freq' -> plot vertical scale in freq
        
        **map_scale** : [ 'deg' | 'm' | 'km' ]
                       Scale of the map coordinates.
                       
                       * 'deg' --> degrees in latitude and longitude
                       
                       * 'm' --> meters for easting and northing
                       
                       * 'km' --> kilometers for easting and northing
             
        **image_dict** : dictionary of image properties
        
                         * 'file' : string
                                   full path to image file name
                                   
                         * 'extent' : tuple (xmin, xmax, ymin, ymax)
                                     coordinates according to map_scale. Must be
                                     input if image file is not None.
         
        **plot_yn** : [ 'y' | 'n' ]
                      *'y' to plot on creating an instance
                      
                      *'n' to not plot on creating an instance           
                       
        **fignum** : int
                     figure number.  *Default* is 1
                     
        **title** : string
                    figure title
                    
        **dpi** : int 
                  dots per inch of the resolution. *default* is 300
                         

        **font_size** : float
                        size of the font that labels the plot, 2 will be added
                        to this number for the axis labels.
                        

        **station_dict** : dictionary
                           font dictionary for station name. Keys can be
                           matplotlib.text properties, common ones are:
                               * 'size'   -> for font size
                               * 'weight' -> for font weight
                               * 'color'  -> for color of font
        

        
        **reference_point** : tuple (x0,y0)
                              reference point estimate relative distance to.  
                              This point will be (0,0) on the map and 
                              everything else is referenced to this point
        
        **rot90** : [ 0 | 90 | 180 | 270 ]
                    rotate the ellipses by this amout, *default* is 90
                    angle of long axis is rot90-azimuth
         
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
        >>> ptmap = mtplot.plot_residual_pt_maps(edilist1, edilist2, freqspot=10,
        >>> ...                                  ellipse_dict={'size':1,
        >>> ...                                              'range':(0,5)})
        >>> 
        >>> #
        >>> #---add an image---
        >>> ptmap.image_file = r"/home/Maps/Basemap.jpg"
        >>> ptmap.image_extent = (0,10,0,10)
        >>> ptmap.redraw_plot()
        >>> #
        >>> #---Save the plot---
        >>> ptmap.save_plot(r"/home/EDIfiles",file_format='pdf')
        >>> 'Saved figure to /home/EDIfile/PTMaps/PTmap_phimin_10.0_Hz.pdf'
        
    :Example: ::
        
        >>> #change the axis label and grid color
        >>> ptmap.ax.set_xlabel('Latitude (deg)')
        >>> ptmap.ax.grid(which='major', color=(.5,1,0))
        >>> ptmap.update_plot()
        
    Attributes:
    -----------
    
        -arrow_color_imag         imaginary induction arrow color
        -arrow_color_real         real induction arrow color
        -arrow_direction          directional convention of arrows
        -arrow_head_length        head length of arrows in relative points
        -arrow_head_width         head width of arrows in relative points
        -arrow_legend_fontdict    font dictionary for arrow legend
        -arrow_legend_fontpad     padding between font of legend and arrow
        -arrow_legend_position    legend position see matplotlib.legend
        -arrow_legend_xborderpad  padding between legend and x axis
        -arrow_legend_yborderpad  padding between legend and y axis
        -arrow_lw                 arrow line width
        -arrow_size               scaling factor to make arrows visible
        -arrow_threshold          threshold for plotting arrows anything above 
                                  this number will not be plotted 
        
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
                              
        -ftol                 tolerance to look for matching freq
        -jj                   index of plot freq

        -map_scale             scale of map
        
        -mt_list               list of mtplot.MTplot instances containing all 
                              the important information for each station
                              
        -plot_freq       freq in Hz to plot
        -plot_reference_point  reference point of map, everything will be 
                               measured relative to this point
        
        -plot_xarr            array of x-coordinates for stations 
        -plot_yarr            array of y-coordinates for stations
        -plot_yn              plot on instance creation
        
        -rot_z                rotates the data by this angle assuming North is
                              0 and angle measures clockwise
                              
        -tickstrfmt           format of tick strings
        -title                title of figure
        -tscale               temporal scale of y-axis ('freq' | 'period')
        -xpad                 padding between furthest station in x-direction 
                              and the axes edge
        -ypad                 padding between furthes station in y-direction 
                              and axes edge
                              
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
        
        self.residual_pt_list = None
        self.plot_data = None
        self.med_filt_kernel = kwargs.pop('med_filt_kernel', None)
        
        
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
        #set some of the properties as attributes much to Lars' discontent
        self.fig_num = kwargs.pop('fig_num', 1)
        self.plot_num = kwargs.pop('plot_num', 1)
        self.plot_style = kwargs.pop('plot_style', 'pseudo')
        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        
        self.tscale = kwargs.pop('tscale', 'period')
        self.fig_size = kwargs.pop('fig_size', [8, 8])
        self.map_scale = kwargs.pop('map_scale', 'deg')
        self.rot90 = kwargs.pop('rot90', 90)
        if self.map_scale == 'deg':        
            self.xpad = kwargs.pop('xpad', .005)
            self.ypad = kwargs.pop('ypad', .005)
            #--> set the ellipse properties -------------------
            self._ellipse_dict = kwargs.pop('ellipse_dict', 
                                            {'range':(0,5), 'cmap':'mt_wh2or',
                                             'size':.005,
                                             'colorby':'geometric_mean'}) 
            self._read_ellipse_dict()
        elif self.map_scale == 'm':        
            self.xpad = kwargs.pop('xpad', 1000)
            self.ypad = kwargs.pop('ypad', 1000)
            #--> set the ellipse properties -------------------
            self._ellipse_dict = kwargs.pop('ellipse_dict', 
                                            {'range':(0,5), 'cmap':'mt_wh2or',
                                             'size':500,
                                             'colorby':'geometric_mean'})
            self._read_ellipse_dict()
        elif self.map_scale == 'km':        
            self.xpad = kwargs.pop('xpad', 1)
            self.ypad = kwargs.pop('ypad', 1)
            #--> set the ellipse properties -------------------
            self._ellipse_dict = kwargs.pop('ellipse_dict', 
                                            {'range':(0,5), 'cmap':'mt_wh2or',
                                             'size':.5,
                                             'colorby':'geometric_mean'})
            self._read_ellipse_dict()
        self.font_size = kwargs.pop('font_size', 7)
        
        #--> set the freq to plot
        self.plot_freq = kwargs.pop('plot_freq', 1.0)
        self.ftol = kwargs.pop('ftol', .1)
        
        #--> set spacing of plot
        self.subplot_wspace = .1
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1
        
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
        
        #--> set background image properties
        image_dict = kwargs.pop('image_dict', None)
        if image_dict!=None:
            # make sure there is a file
            try:
                self.image_file = image_dict['file']
                if not os.path.isfile(image_dict['file']):
                    raise IOError('Image file does not exist')
                    
            except KeyError:
                raise IOError('Need to include filename if plotting an image')
              
            # make sure an extent is given
            try:
                self.image_extent = image_dict['extent']
            except KeyError:
                raise NameError('Need to include the extent of the image as '+\
                                '(left, right, bottom, top)')
                                
        #--> set a central reference point
        self.plot_reference_point = kwargs.pop('reference_point', (0, 0))
            
        #--> set station name properties
        station_dict = kwargs.pop('station_dict', None)
        if station_dict!=None:
            try:
                self.station_id = station_dict['id']
            except KeyError:
                self.station_id = (0,2)
            
            #set spacing of station name and ellipse
            try:
                self.station_pad = station_dict['pad']
            except KeyError:
                self.station_pad = .0005
                
            #set font properties of the station label
            try:
                self.station_font_size = station_dict['font_dict']
            except KeyError:
                self.station_font_dict = {'size':self.font_size,
                                          'weight':'bold'}        
        
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
            
        for ii,mt in enumerate(self.mt_list):
            mt.rot_z = self._rot_z[ii]
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

                    fdict1 = dict([(np.round(ff, 5), ii) 
                                    for ii, ff in enumerate(mt1.freq)])
                    
                    fdict2 = dict([(np.round(ff, 5), jj) 
                                    for jj, ff in enumerate(mt2.freq)])
                    f1_index_list = []
                    f2_index_list = []
                    for ff in sorted(fdict1.keys()):
                        try:
                            f2_index_list.append(fdict2[ff])
                            f1_index_list.append(fdict1[ff])
                        except KeyError:
                            pass

                    try:
                        new_Z1 = mtpl.mtz.Z(z_array=mt1.z[f1_index_list],
                                            zerr_array=mt1.z_err[f1_index_list],
                                            freq=mt1.freq[f1_index_list])
                        new_Z2 = mtpl.mtz.Z(z_array=mt2.z[f2_index_list],
                                            zerr_array=mt2.z_err[f2_index_list],
                                            freq=mt2.freq[f2_index_list])
                       
                        pt1 = mtpt.PhaseTensor(z_object=new_Z1)
                        pt2 = mtpt.PhaseTensor(z_object=new_Z2)

                        rpt = mtpt.ResidualPhaseTensor(pt1, pt2)
                        rpt.compute_residual_pt(pt1, pt2)
                        rpt.station = mt1.station
                        rpt.lat = mt1.lat
                        rpt.lon = mt1.lon
                        rpt.freq = new_Z1.freq
                        self.residual_pt_list.append(rpt)
                        station_find = True
                        break
                    except IndexError:
                        for key in sorted(fdict2.keys()):
                            print key, fdict2[key]
                        
                else:
                    pass
            if station_find == False:
                print 'Did not find {0} from list 1 in list 2'.format(mt1.station) 
                 
    def _get_plot_freq_data(self):
        """
        get the data to plot in the form of arrays
        
        """
        
        if self.residual_pt_list is None:
            self._compute_residual_pt()
            
        data_count = 0
        d_index_list = []
        for r_index, rpt in enumerate(self.residual_pt_list):
            try:
                dd = np.where(rpt.freq == self.plot_freq)[0][0]
                data_count += 1
                d_index_list.append([data_count, r_index, dd])
            except IndexError:
                try:
                    dd = np.where((rpt.freq >= self.plot_freq*(1-self.ftol)) &
                             (rpt.freq <= self.plot_freq*(1+self.ftol)))[0][0]
                    data_count += 1
                    d_index_list.append([data_count-1, r_index, dd])
                except IndexError:
                    pass
            
        if data_count == 0:
            raise mtpl.mtex.MTpyError_value('No data at {0:.5e} Hz'.format(
                                            self.plot_freq))

        self.plot_data = np.zeros(data_count, dtype=[('phimin', np.float),
                                                    ('phimax', np.float),
                                                    ('skew', np.float),
                                                    ('azimuth', np.float),
                                                    ('ellipticity', np.float),
                                                    ('station', '|S20'),
                                                    ('lat', np.float),
                                                    ('lon', np.float),
                                                    ('geometric_mean', np.float)])
        
        for ii, r_index, f_index in d_index_list:
            
            rpt = self.residual_pt_list[r_index]

            self.plot_data[ii]['phimin'] = rpt.residual_pt.phimin[0][f_index]
            self.plot_data[ii]['phimax'] = rpt.residual_pt.phimax[0][f_index]
            self.plot_data[ii]['skew'] = rpt.residual_pt.beta[0][f_index]
            self.plot_data[ii]['azimuth'] = rpt.residual_pt.azimuth[0][f_index]
            self.plot_data[ii]['ellipticity'] = rpt.residual_pt.ellipticity[0][f_index]
            self.plot_data[ii]['station'] = rpt.station
            self.plot_data[ii]['lat'] = rpt.lat
            self.plot_data[ii]['lon'] = rpt.lon
            self.plot_data[ii]['geometric_mean'] = np.sqrt(abs(
                                          rpt.residual_pt.phimin[0][f_index]*\
                                          rpt.residual_pt.phimax[0][f_index]))

            
                     
    def _apply_median_filter(self):
        """
        apply a median filter to the data to remove extreme outliers
        
        kernel is (station, frequency)
        
        """
        
        for key in ['phimin', 'phimax', 'skew', 'azimuth', 'ellipticity']:
            self.plot_data[key] = sps.medfilt(self.plot_data[key], 
                                              self.med_filt_kernel[0])
            
      
    def plot(self):
        """
        plot residual phase tensor
        """                            
        #get residual phase tensor for plotting        
        self._get_plot_freq_data()
        
        #filter data if desired
        if self.med_filt_kernel is not None:
            self._apply_median_filter()
        
        
        #set position properties for the plot
        plt.rcParams['font.size']=self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace
        
        #make figure instance
        if self.tscale == 'period':
            titlefreq = '{0:.5g} (s)'.format(1./self.plot_freq)
        else:
            titlefreq='{0:.5g} (Hz)'.format(self.plot_freq)
            
        self.fig = plt.figure(titlefreq,
                              self.fig_size, dpi=self.fig_dpi)
        
        #clear the figure if there is already one up
        plt.clf()
        
        #make an axes instance
        self.ax = self.fig.add_subplot(1, 1, 1, aspect='equal')
        
        #--> plot the background image if desired-----------------------
        try:
            im = plt.imread(self.image_file)
            self.ax.imshow(im, origin='lower', extent=self.image_extent, 
                           aspect='auto')
        except AttributeError:
            pass
        
        #get the reference point
        refpoint = self.plot_reference_point
            
        #set some local parameters
        es = float(self.ellipse_size)
        cmap = self.ellipse_cmap
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        nseg = float((ckmax-ckmin)/(2*ckstep))
        ck = self.ellipse_colorby


        #--> set the bounds on the segmented colormap
        if cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
            nseg = float((ckmax-ckmin)/(2*ckstep)) 
            
        #set tick parameters depending on the map_scale
        if self.map_scale == 'deg':
            self.tickstrfmt = '%.3f'
            
        elif self.map_scale == 'm' or self.map_scale == 'km':
            self.tickstrfmt = '%.0f'
        
        #make some empty arrays
        elliplist=[]
        self.plot_xarr = np.zeros(len(self.plot_data))
        self.plot_yarr = np.zeros(len(self.plot_data))

        #--> get size of largest ellipse for this frequency for 
        #    normalization to give an indication of the size of 
        #    change.
        emax = self.plot_data[self.ellipse_colorby].max()
        for ii, rpt in enumerate(self.plot_data):
            #if map scale is lat lon set parameters                
            if self.map_scale == 'deg':
                plotx = rpt['lon']-refpoint[0]
                ploty = rpt['lat']-refpoint[1]
            
            #if map scale is in meters easting and northing
            elif self.map_scale == 'm':
                zone, east, north = utm2ll.LLtoUTM(23, rpt['lat'], 
                                                   rpt['lon'])
                
                #set the first point read in as a refernce other points                    
                if ii == 0:
                    zone1 = zone
                    plotx = east-refpoint[0]
                    ploty = north-refpoint[1]
                    
                #read in all the other point
                else:
                    #check to make sure the zone is the same this needs
                    #to be more rigorously done
                    if zone1 != zone:
                        print 'Zone change at station {0}'.format(
                                                            rpt['station'])
                        if zone1[0:2] == zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east += 500000
                        else:
                            east -= -500000
                        plotx = east-refpoint[0]
                        ploty = north-refpoint[1]
                    else:
                        plotx = east-refpoint[0]
                        ploty = north-refpoint[1]
                
            #if map_scale is in km easting and northing
            elif self.map_scale == 'km':
                zone, east, north = utm2ll.LLtoUTM(23, rpt['lat'], 
                                                       rpt['lon'])
                if ii == 0:
                    zone1 = zone
                    plotx = (east-refpoint[0])/1000.
                    ploty = (north-refpoint[1])/1000.
                
                else:
                    if zone1 != zone:
                        print 'Zone change at station {0}'.format(
                                                            rpt['station'])
                        if zone1[0:2] == zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east += 500000
                        else:
                            east -= 500000

                        plotx = (east-refpoint[0])/1000.
                        ploty = (north-refpoint[1])/1000.
                    else:
                        plotx = (east-refpoint[0])/1000.
                        ploty = (north-refpoint[1])/1000.
                
            #put the location of each ellipse into an array in x and y
            self.plot_xarr[ii] = plotx
            self.plot_yarr[ii] = ploty

            #--> get ellipse properties
            #if the ellipse size is not physically correct make it a dot
            if rpt['phimax'] == 0 or rpt['phimax'] > 100 or\
               rpt['phimin'] == 0 or rpt['phimin'] > 100:
                eheight=.0000001*es
                ewidth=.0000001*es
                print 'Bad data at {0}'.format(rpt['station'])
            else:
                scaling = es/emax
                eheight = rpt['phimin']*scaling
                ewidth = rpt['phimax']*scaling
            
            #make an ellipse
            ellipd=patches.Ellipse((plotx,ploty),
                                   width=ewidth,
                                   height=eheight,
                                   angle=self.rot90-rpt['azimuth'])
                                   
            #get ellipse color
            if cmap.find('seg')>0:
                ellipd.set_facecolor(mtcl.get_plot_color(rpt[self.ellipse_colorby],
                                                         self.ellipse_colorby,
                                                         cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
            else:
                ellipd.set_facecolor(mtcl.get_plot_color(rpt[self.ellipse_colorby],
                                                         self.ellipse_colorby,
                                                         cmap,
                                                         ckmin,
                                                         ckmax))
            
            #==> add ellipse to the plot
            elliplist.append(ellipd)
            self.ax.add_artist(ellipd)
                    
            #------------Plot station name------------------------------
            try:
                self.ax.text(plotx,
                        ploty+self.station_pad,
                        rpt['station'][self.station_id[0]:self.station_id[1]],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=self.station_font_dict)
            except AttributeError:
                pass
        
        #--> set axes properties depending on map scale------------------------
        if self.map_scale == 'deg':    
            self.ax.set_xlabel('Longitude',
                               fontsize=self.font_size+2,
                               fontweight='bold')
            self.ax.set_ylabel('Latitude',
                               fontsize=self.font_size+2,
                               fontweight='bold')
            
        elif self.map_scale == 'm':
            self.ax.set_xlabel('Easting (m)',
                               fontsize=self.font_size+2,
                               fontweight='bold')
            self.ax.set_ylabel('Northing (m)',
                               fontsize=self.font_size+2,
                               fontweight='bold')
      
        elif self.map_scale == 'km':
            self.ax.set_xlabel('Easting (km)',
                               fontsize=self.font_size+2,
                               fontweight='bold')
            self.ax.set_ylabel('Northing (km)',
                               fontsize=self.font_size+2,
                               fontweight='bold')

        
        #--> set plot limits
        self.ax.set_xlim(self.plot_xarr.min()-self.xpad,
                             self.plot_xarr.max()+self.xpad)
        self.ax.set_ylim(self.plot_yarr.min()-self.xpad,
                         self.plot_yarr.max()+self.xpad)
                         
        #--> set tick label format
        self.ax.xaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        
        #--> set title in period or freq
        if self.tscale == 'period':
            titlefreq = '{0:.5g} (s)'.format(1./self.plot_freq)
        else:
            titlefreq='{0:.5g} (Hz)'.format(self.plot_freq)
        
        if not self.plot_title:
            self.ax.set_title('Phase Tensor Map for {0}'.format(titlefreq),
                              fontsize=self.font_size+2,fontweight='bold')
        else:
            self.ax.set_title(self.plot_title+titlefreq,
                              fontsize=self.font_size+2,fontweight='bold')
                              
        #make a grid with gray lines
        self.ax.grid(alpha=.25)
        
        #==> make a colorbar with appropriate colors
        if self.cb_position == None:
            self.ax2, kw = mcb.make_axes(self.ax,
                                         orientation=self.cb_orientation,
                                         shrink=.35)
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)
        
        if cmap == 'mt_seg_bl2wh2rd':
            #make a color list
            self.clist = [(cc, cc ,1) 
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
            self.cb=mcb.ColorbarBase(self.ax2,
                                     cmap=mt_seg_bl2wh2rd,
                                     norm=norms,
                                     orientation=self.cb_orientation,
                                     ticks=bounds[1:-1])
        else:
            self.cb=mcb.ColorbarBase(self.ax2,
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
            self.cb.ax.yaxis.set_label_coords(1.25, .5)
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