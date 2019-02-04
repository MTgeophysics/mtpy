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
import mtpy.utils.gis_tools as gis_tools
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
    This will plot residual phase tensors in a map for a single frequency. 
    The data is read in and stored in 2 ways, one as a list ResidualPhaseTensor
    object for each matching station and the other in a structured array with 
    all the important information.  The structured array is the one that is
    used for plotting.  It is computed each time plot() is called so if it is
    manipulated it is reset.  The array is sorted by relative offset, so no
    special order of input is needed for the file names.  However, the 
    station names should be verbatim between surveys, otherwise it will not
    work.  
    
    The residual phase tensor is calculated as I-(Phi_2)^-1 (Phi_1)
    
    The default coloring is by the geometric mean as sqrt(Phi_min*Phi_max), 
    which defines the percent change between measurements.
    
    There are a lot of parameters to change how the plot looks, have a look 
    below if you figure looks a little funny.  The most useful will be 
    ellipse_size
    
    The ellipses are normalized by the largest Phi_max of the survey.
    
     Arguments:
    --------------
        **fn_list1** : list of strings
                        full paths to .edi files for survey 1

        **fn_list2** : list of strings
                        full paths to .edi files for survey 2
                        
        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor 
                          ellipses with keys:
                              * 'size' -> size of ellipse in points 
                                         *default* is 2
                              
                              * 'colorby' : [ 'phimin' | 'phimax' | 'skew' | 
                                              'skew_seg' | 'phidet' | 
                                              'ellipticity' | 'geometric_mean']
                                    
                                    - 'phimin' -> colors by minimum phase
                                    - 'phimax' -> colors by maximum phase
                                    - 'skew' -> colors by beta (skew)
                                    - 'skew_seg' -> colors by beta in 
                                                   discrete segments 
                                                   defined by the range
                                    - 'phidet' -> colors by determinant of
                                                 the phase tensor
                                    - 'ellipticity' -> colors by ellipticity
                                    - 'geometric_mean' -> sqrt(phimin*phimax)
                                    *default* is 'geometric_mean'
                                
                               * 'range' : tuple (min, max, step)
                                     Need to input at least the min and max
                                     and if using 'skew_seg' to plot
                                     discrete values input step as well
                                     *default* depends on 'colorby'
                                     
                          * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 
                                      'mt_wh2bl' | 'mt_rd2bl' | 
                                      'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                                      'mt_rd2gr2bl' | 'mt_wh2or ]
                                      
                                   - 'mt_yl2rd' -> yellow to red
                                   - 'mt_bl2yl2rd' -> blue to yellow to red
                                   - 'mt_wh2bl' -> white to blue
                                   - 'mt_rd2bl' -> red to blue
                                   - 'mt_bl2wh2rd' -> blue to white to red
                                   - 'mt_bl2gr2rd' -> blue to green to red
                                   - 'mt_rd2gr2bl' -> red to green to blue
                                   - 'mt_seg_bl2wh2rd' -> discrete blue to 
                                                         white to red
                                   - 'mt_wh2or' -> white to orange
                                                    *default*
        **ellipse_scale** : float
                            value to which all ellipses are scaled.
                            *default* is None, which means all ellipses
                            will be scaled by the maximum value of Phi_max
                            of the full set of data.
                            
        **rot90** : [ True | False ] True to rotate residual phase tensor
                    by 90 degrees. False to leave as is.


        
        **med_filt_kernel** : tuple(station, period)
                              kernel size for the 2D median filter.  
                              the first is the number of stations to smooth 
                              over. The second number is the number of periods 
                              to smooth over. Both should be odd.
                              *default* is None
        
        
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
                    
                       
        **fig_num** : int
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

        **fig_size** : tuple or list (x, y) in inches
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
         
        **plot_yn** : [ 'y' | 'n' ]
                      *'y' to plot on creating an instance
                      
                      *'n' to not plot on creating an instance           
        
    ==================== ======================================================
      Attributes          Description
    ==================== ======================================================
     ax                   matplotlib.axes instance for the main plot
     ax2                  matplotlib.axes instance for the color bar
     cb                   matplotlib.colors.ColorBar instance for color bar
     cb_orientation       color bar orientation ('vertical' | 'horizontal')
     cb_position          color bar position (x, y, dx, dy)
     ellipse_cmap         ellipse color map, see above for options
     ellipse_colorby      parameter to color ellipse by
     ellipse_range        (min, max, step) values to color ellipses
     ellipse_scale        value to which all ellipses are normalized to.     
     ellipse_size         scaling factor to make ellipses visible
     fig                  matplotlib.figure instance for the figure  
     fig_dpi              dots-per-inch resolution
     fig_num              number of figure being plotted
     fig_size             size of figure in inches
     fn_list1             list of .edi file names for survey 1
     fn_list2             list of .edi file names for survey 2
     font_size            font size of axes tick label, axes labels will be
                          font_size + 2
     freq_list            list of frequencies from all .edi files
     ftol                 tolerance to search fro frequencies *default* is 0.05
     map_scale            [ 'm' | 'km' | 'deg' ] *default* is 'deg'
     med_filt_kernel      (station, frequency) kernel to apply median smoothing
                          to the data.      
     mt_list1             list of mtplot.MTplot instances containing all
                          important information for each station in survey 1
     mt_list2             list of mtplot.MTplot instances containing all
                          important information for each station in survey 2
     plot_freq            frequency to plot in Hz
     plot_freq_index      index in freq_list where plot_freq is
     plot_reference_point point where everything is refrerenced to, i.e. 
                          center point 
     plot_title           title of the plot
     plot_yn              plot the pseudo section on instance creation
     residual_pt_list     list ofmtpy.pt.ResidualPhaseTensor objects
     rot90                [ True | False ] rotates the residual phase tensors 
                          by 90 degrees if set to True
     rot_z                rotates the impedence tensor of each station by 
                          this amount assuming 0 is N and 90 is E. 
     rpt_array            structured array with all the important information.
                          This is the important array from which plotting 
                          occurs.
     station_font_dict    font dictionary for station labels 
     station_id           index [min, max] to reaad station name
     station_list         list of stations plotted
     station_pad          padding between axis and station label
     subplot_bottom       spacing between plot and bottom of figure window
     subplot_hspace       vertical spacing between subplots
     subplot_left         spacing between plot and left of figure window
     subplot_right        spacing between plot and right of figure window
     subplot_top          spacing between plot and top of figure window
     subplot_wspace       horizontal spacing between subplots
     tscale               temporal scale of y-axis ('frequency' | 'period')
     xpad                 padding of edge of plot from extreme ellipses in 
                          map_scale units
     ypad                 padding of edge of plot from extreme ellipses in 
                          map_scale units                  
    ==================== ======================================================        
    
    
    ======================= ===================================================
    Methods                 Description  
    ======================= ===================================================
     plot                   plots the pseudo section
     redraw_plot            on call redraws the plot from scratch
     save_figure            saves figure to a file of given format
     update_plot            updates the plot while still active
     _apply_median_filter   apply a 2D median filter to the data 
     _compute_residual_pt   compute residual pt and fill rpt_array and 
                            fill residaul_pt_list          
     _get_freq_list         get a list of all possible frequencies from .edi's 
     _get_plot_freq_index   get the index from freq_list to plot
     _read_ellipse_dict     read ellipse dictionary and return a ellipse object
    ======================= ===================================================
      
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

    :Example: ::
        
        >>> # plot seismic hypocenters from a file
        >>> lat, lon, depth = np.loadtxt(r"/home/seismic_hypocenter.txt")
        >>> ptmap.ax.scatter(lon, lat, marker='o')
                      
    """
    
    def __init__(self, fn_list1, fn_list2, **kwargs):
        assert len(fn_list1) == len(fn_list2)
        self.fn_list1 = fn_list1
        self.fn_list2 = fn_list2
        
        self.mt_list1 = mtpl.get_mtlist(fn_list=fn_list1)
        self.mt_list2 = mtpl.get_mtlist(fn_list=fn_list2)
        
        self.residual_pt_list = None
        self.rpt_array = None
        self.med_filt_kernel = kwargs.pop('med_filt_kernel', None)
        
        
        #--> set colorbar properties---------------------------------
        #set orientation to horizontal
        cb_dict = kwargs.pop('cb_dict', {})
        self.cb_orientation = cb_dict.pop('orientation', 'vertical')
        self.cb_position = cb_dict.pop('position', None)

        #--> set plot properties ------------------------------
        #set some of the properties as attributes much to Lars' discontent
        self.plot_title = kwargs.pop('plot_title', None)
        self.plot_station_name = kwargs.pop('plot_station_name', False)
        
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        
        self.tscale = kwargs.pop('tscale', 'period')
        self.map_scale = kwargs.pop('map_scale', 'deg')
        self.rot90 = kwargs.pop('rot90', True)
        
        if self.map_scale == 'deg':        
            self.xpad = kwargs.pop('xpad', .005)
            self.ypad = kwargs.pop('ypad', .005)
            #--> set the ellipse properties -------------------
            self._ellipse_dict = kwargs.pop('ellipse_dict', 
                                            {'range':(0,10), 'cmap':'mt_yl2rd',
                                             'size':.005,
                                             'colorby':'geometric_mean'}) 
            self._read_ellipse_dict()
            self.ellipse_scale = kwargs.pop('ellipse_scale', None)
        elif self.map_scale == 'm':        
            self.xpad = kwargs.pop('xpad', 1000)
            self.ypad = kwargs.pop('ypad', 1000)
            #--> set the ellipse properties -------------------
            self._ellipse_dict = kwargs.pop('ellipse_dict', 
                                            {'range':(0,5), 'cmap':'mt_yl2rd',
                                             'size':500,
                                             'colorby':'geometric_mean'})
            self._read_ellipse_dict()
            self.ellipse_scale = kwargs.pop('ellipse_scale', None)

        elif self.map_scale == 'km':        
            self.xpad = kwargs.pop('xpad', 1)
            self.ypad = kwargs.pop('ypad', 1)
            #--> set the ellipse properties -------------------
            self._ellipse_dict = kwargs.pop('ellipse_dict', 
                                            {'range':(0,5), 'cmap':'mt_yl2rd',
                                             'size':.5,
                                             'colorby':'geometric_mean'})
            self._read_ellipse_dict()
            self.ellipse_scale = kwargs.pop('ellipse_scale', None)


        self.font_size = kwargs.pop('font_size', 7)
        
        #--> set the freq to plot
        self.plot_freq = kwargs.pop('plot_freq', 1.0)
        self.ftol = kwargs.pop('ftol', .05)
        self.plot_freq_index = None
        
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
        station_dict = kwargs.pop('station_dict', {})
        self.station_id = station_dict.pop('id', (0, 2))
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
            
        for ii,mt in enumerate(self.mt_list):
            mt.rot_z = self._rot_z[ii]
    def _get_rot_z(self):
        return self._rot_z
        
    rot_z = property(fget=_get_rot_z, fset=_set_rot_z, 
                     doc="""rotation angle(s)""")
     
    #-------------------------------------------------------------------                 
    def _get_freq_list(self):
        """
        get all possible periods to plot
        
        """
        
        freq_list = []
        for mt1 in self.mt_list1:
            freq_list.extend(mt1.freq)
        for mt2 in self.mt_list2:
            freq_list.extend(mt2.freq)
            
        self.freq_list = np.array(sorted(set(freq_list), reverse=True))
                
    #------------------------------------------------------------------ 
    def _compute_residual_pt(self):
        """
        compute residual phase tensor so the result is something useful to 
        plot
        """
        
        self._get_freq_list()
        freq_dict = dict([(np.round(key, 5), value) 
                           for value, key in enumerate(self.freq_list)])
                               
        num_freq = self.freq_list.shape[0]
        num_station = len(self.mt_list1)
        
        #make a structured array to put stuff into for easier manipulation
        self.rpt_array = np.zeros(num_station, 
                                  dtype=[('station', '|S10'),
                                         ('lat', np.float),
                                         ('lon', np.float),
                                         ('plotx', np.float),
                                         ('ploty', np.float),
                                         ('elev', np.float),
                                         ('phimin', (np.float, num_freq)),
                                         ('phimax', (np.float, num_freq)),
                                         ('skew', (np.float, num_freq)),
                                         ('azimuth', (np.float, num_freq)),
                                         ('geometric_mean', (np.float, num_freq))])
                                         
        self.residual_pt_list = []
        for mm, mt1 in enumerate(self.mt_list1):
            station_find = False
            fdict1 = dict([(np.round(ff, 5), ii) 
                            for ii, ff in enumerate(mt1.freq)])
            for mt2 in self.mt_list2:
                if mt2.station == mt1.station:
                    fdict2 = dict([(np.round(ff, 5), ii) 
                                    for ii, ff in enumerate(mt2.freq)])
                    
                    #need to make sure only matched frequencies are compared
                    index_1 = []
                    index_2 = []
                    for key1 in sorted(fdict1.keys()):
                        try:
                            index_2.append(fdict2[key1])
                            index_1.append(fdict1[key1])
                        except KeyError:
                            'Did not find {0:.4e} Hz in {1}'.format(key1, 
                                                                  mt2.fn)
                                                                  
                    #need to sort the index list, otherwise weird things happen
                    index_1.sort()
                    index_2.sort()
                    
                    #create new Z objects that have similar frequencies                                              
                    new_z1 = mtpl.mtz.Z(z_array=mt1.z[index_1],
                                        z_err_array=mt1.z_err[index_1],
                                        freq=mt1.freq[index_1])
                    new_z2 = mtpl.mtz.Z(z_array=mt2.z[index_2],
                                        z_err_array=mt2.z_err[index_2],
                                        freq=mt2.freq[index_2])
                                        
                    #make new phase tensor objects
                    pt1 = mtpt.PhaseTensor(z_object=new_z1)
                    pt2 = mtpt.PhaseTensor(z_object=new_z2)

                    #compute residual phase tensor
                    rpt = mtpt.ResidualPhaseTensor(pt1, pt2)
                    rpt.compute_residual_pt(pt1, pt2)
                    
                    #add some attributes to residual phase tensor object
                    rpt.station = mt1.station
                    rpt.lat = mt1.lat
                    rpt.lon = mt1.lon
                    
                    #append to list for manipulating later
                    self.residual_pt_list.append(rpt)
                    
                    #be sure to tell the program you found the station
                    station_find = True
                   
                    #put stuff into an array because we cannot set values of 
                    #rpt, need this for filtering.
                    st_1, st_2 = self.station_id
                    self.rpt_array[mm]['station'] = mt1.station[st_1:st_2]
                    self.rpt_array[mm]['lat'] = mt1.lat
                    self.rpt_array[mm]['lon'] = mt1.lon
                    self.rpt_array[mm]['elev'] = mt1.elev
                    
                    rpt_fdict = dict([(np.round(key, 5), value)
                                       for value, key in enumerate(rpt.freq)])
                    for f_index, freq in enumerate(rpt.freq):
                        aa = freq_dict[np.round(freq, 5)]
                        try:
                            try:
                                rr = rpt_fdict[np.round(freq, 5)]
                                
                                self.rpt_array[mm]['phimin'][aa] = \
                                            abs(rpt.residual_pt.phimin[0][rr])
                                self.rpt_array[mm]['phimax'][aa] = \
                                            abs(rpt.residual_pt.phimax[0][rr])
                                self.rpt_array[mm]['skew'][aa] = \
                                                    rpt.residual_pt.beta[0][rr]
                                self.rpt_array[mm]['azimuth'][aa] = \
                                                    rpt.residual_pt.azimuth[0][rr]
                                self.rpt_array[mm]['geometric_mean'][aa] = \
                                            np.sqrt(abs(rpt.residual_pt.phimin[0][rr]*\
                                                    rpt.residual_pt.phimax[0][rr]))
                            except IndexError:
                                print('-'*50)
                                print(mt1.station) 
                                print('freq_index for 1:  {0}'.format(f_index))
                                print('freq looking for:  {0}'.format(freq))
                                print('index in big    :  {0}'.format(aa))
                                print('index in 1      :  {0} '.format(rr))
                                print('len_1 = {0}, len_2 = {1}'.format(
                                        len(mt2.freq), len(mt1.freq)))
                                print('len rpt_freq = {0}'.format(len(rpt.freq)))
                        except KeyError:
                            print('Station {0} does not have {1:.5f}Hz'.format(
                                   mt1.station, freq))
                        
                    
                    break
                else:
                    pass
            if station_find == False:
                print('Did not find {0} from list 1 in list 2'.format(mt1.station))
               
        # from the data get the relative offsets and sort the data by them
        self.rpt_array.sort(order=['lon', 'lat'])
        
        # get relative positions for plotting
        self._get_relative_position()
    
    #-------------------------------------------------------------------
    def _apply_median_filter(self, kernel=(3, 3)):
        """
        apply a median filter to the data to remove extreme outliers
        
        kernel is (station, frequency)
        
        """

                                   
        filt_phimin_arr = sps.medfilt2d(self.rpt_array['phimin'], 
                                        kernel_size=kernel)                            
        filt_phimax_arr = sps.medfilt2d(self.rpt_array['phimax'], 
                                        kernel_size=kernel)                            
        filt_skew_arr = sps.medfilt2d(self.rpt_array['skew'],
                                      kernel_size=kernel)                            
        filt_azimuth_arr = sps.medfilt2d(self.rpt_array['azimuth'], 
                                         kernel_size=kernel) 
        
        self.rpt_array['phimin'] = filt_phimin_arr
        self.rpt_array['phimax'] = filt_phimax_arr
        self.rpt_array['skew'] = filt_skew_arr
        self.rpt_array['azimuth'] = filt_azimuth_arr
        self.rpt_array['geometric_mean'] = np.sqrt(abs(filt_phimin_arr*\
                                                   filt_phimax_arr))
        
        print('Applying Median Filter with kernel {0}'.format(kernel))
    
    #-------------------------------------------------------------------------
    def _get_relative_position(self):
        """
        get the relative positions for each station in the plotting 
        coordinates
        """
        #if map scale is lat lon set parameters
        for ii, rpt in enumerate(self.rpt_array):                
            if self.map_scale == 'deg':
                plotx = rpt['lon']-self.plot_reference_point[0]
                ploty = rpt['lat']-self.plot_reference_point[1]
            
            #if map scale is in meters easting and northing
            elif self.map_scale == 'm':
                east, north, zone = gis_tools.project_point_ll2utm(rpt['lat'],
                                                                   rpt['lon'])
                
                #set the first point read in as a refernce other points                    
                if ii == 0:
                    zone1 = zone
                    plotx = east-self.plot_reference_point[0]
                    ploty = north-self.plot_reference_point[1]
                    
                #read in all the other point
                else:
                    #check to make sure the zone is the same this needs
                    #to be more rigorously done
                    if zone1 != zone:
                        print('Zone change at station {0}'.format(
                                                            rpt['station']))
                        if zone1[0:2] == zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east += 500000
                        else:
                            east -= -500000
                        plotx = east-self.plot_reference_point[0]
                        ploty = north-self.plot_reference_point[1]
                    else:
                        plotx = east-self.plot_reference_point[0]
                        ploty = north-self.plot_reference_point[1]
                
            #if map_scale is in km easting and northing
            elif self.map_scale == 'km':
                east, north, zone = gis_tools.project_point_ll2utm(rpt['lat'],
                                                                   rpt['lon'])
                if ii == 0:
                    zone1 = zone
                    plotx = (east-self.plot_reference_point[0])/1000.
                    ploty = (north-self.plot_reference_point[1])/1000.
                
                else:
                    if zone1 != zone:
                        print('Zone change at station {0}'.format(
                                                            rpt['station']))
                        if zone1[0:2] == zone[0:2]:
                            pass
                        elif int(zone1[0:2])<int(zone[0:2]):
                            east += 500000
                        else:
                            east -= 500000

                        plotx = (east-self.plot_reference_point[0])/1000.
                        ploty = (north-self.plot_reference_point[1])/1000.
                    else:
                        plotx = (east-self.plot_reference_point[0])/1000.
                        ploty = (north-self.plot_reference_point[1])/1000.
                
            #put the location of each ellipse into an array in x and y
            rpt['plotx'] = plotx
            rpt['ploty'] = ploty
            
    def _get_plot_freq_index(self):
        """
        get frequency to plot
        """
        ftol_m = 1-self.ftol
        ftol_p = 1+self.ftol
        try:
            self.plot_freq_index = np.where(self.freq_list == 
                                            self.plot_freq)[0][0]
        except IndexError:
            try:
                self.plot_freq_index = np.where((self.freq_list >= 
                                                 self.plot_freq*ftol_m) &
                                                 (self.freq_list <= 
                                                 self.plot_freq*ftol_p))[0][0]
            except IndexError:
                raise ValueError('could not find {0} Hz'.format(
                                                              self.plot_freq))
       
   #------------------------------------------------------------------------  
    def plot(self):
        """
        plot residual phase tensor
        """                            
        #get residual phase tensor for plotting        
        self._compute_residual_pt()
        
        #filter data if desired
        if self.med_filt_kernel is not None:
          try:
            self._apply_median_filter(kernel=self.med_filt_kernel)
          except:
            print('Warning - Could not apply median filter')
            
        #get frequency index
        self._get_plot_freq_index()

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
            title_freq = '{0:.5g} (s)'.format(1./self.plot_freq)
        else:
            title_freq='{0:.5g} (Hz)'.format(self.plot_freq)
            
        self.fig = plt.figure(title_freq,
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
        ckey = self.ellipse_colorby
        
        f_index = self.plot_freq_index


        #--> set the bounds on the segmented colormap
        if cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
            nseg = float((ckmax-ckmin)/(2*ckstep)) 
            
        #set tick parameters depending on the map_scale
        if self.map_scale == 'deg':
            self.tickstrfmt = '%.3f'
            
        elif self.map_scale == 'm' or self.map_scale == 'km':
            self.tickstrfmt = '%.0f'

        #--> get size of largest ellipse for this frequency for 
        #    normalization to give an indication of the size of 
        #    change.
        if self.ellipse_scale is None:
            emax = self.rpt_array['phimax'].max()
        else:
            emax = self.ellipse_scale
        
        #--> plot        
        for ii, rpt in enumerate(self.rpt_array):
            #--> get ellipse properties
            #if the ellipse size is not physically correct make it a dot
            if rpt['phimax'][f_index] == 0 and \
               rpt['phimax'][f_index] == 0: 
                eheight = .0000001*es
                ewidth = .0000001*es
            
            elif rpt['phimax'][f_index] > 100 or \
               rpt['phimax'][f_index] > 100: 
                eheight = .0000001*es
                ewidth = .0000001*es
                print('Bad data at {0}'.format(rpt['station']))
            
            else:
                scaling = es/emax
                eheight = rpt['phimin'][f_index]*scaling
                ewidth = rpt['phimax'][f_index]*scaling
            
            #make an ellipse
            if self.rot90 == True:
                ellipd = patches.Ellipse((rpt['plotx'],rpt['ploty']),
                                         width=ewidth,
                                         height=eheight,
                                         angle=rpt['azimuth'][f_index]-90)
            elif self.rot90 == False:
                ellipd = patches.Ellipse((rpt['plotx'],rpt['ploty']),
                                         width=ewidth,
                                         height=eheight,
                                         angle=rpt['azimuth'][f_index])
                                   
            #get ellipse color
            if cmap.find('seg')>0:
                ellipd.set_facecolor(mtcl.get_plot_color(rpt[ckey][f_index],
                                                         ckey,
                                                         cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
            else:
                ellipd.set_facecolor(mtcl.get_plot_color(rpt[ckey][f_index],
                                                         ckey,
                                                         cmap,
                                                         ckmin,
                                                         ckmax))
            
            #==> add ellipse to the plot
            self.ax.add_artist(ellipd)
                    
            #------------Plot station name------------------------------
            if self.plot_station_name == True:
                self.ax.text(rpt['plotx'], rpt['ploty']+self.station_pad,
                             rpt['station'],
                             horizontalalignment='center',
                             verticalalignment='baseline',
                             fontdict=self.station_font_dict)

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
        self.ax.set_xlim(self.rpt_array['plotx'].min()-self.xpad,
                             self.rpt_array['plotx'].max()+self.xpad)
        self.ax.set_ylim(self.rpt_array['ploty'].min()-self.xpad,
                         self.rpt_array['ploty'].max()+self.xpad)
                         
        #--> set tick label format
        self.ax.xaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        
        #--> set title in period or freq
        if self.tscale == 'period':
            title_freq = '{0:.5g} (s)'.format(1./self.plot_freq)
        else:
            title_freq='{0:.5g} (Hz)'.format(self.plot_freq)
        
        if not self.plot_title:
            self.ax.set_title('Phase Tensor Map for {0}'.format(title_freq),
                              fontsize=self.font_size+2,fontweight='bold')
        else:
            self.ax.set_title(self.plot_title+title_freq,
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
        self.cb.set_label(mtpl.ckdict[ckey],
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
            
        #--> add reference ellipse
        ref_ellip = patches.Ellipse((0, .0),
                                       width=es,
                                       height=es,
                                       angle=0)
        ref_ellip.set_facecolor((0, 0, 0))
        ref_ax_loc = list(self.ax2.get_position().bounds)
        ref_ax_loc[0] *= .95
        ref_ax_loc[1] -= .17
        ref_ax_loc[2] = .1
        ref_ax_loc[3] = .1
        self.ref_ax = self.fig.add_axes(ref_ax_loc, aspect='equal')
        self.ref_ax.add_artist(ref_ellip)
        self.ref_ax.set_xlim(-es/2.*1.05, es/2.*1.05)
        self.ref_ax.set_ylim(-es/2.*1.05, es/2.*1.05)
        plt.setp(self.ref_ax.xaxis.get_ticklabels(), visible=False)
        plt.setp(self.ref_ax.yaxis.get_ticklabels(), visible=False)
        self.ref_ax.set_title(r'$\Delta \Phi$ = 1')
          
        # put the grid lines behind 
#        [line.set_zorder(10000) for line in self.ax.lines]
        self.ax.set_axisbelow(True)
            
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