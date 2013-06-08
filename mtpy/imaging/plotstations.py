# -*- coding: utf-8 -*-
"""
Created on Fri Jun 07 18:20:00 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import mtpy.imaging.mtplottools as mtpt
import mtpy.utils.exceptions as mtex

#==============================================================================

class PlotStations(object):
    """
    plot station locations in map view.
    
    
    """
    
    def __init__(self, **kwargs):
        
        fn_lst = kwargs.pop('fn_lst', None)
        mt_object_lst = kwargs.pop('mt_object_lst', None)
        
        #----set attributes for the class-------------------------
        self.mt_lst = mtpt.MTplot_lst(fn_lst=fn_lst,  
                                      mt_object_lst=mt_object_lst)
        
            
            
        #--> set plot properties
        self.fig_num = kwargs.pop('fig_num', 1)
        self.plot_num = kwargs.pop('plot_num', 1)
        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_size = kwargs.pop('fig_size', None)
        self.font_size = kwargs.pop('font_size', 7)
        self.stationid = kwargs.pop('stationid', [0,4])
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)
        
        self.map_scale = kwargs.pop('map_scale', 'latlon')
        self.marker = kwargs.pop('marker', 'v')
        self.marker_size = kwargs.pop('marker_size', 10)
        self.marker_color = kwargs.pop('marker_color', 'k')
        self.plot_names = kwargs.pop('plot_names', True)

        self.text_size = kwargs.pop('text_size', 7)
        self.text_weight = kwargs.pop('text_weight', 'normal')
        self.text_color = kwargs.pop('text_color', 'k')
        self.text_ha = kwargs.pop('text_ha', 'center')
        self.text_va = kwargs.pop('text_va', 'top')
        self.text_angle = kwargs.pop('text_angle', 0)
        self.text_pad = kwargs.pop('text_pad', None)
        
        self.image_file = kwargs.pop('image_file', None)
        self.image_extent = kwargs.pop('image_extent', None)
        
        if self.image_file is not None:
            if self.image_extent is None:
                raise mtex.MTpyError_inputarguments('Need to input extents '+\
                                                    'of the image as' +\
                                                    '(x0, y0, x1, y1)')
        
        #--> plot if desired
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()
            
    def plot(self):
        """
        plots the station locations
        
        """
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .09
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .09
        plt.rcParams['figure.subplot.top'] = .98
        
        #get station locations
        self.mt_lst.get_station_locations(map_scale=self.map_scale)
        
        text_dict = {'size':self.text_size,
                     'weight':self.text_weight,
                     'rotation':self.text_angle,
                     'color':self.text_color}
                     
        font_dict = {'size':self.font_size+2, 'weight':'bold'}
                     
        if self.xlimits is None:
            if np.sign(self.mt_lst.map_xarr.min()) == -1:
                self.xlimits = (self.mt_lst.map_xarr.min()*1.002, 
                                self.mt_lst.map_xarr.max()*.998)
            else:   
                self.xlimits = (self.mt_lst.map_xarr.min()*.998, 
                                self.mt_lst.map_xarr.max()*1.002)
                            
        if self.ylimits is None:
            if np.sign(self.mt_lst.map_yarr.min()) == -1:
                self.ylimits = (self.mt_lst.map_yarr.min()*1.002, 
                                self.mt_lst.map_yarr.max()*.998)
            else:   
                self.ylimits = (self.mt_lst.map_yarr.min()*.998, 
                                self.mt_lst.map_yarr.max()*1.002)
                            
        if self.map_scale == 'latlon':
            xlabel = 'Longitude (deg)'
            ylabel = 'Latitude (deg)'
            
        elif self.map_scale == 'eastnorth':
            xlabel = 'Easting (m)'
            ylabel = 'Northing (m)'
            
        elif self.map_scale == 'eastnorthkm':
            xlabel = 'Easting (km)'
            ylabel = 'Northing (km)'
                            
        #make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        
        #add and axes
        self.ax = self.fig.add_subplot(1, 1, 1, aspect='equal')
        
        #--> plot the background image if desired-----------------------
        if self.image_file is not None:
            im=plt.imread(self.image_file)
            self.ax.imshow(im, origin='lower', extent=self.image_extent, 
                      aspect='auto')
        
        
        for key in self.mt_lst.map_dict.keys():
            self.ax.scatter(self.mt_lst.map_dict[key][0],
                            self.mt_lst.map_dict[key][1],
                            marker=self.marker,
                            c=self.marker_color,
                            s=self.marker_size)
                            
            if self.plot_names == True:
                if self.text_pad is None:
                    self.text_pad = .0015*self.mt_lst.map_dict[key][1]
                    
                self.ax.text(self.mt_lst.map_dict[key][0],
                             self.mt_lst.map_dict[key][1]+self.text_pad*\
                                 np.sign(self.mt_lst.map_dict[key][1]),
                             key[self.stationid[0]:self.stationid[1]], 
                             verticalalignment=self.text_va,
                             horizontalalignment=self.text_ha,
                             fontdict=text_dict)
                             
        #set axis properties
        self.ax.set_xlabel(xlabel, fontdict=font_dict)
        self.ax.set_ylabel(ylabel, fontdict=font_dict)
        self.ax.grid(alpha=.35, color=(.25, .25, .25))
        self.ax.set_xlim(self.xlimits)
        self.ax.set_ylim(self.ylimits)
        
        plt.show()
                             
    def save_plot(self, save_fn, file_format='pdf', 
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
            >>> p1 = mtplot.PlotPhaseTensorMaps(edilst,freqspot=10)
            >>> p1.save_plot(r'/home/MT', file_format='jpg')
            'Figure saved to /home/MT/PTMaps/PTmap_phimin_10Hz.jpg'
            
        """

        sf='_{0:.6g}'.format(self.plot_freq)
        
        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
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
        return "Plot station locations"              
                             
                        
        
        