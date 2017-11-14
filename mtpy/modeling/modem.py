#!/usr/bin/env python
"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
# ==============================================================================
# Imports
# ==============================================================================
# general packages
import os

import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
# Plotting tools
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# mtpy modules
import mtpy.utils.exceptions as mtex
# vtk tools
from mtpy.modeling.ModEM import Data, Model

try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')


# ==============================================================================


# ==============================================================================
# mesh class
# ==============================================================================
# ==============================================================================
# Residuals
# ==============================================================================
# ==============================================================================
# Control File for inversion
# ==============================================================================
# ==============================================================================
# Control File for inversion
# ==============================================================================
# ==============================================================================
# covariance 
# ==============================================================================
# ==============================================================================
# Write inversion parameters to a config type file
# ==============================================================================
# ==============================================================================
# Manipulate the model to test structures or create a starting model
# ==============================================================================
# ==============================================================================
# plot response       
# ==============================================================================
# ==============================================================================
# plot phase tensors
# ==============================================================================
# ==============================================================================
# plot depth slices
# ==============================================================================
class PlotDepthSlice(object):
    """
    Plots depth slices of resistivity model
    
    :Example: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> sfn = r"/home/MT/ws3dinv/Inv1/WSStationLocations.txt"
        >>> # plot just first layer to check the formating        
        >>> pds = ws.PlotDepthSlice(model_fn=mfn, station_fn=sfn, 
        >>> ...                     depth_index=0, save_plots='n')
        >>> #move color bar up 
        >>> pds.cb_location
        >>> (0.64500000000000002, 0.14999999999999997, 0.3, 0.025)
        >>> pds.cb_location = (.645, .175, .3, .025)
        >>> pds.redraw_plot()
        >>> #looks good now plot all depth slices and save them to a folder
        >>> pds.save_path = r"/home/MT/ws3dinv/Inv1/DepthSlices"
        >>> pds.depth_index = None
        >>> pds.save_plots = 'y'
        >>> pds.redraw_plot()
    
    ======================= ===================================================
    Attributes              Description    
    ======================= ===================================================
    cb_location             location of color bar (x, y, width, height)
                            *default* is None, automatically locates
    cb_orientation          [ 'vertical' | 'horizontal' ] 
                            *default* is horizontal 
    cb_pad                  padding between axes and colorbar
                            *default* is None
    cb_shrink               percentage to shrink colorbar by
                            *default* is None
    climits                 (min, max) of resistivity color on log scale
                            *default* is (0, 4)
    cmap                    name of color map *default* is 'jet_r'
    data_fn                 full path to data file
    depth_index             integer value of depth slice index, shallowest
                            layer is 0
    dscale                  scaling parameter depending on map_scale 
    ew_limits               (min, max) plot limits in e-w direction in 
                            map_scale units. *default* is None, sets viewing
                            area to the station area
    fig_aspect              aspect ratio of plot. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch. *default* is
                            300
    fig_list                list of matplotlib.figure instances for each 
                            depth slice                 
    fig_size                [width, height] in inches of figure size
                            *default* is [6, 6]
    font_size               size of ticklabel font in points, labels are 
                            font_size+2. *default* is 7
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    initial_fn              full path to initial file
    map_scale               [ 'km' | 'm' ] distance units of map. *default* is 
                            km
    mesh_east               np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north              np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                full path to model file
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    ns_limits               (min, max) plot limits in n-s direction in 
                            map_scale units. *default* is None, sets viewing
                            area to the station area
    plot_grid               [ 'y' | 'n' ] 'y' to plot mesh grid lines. 
                            *default* is 'n'
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
    res_model               np.ndarray(n_north, n_east, n_vertical) of 
                            model resistivity values in linear scale
    save_path               path to save figures to
    save_plots              [ 'y' | 'n' ] 'y' to save depth slices to save_path
    station_east            location of stations in east direction in 
                            map_scale units  
    station_fn              full path to station locations file
    station_names           station names
    station_north           location of station in north direction in 
                            map_scale units
    subplot_bottom          distance between axes and bottom of figure window
    subplot_left            distance between axes and left of figure window  
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    title                   titiel of plot *default* is depth of slice
    xminorticks             location of xminorticks
    yminorticks             location of yminorticks
    ======================= ===================================================
    """

    def __init__(self, model_fn=None, data_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.initial_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.initial_fn)

        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

        self.save_plots = kwargs.pop('save_plots', 'y')

        self.depth_index = kwargs.pop('depth_index', None)
        self.map_scale = kwargs.pop('map_scale', 'km')
        # make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)

        self.plot_grid = kwargs.pop('plot_grid', 'n')

        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []

        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)

        self.climits = kwargs.pop('climits', (0, 4))
        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.font_size = kwargs.pop('font_size', 8)

        self.cb_shrink = kwargs.pop('cb_shrink', .8)
        self.cb_pad = kwargs.pop('cb_pad', .01)
        self.cb_orientation = kwargs.pop('cb_orientation', 'horizontal')
        self.cb_location = kwargs.pop('cb_location', None)

        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1

        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None

        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        self.mesh_east = None
        self.mesh_north = None

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def read_files(self):
        """
        read in the files to get appropriate information
        """
        # --> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                md_model = Model()
                md_model.read_model_file(self.model_fn)
                self.res_model = md_model.res_model
                self.grid_east = md_model.grid_east / self.dscale
                self.grid_north = md_model.grid_north / self.dscale
                self.grid_z = md_model.grid_z / self.dscale
                self.nodes_east = md_model.nodes_east / self.dscale
                self.nodes_north = md_model.nodes_north / self.dscale
                self.nodes_z = md_model.nodes_z / self.dscale
            else:
                raise mtex.MTpyError_file_handling(
                    '{0} does not exist, check path'.format(self.model_fn))

        # --> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                md_data = Data()
                md_data.read_data_file(self.data_fn)
                self.station_east = md_data.station_locations.rel_east / self.dscale
                self.station_north = md_data.station_locations.rel_north / self.dscale
                self.station_elev = md_data.station_locations.elev / self.dscale
                self.station_names = md_data.station_locations.station
            else:
                print 'Could not find data file {0}'.format(self.data_fn)

    def plot(self):
        """
        plot depth slices
        """
        # --> get information from files
        self.read_files()

        fdict = {'size': self.font_size + 2, 'weight': 'bold'}

        cblabeldict = {-2: '$10^{-3}$', -1: '$10^{-1}$', 0: '$10^{0}$', 1: '$10^{1}$',
                       2: '$10^{2}$', 3: '$10^{3}$', 4: '$10^{4}$', 5: '$10^{5}$',
                       6: '$10^{6}$', 7: '$10^{7}$', 8: '$10^{8}$'}

        # create an list of depth slices to plot
        if self.depth_index == None:
            zrange = range(self.grid_z.shape[0])
        elif type(self.depth_index) is int:
            zrange = [self.depth_index]
        elif type(self.depth_index) is list or \
                        type(self.depth_index) is np.ndarray:
            zrange = self.depth_index

        # set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                xlimits = (np.floor(self.station_east.min()),
                           np.ceil(self.station_east.max()))
            else:
                xlimits = (self.grid_east[5], self.grid_east[-5])
        else:
            xlimits = self.ew_limits

        if self.ns_limits == None:
            if self.station_north is not None:
                ylimits = (np.floor(self.station_north.min()),
                           np.ceil(self.station_north.max()))
            else:
                ylimits = (self.grid_north[5], self.grid_north[-5])
        else:
            ylimits = self.ns_limits

        # make a mesh grid of north and east
        self.mesh_east, self.mesh_north = np.meshgrid(self.grid_east,
                                                      self.grid_north,
                                                      indexing='ij')

        plt.rcParams['font.size'] = self.font_size

        # --> plot depths into individual figures
        for ii in zrange:
            depth = '{0:.3f} ({1})'.format(self.grid_z[ii],
                                           self.map_scale)
            fig = plt.figure(depth, figsize=self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            ax1 = fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)
            plot_res = np.log10(self.res_model[:, :, ii].T)
            mesh_plot = ax1.pcolormesh(self.mesh_east,
                                       self.mesh_north,
                                       plot_res,
                                       cmap=self.cmap,
                                       vmin=self.climits[0],
                                       vmax=self.climits[1])

            # plot the stations
            if self.station_east is not None:
                for ee, nn in zip(self.station_east, self.station_north):
                    ax1.text(ee, nn, '*',
                             verticalalignment='center',
                             horizontalalignment='center',
                             fontdict={'size': 5, 'weight': 'bold'})

            # set axis properties
            ax1.set_xlim(xlimits)
            ax1.set_ylim(ylimits)
            ax1.xaxis.set_minor_locator(MultipleLocator(self.xminorticks / self.dscale))
            ax1.yaxis.set_minor_locator(MultipleLocator(self.yminorticks / self.dscale))
            ax1.set_ylabel('Northing (' + self.map_scale + ')', fontdict=fdict)
            ax1.set_xlabel('Easting (' + self.map_scale + ')', fontdict=fdict)
            ax1.set_title('Depth = {0}'.format(depth), fontdict=fdict)

            # plot the grid if desired
            if self.plot_grid == 'y':
                east_line_xlist = []
                east_line_ylist = []
                for xx in self.grid_east:
                    east_line_xlist.extend([xx, xx])
                    east_line_xlist.append(None)
                    east_line_ylist.extend([self.grid_north.min(),
                                            self.grid_north.max()])
                    east_line_ylist.append(None)
                ax1.plot(east_line_xlist,
                         east_line_ylist,
                         lw=.25,
                         color='k')

                north_line_xlist = []
                north_line_ylist = []
                for yy in self.grid_north:
                    north_line_xlist.extend([self.grid_east.min(),
                                             self.grid_east.max()])
                    north_line_xlist.append(None)
                    north_line_ylist.extend([yy, yy])
                    north_line_ylist.append(None)
                ax1.plot(north_line_xlist,
                         north_line_ylist,
                         lw=.25,
                         color='k')

            # plot the colorbar
            if self.cb_location is None:
                if self.cb_orientation == 'horizontal':
                    self.cb_location = (ax1.axes.figbox.bounds[3] - .225,
                                        ax1.axes.figbox.bounds[1] + .05, .3, .025)

                elif self.cb_orientation == 'vertical':
                    self.cb_location = ((ax1.axes.figbox.bounds[2] - .15,
                                         ax1.axes.figbox.bounds[3] - .21, .025, .3))

            ax2 = fig.add_axes(self.cb_location)

            cb = mcb.ColorbarBase(ax2,
                                  cmap=self.cmap,
                                  norm=Normalize(vmin=self.climits[0],
                                                 vmax=self.climits[1]),
                                  orientation=self.cb_orientation)

            if self.cb_orientation == 'horizontal':
                cb.ax.xaxis.set_label_position('top')
                cb.ax.xaxis.set_label_coords(.5, 1.3)


            elif self.cb_orientation == 'vertical':
                cb.ax.yaxis.set_label_position('right')
                cb.ax.yaxis.set_label_coords(1.25, .5)
                cb.ax.yaxis.tick_left()
                cb.ax.tick_params(axis='y', direction='in')

            cb.set_label('Resistivity ($\Omega \cdot$m)',
                         fontdict={'size': self.font_size + 1})
            cb.set_ticks(np.arange(self.climits[0], self.climits[1] + 1))
            cb.set_ticklabels([cblabeldict[cc]
                               for cc in np.arange(self.climits[0],
                                                   self.climits[1] + 1)])

            self.fig_list.append(fig)

            # --> save plots to a common folder
            if self.save_plots == 'y':

                fig.savefig(os.path.join(self.save_path,
                                         "Depth_{}_{:.4f}.png".format(ii, self.grid_z[ii])),
                            dpi=self.fig_dpi, bbox_inches='tight')
                fig.clear()
                plt.close()

            else:
                pass

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()

    def update_plot(self, fig):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots depth slices of model from WS3DINV")


# ==============================================================================
# plot slices 
# ==============================================================================
class PlotSlices(object):
    #    """
    #    plot all slices and be able to scroll through the model
    #
    #    :Example: ::
    #
    #        >>> import mtpy.modeling.modem as modem
    #        >>> mfn = r"/home/modem/Inv1/Modular_NLCG_100.rho"
    #        >>> dfn = r"/home/modem/Inv1/ModEM_data.dat"
    #        >>> pds = ws.PlotSlices(model_fn=mfn, data_fn=dfn)
    #
    #    ======================= ===================================================
    #    Buttons                  Description
    #    ======================= ===================================================
    #    'e'                     moves n-s slice east by one model block
    #    'w'                     moves n-s slice west by one model block
    #    'n'                     moves e-w slice north by one model block
    #    'm'                     moves e-w slice south by one model block
    #    'd'                     moves depth slice down by one model block
    #    'u'                     moves depth slice up by one model block
    #    ======================= ===================================================
    #
    #
    #    ======================= ===================================================
    #    Attributes              Description
    #    ======================= ===================================================
    #    ax_en                   matplotlib.axes instance for depth slice  map view
    #    ax_ez                   matplotlib.axes instance for e-w slice
    #    ax_map                  matplotlib.axes instance for location map
    #    ax_nz                   matplotlib.axes instance for n-s slice
    #    climits                 (min , max) color limits on resistivity in log
    #                            scale. *default* is (0, 4)
    #    cmap                    name of color map for resisitiviy.
    #                            *default* is 'jet_r'
    #    data_fn                 full path to data file name
    #    dscale                  scaling parameter depending on map_scale
    #    east_line_xlist         list of line nodes of east grid for faster plotting
    #    east_line_ylist         list of line nodes of east grid for faster plotting
    #    ew_limits               (min, max) limits of e-w in map_scale units
    #                            *default* is None and scales to station area
    #    fig                     matplotlib.figure instance for figure
    #    fig_aspect              aspect ratio of plots. *default* is 1
    #    fig_dpi                 resolution of figure in dots-per-inch
    #                            *default* is 300
    #    fig_num                 figure instance number
    #    fig_size                [width, height] of figure window.
    #                            *default* is [6,6]
    #    font_dict               dictionary of font keywords, internally created
    #    font_size               size of ticklables in points, axes labes are
    #                            font_size+2. *default* is 7
    #    grid_east               relative location of grid nodes in e-w direction
    #                            in map_scale units
    #    grid_north              relative location of grid nodes in n-s direction
    #                            in map_scale units
    #    grid_z                  relative location of grid nodes in z direction
    #                            in map_scale units
    #    index_east              index value of grid_east being plotted
    #    index_north             index value of grid_north being plotted
    #    index_vertical          index value of grid_z being plotted
    #    initial_fn              full path to initial file
    #    key_press               matplotlib.canvas.connect instance
    #    map_scale               [ 'm' | 'km' ] scale of map. *default* is km
    #    mesh_east               np.meshgrid(grid_east, grid_north)[0]
    #    mesh_en_east            np.meshgrid(grid_east, grid_north)[0]
    #    mesh_en_north           np.meshgrid(grid_east, grid_north)[1]
    #    mesh_ez_east            np.meshgrid(grid_east, grid_z)[0]
    #    mesh_ez_vertical        np.meshgrid(grid_east, grid_z)[1]
    #    mesh_north              np.meshgrid(grid_east, grid_north)[1]
    #    mesh_nz_north           np.meshgrid(grid_north, grid_z)[0]
    #    mesh_nz_vertical        np.meshgrid(grid_north, grid_z)[1]
    #    model_fn                full path to model file
    #    ms                      size of station markers in points. *default* is 2
    #    nodes_east              relative distance betwen nodes in e-w direction
    #                            in map_scale units
    #    nodes_north             relative distance betwen nodes in n-s direction
    #                            in map_scale units
    #    nodes_z                 relative distance betwen nodes in z direction
    #                            in map_scale units
    #    north_line_xlist        list of line nodes north grid for faster plotting
    #    north_line_ylist        list of line nodes north grid for faster plotting
    #    ns_limits               (min, max) limits of plots in n-s direction
    #                            *default* is None, set veiwing area to station area
    #    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
    #                            *default* is 'y'
    #    res_model               np.ndarray(n_north, n_east, n_vertical) of
    #                            model resistivity values in linear scale
    #    station_color           color of station marker. *default* is black
    #    station_dict_east       location of stations for each east grid row
    #    station_dict_north      location of stations for each north grid row
    #    station_east            location of stations in east direction
    #    station_fn              full path to station file
    #    station_font_color      color of station label
    #    station_font_pad        padding between station marker and label
    #    station_font_rotation   angle of station label
    #    station_font_size       font size of station label
    #    station_font_weight     weight of font for station label
    #    station_id              [min, max] index values for station labels
    #    station_marker          station marker
    #    station_names           name of stations
    #    station_north           location of stations in north direction
    #    subplot_bottom          distance between axes and bottom of figure window
    #    subplot_hspace          distance between subplots in vertical direction
    #    subplot_left            distance between axes and left of figure window
    #    subplot_right           distance between axes and right of figure window
    #    subplot_top             distance between axes and top of figure window
    #    subplot_wspace          distance between subplots in horizontal direction
    #    title                   title of plot
    #    z_limits                (min, max) limits in vertical direction,
    #    ======================= ===================================================
    #
    #    """

    def __init__(self, model_fn, data_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.font_size = kwargs.pop('font_size', 7)

        self.subplot_wspace = .20
        self.subplot_hspace = .30
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.index_vertical = kwargs.pop('index_vertical', 0)
        self.index_east = kwargs.pop('index_east', 0)
        self.index_north = kwargs.pop('index_north', 0)

        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.climits = kwargs.pop('climits', (0, 4))

        self.map_scale = kwargs.pop('map_scale', 'km')
        # make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)
        self.z_limits = kwargs.pop('z_limits', None)

        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None

        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        self.mesh_east = None
        self.mesh_north = None

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.station_id = kwargs.pop('station_id', None)
        self.station_font_size = kwargs.pop('station_font_size', 8)
        self.station_font_pad = kwargs.pop('station_font_pad', 1.0)
        self.station_font_weight = kwargs.pop('station_font_weight', 'bold')
        self.station_font_rotation = kwargs.pop('station_font_rotation', 60)
        self.station_font_color = kwargs.pop('station_font_color', 'k')
        self.station_marker = kwargs.pop('station_marker',
                                         r"$\blacktriangledown$")
        self.station_color = kwargs.pop('station_color', 'k')
        self.ms = kwargs.pop('ms', 10)

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def read_files(self):
        """
        read in the files to get appropriate information
        """
        # --> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                md_model = Model()
                md_model.read_model_file(self.model_fn)
                self.res_model = md_model.res_model
                self.grid_east = md_model.grid_east / self.dscale
                self.grid_north = md_model.grid_north / self.dscale
                self.grid_z = md_model.grid_z / self.dscale
                self.nodes_east = md_model.nodes_east / self.dscale
                self.nodes_north = md_model.nodes_north / self.dscale
                self.nodes_z = md_model.nodes_z / self.dscale
            else:
                raise mtex.MTpyError_file_handling(
                    '{0} does not exist, check path'.format(self.model_fn))

        # --> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                md_data = Data()
                md_data.read_data_file(self.data_fn)
                self.station_east = md_data.station_locations.rel_east / self.dscale
                self.station_north = md_data.station_locations.rel_north / self.dscale
                self.station_names = md_data.station_locations.station
                self.station_elev = md_data.station_locations.elev / self.dscale
            else:
                print 'Could not find data file {0}'.format(self.data_fn)

    def plot(self):
        """
        plot:
            east vs. vertical,
            north vs. vertical,
            east vs. north
            
        
        """

        self.read_files()

        self.get_station_grid_locations()

        print "=============== ==============================================="
        print "    Buttons                  Description                       "
        print "=============== ==============================================="
        print "     'e'          moves n-s slice east by one model block"
        print "     'w'          moves n-s slice west by one model block"
        print "     'n'          moves e-w slice north by one model block"
        print "     'm'          moves e-w slice south by one model block"
        print "     'd'          moves depth slice down by one model block"
        print "     'u'          moves depth slice up by one model block"
        print "=============== ==============================================="

        self.font_dict = {'size': self.font_size + 2, 'weight': 'bold'}

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size

        # set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                self.ew_limits = (np.floor(self.station_east.min()),
                                  np.ceil(self.station_east.max()))
            else:
                self.ew_limits = (self.grid_east[5], self.grid_east[-5])

        if self.ns_limits == None:
            if self.station_north is not None:
                self.ns_limits = (np.floor(self.station_north.min()),
                                  np.ceil(self.station_north.max()))
            else:
                self.ns_limits = (self.grid_north[5], self.grid_north[-5])

        if self.z_limits == None:
            depth_limit = max([(abs(self.ew_limits[0]) + abs(self.ew_limits[1])),
                               (abs(self.ns_limits[0]) + abs(self.ns_limits[1]))])
            self.z_limits = (-5000 / self.dscale, depth_limit)

        self.fig = plt.figure(self.fig_num, figsize=self.fig_size,
                              dpi=self.fig_dpi)
        plt.clf()
        gs = gridspec.GridSpec(2, 2,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom,
                               right=self.subplot_right,
                               hspace=self.subplot_hspace)

        # make subplots
        self.ax_ez = self.fig.add_subplot(gs[0, 0], aspect=self.fig_aspect)
        self.ax_nz = self.fig.add_subplot(gs[1, 1], aspect=self.fig_aspect)
        self.ax_en = self.fig.add_subplot(gs[1, 0], aspect=self.fig_aspect)
        self.ax_map = self.fig.add_subplot(gs[0, 1])

        # make grid meshes being sure the indexing is correct
        self.mesh_ez_east, self.mesh_ez_vertical = np.meshgrid(self.grid_east,
                                                               self.grid_z,
                                                               indexing='ij')
        self.mesh_nz_north, self.mesh_nz_vertical = np.meshgrid(self.grid_north,
                                                                self.grid_z,
                                                                indexing='ij')
        self.mesh_en_east, self.mesh_en_north = np.meshgrid(self.grid_east,
                                                            self.grid_north,
                                                            indexing='ij')

        # --> plot east vs vertical
        self._update_ax_ez()

        # --> plot north vs vertical
        self._update_ax_nz()

        # --> plot east vs north
        self._update_ax_en()

        # --> plot the grid as a map view
        self._update_map()

        # plot color bar
        cbx = mcb.make_axes(self.ax_map, fraction=.15, shrink=.75, pad=.15)
        cb = mcb.ColorbarBase(cbx[0],
                              cmap=self.cmap,
                              norm=Normalize(vmin=self.climits[0],
                                             vmax=self.climits[1]))

        cb.ax.yaxis.set_label_position('right')
        cb.ax.yaxis.set_label_coords(1.25, .5)
        cb.ax.yaxis.tick_left()
        cb.ax.tick_params(axis='y', direction='in')

        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size': self.font_size + 1})

        cb.set_ticks(np.arange(np.ceil(self.climits[0]),
                               np.floor(self.climits[1] + 1)))
        cblabeldict = {-2: '$10^{-3}$', -1: '$10^{-1}$', 0: '$10^{0}$', 1: '$10^{1}$',
                       2: '$10^{2}$', 3: '$10^{3}$', 4: '$10^{4}$', 5: '$10^{5}$',
                       6: '$10^{6}$', 7: '$10^{7}$', 8: '$10^{8}$'}
        cb.set_ticklabels([cblabeldict[cc]
                           for cc in np.arange(np.ceil(self.climits[0]),
                                               np.floor(self.climits[1] + 1))])

        plt.show()

        self.key_press = self.fig.canvas.mpl_connect('key_press_event',
                                                     self.on_key_press)

    def on_key_press(self, event):
        """
        on a key press change the slices
        
        """

        key_press = event.key

        if key_press == 'n':
            if self.index_north == self.grid_north.size:
                print 'Already at northern most grid cell'
            else:
                self.index_north += 1
                if self.index_north > self.grid_north.size:
                    self.index_north = self.grid_north.size
            self._update_ax_ez()
            self._update_map()

        if key_press == 'm':
            if self.index_north == 0:
                print 'Already at southern most grid cell'
            else:
                self.index_north -= 1
                if self.index_north < 0:
                    self.index_north = 0
            self._update_ax_ez()
            self._update_map()

        if key_press == 'e':
            if self.index_east == self.grid_east.size:
                print 'Already at eastern most grid cell'
            else:
                self.index_east += 1
                if self.index_east > self.grid_east.size:
                    self.index_east = self.grid_east.size
            self._update_ax_nz()
            self._update_map()

        if key_press == 'w':
            if self.index_east == 0:
                print 'Already at western most grid cell'
            else:
                self.index_east -= 1
                if self.index_east < 0:
                    self.index_east = 0
            self._update_ax_nz()
            self._update_map()

        if key_press == 'd':
            if self.index_vertical == self.grid_z.size:
                print 'Already at deepest grid cell'
            else:
                self.index_vertical += 1
                if self.index_vertical > self.grid_z.size:
                    self.index_vertical = self.grid_z.size
            self._update_ax_en()
            print 'Depth = {0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                                 self.map_scale)

        if key_press == 'u':
            if self.index_vertical == 0:
                print 'Already at surface grid cell'
            else:
                self.index_vertical -= 1
                if self.index_vertical < 0:
                    self.index_vertical = 0
            self._update_ax_en()
            print 'Depth = {0:.5gf} ({1})'.format(self.grid_z[self.index_vertical],
                                                  self.map_scale)

    def _update_ax_ez(self):
        """
        update east vs vertical plot
        """
        self.ax_ez.cla()
        plot_ez = np.log10(self.res_model[self.index_north, :, :])
        self.ax_ez.pcolormesh(self.mesh_ez_east,
                              self.mesh_ez_vertical,
                              plot_ez,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        # plot stations
        for sx in self.station_dict_north[self.grid_north[self.index_north]]:
            self.ax_ez.text(sx,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms,
                                      'color': self.station_color})

        self.ax_ez.set_xlim(self.ew_limits)
        self.ax_ez.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_ez.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_ez.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()

    def _update_ax_nz(self):
        """
        update east vs vertical plot
        """
        self.ax_nz.cla()
        plot_nz = np.log10(self.res_model[:, self.index_east, :])
        self.ax_nz.pcolormesh(self.mesh_nz_north,
                              self.mesh_nz_vertical,
                              plot_nz,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        # plot stations
        for sy in self.station_dict_east[self.grid_east[self.index_east]]:
            self.ax_nz.text(sy,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms,
                                      'color': self.station_color})
        self.ax_nz.set_xlim(self.ns_limits)
        self.ax_nz.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_nz.set_xlabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_nz.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()

    def _update_ax_en(self):
        """
        update east vs vertical plot
        """

        self.ax_en.cla()
        plot_en = np.log10(self.res_model[:, :, self.index_vertical].T)
        self.ax_en.pcolormesh(self.mesh_en_east,
                              self.mesh_en_north,
                              plot_en,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        self.ax_en.set_xlim(self.ew_limits)
        self.ax_en.set_ylim(self.ns_limits)
        self.ax_en.set_ylabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_en.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        # --> plot the stations
        if self.station_east is not None:
            for ee, nn, elev, name in zip(self.station_east,
                                          self.station_north,
                                          self.station_elev,
                                          self.station_names):
                if elev <= self.grid_z[self.index_vertical]:
                    self.ax_en.text(ee, nn, '+',
                                    verticalalignment='center',
                                    horizontalalignment='center',
                                    fontdict={'size': 7, 'weight': 'bold',
                                              'color': (.75, 0, 0)})
                    self.ax_en.text(ee, nn, name[2:],
                                    verticalalignment='center',
                                    horizontalalignment='center',
                                    fontdict={'size': 7, 'weight': 'bold',
                                              'color': (.75, 0, 0)})

        self.fig.canvas.draw()
        self._update_map()

    def _update_map(self):
        self.ax_map.cla()
        self.east_line_xlist = []
        self.east_line_ylist = []
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx, xx])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min(),
                                         self.grid_north.max()])
            self.east_line_ylist.append(None)
        self.ax_map.plot(self.east_line_xlist,
                         self.east_line_ylist,
                         lw=.25,
                         color='k')

        self.north_line_xlist = []
        self.north_line_ylist = []
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min(),
                                          self.grid_east.max()])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy, yy])
            self.north_line_ylist.append(None)
        self.ax_map.plot(self.north_line_xlist,
                         self.north_line_ylist,
                         lw=.25,
                         color='k')
        # --> e-w indication line
        self.ax_map.plot([self.grid_east.min(),
                          self.grid_east.max()],
                         [self.grid_north[self.index_north + 1],
                          self.grid_north[self.index_north + 1]],
                         lw=1,
                         color='g')

        # --> e-w indication line
        self.ax_map.plot([self.grid_east[self.index_east + 1],
                          self.grid_east[self.index_east + 1]],
                         [self.grid_north.min(),
                          self.grid_north.max()],
                         lw=1,
                         color='b')
        # --> plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax_map.text(ee, nn, '*',
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size': 5, 'weight': 'bold'})

        self.ax_map.set_xlim(self.ew_limits)
        self.ax_map.set_ylim(self.ns_limits)
        self.ax_map.set_ylabel('Northing ({0})'.format(self.map_scale),
                               fontdict=self.font_dict)
        self.ax_map.set_xlabel('Easting ({0})'.format(self.map_scale),
                               fontdict=self.font_dict)

        # plot stations
        self.ax_map.text(self.ew_limits[0] * .95, self.ns_limits[1] * .95,
                         '{0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                                self.map_scale),
                         horizontalalignment='left',
                         verticalalignment='top',
                         bbox={'facecolor': 'white'},
                         fontdict=self.font_dict)

        self.fig.canvas.draw()

    def get_station_grid_locations(self):
        """
        get the grid line on which a station resides for plotting
        
        """
        self.station_dict_east = dict([(gx, []) for gx in self.grid_east])
        self.station_dict_north = dict([(gy, []) for gy in self.grid_north])
        if self.station_east is not None:
            for ss, sx in enumerate(self.station_east):
                gx = np.where(self.grid_east <= sx)[0][-1]
                self.station_dict_east[self.grid_east[gx]].append(self.station_north[ss])

            for ss, sy in enumerate(self.station_north):
                gy = np.where(self.grid_north <= sy)[0][-1]
                self.station_dict_north[self.grid_north[gy]].append(self.station_east[ss])
        else:
            return

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn=None, fig_dpi=None, file_format='pdf',
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
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
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_E{0}_N{1}_Z{2}.{3}'.format(
                self.index_east, self.index_north,
                self.index_vertical, file_format))
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print 'Saved figure to: ' + self.fig_fn


# ==============================================================================
# plot rms maps
# ==============================================================================
class Plot_RMS_Maps(object):
    """
    plots the RMS as (data-model)/(error) in map view for all components
    of the data file.  Gets this infomration from the .res file output
    by ModEM.
    
    Arguments:
    ------------------
    
        **residual_fn** : string
                          full path to .res file
                          
    =================== =======================================================
    Attributes                   Description    
    =================== =======================================================
    fig                 matplotlib.figure instance for a single plot                       
    fig_dpi             dots-per-inch resolution of figure *default* is 200
    fig_num             number of fig instance *default* is 1
    fig_size            size of figure in inches [width, height] 
                        *default* is [7,6]
    font_size           font size of tick labels, axis labels are +2
                        *default* is 8 
    marker              marker style for station rms, 
                        see matplotlib.line for options,
                        *default* is 's' --> square
    marker_size         size of marker in points. *default* is 10
    pad_x               padding in map units from edge of the axis to stations
                        at the extremeties in longitude. 
                        *default* is 1/2 tick_locator
    pad_y               padding in map units from edge of the axis to stations
                        at the extremeties in latitude. 
                        *default* is 1/2 tick_locator 
    period_index        index of the period you want to plot according to 
                        self.residual.period_list. *default* is 1
    plot_yn             [ 'y' | 'n' ] default is 'y' to plot on instantiation
    plot_z_list         internal variable for plotting
    residual            modem.Data instance that holds all the information 
                        from the residual_fn given
    residual_fn         full path to .res file
    rms_cmap            matplotlib.cm object for coloring the markers
    rms_cmap_dict       dictionary of color values for rms_cmap 
    rms_max             maximum rms to plot. *default* is 5.0
    rms_min             minimum rms to plot. *default* is 1.0
    save_path           path to save figures to. *default* is directory of 
                        residual_fn
    subplot_bottom      spacing from axis to bottom of figure canvas.
                        *default* is .1
    subplot_hspace      horizontal spacing between subplots.
                        *default* is .1
    subplot_left        spacing from axis to left of figure canvas.
                        *default* is .1
    subplot_right       spacing from axis to right of figure canvas.
                        *default* is .9
    subplot_top         spacing from axis to top of figure canvas.
                        *default* is .95
    subplot_vspace      vertical spacing between subplots.
                        *default* is .01
    tick_locator        increment for x and y major ticks. *default* is 
                        limits/5  
    =================== =======================================================
    
    =================== =======================================================
    Methods             Description    
    =================== =======================================================
    plot                plot rms maps for a single period 
    plot_loop           loop over all frequencies and save figures to save_path
    read_residual_fn    read in residual_fn
    redraw_plot         after updating attributes call redraw_plot to 
                        well redraw the plot
    save_figure         save the figure to a file
    =================== =======================================================
    
    
    :Example: ::
    
        >>> import mtpy.modeling.modem as modem
        >>> rms_plot = Plot_RMS_Maps(r"/home/ModEM/Inv1/mb_NLCG_030.res")
        >>> # change some attributes
        >>> rms_plot.fig_size = [6, 4]
        >>> rms_plot.rms_max = 3
        >>> rms_plot.redraw_plot()
        >>> # happy with the look now loop over all periods
        >>> rms_plot.plot_loop()
    """

    def __init__(self, residual_fn, **kwargs):
        self.residual_fn = residual_fn
        self.residual = None
        self.save_path = kwargs.pop('save_path', os.path.dirname(self.residual_fn))

        self.period_index = kwargs.pop('period_index', 0)

        self.subplot_left = kwargs.pop('subplot_left', .1)
        self.subplot_right = kwargs.pop('subplot_right', .9)
        self.subplot_top = kwargs.pop('subplot_top', .95)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .1)
        self.subplot_vspace = kwargs.pop('subplot_vspace', .01)

        self.font_size = kwargs.pop('font_size', 8)

        self.fig_size = kwargs.pop('fig_size', [7.75, 6.75])
        self.fig_dpi = kwargs.pop('fig_dpi', 200)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig = None

        self.marker = kwargs.pop('marker', 's')
        self.marker_size = kwargs.pop('marker_size', 10)

        self.rms_max = kwargs.pop('rms_max', 5)
        self.rms_min = kwargs.pop('rms_min', 0)

        self.tick_locator = kwargs.pop('tick_locator', None)
        self.pad_x = kwargs.pop('pad_x', None)
        self.pad_y = kwargs.pop('pad_y', None)

        self.plot_yn = kwargs.pop('plot_yn', 'y')

        # colormap for rms, goes white to black from 0 to rms max and 
        # red below 1 to show where the data is being over fit
        self.rms_cmap_dict = {'red': ((0.0, 1.0, 1.0),
                                      (0.2, 1.0, 1.0),
                                      (1.0, 0.0, 0.0)),
                              'green': ((0.0, 0.0, 0.0),
                                        (0.2, 1.0, 1.0),
                                        (1.0, 0.0, 0.0)),
                              'blue': ((0.0, 0.0, 0.0),
                                       (0.2, 1.0, 1.0),
                                       (1.0, 0.0, 0.0))}

        self.rms_cmap = colors.LinearSegmentedColormap('rms_cmap',
                                                       self.rms_cmap_dict,
                                                       256)

        self.plot_z_list = [{'label': r'$Z_{xx}$', 'index': (0, 0), 'plot_num': 1},
                            {'label': r'$Z_{xy}$', 'index': (0, 1), 'plot_num': 2},
                            {'label': r'$Z_{yx}$', 'index': (1, 0), 'plot_num': 3},
                            {'label': r'$Z_{yy}$', 'index': (1, 1), 'plot_num': 4},
                            {'label': r'$T_{x}$', 'index': (0, 0), 'plot_num': 5},
                            {'label': r'$T_{y}$', 'index': (0, 1), 'plot_num': 6}]

        if self.plot_yn == 'y':
            self.plot()

    def read_residual_fn(self):
        if self.residual is None:
            self.residual = Data()
            self.residual.read_data_file(self.residual_fn)
        else:
            pass

    def plot(self):
        """
        plot rms in map view
        """

        self.read_residual_fn()

        font_dict = {'size': self.font_size + 2, 'weight': 'bold'}
        rms_1 = 1. / self.rms_max

        if self.tick_locator is None:
            x_locator = np.round((self.residual.data_array['lon'].max() -
                                  self.residual.data_array['lon'].min()) / 5, 2)
            y_locator = np.round((self.residual.data_array['lat'].max() -
                                  self.residual.data_array['lat'].min()) / 5, 2)

            if x_locator > y_locator:
                self.tick_locator = x_locator

            elif x_locator < y_locator:
                self.tick_locator = y_locator

        if self.pad_x is None:
            self.pad_x = self.tick_locator / 2
        if self.pad_y is None:
            self.pad_y = self.tick_locator / 2

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_hspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_vspace
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        for p_dict in self.plot_z_list:
            ax = self.fig.add_subplot(3, 2, p_dict['plot_num'], aspect='equal')

            ii = p_dict['index'][0]
            jj = p_dict['index'][0]

            for r_arr in self.residual.data_array:
                # calulate the rms self.residual/error
                if p_dict['plot_num'] < 5:
                    rms = r_arr['z'][self.period_index, ii, jj].__abs__() / \
                          (r_arr['z_err'][self.period_index, ii, jj].real)

                else:
                    rms = r_arr['tip'][self.period_index, ii, jj].__abs__() / \
                          (r_arr['tip_err'][self.period_index, ii, jj].real)

                # color appropriately
                if np.nan_to_num(rms) == 0.0:
                    marker_color = (1, 1, 1)
                    marker = '.'
                    marker_size = .1
                    marker_edge_color = (1, 1, 1)
                if rms > self.rms_max:
                    marker_color = (0, 0, 0)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)

                elif rms >= 1 and rms <= self.rms_max:
                    r_color = 1 - rms / self.rms_max + rms_1
                    marker_color = (r_color, r_color, r_color)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)

                elif rms < 1:
                    r_color = 1 - rms / self.rms_max
                    marker_color = (1, r_color, r_color)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)

                ax.plot(r_arr['lon'], r_arr['lat'],
                        marker=marker,
                        ms=marker_size,
                        mec=marker_edge_color,
                        mfc=marker_color,
                        zorder=3)

            if p_dict['plot_num'] == 1 or p_dict['plot_num'] == 3:
                ax.set_ylabel('Latitude (deg)', fontdict=font_dict)
                plt.setp(ax.get_xticklabels(), visible=False)

            elif p_dict['plot_num'] == 2 or p_dict['plot_num'] == 4:
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)

            elif p_dict['plot_num'] == 6:
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.set_xlabel('Longitude (deg)', fontdict=font_dict)

            else:
                ax.set_xlabel('Longitude (deg)', fontdict=font_dict)
                ax.set_ylabel('Latitude (deg)', fontdict=font_dict)

            ax.text(self.residual.data_array['lon'].min() + .005 - self.pad_x,
                    self.residual.data_array['lat'].max() - .005 + self.pad_y,
                    p_dict['label'],
                    verticalalignment='top',
                    horizontalalignment='left',
                    bbox={'facecolor': 'white'},
                    zorder=3)

            ax.tick_params(direction='out')
            ax.grid(zorder=0, color=(.75, .75, .75))

            # [line.set_zorder(3) for line in ax.lines]

            ax.set_xlim(self.residual.data_array['lon'].min() - self.pad_x,
                        self.residual.data_array['lon'].max() + self.pad_x)

            ax.set_ylim(self.residual.data_array['lat'].min() - self.pad_y,
                        self.residual.data_array['lat'].max() + self.pad_y)

            ax.xaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.yaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%2.2f'))

        # cb_ax = mcb.make_axes(ax, orientation='vertical', fraction=.1)
        cb_ax = self.fig.add_axes([self.subplot_right + .02, .225, .02, .45])
        color_bar = mcb.ColorbarBase(cb_ax,
                                     cmap=self.rms_cmap,
                                     norm=colors.Normalize(vmin=self.rms_min,
                                                           vmax=self.rms_max),
                                     orientation='vertical')

        color_bar.set_label('RMS', fontdict=font_dict)

        self.fig.suptitle('period = {0:.5g} (s)'.format(self.residual.period_list[self.period_index]),
                          fontdict={'size': self.font_size + 3, 'weight': 'bold'})
        plt.show()

    def redraw_plot(self):
        plt.close('all')
        self.plot()

    def save_figure(self, save_path=None, save_fn_basename=None,
                    save_fig_dpi=None, fig_format='.png', fig_close=True):
        """
        save figure in the desired format
        """
        if save_path is not None:
            self.save_path = save_path

        if save_fn_basename is not None:
            pass
        else:
            save_fn_basename = '{0:02}_RMS_{1:.5g}_s.{2}'.format(self.period_index,
                                                                 self.residual.period_list[self.period_index],
                                                                 fig_format)
        save_fn = os.path.join(self.save_path, save_fn_basename)

        if save_fig_dpi is not None:
            self.fig_dpi = save_fig_dpi

        self.fig.savefig(save_fn, dpi=self.fig_dpi)
        print 'saved file to {0}'.format(save_fn)

        if fig_close == True:
            plt.close('all')

    def plot_loop(self, fig_format='png'):
        """
        loop over all periods and save figures accordingly
        """
        self.read_residual_fn()

        for f_index in range(self.residual.period_list.size):
            self.period_index = f_index
            self.plot()
            self.save_figure(fig_format=fig_format)


# ==============================================================================
# Exceptions
# ==============================================================================
