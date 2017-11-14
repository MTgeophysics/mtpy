"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
import os

import numpy as np
from matplotlib import pyplot as plt, gridspec as gridspec, colorbar as mcb
from matplotlib.colors import Normalize

from mtpy.utils import exceptions as mtex
from .data import Data
from .model import Model

__all__ = ['PlotSlices']


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


