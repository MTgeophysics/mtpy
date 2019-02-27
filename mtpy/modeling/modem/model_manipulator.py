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
from matplotlib import cm as cm, pyplot as plt, colorbar as mcb, colors as colors, widgets as widgets

from mtpy.imaging import mtplottools as mtplottools
from .data import Data
from .model import Model

__all__ = ['ModelManipulator']


class ModelManipulator(Model):
    """
    will plot a model from wsinv3d or init file so the user can manipulate the
    resistivity values relatively easily.  At the moment only plotted
    in map view.


    :Example: ::
        >>> import mtpy.modeling.ws3dinv as ws
        >>> initial_fn = r"/home/MT/ws3dinv/Inv1/WSInitialFile"
        >>> mm = ws.WSModelManipulator(initial_fn=initial_fn)

    =================== =======================================================
    Buttons              Description
    =================== =======================================================
    '='                 increase depth to next vertical node (deeper)
    '-'                 decrease depth to next vertical node (shallower)
    'q'                 quit the plot, rewrites initial file when pressed
    'a'                 copies the above horizontal layer to the present layer
    'b'                 copies the below horizonal layer to present layer
    'u'                 undo previous change
    =================== =======================================================


    =================== =======================================================
    Attributes          Description
    =================== =======================================================
    ax1                 matplotlib.axes instance for mesh plot of the model
    ax2                 matplotlib.axes instance of colorbar
    cb                  matplotlib.colorbar instance for colorbar
    cid_depth           matplotlib.canvas.connect for depth
    cmap                matplotlib.colormap instance
    cmax                maximum value of resistivity for colorbar. (linear)
    cmin                minimum value of resistivity for colorbar (linear)
    data_fn             full path fo data file
    depth_index         integer value of depth slice for plotting
    dpi                 resolution of figure in dots-per-inch
    dscale              depth scaling, computed internally
    east_line_xlist     list of east mesh lines for faster plotting
    east_line_ylist     list of east mesh lines for faster plotting
    fdict               dictionary of font properties
    fig                 matplotlib.figure instance
    fig_num              number of figure instance
    fig_size             size of figure in inches
    font_size           size of font in points
    grid_east           location of east nodes in relative coordinates
    grid_north          location of north nodes in relative coordinates
    grid_z              location of vertical nodes in relative coordinates
    initial_fn          full path to initial file
    m_height            mean height of horizontal cells
    m_width             mean width of horizontal cells
    map_scale            [ 'm' | 'km' ] scale of map
    mesh_east           np.meshgrid of east, north
    mesh_north          np.meshgrid of east, north
    mesh_plot           matplotlib.axes.pcolormesh instance
    model_fn            full path to model file
    new_initial_fn      full path to new initial file
    nodes_east          spacing between east nodes
    nodes_north         spacing between north nodes
    nodes_z             spacing between vertical nodes
    north_line_xlist    list of coordinates of north nodes for faster plotting
    north_line_ylist    list of coordinates of north nodes for faster plotting
    plot_yn             [ 'y' | 'n' ] plot on instantiation
    radio_res           matplotlib.widget.radio instance for change resistivity
    rect_selector       matplotlib.widget.rect_selector
    res                 np.ndarray(nx, ny, nz) for model in linear resistivity
    res_copy            copy of res for undo
    res_dict            dictionary of segmented resistivity values
    res_list            list of resistivity values for model linear scale
    res_model           np.ndarray(nx, ny, nz) of resistivity values from
                        res_list (linear scale)
    res_model_int       np.ndarray(nx, ny, nz) of integer values corresponding
                        to res_list for initial model
    res_value           current resistivty value of radio_res
    save_path           path to save initial file to
    station_east        station locations in east direction
    station_north       station locations in north direction
    xlimits             limits of plot in e-w direction
    ylimits             limits of plot in n-s direction
    =================== =======================================================

    """

    def __init__(self, model_fn=None, data_fn=None, **kwargs):

        # be sure to initialize Model
        Model.__init__(self, model_fn=model_fn, **kwargs)

        self.data_fn = data_fn
        self.model_fn_basename = kwargs.pop('model_fn_basename',
                                            'ModEM_Model_rw.ws')

        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.data_fn is not None:
            self.save_path = os.path.dirname(self.data_fn)
        else:
            self.save_path = os.getcwd()

        # station locations in relative coordinates read from data file
        self.station_east = None
        self.station_north = None

        # --> set map scale
        self.map_scale = kwargs.pop('map_scale', 'km')

        self.m_width = 100
        self.m_height = 100

        # --> scale the map coordinates
        if self.map_scale == 'km':
            self.dscale = 1000.
        if self.map_scale == 'm':
            self.dscale = 1.

        # figure attributes
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.cb = None
        self.east_line_xlist = None
        self.east_line_ylist = None
        self.north_line_xlist = None
        self.north_line_ylist = None

        # make a default resistivity list to change values
        self._res_sea = 0.3
        self._res_air = 1E12
        self.res_dict = None
        self.res_list = kwargs.pop('res_list', None)
        if self.res_list is None:
            self.set_res_list(np.array([self._res_sea, 1, 10, 50, 100, 500,
                                        1000, 5000],
                                       dtype=np.float))

        # set initial resistivity value
        self.res_value = self.res_list[0]
        self.cov_arr = None

        # --> set map limits
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)

        self.font_size = kwargs.pop('font_size', 7)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.cmap = kwargs.pop('cmap', cm.jet_r)
        self.depth_index = kwargs.pop('depth_index', 0)

        self.fdict = {'size': self.font_size + 2, 'weight': 'bold'}

        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .8)
        self.subplot_left = kwargs.pop('subplot_left', .01)
        self.subplot_top = kwargs.pop('subplot_top', .93)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)

        # plot on initialization
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.get_model()
            self.plot()

    def set_res_list(self, res_list):
        """
        on setting res_list also set the res_dict to correspond
        """
        self.res_list = res_list
        # make a dictionary of values to write to file.
        self.res_dict = dict([(res, ii)
                              for ii, res in enumerate(self.res_list, 1)])
        if self.fig is not None:
            plt.close()
            self.plot()

    # ---read files-------------------------------------------------------------
    def get_model(self):
        """
        reads in initial file or model file and set attributes:
            -resmodel
            -northrid
            -eastrid
            -zgrid
            -res_list if initial file

        """
        # --> read in model file
        self.read_model_file()

        self.cov_arr = np.ones_like(self.res_model)

        # --> read in data file if given
        if self.data_fn is not None:
            md_data = Data()
            md_data.read_data_file(self.data_fn)

            # get station locations
            self.station_east = md_data.station_locations.rel_east
            self.station_north = md_data.station_locations.rel_north

        # get cell block sizes
        self.m_height = np.median(self.nodes_north[5:-5]) / self.dscale
        self.m_width = np.median(self.nodes_east[5:-5]) / self.dscale

        # make a copy of original in case there are unwanted changes
        self.res_copy = self.res_model.copy()

    # ---plot model-------------------------------------------------------------
    def plot(self):
        """
        plots the model with:
            -a radio dial for depth slice
            -radio dial for resistivity value

        """
        # set plot properties
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        font_dict = {'size': self.font_size + 2, 'weight': 'bold'}

        # make sure there is a model to plot
        if self.res_model is None:
            self.get_model()

        self.cmin = np.floor(np.log10(min(self.res_list)))
        self.cmax = np.ceil(np.log10(max(self.res_list)))

        # -->Plot properties
        plt.rcParams['font.size'] = self.font_size

        # need to add an extra row and column to east and north to make sure
        # all is plotted see pcolor for details.
        plot_east = self.grid_east / self.dscale
        plot_north = self.grid_north / self.dscale

        # make a mesh grid for plotting
        # the 'ij' makes sure the resulting grid is in east, north
        self.mesh_east, self.mesh_north = np.meshgrid(plot_east,
                                                      plot_north,
                                                      indexing='ij')

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()
        self.ax1 = self.fig.add_subplot(1, 1, 1, aspect='equal')

        # transpose to make x--east and y--north
        plot_res = np.log10(self.res_model[:, :, self.depth_index].T)

        self.mesh_plot = self.ax1.pcolormesh(self.mesh_east,
                                             self.mesh_north,
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)

        # on plus or minus change depth slice
        self.cid_depth = \
            self.mesh_plot.figure.canvas.mpl_connect('key_press_event',
                                                     self._on_key_callback)

        # plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax1.text(ee / self.dscale, nn / self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size': self.font_size - 2,
                                        'weight': 'bold'})

        # set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(xmin=self.grid_east.min() / self.dscale,
                              xmax=self.grid_east.max() / self.dscale)

        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(ymin=self.grid_north.min() / self.dscale,
                              ymax=self.grid_north.max() / self.dscale)

        # self.ax1.xaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        # self.ax1.yaxis.set_minor_locator(MultipleLocator(100*1./dscale))

        self.ax1.set_ylabel('Northing (' + self.map_scale + ')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting (' + self.map_scale + ')',
                            fontdict=self.fdict)

        depth_title = self.grid_z[self.depth_index] / self.dscale

        self.ax1.set_title('Depth = {:.3f} '.format(depth_title) + \
                           '(' + self.map_scale + ')',
                           fontdict=self.fdict)

        # plot the grid if desired
        self.east_line_xlist = []
        self.east_line_ylist = []
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx / self.dscale, xx / self.dscale])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min() / self.dscale,
                                         self.grid_north.max() / self.dscale])
            self.east_line_ylist.append(None)
        self.ax1.plot(self.east_line_xlist,
                      self.east_line_ylist,
                      lw=.25,
                      color='k')

        self.north_line_xlist = []
        self.north_line_ylist = []
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min() / self.dscale,
                                          self.grid_east.max() / self.dscale])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy / self.dscale, yy / self.dscale])
            self.north_line_ylist.append(None)
        self.ax1.plot(self.north_line_xlist,
                      self.north_line_ylist,
                      lw=.25,
                      color='k')

        # plot the colorbar
        #        self.ax2 = mcb.make_axes(self.ax1, orientation='vertical', shrink=.35)
        self.ax2 = self.fig.add_axes([.81, .45, .16, .03])
        self.ax2.xaxis.set_ticks_position('top')
        # seg_cmap = ws.cmap_discretize(self.cmap, len(self.res_list))
        self.cb = mcb.ColorbarBase(self.ax2, cmap=self.cmap,
                                   norm=colors.Normalize(vmin=self.cmin,
                                                         vmax=self.cmax),
                                   orientation='horizontal')

        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
                          fontdict={'size': self.font_size})
        self.cb.set_ticks(np.arange(self.cmin, self.cmax + 1))
        self.cb.set_ticklabels([mtplottools.labeldict[cc]
                                for cc in np.arange(self.cmin, self.cmax + 1)])

        # make a resistivity radio button
        # resrb = self.fig.add_axes([.85,.1,.1,.2])
        # reslabels = ['{0:.4g}'.format(res) for res in self.res_list]
        # self.radio_res = widgets.RadioButtons(resrb, reslabels,
        #                                active=self.res_dict[self.res_value])

        #        slider_ax_bounds = list(self.cb.ax.get_position().bounds)
        #        slider_ax_bounds[0] += .1
        slider_ax = self.fig.add_axes([.81, .5, .16, .03])
        self.slider_res = widgets.Slider(slider_ax, 'Resistivity',
                                         self.cmin, self.cmax,
                                         valinit=2)

        # make a rectangular selector
        self.rect_selector = widgets.RectangleSelector(self.ax1,
                                                       self.rect_onselect,
                                                       drawtype='box',
                                                       useblit=True)

        plt.show()

        # needs to go after show()
        self.slider_res.on_changed(self.set_res_value)
        # self.radio_res.on_clicked(self.set_res_value)

    def redraw_plot(self):
        """
        redraws the plot
        """

        current_xlimits = self.ax1.get_xlim()
        current_ylimits = self.ax1.get_ylim()

        self.ax1.cla()

        plot_res = np.log10(self.res_model[:, :, self.depth_index].T)

        self.mesh_plot = self.ax1.pcolormesh(self.mesh_east,
                                             self.mesh_north,
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)

        # plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax1.text(ee / self.dscale, nn / self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size': self.font_size - 2,
                                        'weight': 'bold'})

        # set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(current_xlimits)

        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(current_ylimits)

        self.ax1.set_ylabel('Northing (' + self.map_scale + ')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting (' + self.map_scale + ')',
                            fontdict=self.fdict)

        depth_title = self.grid_z[self.depth_index] / self.dscale

        self.ax1.set_title('Depth = {:.3f} '.format(depth_title) + \
                           '(' + self.map_scale + ')',
                           fontdict=self.fdict)

        # plot finite element mesh
        self.ax1.plot(self.east_line_xlist,
                      self.east_line_ylist,
                      lw=.25,
                      color='k')

        self.ax1.plot(self.north_line_xlist,
                      self.north_line_ylist,
                      lw=.25,
                      color='k')

        # be sure to redraw the canvas
        self.fig.canvas.draw()

    #    def set_res_value(self, label):
    #        self.res_value = float(label)
    #        print 'set resistivity to ', label
    #        print self.res_value
    def set_res_value(self, val):
        self.res_value = 10 ** val
        print('set resistivity to ', self.res_value)

    def _on_key_callback(self, event):
        """
        on pressing a key do something

        """

        self.event_change_depth = event

        # go down a layer on push of +/= keys
        if self.event_change_depth.key == '=':
            self.depth_index += 1

            if self.depth_index > len(self.grid_z) - 1:
                self.depth_index = len(self.grid_z) - 1
                print('already at deepest depth')

            print('Plotting Depth {0:.3f}'.format(self.grid_z[self.depth_index] / \
                                                  self.dscale) + '(' + self.map_scale + ')')

            self.redraw_plot()
        # go up a layer on push of - key
        elif self.event_change_depth.key == '-':
            self.depth_index -= 1

            if self.depth_index < 0:
                self.depth_index = 0

            print('Plotting Depth {0:.3f} '.format(self.grid_z[self.depth_index] / \
                                                   self.dscale) + '(' + self.map_scale + ')')

            self.redraw_plot()

        # exit plot on press of q
        elif self.event_change_depth.key == 'q':
            self.event_change_depth.canvas.mpl_disconnect(self.cid_depth)
            plt.close(self.event_change_depth.canvas.figure)
            self.rewrite_model_file()

        # copy the layer above
        elif self.event_change_depth.key == 'a':
            try:
                if self.depth_index == 0:
                    print('No layers above')
                else:
                    self.res_model[:, :, self.depth_index] = \
                        self.res_model[:, :, self.depth_index - 1]
            except IndexError:
                print('No layers above')

            self.redraw_plot()

        # copy the layer below
        elif self.event_change_depth.key == 'b':
            try:
                self.res_model[:, :, self.depth_index] = \
                    self.res_model[:, :, self.depth_index + 1]
            except IndexError:
                print('No more layers below')

            self.redraw_plot()

            # undo
        elif self.event_change_depth.key == 'u':
            if type(self.xchange) is int and type(self.ychange) is int:
                self.res_model[self.ychange, self.xchange, self.depth_index] = \
                    self.res_copy[self.ychange, self.xchange, self.depth_index]
            else:
                for xx in self.xchange:
                    for yy in self.ychange:
                        self.res_model[yy, xx, self.depth_index] = \
                            self.res_copy[yy, xx, self.depth_index]

            self.redraw_plot()

    def change_model_res(self, xchange, ychange):
        """
        change resistivity values of resistivity model

        """
        if type(xchange) is int and type(ychange) is int:
            self.res_model[ychange, xchange, self.depth_index] = self.res_value
        else:
            for xx in xchange:
                for yy in ychange:
                    self.res_model[yy, xx, self.depth_index] = self.res_value

        self.redraw_plot()

    def rect_onselect(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        self.xchange = self._get_east_index(x1, x2)
        self.ychange = self._get_north_index(y1, y2)

        # reset values of resistivity
        self.change_model_res(self.xchange, self.ychange)

    def _get_east_index(self, x1, x2):
        """
        get the index value of the points to be changed

        """
        if x1 < x2:
            xchange = np.where((self.grid_east / self.dscale >= x1) & \
                               (self.grid_east / self.dscale <= x2))[0]
            if len(xchange) == 0:
                xchange = np.where(self.grid_east / self.dscale >= x1)[0][0] - 1
                return [xchange]

        if x1 > x2:
            xchange = np.where((self.grid_east / self.dscale <= x1) & \
                               (self.grid_east / self.dscale >= x2))[0]
            if len(xchange) == 0:
                xchange = np.where(self.grid_east / self.dscale >= x2)[0][0] - 1
                return [xchange]

        # check the edges to see if the selection should include the square
        xchange = np.append(xchange, xchange[0] - 1)
        xchange.sort()

        return xchange

    def _get_north_index(self, y1, y2):
        """
        get the index value of the points to be changed in north direction

        need to flip the index because the plot is flipped

        """

        if y1 < y2:
            ychange = np.where((self.grid_north / self.dscale > y1) & \
                               (self.grid_north / self.dscale < y2))[0]
            if len(ychange) == 0:
                ychange = np.where(self.grid_north / self.dscale >= y1)[0][0] - 1
                return [ychange]

        elif y1 > y2:
            ychange = np.where((self.grid_north / self.dscale < y1) & \
                               (self.grid_north / self.dscale > y2))[0]
            if len(ychange) == 0:
                ychange = np.where(self.grid_north / self.dscale >= y2)[0][0] - 1
                return [ychange]

        ychange -= 1
        ychange = np.append(ychange, ychange[-1] + 1)

        return ychange

    def rewrite_model_file(self, model_fn=None, save_path=None,
                           model_fn_basename=None):
        """
        write an initial file for wsinv3d from the model created.
        """
        if save_path is not None:
            self.save_path = save_path

        self.model_fn = model_fn

        if model_fn_basename is not None:
            self.model_fn_basename = model_fn_basename

        self.write_model_file()


