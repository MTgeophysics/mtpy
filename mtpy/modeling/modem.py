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

import matplotlib.cm as cm
import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
# Plotting tools
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import mtpy.analysis.pt as mtpt
# mtpy modules
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtplottools
import mtpy.modeling.ws3dinv as ws
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
        print 'set resistivity to ', self.res_value

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
                print 'already at deepest depth'

            print 'Plotting Depth {0:.3f}'.format(self.grid_z[self.depth_index] / \
                                                  self.dscale) + '(' + self.map_scale + ')'

            self.redraw_plot()
        # go up a layer on push of - key
        elif self.event_change_depth.key == '-':
            self.depth_index -= 1

            if self.depth_index < 0:
                self.depth_index = 0

            print 'Plotting Depth {0:.3f} '.format(self.grid_z[self.depth_index] / \
                                                   self.dscale) + '(' + self.map_scale + ')'

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
                    print 'No layers above'
                else:
                    self.res_model[:, :, self.depth_index] = \
                        self.res_model[:, :, self.depth_index - 1]
            except IndexError:
                print 'No layers above'

            self.redraw_plot()

        # copy the layer below
        elif self.event_change_depth.key == 'b':
            try:
                self.res_model[:, :, self.depth_index] = \
                    self.res_model[:, :, self.depth_index + 1]
            except IndexError:
                print 'No more layers below'

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


# ==============================================================================
# plot response       
# ==============================================================================
class PlotResponse(object):
    """
    plot data and response 
    
    Plots the real and imaginary impedance and induction vector if present.
    
    :Example: ::
    
        >>> import mtpy.modeling.new_modem as modem
        >>> dfn = r"/home/MT/ModEM/Inv1/DataFile.dat"
        >>> rfn = r"/home/MT/ModEM/Inv1/Test_resp_000.dat"
        >>> mrp = modem.PlotResponse(data_fn=dfn, resp_fn=rfn)
        >>> # plot only the TE and TM modes
        >>> mrp.plot_component = 2
        >>> mrp.redraw_plot()
    
    ======================== ==================================================
    Attributes               Description    
    ======================== ==================================================
    color_mode               [ 'color' | 'bw' ] color or black and white plots
    cted                     color for data TE mode
    ctem                     color for data TM mode
    ctmd                     color for model TE mode
    ctmm                     color for model TM mode
    data_fn                  full path to data file
    data_object              WSResponse instance
    e_capsize                cap size of error bars in points (*default* is .5)
    e_capthick               cap thickness of error bars in points (*default*
                             is 1)
    fig_dpi                  resolution of figure in dots-per-inch (300)
    fig_list                 list of matplotlib.figure instances for plots
    fig_size                 size of figure in inches (*default* is [6, 6])
    font_size                size of font for tick labels, axes labels are
                             font_size+2 (*default* is 7)
    legend_border_axes_pad   padding between legend box and axes 
    legend_border_pad        padding between border of legend and symbols
    legend_handle_text_pad   padding between text labels and symbols of legend
    legend_label_spacing     padding between labels
    legend_loc               location of legend 
    legend_marker_scale      scale of symbols in legend
    lw                       line width response curves (*default* is .5)
    ms                       size of markers (*default* is 1.5)
    mted                     marker for data TE mode
    mtem                     marker for data TM mode
    mtmd                     marker for model TE mode
    mtmm                     marker for model TM mode 
    phase_limits             limits of phase
    plot_component           [ 2 | 4 ] 2 for TE and TM or 4 for all components
    plot_style               [ 1 | 2 ] 1 to plot each mode in a seperate
                             subplot and 2 to plot xx, xy and yx, yy in same 
                             plots
    plot_type                [ '1' | list of station name ] '1' to plot all 
                             stations in data file or input a list of station
                             names to plot if station_fn is input, otherwise
                             input a list of integers associated with the 
                             index with in the data file, ie 2 for 2nd station
    plot_z                   [ True | False ] *default* is True to plot 
                             impedance, False for plotting resistivity and 
                             phase
    plot_yn                  [ 'n' | 'y' ] to plot on instantiation
    res_limits               limits of resistivity in linear scale
    resp_fn                  full path to response file
    resp_object              WSResponse object for resp_fn, or list of 
                             WSResponse objects if resp_fn is a list of
                             response files
    station_fn               full path to station file written by WSStation
    subplot_bottom           space between axes and bottom of figure
    subplot_hspace           space between subplots in vertical direction
    subplot_left             space between axes and left of figure
    subplot_right            space between axes and right of figure
    subplot_top              space between axes and top of figure
    subplot_wspace           space between subplots in horizontal direction    
    ======================== ==================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, **kwargs):
        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.data_object = None
        self.resp_object = []

        self.color_mode = kwargs.pop('color_mode', 'color')

        self.ms = kwargs.pop('ms', 1.5)
        self.ms_r = kwargs.pop('ms_r', 3)
        self.lw = kwargs.pop('lw', .5)
        self.lw_r = kwargs.pop('lw_r', 1.0)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')

        # black and white mode
        elif self.color_mode == 'bw':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')

        self.phase_limits_d = kwargs.pop('phase_limits_d', None)
        self.phase_limits_od = kwargs.pop('phase_limits_od', None)
        self.res_limits_d = kwargs.pop('res_limits_d', None)
        self.res_limits_od = kwargs.pop('res_limits_od', None)
        self.tipper_limits = kwargs.pop('tipper_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .85)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)

        self.legend_loc = 'upper center'
        self.legend_pos = (.5, 1.18)
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_style = kwargs.pop('plot_style', 1)
        self.plot_component = kwargs.pop('plot_component', 4)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.plot_z = kwargs.pop('plot_z', True)
        self.ylabel_pad = kwargs.pop('ylabel_pad', 1.25)

        self.fig_list = []

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plot
        """

        self.data_object = Data()
        self.data_object.read_data_file(self.data_fn)

        # get shape of impedance tensors
        ns = len(self.data_object.mt_dict.keys())

        # read in response files
        if self.resp_fn != None:
            self.resp_object = []
            if type(self.resp_fn) is not list:
                resp_obj = Data()
                resp_obj.read_data_file(self.resp_fn)
                self.resp_object = [resp_obj]
            else:
                for rfile in self.resp_fn:
                    resp_obj = Data()
                    resp_obj.read_data_file(rfile)
                    self.resp_object.append(resp_obj)

        # get number of response files
        nr = len(self.resp_object)

        if type(self.plot_type) is list:
            ns = len(self.plot_type)

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size': self.font_size + 2, 'weight': 'bold'}
        if self.plot_z == True:
            h_ratio = [1, 1, .5]
        elif self.plot_z == False:
            h_ratio = [1.5, 1, .5]

        ax_list = []
        line_list = []
        label_list = []

        # --> make key word dictionaries for plotting
        kw_xx = {'color': self.cted,
                 'marker': self.mted,
                 'ms': self.ms,
                 'ls': ':',
                 'lw': self.lw,
                 'e_capsize': self.e_capsize,
                 'e_capthick': self.e_capthick}

        kw_yy = {'color': self.ctmd,
                 'marker': self.mtmd,
                 'ms': self.ms,
                 'ls': ':',
                 'lw': self.lw,
                 'e_capsize': self.e_capsize,
                 'e_capthick': self.e_capthick}

        if self.plot_type != '1':
            pstation_list = []
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]
            for ii, station in enumerate(self.data_object.mt_dict.keys()):
                if type(station) is not int:
                    for pstation in self.plot_type:
                        if station.find(str(pstation)) >= 0:
                            pstation_list.append(station)
                else:
                    for pstation in self.plot_type:
                        if station == int(pstation):
                            pstation_list.append(ii)
        else:
            pstation_list = self.data_object.mt_dict.keys()

        for jj, station in enumerate(pstation_list):
            z_obj = self.data_object.mt_dict[station].Z
            t_obj = self.data_object.mt_dict[station].Tipper
            period = self.data_object.period_list
            print 'Plotting: {0}'.format(station)

            # convert to apparent resistivity and phase
            z_obj.compute_resistivity_phase()

            # find locations where points have been masked
            nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
            nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
            nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
            nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
            ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
            nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]

            # convert to apparent resistivity and phase
            if self.plot_z == True:
                scaling = np.zeros_like(z_obj.z)
                for ii in range(2):
                    for jj in range(2):
                        scaling[:, ii, jj] = 1. / np.sqrt(z_obj.freq)
                plot_res = abs(z_obj.z.real * scaling)
                plot_res_err = abs(z_obj.z_err * scaling)
                plot_phase = abs(z_obj.z.imag * scaling)
                plot_phase_err = abs(z_obj.z_err * scaling)
                h_ratio = [1, 1, .5]

            elif self.plot_z == False:
                plot_res = z_obj.resistivity
                plot_res_err = z_obj.resistivity_err
                plot_phase = z_obj.phase
                plot_phase_err = z_obj.phase_err
                h_ratio = [1.5, 1, .5]

                try:
                    self.res_limits_d = (10 ** (np.floor(np.log10(min([plot_res[nzxx, 0, 0].min(),
                                                                       plot_res[nzyy, 1, 1].min()])))),
                                         10 ** (np.ceil(np.log10(max([plot_res[nzxx, 0, 0].max(),
                                                                      plot_res[nzyy, 1, 1].max()])))))
                except ValueError:
                    self.res_limits_d = None
                try:
                    self.res_limits_od = (10 ** (np.floor(np.log10(min([plot_res[nzxy, 0, 1].min(),
                                                                        plot_res[nzyx, 1, 0].min()])))),
                                          10 ** (np.ceil(np.log10(max([plot_res[nzxy, 0, 1].max(),
                                                                       plot_res[nzyx, 1, 0].max()])))))
                except ValueError:
                    self.res_limits_od = None

            # make figure
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            fig.suptitle(str(station), fontdict=fontdict)

            # set the grid of subplots
            if np.all(t_obj.tipper == 0.0) == True:
                self.plot_tipper = False
            else:
                self.plot_tipper = True
                self.tipper_limits = (np.round(min([t_obj.tipper[ntx, 0, 0].real.min(),
                                                    t_obj.tipper[nty, 0, 1].real.min(),
                                                    t_obj.tipper[ntx, 0, 0].imag.min(),
                                                    t_obj.tipper[nty, 0, 1].imag.min()]),
                                               1),
                                      np.round(max([t_obj.tipper[ntx, 0, 0].real.max(),
                                                    t_obj.tipper[nty, 0, 1].real.max(),
                                                    t_obj.tipper[ntx, 0, 0].imag.max(),
                                                    t_obj.tipper[nty, 0, 1].imag.max()]),
                                               1))

            gs = gridspec.GridSpec(3, 4,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=h_ratio)

            axrxx = fig.add_subplot(gs[0, 0])
            axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
            axryx = fig.add_subplot(gs[0, 2], sharex=axrxx, sharey=axrxy)
            axryy = fig.add_subplot(gs[0, 3], sharex=axrxx, sharey=axrxx)

            axpxx = fig.add_subplot(gs[1, 0])
            axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
            axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
            axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)

            axtxr = fig.add_subplot(gs[2, 0], sharex=axrxx)
            axtxi = fig.add_subplot(gs[2, 1], sharex=axrxx, sharey=axtxr)
            axtyr = fig.add_subplot(gs[2, 2], sharex=axrxx)
            axtyi = fig.add_subplot(gs[2, 3], sharex=axrxx, sharey=axtyr)

            self.ax_list = [axrxx, axrxy, axryx, axryy,
                            axpxx, axpxy, axpyx, axpyy,
                            axtxr, axtxi, axtyr, axtyi]

            # ---------plot the apparent resistivity-----------------------------------
            # plot each component in its own subplot
            # plot data response
            erxx = mtplottools.plot_errorbar(axrxx,
                                             period[nzxx],
                                             plot_res[nzxx, 0, 0],
                                             plot_res_err[nzxx, 0, 0],
                                             **kw_xx)
            erxy = mtplottools.plot_errorbar(axrxy,
                                             period[nzxy],
                                             plot_res[nzxy, 0, 1],
                                             plot_res_err[nzxy, 0, 1],
                                             **kw_xx)
            eryx = mtplottools.plot_errorbar(axryx,
                                             period[nzyx],
                                             plot_res[nzyx, 1, 0],
                                             plot_res_err[nzyx, 1, 0],
                                             **kw_yy)
            eryy = mtplottools.plot_errorbar(axryy,
                                             period[nzyy],
                                             plot_res[nzyy, 1, 1],
                                             plot_res_err[nzyy, 1, 1],
                                             **kw_yy)
            # plot phase
            epxx = mtplottools.plot_errorbar(axpxx,
                                             period[nzxx],
                                             plot_phase[nzxx, 0, 0],
                                             plot_phase_err[nzxx, 0, 0],
                                             **kw_xx)
            epxy = mtplottools.plot_errorbar(axpxy,
                                             period[nzxy],
                                             plot_phase[nzxy, 0, 1],
                                             plot_phase_err[nzxy, 0, 1],
                                             **kw_xx)
            epyx = mtplottools.plot_errorbar(axpyx,
                                             period[nzyx],
                                             plot_phase[nzyx, 1, 0],
                                             plot_phase_err[nzyx, 1, 0],
                                             **kw_yy)
            epyy = mtplottools.plot_errorbar(axpyy,
                                             period[nzyy],
                                             plot_phase[nzyy, 1, 1],
                                             plot_phase_err[nzyy, 1, 1],
                                             **kw_yy)

            # plot tipper
            if self.plot_tipper == True:
                ertx = mtplottools.plot_errorbar(axtxr,
                                                 period[ntx],
                                                 t_obj.tipper[ntx, 0, 0].real,
                                                 t_obj.tipper_err[ntx, 0, 0],
                                                 **kw_xx)
                erty = mtplottools.plot_errorbar(axtyr,
                                                 period[nty],
                                                 t_obj.tipper[nty, 0, 1].real,
                                                 t_obj.tipper_err[nty, 0, 1],
                                                 **kw_yy)

                eptx = mtplottools.plot_errorbar(axtxi,
                                                 period[ntx],
                                                 t_obj.tipper[ntx, 0, 0].imag,
                                                 t_obj.tipper_err[ntx, 0, 0],
                                                 **kw_xx)
                epty = mtplottools.plot_errorbar(axtyi,
                                                 period[nty],
                                                 t_obj.tipper[nty, 0, 1].imag,
                                                 t_obj.tipper_err[nty, 0, 1],
                                                 **kw_yy)

            # ----------------------------------------------
            # get error bar list for editing later        
            if self.plot_tipper == False:
                try:
                    self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                      [erxy[1][0], erxy[1][1], erxy[2][0]],
                                      [eryx[1][0], eryx[1][1], eryx[2][0]],
                                      [eryy[1][0], eryy[1][1], eryy[2][0]]]
                    line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
                except IndexError:
                    print 'Found no Z components for {0}'.format(self.station)
                    line_list = [[None], [None],
                                 [None], [None]]

                    self._err_list = [[None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None]]

            else:
                try:
                    line_list = [[erxx[0]], [erxy[0]],
                                 [eryx[0]], [eryy[0]],
                                 [ertx[0]], [erty[0]]]

                    self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                      [erxy[1][0], erxy[1][1], erxy[2][0]],
                                      [eryx[1][0], eryx[1][1], eryx[2][0]],
                                      [eryy[1][0], eryy[1][1], eryy[2][0]],
                                      [ertx[1][0], ertx[1][1], ertx[2][0]],
                                      [erty[1][0], erty[1][1], erty[2][0]]]
                except IndexError:
                    print 'Found no Z components for {0}'.format(station)
                    line_list = [[None], [None],
                                 [None], [None],
                                 [None], [None]]

                    self._err_list = [[None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None]]
            # ------------------------------------------
            # make things look nice        
            # set titles of the Z components
            label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                          ['$Z_{yx}$'], ['$Z_{yy}$']]
            for ax, label in zip(self.ax_list[0:4], label_list):
                ax.set_title(label[0], fontdict={'size': self.font_size + 2,
                                                 'weight': 'bold'})

                # set legends for tipper components
            # fake a line
            l1 = plt.Line2D([0], [0], linewidth=0, color='w', linestyle='None',
                            marker='.')
            t_label_list = ['Re{$T_x$}', 'Im{$T_x$}', 'Re{$T_y$}', 'Im{$T_y$}']
            label_list += [['$T_{x}$'], ['$T_{y}$']]
            for ax, label in zip(self.ax_list[-4:], t_label_list):
                ax.legend([l1], [label], loc='upper left',
                          markerscale=.01,
                          borderaxespad=.05,
                          labelspacing=.01,
                          handletextpad=.05,
                          borderpad=.05,
                          prop={'size': max([self.font_size, 6])})



                # set axis properties
            for aa, ax in enumerate(self.ax_list):
                ax.tick_params(axis='y', pad=self.ylabel_pad)

                if aa < 8:
                    #                    ylabels[-1] = ''
                    #                    ylabels[0] = ''
                    #                    ax.set_yticklabels(ylabels)
                    #                    plt.setp(ax.get_xticklabels(), visible=False)
                    if self.plot_z == True:
                        ax.set_yscale('log', nonposy='clip')

                else:
                    ax.set_xlabel('Period (s)', fontdict=fontdict)

                if aa < 4 and self.plot_z is False:
                    ax.set_yscale('log', nonposy='clip')
                    if aa == 0 or aa == 3:
                        ax.set_ylim(self.res_limits_d)
                    elif aa == 1 or aa == 2:
                        ax.set_ylim(self.res_limits_od)

                if aa > 3 and aa < 8 and self.plot_z is False:
                    ax.yaxis.set_major_formatter(MultipleLocator(10))
                    if self.phase_limits_d is not None:
                        ax.set_ylim(self.phase_limits_d)
                # set axes labels
                if aa == 0:
                    if self.plot_z == False:
                        ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                      fontdict=fontdict)
                    elif self.plot_z == True:
                        ax.set_ylabel('Re[Z (mV/km nT)]',
                                      fontdict=fontdict)
                elif aa == 4:
                    if self.plot_z == False:
                        ax.set_ylabel('Phase (deg)',
                                      fontdict=fontdict)
                    elif self.plot_z == True:
                        ax.set_ylabel('Im[Z (mV/km nT)]',
                                      fontdict=fontdict)
                elif aa == 8:
                    ax.set_ylabel('Tipper',
                                  fontdict=fontdict)

                if aa > 7:
                    ax.yaxis.set_major_locator(MultipleLocator(.1))
                    if self.tipper_limits is not None:
                        ax.set_ylim(self.tipper_limits)
                    else:
                        pass

                ax.set_xscale('log', nonposx='clip')
                ax.set_xlim(xmin=10 ** (np.floor(np.log10(period[0]))) * 1.01,
                            xmax=10 ** (np.ceil(np.log10(period[-1]))) * .99)
                ax.grid(True, alpha=.25)

                ylabels = ax.get_yticks().tolist()
                if aa < 8:
                    ylabels[-1] = ''
                    ylabels[0] = ''
                    ax.set_yticklabels(ylabels)
                    plt.setp(ax.get_xticklabels(), visible=False)


                    ##----------------------------------------------
            # plot model response
            if self.resp_object is not None:
                for resp_obj in self.resp_object:
                    resp_z_obj = resp_obj.mt_dict[station].Z
                    resp_z_err = np.nan_to_num((z_obj.z - resp_z_obj.z) / z_obj.z_err)
                    resp_z_obj.compute_resistivity_phase()

                    resp_t_obj = resp_obj.mt_dict[station].Tipper
                    resp_t_err = np.nan_to_num((t_obj.tipper - resp_t_obj.tipper) / t_obj.tipper_err)

                    # convert to apparent resistivity and phase
                    if self.plot_z == True:
                        scaling = np.zeros_like(resp_z_obj.z)
                        for ii in range(2):
                            for jj in range(2):
                                scaling[:, ii, jj] = 1. / np.sqrt(resp_z_obj.freq)
                        r_plot_res = abs(resp_z_obj.z.real * scaling)
                        r_plot_phase = abs(resp_z_obj.z.imag * scaling)

                    elif self.plot_z == False:
                        r_plot_res = resp_z_obj.resistivity
                        r_plot_phase = resp_z_obj.phase

                    rms_xx = resp_z_err[:, 0, 0].std()
                    rms_xy = resp_z_err[:, 0, 1].std()
                    rms_yx = resp_z_err[:, 1, 0].std()
                    rms_yy = resp_z_err[:, 1, 1].std()

                    # --> make key word dictionaries for plotting
                    kw_xx = {'color': self.ctem,
                             'marker': self.mtem,
                             'ms': self.ms_r,
                             'ls': ':',
                             'lw': self.lw_r,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    kw_yy = {'color': self.ctmm,
                             'marker': self.mtmm,
                             'ms': self.ms_r,
                             'ls': ':',
                             'lw': self.lw_r,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    # plot data response
                    rerxx = mtplottools.plot_errorbar(axrxx,
                                                      period[nzxx],
                                                      r_plot_res[nzxx, 0, 0],
                                                      None,
                                                      **kw_xx)
                    rerxy = mtplottools.plot_errorbar(axrxy,
                                                      period[nzxy],
                                                      r_plot_res[nzxy, 0, 1],
                                                      None,
                                                      **kw_xx)
                    reryx = mtplottools.plot_errorbar(axryx,
                                                      period[nzyx],
                                                      r_plot_res[nzyx, 1, 0],
                                                      None,
                                                      **kw_yy)
                    reryy = mtplottools.plot_errorbar(axryy,
                                                      period[nzyy],
                                                      r_plot_res[nzyy, 1, 1],
                                                      None,
                                                      **kw_yy)
                    # plot phase
                    repxx = mtplottools.plot_errorbar(axpxx,
                                                      period[nzxx],
                                                      r_plot_phase[nzxx, 0, 0],
                                                      None,
                                                      **kw_xx)
                    repxy = mtplottools.plot_errorbar(axpxy,
                                                      period[nzxy],
                                                      r_plot_phase[nzxy, 0, 1],
                                                      None,
                                                      **kw_xx)
                    repyx = mtplottools.plot_errorbar(axpyx,
                                                      period[nzyx],
                                                      r_plot_phase[nzyx, 1, 0],
                                                      None,
                                                      **kw_yy)
                    repyy = mtplottools.plot_errorbar(axpyy,
                                                      period[nzyy],
                                                      r_plot_phase[nzyy, 1, 1],
                                                      None,
                                                      **kw_yy)

                    # plot tipper
                    if self.plot_tipper == True:
                        rertx = mtplottools.plot_errorbar(axtxr,
                                                          period[ntx],
                                                          resp_t_obj.tipper[ntx, 0, 0].real,
                                                          None,
                                                          **kw_xx)
                        rerty = mtplottools.plot_errorbar(axtyr,
                                                          period[nty],
                                                          resp_t_obj.tipper[nty, 0, 1].real,
                                                          None,
                                                          **kw_yy)

                        reptx = mtplottools.plot_errorbar(axtxi,
                                                          period[ntx],
                                                          resp_t_obj.tipper[ntx, 0, 0].imag,
                                                          None,
                                                          **kw_xx)
                        repty = mtplottools.plot_errorbar(axtyi,
                                                          period[nty],
                                                          resp_t_obj.tipper[nty, 0, 1].imag,
                                                          None,
                                                          **kw_yy)

                    if self.plot_tipper == False:
                        line_list[0] += [rerxx[0]]
                        line_list[1] += [rerxy[0]]
                        line_list[2] += [reryx[0]]
                        line_list[3] += [reryy[0]]
                        label_list[0] += ['$Z^m_{xx}$ ' +
                                          'rms={0:.2f}'.format(rms_xx)]
                        label_list[1] += ['$Z^m_{xy}$ ' +
                                          'rms={0:.2f}'.format(rms_xy)]
                        label_list[2] += ['$Z^m_{yx}$ ' +
                                          'rms={0:.2f}'.format(rms_yx)]
                        label_list[3] += ['$Z^m_{yy}$ ' +
                                          'rms={0:.2f}'.format(rms_yy)]
                    else:
                        line_list[0] += [rerxx[0]]
                        line_list[1] += [rerxy[0]]
                        line_list[2] += [reryx[0]]
                        line_list[3] += [reryy[0]]
                        line_list[4] += [rertx[0]]
                        line_list[5] += [rerty[0]]
                        label_list[0] += ['$Z^m_{xx}$ ' +
                                          'rms={0:.2f}'.format(rms_xx)]
                        label_list[1] += ['$Z^m_{xy}$ ' +
                                          'rms={0:.2f}'.format(rms_xy)]
                        label_list[2] += ['$Z^m_{yx}$ ' +
                                          'rms={0:.2f}'.format(rms_yx)]
                        label_list[3] += ['$Z^m_{yy}$ ' +
                                          'rms={0:.2f}'.format(rms_yy)]
                        label_list[4] += ['$T^m_{x}$ ' +
                                          'rms={0:.2f}'.format(resp_t_err[:, 0, 0].std())]
                        label_list[5] += ['$T^m_{y}$' +
                                          'rms={0:.2f}'.format(resp_t_err[:, 0, 1].std())]

                legend_ax_list = self.ax_list[0:4]
                #                if self.plot_tipper == True:
                #                    legend_ax_list += [self.ax_list[-4], self.ax_list[-2]]

                for aa, ax in enumerate(legend_ax_list):
                    ax.legend(line_list[aa],
                              label_list[aa],
                              loc=self.legend_loc,
                              bbox_to_anchor=self.legend_pos,
                              markerscale=self.legend_marker_scale,
                              borderaxespad=self.legend_border_axes_pad,
                              labelspacing=self.legend_label_spacing,
                              handletextpad=self.legend_handle_text_pad,
                              borderpad=self.legend_border_pad,
                              prop={'size': max([self.font_size, 5])})

            plt.show()

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

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                    fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.
        
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

        fig = plt.gcf()
        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_L2.' +
                                   file_format)
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(fig)

        else:
            pass

        self.fig_fn = save_fn
        print 'Saved figure to: ' + self.fig_fn

    def update_plot(self):
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

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots data vs model response computed by WS3DINV")


# ==============================================================================
# plot phase tensors
# ==============================================================================
class PlotPTMaps(mtplottools.MTEllipse):
    """
    Plot phase tensor maps including residual pt if response file is input.
    
    :Plot only data for one period: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, plot_period_list=[0])
        
    :Plot data and model response: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> rfn = r"/home/MT/ws3dinv/Inv1/Test_resp.00"
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, resp_fn=rfn, model_fn=mfn,
        >>> ...                 plot_period_list=[0])
        >>> # adjust colorbar
        >>> ptm.cb_res_pad = 1.25
        >>> ptm.redraw_plot()
    
    
    ========================== ================================================
    Attributes                 Description    
    ========================== ================================================
    cb_pt_pad                  percentage from top of axes to place pt 
                               color bar. *default* is .90
    cb_res_pad                 percentage from bottom of axes to place
                               resistivity color bar. *default* is 1.2
    cb_residual_tick_step      tick step for residual pt. *default* is 3
    cb_tick_step               tick step for phase tensor color bar, 
                               *default* is 45
    data                       np.ndarray(n_station, n_periods, 2, 2)
                               impedance tensors for station data                     
    data_fn                    full path to data fle               
    dscale                     scaling parameter depending on map_scale
    ellipse_cmap               color map for pt ellipses. *default* is
                               mt_bl2gr2rd
    ellipse_colorby            [ 'skew' | 'skew_seg' | 'phimin' | 'phimax'|
                                 'phidet' | 'ellipticity' ] parameter to color
                                 ellipses by. *default* is 'phimin'
    ellipse_range              (min, max, step) min and max of colormap, need
                               to input step if plotting skew_seg
    ellipse_size               relative size of ellipses in map_scale
    ew_limits                  limits of plot in e-w direction in map_scale
                               units.  *default* is None, scales to station 
                               area
    fig_aspect                 aspect of figure. *default* is 1
    fig_dpi                    resolution in dots-per-inch. *default* is 300
    fig_list                   list of matplotlib.figure instances for each
                               figure plotted.
    fig_size                   [width, height] in inches of figure window
                               *default* is [6, 6]
    font_size                  font size of ticklabels, axes labels are 
                               font_size+2. *default* is 7
    grid_east                  relative location of grid nodes in e-w direction
                               in map_scale units
    grid_north                 relative location of grid nodes in n-s direction
                               in map_scale units
    grid_z                     relative location of grid nodes in z direction
                               in map_scale units
    model_fn                 full path to initial file
    map_scale                  [ 'km' | 'm' ] distance units of map. 
                               *default* is km
    mesh_east                  np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north                 np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                   full path to model file
    nodes_east                 relative distance betwen nodes in e-w direction
                               in map_scale units
    nodes_north                relative distance betwen nodes in n-s direction
                               in map_scale units
    nodes_z                    relative distance betwen nodes in z direction
                               in map_scale units
    ns_limits                  (min, max) limits of plot in n-s direction
                               *default* is None, viewing area is station area
    pad_east                   padding from extreme stations in east direction
    pad_north                  padding from extreme stations in north direction
    period_list                list of periods from data
    plot_grid                  [ 'y' | 'n' ] 'y' to plot grid lines
                               *default* is 'n'
    plot_period_list           list of period index values to plot
                               *default* is None 
    plot_yn                    ['y' | 'n' ] 'y' to plot on instantiation
                               *default* is 'y'
    res_cmap                   colormap for resisitivity values. 
                               *default* is 'jet_r'
    res_limits                 (min, max) resistivity limits in log scale
                               *default* is (0, 4)
    res_model                  np.ndarray(n_north, n_east, n_vertical) of 
                               model resistivity values in linear scale
    residual_cmap              color map for pt residuals. 
                               *default* is 'mt_wh2or' 
    resp                       np.ndarray(n_stations, n_periods, 2, 2)
                               impedance tensors for model response  
    resp_fn                    full path to response file
    save_path                  directory to save figures to
    save_plots                 [ 'y' | 'n' ] 'y' to save plots to save_path
    station_east               location of stations in east direction in 
                               map_scale units  
    station_fn                 full path to station locations file
    station_names              station names
    station_north              location of station in north direction in 
                               map_scale units
    subplot_bottom             distance between axes and bottom of figure window
    subplot_left               distance between axes and left of figure window  
    subplot_right              distance between axes and right of figure window
    subplot_top                distance between axes and top of figure window
    title                      titiel of plot *default* is depth of slice
    xminorticks                location of xminorticks
    yminorticks                location of yminorticks
    ========================== ================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, model_fn=None, **kwargs):

        mtplottools.MTEllipse.__init__(self, **kwargs)

        self.model_fn = model_fn
        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)

        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

        self.save_plots = kwargs.pop('save_plots', 'y')
        self.plot_period_list = kwargs.pop('plot_period_list', None)
        self.period_dict = None

        self.map_scale = kwargs.pop('map_scale', 'km')
        # make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)

        self.pad_east = kwargs.pop('pad_east', 2000)
        self.pad_north = kwargs.pop('pad_north', 2000)

        self.plot_grid = kwargs.pop('plot_grid', 'n')

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []

        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)

        self.residual_cmap = kwargs.pop('residual_cmap', 'mt_wh2or')
        self.font_size = kwargs.pop('font_size', 7)

        self.cb_tick_step = kwargs.pop('cb_tick_step', 45)
        self.cb_residual_tick_step = kwargs.pop('cb_residual_tick_step', 3)
        self.cb_pt_pad = kwargs.pop('cb_pt_pad', 1.2)
        self.cb_res_pad = kwargs.pop('cb_res_pad', .5)

        self.res_limits = kwargs.pop('res_limits', (0, 4))
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')

        # --> set the ellipse properties -------------------
        self._ellipse_dict = kwargs.pop('ellipse_dict', {'size': 2})
        self._read_ellipse_dict(self._ellipse_dict)

        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .2
        self.subplot_wspace = .05

        self.data_obj = None
        self.resp_obj = None
        self.model_obj = None
        self.period_list = None

        self.pt_data_arr = None
        self.pt_resp_arr = None
        self.pt_resid_arr = None

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def _read_files(self):
        """
        get information from files
        """

        # --> read in data file
        self.data_obj = Data()
        self.data_obj.read_data_file(self.data_fn)

        # --> read response file
        if self.resp_fn is not None:
            self.resp_obj = Data()
            self.resp_obj.read_data_file(self.resp_fn)

        # --> read mode file
        if self.model_fn is not None:
            self.model_obj = Model()
            self.model_obj.read_model_file(self.model_fn)

        self._get_plot_period_list()
        self._get_pt()

    def _get_plot_period_list(self):
        """
        get periods to plot from input or data file
        """
        # --> get period list to plot
        if self.plot_period_list is None:
            self.plot_period_list = self.data_obj.period_list
        else:
            if type(self.plot_period_list) is list:
                # check if entries are index values or actual periods
                if type(self.plot_period_list[0]) is int:
                    self.plot_period_list = [self.data_obj.period_list[ii]
                                             for ii in self.plot_period_list]
                else:
                    pass
            elif type(self.plot_period_list) is int:
                self.plot_period_list = self.data_obj.period_list[self.plot_period_list]
            elif type(self.plot_period_list) is float:
                self.plot_period_list = [self.plot_period_list]

        self.period_dict = dict([(key, value) for value, key in
                                 enumerate(self.data_obj.period_list)])

    def _get_pt(self):
        """
        put pt parameters into something useful for plotting
        """

        ns = len(self.data_obj.mt_dict.keys())
        nf = len(self.data_obj.period_list)

        data_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                ('phimax', np.float),
                                                ('skew', np.float),
                                                ('azimuth', np.float),
                                                ('east', np.float),
                                                ('north', np.float)])
        if self.resp_fn is not None:
            model_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('skew', np.float),
                                                     ('azimuth', np.float),
                                                     ('east', np.float),
                                                     ('north', np.float)])

            res_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                   ('phimax', np.float),
                                                   ('skew', np.float),
                                                   ('azimuth', np.float),
                                                   ('east', np.float),
                                                   ('north', np.float),
                                                   ('geometric_mean', np.float)])

        for ii, key in enumerate(self.data_obj.mt_dict.keys()):
            east = self.data_obj.mt_dict[key].grid_east / self.dscale
            north = self.data_obj.mt_dict[key].grid_north / self.dscale
            dpt = self.data_obj.mt_dict[key].pt
            data_pt_arr[:, ii]['east'] = east
            data_pt_arr[:, ii]['north'] = north
            data_pt_arr[:, ii]['phimin'] = dpt.phimin[0]
            data_pt_arr[:, ii]['phimax'] = dpt.phimax[0]
            data_pt_arr[:, ii]['azimuth'] = dpt.azimuth[0]
            data_pt_arr[:, ii]['skew'] = dpt.beta[0]
            if self.resp_fn is not None:
                mpt = self.resp_obj.mt_dict[key].pt
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt,
                                                   pt_object2=mpt)
                    rpt = rpt.residual_pt
                    res_pt_arr[:, ii]['east'] = east
                    res_pt_arr[:, ii]['north'] = north
                    res_pt_arr[:, ii]['phimin'] = rpt.phimin[0]
                    res_pt_arr[:, ii]['phimax'] = rpt.phimax[0]
                    res_pt_arr[:, ii]['azimuth'] = rpt.azimuth[0]
                    res_pt_arr[:, ii]['skew'] = rpt.beta[0]
                    res_pt_arr[:, ii]['geometric_mean'] = np.sqrt(abs(rpt.phimin[0] * \
                                                                      rpt.phimax[0]))
                except mtex.MTpyError_PT:
                    print key, dpt.pt.shape, mpt.pt.shape

                model_pt_arr[:, ii]['east'] = east
                model_pt_arr[:, ii]['north'] = north
                model_pt_arr[:, ii]['phimin'] = mpt.phimin[0]
                model_pt_arr[:, ii]['phimax'] = mpt.phimax[0]
                model_pt_arr[:, ii]['azimuth'] = mpt.azimuth[0]
                model_pt_arr[:, ii]['skew'] = mpt.beta[0]

        # make these attributes
        self.pt_data_arr = data_pt_arr
        if self.resp_fn is not None:
            self.pt_resp_arr = model_pt_arr
            self.pt_resid_arr = res_pt_arr

    def plot(self):
        """
        plot phase tensor maps for data and or response, each figure is of a
        different period.  If response is input a third column is added which is 
        the residual phase tensor showing where the model is not fitting the data 
        well.  The data is plotted in km.
        
        """
        # --> read in data first
        if self.data_obj is None:
            self._read_files()

        # set plot properties
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        font_dict = {'size': self.font_size + 2, 'weight': 'bold'}

        # make a grid of subplots 
        gs = gridspec.GridSpec(1, 3, hspace=self.subplot_hspace,
                               wspace=self.subplot_wspace)

        # set some parameters for the colorbar
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        bounds = np.arange(ckmin, ckmax + ckstep, ckstep)

        # set plot limits to be the station area
        if self.ew_limits == None:
            east_min = self.data_obj.data_array['rel_east'].min() - \
                       self.pad_east
            east_max = self.data_obj.data_array['rel_east'].max() + \
                       self.pad_east
            self.ew_limits = (east_min / self.dscale, east_max / self.dscale)

        if self.ns_limits == None:
            north_min = self.data_obj.data_array['rel_north'].min() - \
                        self.pad_north
            north_max = self.data_obj.data_array['rel_north'].max() + \
                        self.pad_north
            self.ns_limits = (north_min / self.dscale, north_max / self.dscale)

        # -------------plot phase tensors------------------------------------
        for ff, per in enumerate(self.plot_period_list):
            data_ii = self.period_dict[per]

            print 'Plotting Period: {0:.5g}'.format(per)
            fig = plt.figure('{0:.5g}'.format(per), figsize=self.fig_size,
                             dpi=self.fig_dpi)
            fig.clf()

            if self.resp_fn is not None:
                axd = fig.add_subplot(gs[0, 0], aspect='equal')
                axm = fig.add_subplot(gs[0, 1], aspect='equal')
                axr = fig.add_subplot(gs[0, 2], aspect='equal')
                ax_list = [axd, axm, axr]

            else:
                axd = fig.add_subplot(gs[0, :], aspect='equal')
                ax_list = [axd]

            # plot model below the phase tensors
            if self.model_fn is not None:
                approx_depth, d_index = ws.estimate_skin_depth(self.model_obj.res_model.copy(),
                                                               self.model_obj.grid_z.copy() / self.dscale,
                                                               per,
                                                               dscale=self.dscale)
                # need to add an extra row and column to east and north to make sure
                # all is plotted see pcolor for details.
                plot_east = np.append(self.model_obj.grid_east,
                                      self.model_obj.grid_east[-1] * 1.25) / \
                            self.dscale
                plot_north = np.append(self.model_obj.grid_north,
                                       self.model_obj.grid_north[-1] * 1.25) / \
                             self.dscale

                # make a mesh grid for plotting
                # the 'ij' makes sure the resulting grid is in east, north
                self.mesh_east, self.mesh_north = np.meshgrid(plot_east,
                                                              plot_north,
                                                              indexing='ij')

                for ax in ax_list:
                    plot_res = np.log10(self.model_obj.res_model[:, :, d_index].T)
                    ax.pcolormesh(self.mesh_east,
                                  self.mesh_north,
                                  plot_res,
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits[0],
                                  vmax=self.res_limits[1])

            # --> plot data phase tensors
            for pt in self.pt_data_arr[data_ii]:
                eheight = pt['phimin'] / \
                          self.pt_data_arr[data_ii]['phimax'].max() * \
                          self.ellipse_size
                ewidth = pt['phimax'] / \
                         self.pt_data_arr[data_ii]['phimax'].max() * \
                         self.ellipse_size

                ellipse = Ellipse((pt['east'],
                                   pt['north']),
                                  width=ewidth,
                                  height=eheight,
                                  angle=90 - pt['azimuth'])

                # get ellipse color
                if self.ellipse_cmap.find('seg') > 0:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                              self.ellipse_colorby,
                                                              self.ellipse_cmap,
                                                              ckmin,
                                                              ckmax,
                                                              bounds=bounds))
                else:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                              self.ellipse_colorby,
                                                              self.ellipse_cmap,
                                                              ckmin,
                                                              ckmax))

                axd.add_artist(ellipse)

            # -----------plot response phase tensors---------------
            if self.resp_fn is not None:
                rcmin = np.floor(self.pt_resid_arr['geometric_mean'].min())
                rcmax = np.floor(self.pt_resid_arr['geometric_mean'].max())
                for mpt, rpt in zip(self.pt_resp_arr[data_ii],
                                    self.pt_resid_arr[data_ii]):
                    eheight = mpt['phimin'] / \
                              self.pt_resp_arr[data_ii]['phimax'].max() * \
                              self.ellipse_size
                    ewidth = mpt['phimax'] / \
                             self.pt_resp_arr[data_ii]['phimax'].max() * \
                             self.ellipse_size

                    ellipsem = Ellipse((mpt['east'],
                                        mpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=90 - mpt['azimuth'])

                    # get ellipse color
                    if self.ellipse_cmap.find('seg') > 0:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                                   self.ellipse_colorby,
                                                                   self.ellipse_cmap,
                                                                   ckmin,
                                                                   ckmax,
                                                                   bounds=bounds))
                    else:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                                   self.ellipse_colorby,
                                                                   self.ellipse_cmap,
                                                                   ckmin,
                                                                   ckmax))

                    axm.add_artist(ellipsem)

                    # -----------plot residual phase tensors---------------
                    eheight = rpt['phimin'] / \
                              self.pt_resid_arr[data_ii]['phimax'].max() * \
                              self.ellipse_size
                    ewidth = rpt['phimax'] / \
                             self.pt_resid_arr[data_ii]['phimax'].max() * \
                             self.ellipse_size

                    ellipser = Ellipse((rpt['east'],
                                        rpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=rpt['azimuth'])

                    # get ellipse color
                    rpt_color = np.sqrt(abs(rpt['phimin'] * rpt['phimax']))
                    if self.ellipse_cmap.find('seg') > 0:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                                   'geometric_mean',
                                                                   self.residual_cmap,
                                                                   ckmin,
                                                                   ckmax,
                                                                   bounds=bounds))
                    else:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                                   'geometric_mean',
                                                                   self.residual_cmap,
                                                                   ckmin,
                                                                   ckmax))

                    axr.add_artist(ellipser)

            # --> set axes properties
            # data
            axd.set_xlim(self.ew_limits)
            axd.set_ylim(self.ns_limits)
            axd.set_xlabel('Easting ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            axd.set_ylabel('Northing ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            # make a colorbar for phase tensors
            # bb = axd.axes.get_position().bounds
            bb = axd.get_position().bounds
            y1 = .25 * (2 + (self.ns_limits[1] - self.ns_limits[0]) /
                        (self.ew_limits[1] - self.ew_limits[0]))
            cb_location = (3.35 * bb[2] / 5 + bb[0],
                           y1 * self.cb_pt_pad, .295 * bb[2], .02)
            cbaxd = fig.add_axes(cb_location)
            cbd = mcb.ColorbarBase(cbaxd,
                                   cmap=mtcl.cmapdict[self.ellipse_cmap],
                                   norm=Normalize(vmin=ckmin,
                                                  vmax=ckmax),
                                   orientation='horizontal')
            cbd.ax.xaxis.set_label_position('top')
            cbd.ax.xaxis.set_label_coords(.5, 1.75)
            cbd.set_label(mtplottools.ckdict[self.ellipse_colorby])
            cbd.set_ticks(np.arange(ckmin, ckmax + self.cb_tick_step,
                                    self.cb_tick_step))

            axd.text(self.ew_limits[0] * .95,
                     self.ns_limits[1] * .95,
                     'Data',
                     horizontalalignment='left',
                     verticalalignment='top',
                     bbox={'facecolor': 'white'},
                     fontdict={'size': self.font_size + 1})

            # Model and residual
            if self.resp_fn is not None:
                for aa, ax in enumerate([axm, axr]):
                    ax.set_xlim(self.ew_limits)
                    ax.set_ylim(self.ns_limits)
                    ax.set_xlabel('Easting ({0})'.format(self.map_scale),
                                  fontdict=font_dict)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    # make a colorbar ontop of axis
                    bb = ax.axes.get_position().bounds
                    y1 = .25 * (2 + (self.ns_limits[1] - self.ns_limits[0]) /
                                (self.ew_limits[1] - self.ew_limits[0]))
                    cb_location = (3.35 * bb[2] / 5 + bb[0],
                                   y1 * self.cb_pt_pad, .295 * bb[2], .02)
                    cbax = fig.add_axes(cb_location)
                    if aa == 0:
                        cb = mcb.ColorbarBase(cbax,
                                              cmap=mtcl.cmapdict[self.ellipse_cmap],
                                              norm=Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                              orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(mtplottools.ckdict[self.ellipse_colorby])
                        cb.set_ticks(np.arange(ckmin, ckmax + self.cb_tick_step,
                                               self.cb_tick_step))
                        ax.text(self.ew_limits[0] * .95,
                                self.ns_limits[1] * .95,
                                'Model',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor': 'white'},
                                fontdict={'size': self.font_size + 1})
                    else:
                        cb = mcb.ColorbarBase(cbax,
                                              cmap=mtcl.cmapdict[self.residual_cmap],
                                              norm=Normalize(vmin=rcmin,
                                                             vmax=rcmax),
                                              orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(r"$\sqrt{\Phi_{min} \Phi_{max}}$")
                        cb_ticks = [rcmin, (rcmax - rcmin) / 2, rcmax]
                        cb.set_ticks(cb_ticks)
                        ax.text(self.ew_limits[0] * .95,
                                self.ns_limits[1] * .95,
                                'Residual',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor': 'white'},
                                fontdict={'size': self.font_size + 1})

            if self.model_fn is not None:
                for ax in ax_list:
                    ax.tick_params(direction='out')
                    bb = ax.axes.get_position().bounds
                    y1 = .25 * (2 - (self.ns_limits[1] - self.ns_limits[0]) /
                                (self.ew_limits[1] - self.ew_limits[0]))
                    cb_position = (3.0 * bb[2] / 5 + bb[0],
                                   y1 * self.cb_res_pad, .35 * bb[2], .02)
                    cbax = fig.add_axes(cb_position)
                    cb = mcb.ColorbarBase(cbax,
                                          cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits[0],
                                                         vmax=self.res_limits[1]),
                                          orientation='horizontal')
                    cb.ax.xaxis.set_label_position('top')
                    cb.ax.xaxis.set_label_coords(.5, 1.5)
                    cb.set_label('Resistivity ($\Omega \cdot$m)')
                    cb_ticks = np.arange(np.floor(self.res_limits[0]),
                                         np.ceil(self.res_limits[1] + 1), 1)
                    cb.set_ticks(cb_ticks)
                    cb.set_ticklabels([mtplottools.labeldict[ctk] for ctk in cb_ticks])

            plt.show()
            self.fig_list.append(fig)

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

    def save_figure(self, save_path=None, fig_dpi=None, file_format='pdf',
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

        if os.path.isdir(save_path) == False:
            try:
                os.mkdir(save_path)
            except:
                raise IOError('Need to input a correct directory path')

        for fig in self.fig_list:
            per = fig.canvas.get_window_title()
            save_fn = os.path.join(save_path, 'PT_DepthSlice_{0}s.{1}'.format(
                per, file_format))
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

            if close_fig == 'y':
                plt.close(fig)

            else:
                pass

            self.fig_fn = save_fn
            print 'Saved figure to: ' + self.fig_fn


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
