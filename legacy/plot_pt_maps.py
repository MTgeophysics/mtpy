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
from matplotlib.patches import Ellipse

from mtpy.analysis import pt as mtpt
from mtpy.imaging import mtplottools as mtplottools, mtcolors as mtcl
from mtpy.imaging.mtplottools import MTEllipse
from mtpy.modeling import ws3dinv as ws
from mtpy.utils import exceptions as mtex
from mtpy.modeling.modem import Data
from mtpy.modeling.modem import Model

__all__ = ["PlotPTMaps"]


class PlotPTMaps(MTEllipse):
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
        # MTEllipse.__init__(self, **kwargs)
        super(PlotPTMaps, self).__init__(**kwargs)

        self.model_fn = model_fn
        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.save_path = kwargs.pop("save_path", None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)

        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

        self.save_plots = kwargs.pop("save_plots", "y")
        self.plot_period_list = kwargs.pop("plot_period_list", None)
        self.period_dict = None

        self.map_scale = kwargs.pop("map_scale", "km")
        # make map scale
        if self.map_scale == "km":
            self.dscale = 1000.0
        elif self.map_scale == "m":
            self.dscale = 1.0
        self.ew_limits = kwargs.pop("ew_limits", None)
        self.ns_limits = kwargs.pop("ns_limits", None)

        self.pad_east = kwargs.pop("pad_east", 2000)
        self.pad_north = kwargs.pop("pad_north", 2000)

        self.plot_grid = kwargs.pop("plot_grid", "n")

        self.fig_num = kwargs.pop("fig_num", 1)
        self.fig_size = kwargs.pop("fig_size", [6, 6])
        self.fig_dpi = kwargs.pop("dpi", 300)
        self.fig_aspect = kwargs.pop("fig_aspect", 1)
        self.title = kwargs.pop("title", "on")
        self.fig_list = []

        self.xminorticks = kwargs.pop("xminorticks", 1000)
        self.yminorticks = kwargs.pop("yminorticks", 1000)

        self.residual_cmap = kwargs.pop("residual_cmap", "mt_wh2or")
        self.font_size = kwargs.pop("font_size", 7)

        self.cb_tick_step = kwargs.pop("cb_tick_step", 45)
        self.cb_residual_tick_step = kwargs.pop("cb_residual_tick_step", 3)
        self.cb_pt_pad = kwargs.pop("cb_pt_pad", 1.2)
        self.cb_res_pad = kwargs.pop("cb_res_pad", 0.5)

        self.res_limits = kwargs.pop("res_limits", (0, 4))
        self.res_cmap = kwargs.pop("res_cmap", "jet_r")

        # --> set the ellipse properties -------------------
        self._ellipse_dict = kwargs.pop("ellipse_dict", {"size": 2})
        self._read_ellipse_dict(self._ellipse_dict)

        self.subplot_right = 0.99
        self.subplot_left = 0.085
        self.subplot_top = 0.92
        self.subplot_bottom = 0.1
        self.subplot_hspace = 0.2
        self.subplot_wspace = 0.05

        self.data_obj = None
        self.resp_obj = None
        self.model_obj = None
        self.period_list = None

        self.pt_data_arr = None
        self.pt_resp_arr = None
        self.pt_resid_arr = None

        self.plot_yn = kwargs.pop("plot_yn", "y")
        if self.plot_yn == "y":
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
                    self.plot_period_list = [
                        self.data_obj.period_list[ii] for ii in self.plot_period_list
                    ]
                else:
                    pass
            elif type(self.plot_period_list) is int:
                self.plot_period_list = self.data_obj.period_list[self.plot_period_list]
            elif type(self.plot_period_list) is float:
                self.plot_period_list = [self.plot_period_list]

        self.period_dict = dict(
            [(key, value) for value, key in enumerate(self.data_obj.period_list)]
        )

    def _get_pt(self):
        """
        put pt parameters into something useful for plotting
        """

        ns = len(self.data_obj.mt_dict.keys())
        nf = len(self.data_obj.period_list)

        data_pt_arr = np.zeros(
            (nf, ns),
            dtype=[
                ("phimin", np.float),
                ("phimax", np.float),
                ("skew", np.float),
                ("azimuth", np.float),
                ("east", np.float),
                ("north", np.float),
            ],
        )
        if self.resp_fn is not None:
            model_pt_arr = np.zeros(
                (nf, ns),
                dtype=[
                    ("phimin", np.float),
                    ("phimax", np.float),
                    ("skew", np.float),
                    ("azimuth", np.float),
                    ("east", np.float),
                    ("north", np.float),
                ],
            )

            res_pt_arr = np.zeros(
                (nf, ns),
                dtype=[
                    ("phimin", np.float),
                    ("phimax", np.float),
                    ("skew", np.float),
                    ("azimuth", np.float),
                    ("east", np.float),
                    ("north", np.float),
                    ("geometric_mean", np.float),
                ],
            )

        for ii, key in enumerate(self.data_obj.mt_dict.keys()):
            east = self.data_obj.mt_dict[key].grid_east / self.dscale
            north = self.data_obj.mt_dict[key].grid_north / self.dscale
            dpt = self.data_obj.mt_dict[key].pt
            data_pt_arr[:, ii]["east"] = east
            data_pt_arr[:, ii]["north"] = north
            data_pt_arr[:, ii]["phimin"] = dpt.phimin[0]
            data_pt_arr[:, ii]["phimax"] = dpt.phimax[0]
            data_pt_arr[:, ii]["azimuth"] = dpt.azimuth[0]
            data_pt_arr[:, ii]["skew"] = dpt.beta[0]
            if self.resp_fn is not None:
                mpt = self.resp_obj.mt_dict[key].pt
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt, pt_object2=mpt)
                    rpt = rpt.residual_pt
                    res_pt_arr[:, ii]["east"] = east
                    res_pt_arr[:, ii]["north"] = north
                    res_pt_arr[:, ii]["phimin"] = rpt.phimin[0]
                    res_pt_arr[:, ii]["phimax"] = rpt.phimax[0]
                    res_pt_arr[:, ii]["azimuth"] = rpt.azimuth[0]
                    res_pt_arr[:, ii]["skew"] = rpt.beta[0]
                    res_pt_arr[:, ii]["geometric_mean"] = np.sqrt(
                        abs(rpt.phimin[0] * rpt.phimax[0])
                    )
                except mtex.MTpyError_PT:
                    print key, dpt.pt.shape, mpt.pt.shape

                model_pt_arr[:, ii]["east"] = east
                model_pt_arr[:, ii]["north"] = north
                model_pt_arr[:, ii]["phimin"] = mpt.phimin[0]
                model_pt_arr[:, ii]["phimax"] = mpt.phimax[0]
                model_pt_arr[:, ii]["azimuth"] = mpt.azimuth[0]
                model_pt_arr[:, ii]["skew"] = mpt.beta[0]

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
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.left"] = self.subplot_left
        plt.rcParams["figure.subplot.right"] = self.subplot_right
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top
        font_dict = {"size": self.font_size + 2, "weight": "bold"}

        # make a grid of subplots
        gs = gridspec.GridSpec(
            1, 3, hspace=self.subplot_hspace, wspace=self.subplot_wspace
        )

        # set some parameters for the colorbar
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == "mt_seg_bl2wh2rd":
                raise ValueError("Need to input range as (min, max, step)")
            else:
                ckstep = 3
        bounds = np.arange(ckmin, ckmax + ckstep, ckstep)

        # set plot limits to be the station area
        if self.ew_limits == None:
            east_min = self.data_obj.data_array["rel_east"].min() - self.pad_east
            east_max = self.data_obj.data_array["rel_east"].max() + self.pad_east
            self.ew_limits = (east_min / self.dscale, east_max / self.dscale)

        if self.ns_limits == None:
            north_min = self.data_obj.data_array["rel_north"].min() - self.pad_north
            north_max = self.data_obj.data_array["rel_north"].max() + self.pad_north
            self.ns_limits = (north_min / self.dscale, north_max / self.dscale)

        # -------------plot phase tensors------------------------------------
        for ff, per in enumerate(self.plot_period_list):
            data_ii = self.period_dict[per]

            print "Plotting Period: {0:.5g}".format(per)
            fig = plt.figure(
                "{0:.5g}".format(per), figsize=self.fig_size, dpi=self.fig_dpi
            )
            fig.clf()

            if self.resp_fn is not None:
                axd = fig.add_subplot(gs[0, 0], aspect="equal")
                axm = fig.add_subplot(gs[0, 1], aspect="equal")
                axr = fig.add_subplot(gs[0, 2], aspect="equal")
                ax_list = [axd, axm, axr]

            else:
                axd = fig.add_subplot(gs[0, :], aspect="equal")
                ax_list = [axd]

            # plot model below the phase tensors
            if self.model_fn is not None:
                approx_depth, d_index = ws.estimate_skin_depth(
                    self.model_obj.res_model.copy(),
                    self.model_obj.grid_z.copy() / self.dscale,
                    per,
                    dscale=self.dscale,
                )
                # need to add an extra row and column to east and north to make sure
                # all is plotted see pcolor for details.
                plot_east = (
                    np.append(
                        self.model_obj.grid_east, self.model_obj.grid_east[-1] * 1.25
                    )
                    / self.dscale
                )
                plot_north = (
                    np.append(
                        self.model_obj.grid_north, self.model_obj.grid_north[-1] * 1.25
                    )
                    / self.dscale
                )

                # make a mesh grid for plotting
                # the 'ij' makes sure the resulting grid is in east, north
                self.mesh_east, self.mesh_north = np.meshgrid(
                    plot_east, plot_north, indexing="ij"
                )

                for ax in ax_list:
                    plot_res = np.log10(self.model_obj.res_model[:, :, d_index].T)
                    ax.pcolormesh(
                        self.mesh_east,
                        self.mesh_north,
                        plot_res,
                        cmap=self.res_cmap,
                        vmin=self.res_limits[0],
                        vmax=self.res_limits[1],
                    )

            # --> plot data phase tensors
            for pt in self.pt_data_arr[data_ii]:
                eheight = (
                    pt["phimin"]
                    / self.pt_data_arr[data_ii]["phimax"].max()
                    * self.ellipse_size
                )
                ewidth = (
                    pt["phimax"]
                    / self.pt_data_arr[data_ii]["phimax"].max()
                    * self.ellipse_size
                )

                ellipse = Ellipse(
                    (pt["east"], pt["north"]),
                    width=ewidth,
                    height=eheight,
                    angle=90 - pt["azimuth"],
                )

                # get ellipse color
                if self.ellipse_cmap.find("seg") > 0:
                    ellipse.set_facecolor(
                        mtcl.get_plot_color(
                            pt[self.ellipse_colorby],
                            self.ellipse_colorby,
                            self.ellipse_cmap,
                            ckmin,
                            ckmax,
                            bounds=bounds,
                        )
                    )
                else:
                    ellipse.set_facecolor(
                        mtcl.get_plot_color(
                            pt[self.ellipse_colorby],
                            self.ellipse_colorby,
                            self.ellipse_cmap,
                            ckmin,
                            ckmax,
                        )
                    )

                axd.add_artist(ellipse)

            # -----------plot response phase tensors---------------
            if self.resp_fn is not None:
                rcmin = np.floor(self.pt_resid_arr["geometric_mean"].min())
                rcmax = np.floor(self.pt_resid_arr["geometric_mean"].max())
                for mpt, rpt in zip(
                    self.pt_resp_arr[data_ii], self.pt_resid_arr[data_ii]
                ):
                    eheight = (
                        mpt["phimin"]
                        / self.pt_resp_arr[data_ii]["phimax"].max()
                        * self.ellipse_size
                    )
                    ewidth = (
                        mpt["phimax"]
                        / self.pt_resp_arr[data_ii]["phimax"].max()
                        * self.ellipse_size
                    )

                    ellipsem = Ellipse(
                        (mpt["east"], mpt["north"]),
                        width=ewidth,
                        height=eheight,
                        angle=90 - mpt["azimuth"],
                    )

                    # get ellipse color
                    if self.ellipse_cmap.find("seg") > 0:
                        ellipsem.set_facecolor(
                            mtcl.get_plot_color(
                                mpt[self.ellipse_colorby],
                                self.ellipse_colorby,
                                self.ellipse_cmap,
                                ckmin,
                                ckmax,
                                bounds=bounds,
                            )
                        )
                    else:
                        ellipsem.set_facecolor(
                            mtcl.get_plot_color(
                                mpt[self.ellipse_colorby],
                                self.ellipse_colorby,
                                self.ellipse_cmap,
                                ckmin,
                                ckmax,
                            )
                        )

                    axm.add_artist(ellipsem)

                    # -----------plot residual phase tensors---------------
                    eheight = (
                        rpt["phimin"]
                        / self.pt_resid_arr[data_ii]["phimax"].max()
                        * self.ellipse_size
                    )
                    ewidth = (
                        rpt["phimax"]
                        / self.pt_resid_arr[data_ii]["phimax"].max()
                        * self.ellipse_size
                    )

                    ellipser = Ellipse(
                        (rpt["east"], rpt["north"]),
                        width=ewidth,
                        height=eheight,
                        angle=rpt["azimuth"],
                    )

                    # get ellipse color
                    rpt_color = np.sqrt(abs(rpt["phimin"] * rpt["phimax"]))
                    if self.ellipse_cmap.find("seg") > 0:
                        ellipser.set_facecolor(
                            mtcl.get_plot_color(
                                rpt_color,
                                "geometric_mean",
                                self.residual_cmap,
                                ckmin,
                                ckmax,
                                bounds=bounds,
                            )
                        )
                    else:
                        ellipser.set_facecolor(
                            mtcl.get_plot_color(
                                rpt_color,
                                "geometric_mean",
                                self.residual_cmap,
                                ckmin,
                                ckmax,
                            )
                        )

                    axr.add_artist(ellipser)

            # --> set axes properties
            # data
            axd.set_xlim(self.ew_limits)
            axd.set_ylim(self.ns_limits)
            axd.set_xlabel("Easting ({0})".format(self.map_scale), fontdict=font_dict)
            axd.set_ylabel("Northing ({0})".format(self.map_scale), fontdict=font_dict)
            # make a colorbar for phase tensors
            # bb = axd.axes.get_position().bounds
            bb = axd.get_position().bounds
            y1 = 0.25 * (
                2
                + (self.ns_limits[1] - self.ns_limits[0])
                / (self.ew_limits[1] - self.ew_limits[0])
            )
            cb_location = (
                3.35 * bb[2] / 5 + bb[0],
                y1 * self.cb_pt_pad,
                0.295 * bb[2],
                0.02,
            )
            cbaxd = fig.add_axes(cb_location)
            cbd = mcb.ColorbarBase(
                cbaxd,
                cmap=mtcl.cmapdict[self.ellipse_cmap],
                norm=Normalize(vmin=ckmin, vmax=ckmax),
                orientation="horizontal",
            )
            cbd.ax.xaxis.set_label_position("top")
            cbd.ax.xaxis.set_label_coords(0.5, 1.75)
            cbd.set_label(mtplottools.ckdict[self.ellipse_colorby])
            cbd.set_ticks(
                np.arange(ckmin, ckmax + self.cb_tick_step, self.cb_tick_step)
            )

            axd.text(
                self.ew_limits[0] * 0.95,
                self.ns_limits[1] * 0.95,
                "Data",
                horizontalalignment="left",
                verticalalignment="top",
                bbox={"facecolor": "white"},
                fontdict={"size": self.font_size + 1},
            )

            # Model and residual
            if self.resp_fn is not None:
                for aa, ax in enumerate([axm, axr]):
                    ax.set_xlim(self.ew_limits)
                    ax.set_ylim(self.ns_limits)
                    ax.set_xlabel(
                        "Easting ({0})".format(self.map_scale), fontdict=font_dict
                    )
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    # make a colorbar ontop of axis
                    bb = ax.axes.get_position().bounds
                    y1 = 0.25 * (
                        2
                        + (self.ns_limits[1] - self.ns_limits[0])
                        / (self.ew_limits[1] - self.ew_limits[0])
                    )
                    cb_location = (
                        3.35 * bb[2] / 5 + bb[0],
                        y1 * self.cb_pt_pad,
                        0.295 * bb[2],
                        0.02,
                    )
                    cbax = fig.add_axes(cb_location)
                    if aa == 0:
                        cb = mcb.ColorbarBase(
                            cbax,
                            cmap=mtcl.cmapdict[self.ellipse_cmap],
                            norm=Normalize(vmin=ckmin, vmax=ckmax),
                            orientation="horizontal",
                        )
                        cb.ax.xaxis.set_label_position("top")
                        cb.ax.xaxis.set_label_coords(0.5, 1.75)
                        cb.set_label(mtplottools.ckdict[self.ellipse_colorby])
                        cb.set_ticks(
                            np.arange(
                                ckmin, ckmax + self.cb_tick_step, self.cb_tick_step
                            )
                        )
                        ax.text(
                            self.ew_limits[0] * 0.95,
                            self.ns_limits[1] * 0.95,
                            "Model",
                            horizontalalignment="left",
                            verticalalignment="top",
                            bbox={"facecolor": "white"},
                            fontdict={"size": self.font_size + 1},
                        )
                    else:
                        cb = mcb.ColorbarBase(
                            cbax,
                            cmap=mtcl.cmapdict[self.residual_cmap],
                            norm=Normalize(vmin=rcmin, vmax=rcmax),
                            orientation="horizontal",
                        )
                        cb.ax.xaxis.set_label_position("top")
                        cb.ax.xaxis.set_label_coords(0.5, 1.75)
                        cb.set_label(r"$\sqrt{\Phi_{min} \Phi_{max}}$")
                        cb_ticks = [rcmin, (rcmax - rcmin) / 2, rcmax]
                        cb.set_ticks(cb_ticks)
                        ax.text(
                            self.ew_limits[0] * 0.95,
                            self.ns_limits[1] * 0.95,
                            "Residual",
                            horizontalalignment="left",
                            verticalalignment="top",
                            bbox={"facecolor": "white"},
                            fontdict={"size": self.font_size + 1},
                        )

            if self.model_fn is not None:
                for ax in ax_list:
                    ax.tick_params(direction="out")
                    bb = ax.axes.get_position().bounds
                    y1 = 0.25 * (
                        2
                        - (self.ns_limits[1] - self.ns_limits[0])
                        / (self.ew_limits[1] - self.ew_limits[0])
                    )
                    cb_position = (
                        3.0 * bb[2] / 5 + bb[0],
                        y1 * self.cb_res_pad,
                        0.35 * bb[2],
                        0.02,
                    )
                    cbax = fig.add_axes(cb_position)
                    cb = mcb.ColorbarBase(
                        cbax,
                        cmap=self.res_cmap,
                        norm=Normalize(
                            vmin=self.res_limits[0], vmax=self.res_limits[1]
                        ),
                        orientation="horizontal",
                    )
                    cb.ax.xaxis.set_label_position("top")
                    cb.ax.xaxis.set_label_coords(0.5, 1.5)
                    cb.set_label("Resistivity ($\Omega \cdot$m)")
                    cb_ticks = np.arange(
                        np.floor(self.res_limits[0]), np.ceil(self.res_limits[1] + 1), 1
                    )
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

    def save_figure(
        self,
        save_path=None,
        fig_dpi=None,
        file_format="pdf",
        orientation="landscape",
        close_fig="y",
    ):
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

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if not os.path.isdir(save_path):
            try:
                os.mkdir(save_path)
            except:
                raise IOError("Need to input a correct directory path")

        for fig in self.fig_list:
            per = fig.canvas.get_window_title()
            save_fn = os.path.join(
                save_path, "PT_DepthSlice_{0}s.{1}".format(per, file_format)
            )
            fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

            if close_fig == "y":
                plt.close(fig)

            else:
                pass

            self.fig_fn = save_fn
            print ("Saved figure to: " + self.fig_fn)


# ==================================================================================
# FZ: add example usage code
# Justdo>  python mtpy/modeling/modem/plot_pt_maps.py
# ==================================================================================
if __name__ == "__main__":
    from mtpy.mtpy_globals import *

    # directory where files are located
    wd = os.path.join(SAMPLE_DIR, "ModEM")

    # file stem for inversion result
    filestem = "Modular_MPI_NLCG_004"

    filestem = "Modular_MPI_NLCG_004"
    datafn = "ModEM_Data.dat"
    PlotPTMaps(
        data_fn=os.path.join(wd, datafn),
        resp_fn=os.path.join(wd, filestem + ".dat"),
        ellipse_size=20,
    )
