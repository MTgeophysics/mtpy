# -*- coding: utf-8 -*-
"""
PlotResidualPhaseTensor
=======================

    *plots the residual phase tensor for two different sets of measurments
    
    
Created on Wed Oct 16 14:56:04 2013

@author: jpeacock-pr
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import scipy.signal as sps

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as patches
import matplotlib.colorbar as mcb

import mtpy.utils.gis_tools as gis_tools
import mtpy.imaging.mtcolors as mtcl
from mtpy.imaging.mtplot_tools import PlotBase
from mtpy.analysis.residual_phase_tensor import ResidualPhaseTensor

try:
    import contextily as cx

    has_cx = True
except ModuleNotFoundError:
    has_cx = False

# ==============================================================================


class PlotResidualPTMaps(PlotBase):
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

    def __init__(
        self,
        mt_data_01,
        mt_data_02,
        frequencies=np.logspace(-3, 3, 40),
        **kwargs,
    ):

        super().__init__(**kwargs)

        self.map_epsg = 4326

        self.image_file = None
        self.image_extent = None

        self.freq_list = frequencies
        self.mt_data_01 = mt_data_01
        self.mt_data_02 = mt_data_02

        self.ellipse_range = (0, 25)
        self.ellipse_cmap = "mt_yl2rd"
        self.ellipse_colorby = "geometric_mean"
        self.ellipse_scale = None

        self.cx_source = None
        self.cx_zoom = None
        if has_cx:
            self.cx_source = cx.providers.USGS.USTopo

        self.residual_pt_list = None
        self.rpt_array = None
        self.med_filt_kernel = None
        self._filt_applied = False

        self._rotation_angle = 0

        # --> set plot properties ------------------------------
        # set some of the properties as attributes much to Lars' discontent
        self.plot_station_name = False

        self.tscale = "period"
        self.map_scale = "deg"
        self.map_epsg = 4326
        self.map_utm_zone = None
        self.rot90 = True

        # --> set the frequency to plot
        self.plot_freq = kwargs.pop("plot_freq", 1.0)
        self.ftol = kwargs.pop("ftol", 0.05)
        self.plot_freq_index = None

        # --> set background image properties
        if self.image_file is not None:
            if self.image_extent is None:
                raise ValueError(
                    "Need to input extents of the image as (x0, y0, x1, y1)"
                )
        # --> set a central reference point
        self.plot_reference_point = (0, 0)

        # --> set station name properties
        self.station_id = (None, None)

        for key, value in kwargs.items():
            setattr(self, key, value)

        # --> plot if desired ------------------------
        if self.show_plot:
            self.plot()

    @property
    def map_scale(self):
        return self._map_scale

    @map_scale.setter
    def map_scale(self, map_scale):
        self._map_scale = map_scale

        if self._map_scale == "deg":
            self.xpad = 0.005
            self.ypad = 0.005
            self.ellipse_size = 0.005
            self.tickstrfmt = "%.3f"
            self.y_label = "Latitude (deg)"
            self.x_label = "Longitude (deg)"

        elif self._map_scale == "m":
            self.xpad = 1000
            self.ypad = 1000
            self.ellipse_size = 500
            self.tickstrfmt = "%.0f"
            self.x_label = "Easting (m)"
            self.y_label = "Northing (m)"

        elif self._map_scale == "km":
            self.xpad = 1
            self.ypad = 1
            self.ellipse_size = 0.500
            self.tickstrfmt = "%.0f"
            self.x_label = "Easting (km)"
            self.y_label = "Northing (km)"

        else:
            raise ValueError(f"map scale {map_scale} is not supported.")

    @property
    def rotation_angle(self):
        return self._rotation_angle

    # ---need to rotate data on setting rotz
    @rotation_angle.setter
    def rotation_angle(self, rot_z):
        """
        need to rotate data when setting z
        """

        for tf_obj in self.mt_data_01 + self.mt_data_02:
            tf_obj.rotation_angle = rot_z

    def _match_lists(self, one, two):
        """
        Match the input lists by transfer function id

        :return: DESCRIPTION
        :rtype: TYPE

        """

        matches = []
        two_found = []
        for key, mt1 in one.items():
            station_find = False
            station_id = key.split(".")[1]
            try:
                mt2 = two.get_station(station_id=station_id)
                matches.append([mt1, mt2])
                station_find = True
                two_found.append(mt2.tf_id)

            except KeyError:
                for mt2 in two.values():
                    if mt2.tf_id in two_found:
                        continue
                    if mt1.tf_id == mt2.tf_id:
                        if (
                            abs(mt1.latitude - mt2.latitude) < 0.001
                            and abs(mt1.longitude - mt2.longitude) < 0.001
                        ):
                            matches.append([mt1, mt2])
                            station_find = True
                            two_found.append(mt2.tf_id)
                        break
            if not station_find:
                self.logger.warning(
                    f"Could not find tf {mt1.tf_id} in second list"
                )

        return matches

    # ------------------------------------------------------------------
    def _compute_residual_pt(self):
        """
        compute residual phase tensor so the result is something useful to
        plot
        """

        matches = self._match_lists(self.mt_data_01, self.mt_data_02)
        num_freq = self.freq_list.shape[0]
        num_station = len(matches)

        # make a structured array to put stuff into for easier manipulation
        self.rpt_array = np.zeros(
            num_station,
            dtype=[
                ("station", "|U20"),
                ("lat", float),
                ("lon", float),
                ("plotx", float),
                ("ploty", float),
                ("elev", float),
                ("phimin", (float, num_freq)),
                ("phimax", (float, num_freq)),
                ("skew", (float, num_freq)),
                ("azimuth", (float, num_freq)),
                ("geometric_mean", (float, num_freq)),
            ],
        )

        self.residual_pt_list = []

        for mm, match in enumerate(matches):
            mt1 = match[0]
            mt2 = match[1]

            new_1 = mt1.interpolate(self.freq_list, bounds_error=False)
            new_2 = mt2.interpolate(self.freq_list, bounds_error=False)

            # compute residual phase tensor
            with np.errstate(divide="ignore", invalid="ignore"):
                rpt = ResidualPhaseTensor(new_1.pt, new_2.pt)

            # add some attributes to residual phase tensor object
            rpt.station = mt1.station
            rpt.lat = mt1.latitude
            rpt.lon = mt1.longitude

            # append to list for manipulating later
            self.residual_pt_list.append(rpt)

            # put stuff into an array because we cannot set values of
            # rpt, need this for filtering.
            st_1, st_2 = self.station_id
            self.rpt_array[mm]["station"] = mt1.station[st_1:st_2]
            self.rpt_array[mm]["lat"] = mt1.latitude
            self.rpt_array[mm]["lon"] = mt1.longitude
            self.rpt_array[mm]["elev"] = mt1.elevation

            self.rpt_array[mm]["phimin"][:] = np.abs(rpt.residual_pt.phimin)
            self.rpt_array[mm]["phimax"][:] = np.abs(rpt.residual_pt.phimax)
            self.rpt_array[mm]["skew"][:] = rpt.residual_pt.beta
            self.rpt_array[mm]["azimuth"][:] = rpt.residual_pt.azimuth
            self.rpt_array[mm]["geometric_mean"][:] = np.sqrt(
                np.abs(rpt.residual_pt.phimin * rpt.residual_pt.phimax)
            )

        # from the data get the relative offsets and sort the data by them
        self.rpt_array.sort(order=["lon", "lat"])

        # get relative positions for plotting
        self._get_relative_position()

        if self.med_filt_kernel is not None:
            try:
                self._apply_median_filter(kernel=self.med_filt_kernel)
            except:
                self.logger.warning("Could not apply median filter")

    # -------------------------------------------------------------------
    def _apply_median_filter(self, kernel=(3, 3)):
        """
        apply a median filter to the data to remove extreme outliers

        kernel is (station, frequency)

        """

        filt_phimin_arr = sps.medfilt2d(
            self.rpt_array["phimin"], kernel_size=kernel
        )
        filt_phimax_arr = sps.medfilt2d(
            self.rpt_array["phimax"], kernel_size=kernel
        )
        filt_skew_arr = sps.medfilt2d(
            self.rpt_array["skew"], kernel_size=kernel
        )
        filt_azimuth_arr = sps.medfilt2d(
            self.rpt_array["azimuth"], kernel_size=kernel
        )

        self.rpt_array["phimin"] = filt_phimin_arr
        self.rpt_array["phimax"] = filt_phimax_arr
        self.rpt_array["skew"] = filt_skew_arr
        self.rpt_array["azimuth"] = filt_azimuth_arr
        self.rpt_array["geometric_mean"] = np.sqrt(
            abs(filt_phimin_arr * filt_phimax_arr)
        )

        self.logger.info(f"Applying Median Filter with kernel {kernel}")

    # -------------------------------------------------------------------------
    def _get_relative_position(self):
        """
        get the relative positions for each station in the plotting
        coordinates
        """
        # if map scale is lat lon set parameters
        for ii, rpt in enumerate(self.rpt_array):
            if self.map_scale == "deg":
                plotx = rpt["lon"] - self.plot_reference_point[0]
                ploty = rpt["lat"] - self.plot_reference_point[1]
            # if map scale is in meters easting and northing
            elif self.map_scale in ["m", "km"]:
                if self.map_utm_zone is None:
                    raise ValueError(
                        "map_utm_zone must be specified for units of m or km"
                    )
                east, north, zone = gis_tools.project_point_ll2utm(
                    rpt["lat"],
                    rpt["lon"],
                    utm_zone=self.map_utm_zone,
                    epsg=self.map_epsg,
                )

                plotx = east - self.plot_reference_point[0]
                ploty = north - self.plot_reference_point[1]

                if self.map_scale in ["km"]:
                    plotx /= 1000
                    ploty /= 1000

            else:
                raise ValueError(f"map scale {self.map_scale} not supported")
            # put the location of each ellipse into an array in x and y
            rpt["plotx"] = plotx
            rpt["ploty"] = ploty

    def _get_plot_freq_index(self):
        """
        get frequency to plot
        """

        array = np.asarray(self.freq_list)
        self.plot_freq_index = np.abs(array - self.plot_freq).argmin()
        print(f"--> Plotting {self.freq_list[self.plot_freq_index]:.6g} Hz")

    def _get_title(self):
        if self.tscale == "period":
            return f"{(1.0 / self.plot_freq):.5g} (s)"
        else:
            return f"{self.plot_freq:.5g} (Hz)"

    def _get_ellipse_max(self):
        if self.ellipse_scale is None:
            return np.nanmax(self.rpt_array["phimax"])
        else:
            return self.ellipse_scale

    def _make_reference_ellipse(self):
        ref_ellip = patches.Ellipse(
            (0, 0.0),
            width=self.ellipse_size,
            height=self.ellipse_size,
            angle=0,
        )
        ref_ellip.set_facecolor((0, 0, 0))
        ref_ax_loc = list(self.ax2.get_position().bounds)
        ref_ax_loc[0] *= 0.95
        ref_ax_loc[1] -= 0.17
        ref_ax_loc[2] = 0.1
        ref_ax_loc[3] = 0.1
        self.ref_ax = self.fig.add_axes(ref_ax_loc, aspect="equal")
        self.ref_ax.add_artist(ref_ellip)
        self.ref_ax.set_xlim(
            -self.ellipse_size / 2.0 * 1.05, self.ellipse_size / 2.0 * 1.05
        )
        self.ref_ax.set_ylim(
            -self.ellipse_size / 2.0 * 1.05, self.ellipse_size / 2.0 * 1.05
        )
        plt.setp(self.ref_ax.xaxis.get_ticklabels(), visible=False)
        plt.setp(self.ref_ax.yaxis.get_ticklabels(), visible=False)
        self.ref_ax.set_title(r"$\Delta \Phi$ = 1")

    # ------------------------------------------------------------------------
    def plot(self):
        """
        plot residual phase tensor
        """
        self._set_subplot_params()

        # get residual phase tensor for plotting
        if self.rpt_array is None:
            self._compute_residual_pt()
        # get frequency index
        self._get_plot_freq_index()

        # make figure instance

        self.fig = plt.figure(
            self._get_title(), self.fig_size, dpi=self.fig_dpi
        )

        # clear the figure if there is already one up
        plt.clf()

        # make an axes instance
        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        # --> plot the background image if desired-----------------------
        try:
            im = plt.imread(self.image_file)
            self.ax.imshow(
                im, origin="lower", extent=self.image_extent, aspect="auto"
            )
        except AttributeError:
            pass

        f_index = self.plot_freq_index

        # --> get size of largest ellipse for this frequency for
        #    normalization to give an indication of the size of
        #    change.
        emax = self._get_ellipse_max()

        # --> plot
        for ii, rpt in enumerate(self.rpt_array):
            # --> get ellipse properties
            # if the ellipse size is not physically correct make it a dot
            if rpt["phimax"][f_index] == 0 and rpt["phimin"][f_index] == 0:
                continue
            elif rpt["phimax"][f_index] > 100 or rpt["phimin"][f_index] > 100:
                continue
                self.logger.warning(
                    "Bad data at {rpt['station']} for frequency {self.plot_freq}"
                )
            else:
                scaling = self.ellipse_size / emax
                e_height = rpt["phimin"][f_index] * scaling
                e_width = rpt["phimax"][f_index] * scaling
            # make an ellipse
            if self.rot90:
                e_angle = rpt["azimuth"][f_index] - 90
            else:
                e_angle = rpt["azimuth"][f_index]

            ellipd = patches.Ellipse(
                (rpt["plotx"], rpt["ploty"]),
                width=e_width,
                height=e_height,
                angle=e_angle,
            )
            # get ellipse color
            if self.ellipse_cmap.find("seg") > 0:
                ellipd.set_facecolor(
                    mtcl.get_plot_color(
                        rpt[self.ellipse_colorby][f_index],
                        self.ellipse_colorby,
                        self.ellipse_cmap,
                        self.ellipse_range[0],
                        self.ellipse_range[1],
                        bounds=self.ellipse_cmap_bounds,
                    )
                )
            else:
                ellipd.set_facecolor(
                    mtcl.get_plot_color(
                        rpt[self.ellipse_colorby][f_index],
                        self.ellipse_colorby,
                        self.ellipse_cmap,
                        self.ellipse_range[0],
                        self.ellipse_range[1],
                    )
                )
            # ==> add ellipse to the plot
            self.ax.add_artist(ellipd)

            # ------------Plot station name------------------------------
            if self.plot_station_name == True:
                self.ax.annotate(
                    rpt["station"],
                    xy=(rpt["plotx"], rpt["ploty"]),
                    ha="center",
                    va="baseline",
                    xytext=(rpt["plotx"], rpt["ploty"] + self.text_y_pad),
                    rotation=self.text_angle,
                    color=self.text_color,
                    fontsize=self.text_size,
                    fontweight=self.text_weight,
                )

        if self.image_file is None:
            if has_cx:
                try:
                    cx_kwargs = {
                        "crs": f"epsg:{self.map_epsg}",
                        "source": self.cx_source,
                    }
                    if self.cx_zoom is not None:
                        cx_kwargs["zoom"] = self.cx_zoom
                    cx.add_basemap(
                        self.ax,
                        **cx_kwargs,
                    )
                except Exception as error:
                    self.logger.warning(
                        f"Could not add base map because {error}"
                    )

        # --> set axes properties depending on map scale------------------------
        self.ax.set_xlabel(self.x_label, fontdict=self.font_dict)
        self.ax.set_ylabel(self.y_label, fontdict=self.font_dict)

        # --> set plot limits
        self.ax.set_xlim(
            self.rpt_array["plotx"].min() - self.xpad,
            self.rpt_array["plotx"].max() + self.xpad,
        )
        self.ax.set_ylim(
            self.rpt_array["ploty"].min() - self.xpad,
            self.rpt_array["ploty"].max() + self.xpad,
        )

        # --> set tick label format
        self.ax.xaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))

        # --> set title in period or frequency
        if not self.plot_title:
            self.ax.set_title(
                f"Phase Tensor Map for {self._get_title()}",
                fontsize=self.font_size + 2,
                fontweight="bold",
            )
        else:
            self.ax.set_title(
                self.plot_title + self._get_title(),
                fontsize=self.font_size + 2,
                fontweight="bold",
            )
        # make a grid with gray lines
        self.ax.grid(alpha=0.25)

        # ==> make a colorbar with appropriate colors
        if self.cb_position == None:
            self.ax2, kw = mcb.make_axes(
                self.ax, orientation=self.cb_orientation, shrink=0.35
            )
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)
        self.cb = self.make_pt_cb(self.ax2)

        # --> add reference ellipse
        self._make_reference_ellipse()

        # put the grid lines behind
        self.ax.set_axisbelow(True)

        plt.show()
