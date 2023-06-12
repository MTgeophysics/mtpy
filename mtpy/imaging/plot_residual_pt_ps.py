# -*- coding: utf-8 -*-
"""
PlotResidualPhaseTensorPseudoSection
=======================

    *plots the residual phase tensor for two different sets of measurments


Created on Wed Oct 16 14:56:04 2013

@author: jpeacock-pr
"""
# ==============================================================================
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
from mtpy.imaging.mtplot_tools import PlotBaseProfile
from mtpy.analysis.residual_phase_tensor import ResidualPhaseTensor
import scipy.signal as sps

# ==============================================================================


class PlotResidualPTPseudoSection(PlotBaseProfile):
    """
    This will plot residual phase tensors in a pseudo section.  The data is
    read in and stored in 2 ways, one as a list ResidualPhaseTensor object for
    each matching station and the other in a structured array with all the
    important information.  The structured array is the one that is used for
    plotting.  It is computed each time plot() is called so if it is
    manipulated it is reset.  The array is sorted by relative offset, so no
    special order of input is needed for the file names.  However, the
    station names should be verbatim between surveys, otherwise it will not
    work.

    The residual phase tensor is calculated as I-(Phi_2)^-1 (Phi_1)

    The default coloring is by the geometric mean as sqrt(Phi_min*Phi_max),
    which defines the percent change between measurements.

    There are a lot of parameters to change how the plot looks, have a look
    below if you figure looks a little funny.  The most useful will be
    xstretch, ystretch and ellipse_size

    The ellipses are normalized by the largest Phi_max of the survey.


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

    * If you want to plot geometric mean colored from white to orange in a
      range of 0 to 10 percent you can do it one of two ways-->

    1)
    :Example: ::

        >>> edict = {'range':(0,10), 'cmap':'mt_wh2or', \
                     'colorby':'geometric_mean', 'size':10}
        >>> pt1 = mtplot.residual_pt_ps(edilist1, edilst2, ellipse_dict=edict)

    2)
    :Example: ::

        >>> pt1 = mtplot.residual_pt_ps(edilist1, edilst2, ellipse_dict=edict,\
                                        plot_yn='n')
        >>> pt1.ellipse_colorby = 'geometric_mean'
        >>> pt1.ellipse_cmap = 'mt_wh2or'
        >>> pt1.ellipse_range = (0, 10)
        >>> pt1.ellipse_size = 10
        >>> pt1.plot()


    * If you want to save the plot as a pdf with a generic name -->
    :Example: ::
        >>> pt1.save_figure(r"/home/PTFigures", file_format='pdf', dpi=300)
        File saved to '/home/PTFigures/PTPseudoSection.pdf'

    """

    def __init__(
        self,
        mt_data_01,
        mt_data_02,
        frequencies=np.logspace(-3, 3, 40),
        **kwargs,
    ):
        assert len(mt_data_01) == len(mt_data_02)

        super().__init__(None, **kwargs)

        self.freq_list = frequencies
        self.mt_data_01 = mt_data_01
        self.mt_data_02 = mt_data_02
        self.mt_data = mt_data_01

        self.ellipse_range = (0, 25)
        self.ellipse_cmap = "mt_yl2rd"
        self.ellipse_colorby = "geometric_mean"
        self.ellipse_scale = None

        self.residual_pt_list = None
        self.rpt_array = None
        self.med_filt_kernel = None
        self._filt_applied = False
        self.rot90 = True

        self._rotation_angle = 0
        # --> set station name properties
        self.station_y_pad = 0.0005

        # --> set station name properties
        self.station_id = (None, None)

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.mt_data_01 is not None and self.mt_data_02 is not None:
            self._compute_residual_pt()

        # --> plot if desired ------------------------
        if self.show_plot:
            self.plot()

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
            try:
                mt2 = two[key]
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
                ("station", "U20"),
                ("lat", float),
                ("lon", float),
                ("elev", float),
                ("offset", float),
                ("phimin", (float, num_freq)),
                ("phimax", (float, num_freq)),
                ("skew", (float, num_freq)),
                ("azimuth", (float, num_freq)),
                ("geometric_mean", (float, num_freq)),
            ],
        )

        self.residual_pt_list = []
        freq_dict = dict(
            [(np.round(ff, 5), ii) for ii, ff in enumerate(self.freq_list)]
        )
        for mm, match in enumerate(matches):
            mt1 = match[0]
            mt2 = match[1]

            new_mt1 = mt1.interpolate(self.freq_list, bounds_error=False)
            new_mt2 = mt2.interpolate(self.freq_list, bounds_error=False)

            # compute residual phase tensor
            rpt = ResidualPhaseTensor(new_mt1.pt, new_mt2.pt)

            # add some attributes to residual phase tensor object
            rpt.station = mt1.station
            rpt.lat = mt1.latitude
            rpt.lon = mt1.longitude

            # append to list for manipulating later
            self.residual_pt_list.append(rpt)

            # put stuff into an array because we cannot set values of
            # rpt, need this for filtering.
            self.rpt_array[mm]["station"] = mt1.station
            self.rpt_array[mm]["lat"] = mt1.latitude
            self.rpt_array[mm]["lon"] = mt1.longitude
            self.rpt_array[mm]["elev"] = mt1.elevation

            rpt_fdict = dict(
                [
                    (np.round(key, 5), value)
                    for value, key in enumerate(rpt.frequency)
                ]
            )
            for f_index, frequency in enumerate(rpt.frequency):
                aa = freq_dict[np.round(frequency, 5)]
                try:
                    try:
                        rr = rpt_fdict[np.round(frequency, 5)]

                        self.rpt_array[mm]["phimin"][aa] = abs(
                            rpt.residual_pt.phimin[rr]
                        )
                        self.rpt_array[mm]["phimax"][aa] = abs(
                            rpt.residual_pt.phimax[rr]
                        )
                        self.rpt_array[mm]["skew"][aa] = rpt.residual_pt.beta[
                            rr
                        ]
                        self.rpt_array[mm]["azimuth"][
                            aa
                        ] = rpt.residual_pt.azimuth[rr]
                        self.rpt_array[mm]["geometric_mean"][aa] = np.sqrt(
                            abs(
                                rpt.residual_pt.phimin[rr]
                                * rpt.residual_pt.phimax[rr]
                            )
                        )
                    except IndexError:
                        self.logger.info("-" * 50)
                        self.logger.info(mt1.station)
                        self.logger.info(f"freq_index for 1:  {f_index}")
                        self.logger.info(f"frequency looking for:  {frequency}")
                        self.logger.info(f"index in big    :  {aa}")
                        self.logger.info(f"index in 1      :  {rr}")
                        self.logger.info(
                            f"len_1 = {len(mt1.frequency)}, len_2 = {len(mt2.frequency)}"
                        )
                        self.logger.info(f"len rpt_freq = {len(rpt.frequency)}")
                except KeyError:
                    self.logger.info(
                        f"Station {mt1.station} does not have {frequency:.5f}Hz"
                    )

        # get profile line
        self._get_profile_line()

        # get offsets
        for rpt, mt_obj in zip(self.rpt_array, self.mt_data.values()):
            rpt["offset"] = self._get_offset(mt_obj)

        # from the data get the relative offsets and sort the data by them
        self.rpt_array.sort(order=["offset"])

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

    def get_pt_color_array(self, rpt_array):
        """
        Get the appropriat color by array
        """

        # get the properties to color the ellipses by
        if (
            self.ellipse_colorby == "phiminang"
            or self.ellipse_colorby == "phimin"
        ):
            color_array = rpt_array["phimin"]
        elif (
            self.ellipse_colorby == "phimaxang"
            or self.ellipse_colorby == "phimax"
        ):
            color_array = rpt_array["phimax"]

        elif (
            self.ellipse_colorby == "skew" or self.ellipse_colorby == "skew_seg"
        ):
            color_array = rpt_array["skew"]

        elif self.ellipse_colorby in ["strike", "azimuth"]:
            color_array = rpt_array["azimuth"] % 180
            color_array[np.where(color_array > 90)] -= 180

        elif self.ellipse_colorby in ["geometric_mean"]:
            color_array = rpt_array["geometric_mean"]
        else:
            raise NameError(self.ellipse_colorby + " is not supported")
        return color_array

    def _get_ellipse_max(self):
        if self.ellipse_scale is None:
            return np.nanmax(self.rpt_array["phimax"])
        else:
            return self.ellipse_scale

    def _get_patch(self, rpt):
        """
        Get ellipse patch

        :param tf: DESCRIPTION
        :type tf: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        # --> get size of largest ellipse for this frequency for
        #    normalization to give an indication of the size of
        #    change.
        emax = self._get_ellipse_max()

        for f_index, ff in enumerate(self.freq_list):
            if self.y_scale == "period":
                plot_y = np.log10(1.0 / ff) * self.y_stretch
            else:
                plot_y = np.log10(ff) * self.y_stretch

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
                (rpt["offset"] * self.x_stretch, plot_y),
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

    def _add_colorbar(self):
        """
        Add phase tensor color bar

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.cb_position is None:
            self.ax2, kw = mcb.make_axes(
                self.ax, orientation=self.cb_orientation, shrink=0.35
            )

        else:
            self.ax2 = self.fig.add_axes(self.cb_position)

        # make the colorbar
        if self.ellipse_cmap in list(mtcl.cmapdict.keys()):
            cmap_input = mtcl.cmapdict[self.ellipse_cmap]
        else:
            cmap_input = mtcl.cm.get_cmap(self.ellipse_cmap)

        if "seg" in self.ellipse_cmap:
            norms = colors.BoundaryNorm(self.ellipse_cmap_bounds, cmap_input.N)
            self.cb = mcb.ColorbarBase(
                self.ax2,
                cmap=cmap_input,
                norm=norms,
                orientation=self.cb_orientation,
                ticks=self.ellipse_cmap_bounds,
            )
        else:
            self.cb = mcb.ColorbarBase(
                self.ax2,
                cmap=cmap_input,
                norm=colors.Normalize(
                    vmin=self.ellipse_range[0], vmax=self.ellipse_range[1]
                ),
                orientation=self.cb_orientation,
            )

        # label the color bar accordingly
        self.cb.set_label(
            self.cb_label_dict[self.ellipse_colorby],
            fontdict={"size": self.font_size, "weight": "bold"},
        )

        # place the label in the correct location
        if self.cb_orientation == "horizontal":
            self.cb.ax.xaxis.set_label_position("top")
            self.cb.ax.xaxis.set_label_coords(0.5, 1.3)
        elif self.cb_orientation == "vertical":
            self.cb.ax.yaxis.set_label_position("right")
            self.cb.ax.yaxis.set_label_coords(1.25, 0.5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis="y", direction="in")

    # --------------------------------------------------------------------------
    def plot(self):
        """
        plot residual phase tensor
        """

        self._set_subplot_params()

        if self.med_filt_kernel is not None:
            if self._filt_applied:
                self._compute_residual_pt()
                self._apply_median_filter(self.med_filt_kernel)
            else:
                self._apply_median_filter(self.med_filt_kernel)

            self._filt_applied = True

        # create a plot instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        self.fig.clf()
        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        y_min = 1
        y_max = 1
        station_list = np.zeros(
            self.rpt_array.size,
            dtype=[("offset", np.float), ("station", "U10")],
        )

        for ii, rpt in enumerate(self.rpt_array):
            self._get_patch(rpt)
            station_list[ii]["station"] = rpt["station"][
                self.station_id[0] : self.station_id[1]
            ]
            station_list[ii]["offset"] = rpt["offset"] * self.x_stretch

        y_min = np.floor(np.log10(self.freq_list.min())) * self.y_stretch
        y_max = np.ceil(np.log10(self.freq_list.max())) * self.y_stretch

        # --> Set plot parameters
        self.station_list = np.sort(station_list, order="offset")

        # set y-ticklabels
        if self.y_scale == "period":
            y_label = "Period (s)"
            p_min = float(y_min)
            p_max = float(y_max)
            y_min = -1 * p_min
            y_max = -1 * p_max

        else:
            y_label = "Frequency (Hz)"

        self.ax.set_ylabel(y_label, fontdict=self.font_dict)

        # set y-axis tick labels
        self.ax.yaxis.set_ticks(
            np.arange(y_min, (y_max), self.y_stretch * np.sign(y_max))
        )

        y_tick_labels = []

        for tk in self.ax.get_yticks():
            try:
                y_tick_labels.append(
                    self.period_label_dict[int(tk / self.y_stretch)]
                )
            except KeyError:
                y_tick_labels.append("")

        self.ax.set_yticklabels(y_tick_labels)

        # --> set tick locations and labels

        # set x-axis ticks
        self.ax.set_xticks(self.station_list["offset"])

        # set x-axis tick labels as station names
        self.ax.set_xticklabels(self.station_list["station"])

        # set x-axis label
        self.ax.set_xlabel(
            "Station", fontsize=self.font_size + 2, fontweight="bold"
        )

        # --> set x-limits
        if self.x_limits is None:
            self.ax.set_xlim(
                np.floor(self.station_list["offset"].min())
                - self.ellipse_size / 2,
                np.ceil(self.station_list["offset"].max())
                + self.ellipse_size / 2,
            )
        else:
            self.ax.set_xlim(self.x_limits)
        # --> set y-limits
        if self.y_limits is None:
            self.ax.set_ylim(y_min, y_max)
        else:
            pmin = np.log10(self.y_limits[0]) * self.y_stretch
            pmax = np.log10(self.y_limits[1]) * self.y_stretch
            self.ax.set_ylim(np.floor(pmin), np.ceil(pmax))

        # --> set title of the plot
        if self.plot_title is None:
            pass
        else:
            self.ax.set_title(self.plot_title, fontsize=self.font_size + 2)

        # put a grid on the plot
        self.ax.grid(alpha=0.25, which="both", color=(0.25, 0.25, 0.25))

        # ==> make a colorbar with appropriate colors
        self._add_colorbar()

        self.ax.set_axisbelow(True)

        plt.show()
