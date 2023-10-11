# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:10:55 2013

@author: jpeacock-pr
"""

# ==============================================================================

import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from mtpy.imaging.mtcolors import get_plot_color
from mtpy.imaging.mtplot_tools import (
    PlotBaseProfile,
    period_label_dict,
)


# ==============================================================================


class PlotPhaseTensorPseudoSection(PlotBaseProfile):
    """
    PlotPhaseTensorPseudoSection will plot the phase tensor ellipses in a
    pseudo section format


    To get a list of .edi files that you want to plot -->
    :Example: ::

        >>> import mtpy.imaging.mtplot as mtplot
        >>> import os
        >>> edipath = r"/home/EDIfiles"
        >>> edilist = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')>0]

    * If you want to plot minimum phase colored from blue to red in a range of
     20 to 70 degrees you can do it one of two ways-->

    1)
    :Example: ::

        >>> edict = {'range':(20,70), 'cmap':'mt_bl2gr2rd','colorby':'phimin'}
        >>> pt1 = mtplot.plot_pt_pseudosection(fn_list=edilist,
                                               ellipse_dict=edict)

    2)
    :Example: ::

        >>> pt1 = mtplot.plot_pt_pseudosection(fn_list=edilist, plot_yn='n')
        >>> pt1.ellipse_colorby = 'phimin'
        >>> pt1.ellipse_cmap = 'mt_bl2gr2rd'
        >>> pt1.ellipse_range = (20,70)
        >>> pt1.plot()

    * If you want to add real induction arrows that are scaled by 10 and point
     away from a conductor -->
    :Example: ::

        >>> pt1.plot_tipper = 'yr'
        >>> pt1.arrow_size = 10
        >>> pt1.arrow_direction = -1
        >>> pt1.redraw_plot()

    * If you want to save the plot as a pdf with a generic name -->
    :Example: ::
        >>> pt1.save_figure(r"/home/PTFigures", file_format='pdf', dpi=300)
        File saved to '/home/PTFigures/PTPseudoSection.pdf'

    """

    def __init__(self, mt_data, **kwargs):
        super().__init__(mt_data, **kwargs)

        self._rotation_angle = 0

        self.x_stretch = 1
        self.y_stretch = 10000
        self.y_scale = "period"
        # --> set plot properties

        self.station_id = [0, None]
        self.profile_vector = None
        self.profile_angle = None
        self.profile_line = None
        self.profile_reverse = False

        self.ellipse_size = 2000
        self.arrow_size = 4000
        self.arrow_head_width = 250
        self.arrow_head_length = 400

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.show_plot:
            self.plot()

    def _get_patch(self, tf):
        """
        Get ellipse patch

        :param tf: DESCRIPTION
        :type tf: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        plot_x = self._get_offset(tf)
        # --> set local variables

        pt_obj = tf.pt

        color_array = self.get_pt_color_array(pt_obj)
        for index, ff in enumerate(pt_obj.frequency):
            if self.y_scale == "period":
                plot_y = np.log10(1.0 / ff) * self.y_stretch
            else:
                plot_y = np.log10(ff) * self.y_stretch

            # --> get ellipse properties
            # if the ellipse size is not physically correct make it a dot
            if (
                pt_obj.phimax[index] == 0
                or pt_obj.phimax[index] > 100
                or pt_obj.phimin[index] == 0
                or pt_obj.phimin[index] > 100
            ):
                continue
            else:
                scaling = self.ellipse_size / pt_obj.phimax.max()
                eheight = pt_obj.phimin[index] * scaling
                ewidth = pt_obj.phimax[index] * scaling
                azimuth = 90 - pt_obj.azimuth[index]
                if self.y_scale == "period":
                    azimuth = 90 + pt_obj.azimuth[index]
            # make an ellipse
            ellipd = patches.Ellipse(
                (plot_x, plot_y),
                width=ewidth,
                height=eheight,
                angle=azimuth,
                lw=self.lw,
            )

            # get ellipse color
            ellipd.set_facecolor(
                get_plot_color(
                    color_array[index],
                    self.ellipse_colorby,
                    self.ellipse_cmap,
                    self.ellipse_range[0],
                    self.ellipse_range[1],
                    bounds=self.ellipse_cmap_bounds,
                )
            )

            self.ax.add_artist(ellipd)

            if tf.Tipper is not None:
                t_obj = tf.Tipper
                if "r" in self.plot_tipper:

                    if t_obj.mag_real[index] <= self.arrow_threshold:
                        if self.y_scale == "period":
                            txr = (
                                t_obj.mag_real[index]
                                * self.arrow_size
                                * np.sin(
                                    np.deg2rad(-t_obj.angle_real[index] + 180)
                                    + self.arrow_direction * np.pi
                                )
                            )
                            tyr = (
                                t_obj.mag_real[index]
                                * self.arrow_size
                                * np.cos(
                                    np.deg2rad(-t_obj.angle_real[index] + 180)
                                    + self.arrow_direction * np.pi
                                )
                            )
                        else:
                            txr = (
                                t_obj.mag_real[index]
                                * self.arrow_size
                                * np.sin(
                                    np.deg2rad(t_obj.angle_real[index])
                                    + self.arrow_direction * np.pi
                                )
                            )
                            tyr = (
                                t_obj.mag_real[index]
                                * self.arrow_size
                                * np.cos(
                                    np.deg2rad(t_obj.angle_real[index])
                                    + self.arrow_direction * np.pi
                                )
                            )

                        self.ax.arrow(
                            plot_x,
                            plot_y,
                            txr,
                            tyr,
                            width=self.arrow_lw,
                            facecolor=self.arrow_color_real,
                            edgecolor=self.arrow_color_real,
                            length_includes_head=False,
                            head_width=self.arrow_head_width,
                            head_length=self.arrow_head_length,
                        )
                    else:
                        pass
                # plot imaginary tipper
                if "i" in self.plot_tipper:
                    if t_obj.mag_imag[index] <= self.arrow_threshold:
                        if self.y_scale == "period":
                            txi = (
                                t_obj.mag_imag[index]
                                * self.arrow_size
                                * np.sin(
                                    np.deg2rad(-t_obj.angle_imag[index] + 180)
                                    + self.arrow_direction * np.pi
                                )
                            )
                            tyi = (
                                t_obj.mag_imag[index]
                                * self.arrow_size
                                * np.cos(
                                    np.deg2rad(-t_obj.angle_imag[index] + 180)
                                    + self.arrow_direction * np.pi
                                )
                            )
                        else:
                            txi = (
                                t_obj.mag_imag[index]
                                * self.arrow_size
                                * np.sin(
                                    np.deg2rad(t_obj.angle_imag[index])
                                    + self.arrow_direction * np.pi
                                )
                            )
                            tyi = (
                                t_obj.mag_imag[index]
                                * self.arrow_size
                                * np.cos(
                                    np.deg2rad(t_obj.angle_imag[index])
                                    + self.arrow_direction * np.pi
                                )
                            )

                        self.ax.arrow(
                            plot_x,
                            plot_y,
                            txi,
                            tyi,
                            width=self.arrow_lw,
                            facecolor=self.arrow_color_imag,
                            edgecolor=self.arrow_color_imag,
                            length_includes_head=False,
                            head_width=self.arrow_head_width,
                            head_length=self.arrow_head_length,
                        )

        return plot_x, tf.tf_id[self.station_id[0] : self.station_id[1]]

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
        cmap_input = plt.get_cmap(self.ellipse_cmap)

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

    def plot(self):
        """
        plots the phase tensor pseudo section.  See class doc string for
        more details.
        """

        self._set_subplot_params()

        # create a plot instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        self.fig.clf()
        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        self._get_profile_line()

        y_min = 1
        y_max = 1
        station_list = np.zeros(
            self.mt_data.n_stations,
            dtype=[("offset", np.float), ("station", "U10")],
        )

        for ii, tf in enumerate(self.mt_data.values()):
            offset, station = self._get_patch(tf)
            station_list[ii]["station"] = station
            station_list[ii]["offset"] = offset

            if np.log10(tf.Z.frequency.min()) < y_min:
                y_min = np.log10(tf.Z.frequency.min()) * self.y_stretch
            if np.log10(tf.Z.frequency.max()) > y_max:
                y_max = np.log10(tf.Z.frequency.max()) * self.y_stretch

        y_min = np.floor(y_min / self.y_stretch) * self.y_stretch
        y_max = np.ceil(y_max / self.y_stretch) * self.y_stretch

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
                    period_label_dict[int(tk / self.y_stretch)]
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
