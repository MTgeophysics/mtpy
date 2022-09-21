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
from scipy import stats

import mtpy.imaging.mtcolors as mtcl
from mtpy.imaging.mtplot_tools import PlotBase


# ==============================================================================


class PlotPhaseTensorPseudoSection(PlotBase):
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

    def __init__(self, tf_list, **kwargs):
        super().__init__(**kwargs)

        self._rotation_angle = 0
        self.tf_list = tf_list

        self.x_stretch = 200
        self.y_stretch = 10
        # --> set plot properties

        self.linedir = "ew"
        self.station_id = [0, 4]
        self.y_step = 4
        self.x_step = 1
        self.profile_vector = None
        self.profile_angle = None
        self.profile_line = None

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.show_plot:
            self.plot()

    # ---need to rotate data on setting rotz
    @property
    def rotation_angle(self):
        return self._rotation_angle

    @rotation_angle.setter
    def rotation_angle(self, value):
        """
        only a single value is allowed
        """
        for tf in self.tf_list:
            tf.rotation_angle = value
        self._rotation_angle = value

    def _get_profile_line(self):
        east = np.zeros(len(self.tf_list))
        north = np.zeros(len(self.tf_list))

        for ii, tf in enumerate(self.tf_list):
            tf.project_to_utm()
            east.append(tf.east)
            north.append(tf.north)

        # check regression for 2 profile orientations:
        # horizontal (N=N(E)) or vertical(E=E(N))
        # use the one with the lower standard deviation
        profile1 = stats.linregress(east, north)
        profile2 = stats.linregress(north, east)
        # if the profile is rather E=E(N), the parameters have to converted
        # into N=N(E) form:
        if profile2.stderr < profile1.stderr:
            self.profile_line = (
                1.0 / profile2.slope,
                -profile2.intercept / profile2.slope,
            )
        else:
            self.profile_line = profile1[:2]

        self.profile_angle = (
            90 - np.rad2deg(np.arctan(self.profile_line.slope))
        ) % 180

        self.profile_vector = np.array([1, self.profile_line.slope])
        self.profile_vector /= np.linalg.norm(self.profile_vector)

    def _get_patch(self, tf):
        """
        Get ellipse patch

        :param tf: DESCRIPTION
        :type tf: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        pt_obj = tf.pt
        t_obj = tf.Tipper

        station_vector = np.array([tf.east, tf.north - self.profile_line[1]])
        plot_x = (
            np.linalg.norm(
                np.dot(self.profile_vector, station_vector)
                * self.profile_vector
            )
            * self.x_stretch
        )
        # --> set local variables

        color_array = self.get_pt_color_array(pt_obj)
        for index, ff in enumerate(pt_obj.freq):
            plot_y = np.log10(ff) * self.y_stretch
            # --> get ellipse properties
            # if the ellipse size is not physically correct make it a dot
            if (
                pt_obj.phimax[index] == 0
                or pt_obj.phimax[index] > 100
                or pt_obj.phimin[index] == 0
                or pt_obj.phimin[index] > 100
            ):
                eheight = 0.0000001
                ewidth = 0.0000001
            else:
                scaling = self.ellipse_size / self.pt_obj.phimax.max()
                eheight = self.pt_obj.phimin[index] * scaling
                ewidth = self.pt_obj.phimax[index] * scaling
            # make an ellipse
            ellipd = patches.Ellipse(
                (plot_x, plot_y),
                width=ewidth,
                height=eheight,
                angle=90 - self.pt_obj.azimuth[index],
                lw=self.lw,
            )

            # get ellipse color
            ellipd.set_facecolor(
                mtcl.get_plot_color(
                    color_array[0],
                    self.ellipse_colorby,
                    self.ellipse_cmap,
                    self.ellipse_range[0],
                    self.ellipse_range[1],
                    bounds=self.ellipse_cmap_bounds,
                )
            )

            self.ax.add_artist(ellipd)

            if t_obj is not None:
                if "r" in self.plot_tipper == "yri":

                    if t_obj.mag_real[0] <= self.arrow_threshold:
                        txr = (
                            t_obj.mag_real[0]
                            * self.arrow_size
                            * np.sin(
                                (t_obj.angle_real[0]) * np.pi / 180
                                + self.arrow_direction * np.pi
                            )
                        )
                        tyr = (
                            t_obj.mag_real[0]
                            * self.arrow_size
                            * np.cos(
                                (t_obj.angle_real[0]) * np.pi / 180
                                + self.arrow_direction * np.pi
                            )
                        )

                        self.ax.arrow(
                            plot_x * self.x_stretch,
                            plot_y * self.y_stretch,
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
                    if t_obj.mag_imag[0] <= self.arrow_threshold:
                        txi = (
                            t_obj.mag_imag[0]
                            * self.arrow_size
                            * np.sin(
                                (t_obj.angle_imag[0]) * np.pi / 180
                                + self.arrow_direction * np.pi
                            )
                        )
                        tyi = (
                            t_obj.mag_imag[0]
                            * self.arrow_size
                            * np.cos(
                                (t_obj.angle_imag[0]) * np.pi / 180
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

        for tf in self.tf_list:
            station, offset = self._get_path(tf)

        # --> Set plot parameters
        self._plot_period_list = plot_period_list
        n = len(plot_period_list)

        pmin = int(np.floor(np.log10(plot_period_list.min())))
        pmax = int(np.ceil(np.log10(plot_period_list.max())))

        # need to sort the offsets and station labels so they plot correctly
        sdtype = [("offset", np.float), ("station", "U10")]
        slist = np.array(
            [(oo, ss) for oo, ss in zip(self.offset_list, self.station_list)],
            dtype=sdtype,
        )
        offset_sort = np.sort(slist, order="offset")

        self.offset_list = offset_sort["offset"]
        self.station_list = offset_sort["station"]

        # set y-ticklabels
        if self.tscale == "period":
            yticklabels = [
                mtpl.labeldict[ii] for ii in range(pmin, pmax + 1, 1)
            ]

            self.ax.set_ylabel(
                "Period (s)", fontsize=self.font_size + 2, fontweight="bold"
            )
        elif self.tscale == "frequency":
            yticklabels = [
                mtpl.labeldict[-ii] for ii in range(pmin, pmax + 1, 1)
            ]
            self.ax.set_ylabel(
                "Frequency (Hz)",
                fontsize=self.font_size + 2,
                fontweight="bold",
            )
        # set x-axis label
        self.ax.set_xlabel(
            "Station", fontsize=self.font_size + 2, fontweight="bold"
        )

        # --> set tick locations and labels
        # set y-axis major ticks
        self.ax.yaxis.set_ticks(
            np.arange(
                pmin * self.ystretch, (pmax + 1) * self.ystretch, self.ystretch
            )
        )

        # set y-axis tick labels
        self.ax.set_yticklabels(yticklabels)

        # set x-axis ticks
        self.ax.set_xticks(self.offset_list * self.xstretch)

        # set x-axis tick labels as station names
        xticklabels = self.station_list
        if self.xstep != 1:
            xticklabels = np.zeros(
                len(self.station_list), dtype=self.station_list.dtype
            )
            for xx in range(0, len(self.station_list), self.xstep):
                xticklabels[xx] = self.station_list[xx]
        self.ax.set_xticklabels(xticklabels)

        # --> set x-limits
        if self.xlimits is None:
            self.ax.set_xlim(
                self.offset_list.min() * self.xstretch - es * 2,
                self.offset_list.max() * self.xstretch + es * 2,
            )
        else:
            self.ax.set_xlim(self.xlimits)
        # --> set y-limits
        if self.ylimits is None:
            #            self.ax.set_ylim(pmax+es*2, pmin-es*2)
            self.ax.set_ylim(pmax * self.ystretch, pmin * self.ystretch)
        else:
            pmin = np.log10(self.ylimits[0]) * self.ystretch
            pmax = np.log10(self.ylimits[1]) * self.ystretch
            self.ax.set_ylim(pmax + es * 2, pmin - es * 2)
        #            self.ax.set_ylim(pmax, pmin)

        # --> set title of the plot
        if self.plot_title is None:
            pass
        else:
            self.ax.set_title(self.plot_title, fontsize=self.font_size + 2)
        # make a legend for the induction arrows
        if self.plot_tipper.find("y") == 0:
            if self.plot_tipper == "yri":
                treal = self.ax.plot(
                    np.arange(10) * 0.000005,
                    np.arange(10) * 0.00005,
                    color=self.arrow_color_real,
                )
                timag = self.ax.plot(
                    np.arange(10) * 0.000005,
                    np.arange(10) * 0.00005,
                    color=self.arrow_color_imag,
                )
                self.ax.legend(
                    [treal[0], timag[0]],
                    ["Tipper_real", "Tipper_imag"],
                    loc="lower right",
                    prop={"size": self.font_size - 1, "weight": "bold"},
                    ncol=2,
                    markerscale=0.5,
                    borderaxespad=0.005,
                    borderpad=0.25,
                )
            elif self.plot_tipper == "yr":
                treal = self.ax.plot(
                    np.arange(10) * 0.000005,
                    np.arange(10) * 0.00005,
                    color=self.arrow_color_real,
                )
                self.ax.legend(
                    [treal[0]],
                    ["Tipper_real"],
                    loc="lower right",
                    prop={"size": self.font_size - 1, "weight": "bold"},
                    ncol=2,
                    markerscale=0.5,
                    borderaxespad=0.005,
                    borderpad=0.25,
                )
            elif self.plot_tipper == "yi":
                timag = self.ax.plot(
                    np.arange(10) * 0.000005,
                    np.arange(10) * 0.00005,
                    color=self.arrow_color_imag,
                )
                self.ax.legend(
                    [timag[0]],
                    ["Tipper_imag"],
                    loc="lower right",
                    prop={"size": self.font_size - 1, "weight": "bold"},
                    ncol=2,
                    markerscale=0.5,
                    borderaxespad=0.005,
                    borderpad=0.25,
                )
            # make a scale arrow
            if self.scale_arrow:
                print(
                    (
                        np.log10(
                            self.ylimits[1]
                            - self.scale_arrow_dict["text_offset_y"]
                        )
                    )
                    * self.ystretch
                )
                txrl = self.scale_arrow_dict["size"]
                self.ax.arrow(
                    min(self.offset_list) * self.xstretch,
                    np.log10(self.ylimits[1]) * self.ystretch,
                    txrl * self.arrow_size,
                    0.0,
                    lw=alw,
                    facecolor=self.arrow_color_real,
                    edgecolor=self.arrow_color_real,
                    length_includes_head=False,
                    head_width=awidth,
                    head_length=aheight,
                )
                self.ax.text(
                    min(self.offset_list) * self.xstretch,
                    (
                        np.log10(
                            self.ylimits[1]
                            - self.scale_arrow_dict["text_offset_y"]
                        )
                    )
                    * self.ystretch,
                    "|T| = %3.1f" % txrl,
                )
        # put a grid on the plot
        self.ax.grid(alpha=0.25, which="both", color=(0.25, 0.25, 0.25))

        # print out the min an max of the parameter plotted
        print("-" * 25)
        print(ck + " min = {0:.2f}".format(min(minlist)))
        print(ck + " max = {0:.2f}".format(max(maxlist)))
        print("-" * 25)

        # ==> make a colorbar with appropriate colors
        if self.cb_position is None:
            self.ax2, kw = mcb.make_axes(
                self.ax, orientation=self.cb_orientation, shrink=0.35
            )
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)
        if cmap == "mt_seg_bl2wh2rd":
            # make a color list
            self.clist = [
                (cc, cc, 1)
                for cc in np.arange(0, 1 + 1.0 / (nseg), 1.0 / (nseg))
            ] + [
                (1, cc, cc)
                for cc in np.arange(1, -1.0 / (nseg), -1.0 / (nseg))
            ]

            # make segmented colormap
            mt_seg_bl2wh2rd = colors.ListedColormap(self.clist)

            # make bounds so that the middle is white
            bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

            # normalize the colors
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

            # make the colorbar
            self.cb = mcb.ColorbarBase(
                self.ax2,
                cmap=mt_seg_bl2wh2rd,
                norm=norms,
                orientation=self.cb_orientation,
                ticks=bounds[1:-1],
            )
        else:
            self.cb = mcb.ColorbarBase(
                self.ax2,
                cmap=mtcl.cmapdict[cmap],
                norm=colors.Normalize(vmin=ckmin, vmax=ckmax),
                orientation=self.cb_orientation,
            )
        # label the color bar accordingly
        self.cb.set_label(
            mtpl.ckdict[ck],
            fontdict={"size": self.font_size, "weight": "bold"},
        )

        # place the label in the correct location
        if self.cb_orientation == "horizontal":
            self.cb.ax.xaxis.set_label_position("top")
            self.cb.ax.xaxis.set_label_coords(0.5, 1.3)
        elif self.cb_orientation == "vertical":
            self.cb.ax.yaxis.set_label_position("right")
            self.cb.ax.yaxis.set_label_coords(1.5, 0.5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis="y", direction="in")
        # --> add reference ellipse
        ref_ellip = patches.Ellipse((0, 0.0), width=es, height=es, angle=0)
        ref_ellip.set_facecolor((0, 0, 0))
        ref_ax_loc = list(self.ax2.get_position().bounds)
        ref_ax_loc[0] *= 0.95
        ref_ax_loc[1] -= 0.17
        ref_ax_loc[2] = 0.1
        ref_ax_loc[3] = 0.1
        self.ref_ax = self.fig.add_axes(ref_ax_loc, aspect="equal")
        self.ref_ax.add_artist(ref_ellip)
        self.ref_ax.set_xlim(-es / 2.0 * 1.05, es / 2.0 * 1.05)
        self.ref_ax.set_ylim(-es / 2.0 * 1.05, es / 2.0 * 1.05)
        plt.setp(self.ref_ax.xaxis.get_ticklabels(), visible=False)
        plt.setp(self.ref_ax.yaxis.get_ticklabels(), visible=False)
        self.ref_ax.set_title(r"$\Phi$ = 1")

        self.ax.set_axisbelow(True)
