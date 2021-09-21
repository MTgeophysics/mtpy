# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:47:37 2021

:copyright:
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""

#
# ==============================================================================
# Imports
# ==============================================================================
# standard imports
import os

try:
    from PyQt5 import QtCore, QtWidgets, QtGui
except ImportError:
    raise ImportError("This version needs PyQt5")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import matplotlib.widgets as mplwidgets

import mtpy.imaging.mtplottools as mtplottools
import mtpy.modeling.modem as modem
from .response_plot_settings import PlotSettings

# ==============================================================================
# plot part
# ==============================================================================


class PlotResponses(QtWidgets.QWidget):
    """
    the plot and list of stations
    """

    def __init__(self, data_fn=None, resp_fn=None):
        super(PlotResponses, self).__init__()

        self.file_watcher_dfn = QtCore.QFileSystemWatcher()
        self.file_watcher_dfn.fileChanged.connect(self.file_changed_dfn)

        self.modem_data = None
        self.modem_resp = None

        self.station = None

        self._modem_data_copy = None

        self._plot_z = False
        self.plot_settings = PlotSettings()

        self._ax = None
        self._ax2 = None
        self._key = "z"
        self._ax_index = 0
        self.ax_list = None
        self.phase_flip_comp = None
        self.add_error_comp = None
        self.ss_comp = None
        self.add_t_error = 0.02
        self.add_z_error = 5.0
        self.static_shift = 1.0

        self.setup_ui()

        self._data_fn = data_fn
        self._resp_fn = resp_fn

    # ------------------------------------------------
    # make the data_fn and resp_fn properties so that if they are reset
    # they will read in the data to a new modem.Data object
    # trying to use decorators for syntactical sugar
    @property
    def data_fn(self):
        return self._data_fn

    @data_fn.setter
    def data_fn(self, data_fn):
        self._data_fn = os.path.abspath(data_fn)
        self.file_watcher_dfn.addPath(self._data_fn)

        # create new modem data object
        self.modem_data = modem.Data()
        self.modem_data.read_data_file(self._data_fn)

        # make a back up copy that will be unchanged
        # that way we can revert back
        self._modem_data_copy = modem.Data()
        self._modem_data_copy.read_data_file(self._data_fn)

        self.dirpath = os.path.dirname(self._data_fn)

        # fill list of stations
        station_list = np.array(sorted(self.modem_data.mt_dict.keys()))
        self.list_widget.clear()
        for station in station_list:
            self.list_widget.addItem(station)

        if self.station is None:
            self.station = station_list[0]

        self.plot()

    @property
    def resp_fn(self):
        return self._resp_fn

    @resp_fn.setter
    def resp_fn(self, resp_fn):
        self._resp_fn = os.path.abspath(resp_fn)
        self.modem_resp = modem.Data()

        self.modem_resp.read_data_file(self._resp_fn)
        self.plot()

    @property
    def plot_z(self):
        return self._plot_z

    @plot_z.setter
    def plot_z(self, value):
        self._plot_z = value
        self.plot()

    @staticmethod
    def fmt_button(color):
        return "QPushButton {background-color: " + f"{color}" + "; font-weight: bold}"

    # ----------------------------
    def setup_ui(self):
        """
        setup the user interface with list of stations on the left and the
        plot on the right.  There will be a button for save edits.
        """

        # make a widget that will be the station list
        self.list_widget = QtWidgets.QListWidget()
        self.list_widget.itemClicked.connect(self.get_station)
        self.list_widget.currentItemChanged.connect(self.get_station)
        self.list_widget.setMaximumWidth(150)

        self.save_edits_button = QtWidgets.QPushButton()
        self.save_edits_button.setText("Save Edits")
        self.save_edits_button.setStyleSheet(self.fmt_button("#FF9E9E"))
        self.save_edits_button.pressed.connect(self.save_edits)

        self.apply_edits_button = QtWidgets.QPushButton()
        self.apply_edits_button.setText("Apply Edits")
        self.apply_edits_button.setStyleSheet(self.fmt_button("#ffab2e"))
        self.apply_edits_button.pressed.connect(self.apply_edits)

        self.interpolate_button = QtWidgets.QPushButton()
        self.interpolate_button.setText("Interpolate")
        self.interpolate_button.setStyleSheet(self.fmt_button("#ffff30"))
        self.interpolate_button.pressed.connect(self.apply_interpolation)

        self.flip_phase_button = QtWidgets.QPushButton()
        self.flip_phase_button.setText("Flip Phase")
        self.flip_phase_button.setStyleSheet(self.fmt_button("#A3FF8C"))
        self.flip_phase_button.pressed.connect(self.apply_flip_phase)

        self.flip_phase_combo = QtWidgets.QComboBox()
        self.flip_phase_combo.addItems(["", "Zxx", "Zxy", "Zyx", "Zyy", "Tx", "Ty"])
        self.flip_phase_combo.currentIndexChanged.connect(self.set_phase_flip_comp)
        flip_phase_layout = QtWidgets.QHBoxLayout()
        flip_phase_layout.addWidget(self.flip_phase_button)
        flip_phase_layout.addWidget(self.flip_phase_combo)

        self.add_error_button = QtWidgets.QPushButton()
        self.add_error_button.setText("Add Error")
        self.add_error_button.setStyleSheet(self.fmt_button("#8FFFF0"))
        self.add_error_button.pressed.connect(self.apply_add_error)

        self.add_error_combo = QtWidgets.QComboBox()
        self.add_error_combo.addItems(["", "Zxx", "Zxy", "Zyx", "Zyy", "Tx", "Ty"])
        self.add_error_combo.currentIndexChanged.connect(self.set_error_comp)
        add_error_layout = QtWidgets.QHBoxLayout()
        add_error_layout.addWidget(self.add_error_button)
        add_error_layout.addWidget(self.add_error_combo)

        self.add_z_error_text = QtWidgets.QLineEdit(f"{self.add_z_error:.2f}")
        self.add_z_error_text.setValidator(QtGui.QDoubleValidator(0, 100, 2))
        self.add_z_error_text.setMaximumWidth(70)
        self.add_z_error_text.editingFinished.connect(self.set_z_error_value)
        self.add_z_error_label = QtWidgets.QLabel("Z (%)")
        self.add_z_error_label.setMaximumWidth(50)

        self.add_t_error_text = QtWidgets.QLineEdit(f"{self.add_t_error:.2f}")
        self.add_t_error_text.setValidator(QtGui.QDoubleValidator(0, 1, 2))
        self.add_t_error_text.setMaximumWidth(70)
        self.add_t_error_text.editingFinished.connect(self.set_t_error_value)
        self.add_t_error_label = QtWidgets.QLabel("T (abs)")
        self.add_t_error_label.setMaximumWidth(50)
        add_z_error_layout = QtWidgets.QHBoxLayout()
        add_z_error_layout.addWidget(self.add_z_error_label)
        add_z_error_layout.addWidget(self.add_z_error_text)
        add_t_error_layout = QtWidgets.QHBoxLayout()
        add_t_error_layout.addWidget(self.add_t_error_label)
        add_t_error_layout.addWidget(self.add_t_error_text)

        self.static_shift_button = QtWidgets.QPushButton()
        self.static_shift_button.setText("Static Shift")
        self.static_shift_button.setStyleSheet(self.fmt_button("#a7d7cd"))
        self.static_shift_button.pressed.connect(self.apply_static_shift)
        self.static_shift_button.setMaximumWidth(80)
        self.static_shift_combo = QtWidgets.QComboBox()
        self.static_shift_combo.addItems(["", "Zx", "Zy"])
        self.static_shift_combo.currentIndexChanged.connect(self.set_ss_comp)
        self.static_shift_combo.setMaximumWidth(35)
        self.ss_text = QtWidgets.QLineEdit(f"{self.static_shift:.2f}")
        self.ss_text.setValidator(QtGui.QDoubleValidator(-100, 100, 2))
        self.ss_text.editingFinished.connect(self.set_ss_value)
        self.ss_text.setMaximumWidth(35)
        static_shift_layout = QtWidgets.QHBoxLayout()

        static_shift_layout.addWidget(self.static_shift_button)
        static_shift_layout.addWidget(self.static_shift_combo)
        static_shift_layout.addWidget(self.ss_text)

        self.undo_button = QtWidgets.QPushButton()
        self.undo_button.setText("Undo")
        self.undo_button.setStyleSheet(self.fmt_button("#9C9CFF"))
        self.undo_button.pressed.connect(self.apply_undo)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        self.mpl_widget.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.mpl_widget.setFocus()

        # be able to edit the data
        self.mpl_widget.mpl_connect("pick_event", self.on_pick)
        self.mpl_widget.mpl_connect("axes_enter_event", self.in_axes)
        self.mpl_widget.mpl_connect("button_press_event", self.on_pick)

        # make sure the figure takes up the entire plottable space
        self.mpl_widget.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, self)

        # set the layout for the plot
        mpl_vbox = QtWidgets.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)

        left_layout = QtWidgets.QVBoxLayout()
        left_layout.addWidget(self.list_widget)
        left_layout.addWidget(self.apply_edits_button)
        left_layout.addWidget(self.interpolate_button)
        left_layout.addLayout(flip_phase_layout)
        left_layout.addLayout(add_error_layout)
        left_layout.addLayout(add_z_error_layout)
        left_layout.addLayout(add_t_error_layout)
        left_layout.addLayout(static_shift_layout)
        left_layout.addWidget(self.undo_button)
        left_layout.addWidget(self.save_edits_button)

        # set the layout the main window
        layout = QtWidgets.QHBoxLayout()
        layout.addLayout(left_layout)
        layout.addLayout(mpl_vbox)

        self.setLayout(layout)
        self.mpl_widget.updateGeometry()

    def get_station(self, widget_item):
        """
        get the station name from the clicked station
        """
        try:
            self.station = str(widget_item.text())
        except AttributeError:
            self.station = self.list_widget.item(0).text()
            print(f"Station selected does not exist, setting to {self.station}")
        self.plot()

    def file_changed_dfn(self):
        """
        data file changed outside the program reload it
        """

        print("{0} changed".format(self.data_fn))
        self.data_fn = os.path.abspath(self._data_fn)

    def save_edits(self):
        """
        save edits to another file
        """
        fn_dialog = QtWidgets.QFileDialog()
        save_fn = str(
            fn_dialog.getSaveFileName(caption="Choose File to save", filter="*.dat")[0]
        )

        self.modem_data.write_data_file(
            save_path=os.path.dirname(save_fn),
            fn_basename=os.path.basename(save_fn),
            compute_error=False,
            fill=True,
            elevation=self.modem_data.topography,
        )

    def apply_edits(self):
        self.plot()

    def apply_interpolation(self):
        mt_obj = self.modem_data.mt_dict[self.station]
        mt_obj.Z, mt_obj.Tipper = mt_obj.interpolate(mt_obj.Z.freq)

        self.plot()

    def apply_undo(self):
        self.modem_data.mt_dict[self.station] = self._modem_data_copy.mt_dict[
            self.station
        ].copy()
        self.plot()

    def set_phase_flip_comp(self):
        self.phase_flip_comp = str(self.flip_phase_combo.currentText()).lower()

    def apply_flip_phase(self):
        (
            self.modem_data.data_array,
            self.modem_data.mt_dict,
        ) = self.modem_data.flip_phase(self.station, **{self.phase_flip_comp: True})
        self.plot()

    def set_error_comp(self):
        self.add_error_comp = str(self.add_error_combo.currentText()).lower()

    def set_z_error_value(self):
        try:
            self.add_z_error = float(self.add_z_error_text.text().strip())
        except ValueError:
            self.add_z_error = 1.0
        self.add_z_error_text.setText(f"{self.add_z_error:.2f}")

    def set_t_error_value(self):
        try:
            self.add_t_error = float(self.add_t_error_text.text().strip())
        except ValueError:
            self.add_t_error = 0.0
        self.add_t_error_text.setText(f"{self.add_t_error:.2f}")

    def apply_add_error(self):
        self.modem_data.data_array, self.modem_data.mt_dict = self.modem_data.add_error(
            self.station,
            [self.add_error_comp],
            z_value=self.add_z_error,
            t_value=self.add_t_error,
        )
        self.plot()

    def set_ss_comp(self):
        self.ss_comp = str(self.static_shift_combo.currentText()).lower()

    def set_ss_value(self):
        try:
            self.static_shift = float(self.ss_text.text().strip())
        except ValueError:
            self.static_shift = 1.0
        self.ss_text.setText(f"{self.static_shift:.2f}")

    def apply_static_shift(self):
        """
        Remove static shift
        """

        # be sure to apply the static shift to the original data
        if self.ss_comp.lower() == "zx":
            new_z_obj = self.modem_data.mt_dict[self.station].remove_static_shift(
                ss_x=self.static_shift
            )
            self.modem_data.mt_dict[self.station].Z = new_z_obj
            self.plot()
        elif self.ss_comp.lower() == "zy":
            new_z_obj = self.modem_data.mt_dict[self.station].remove_static_shift(
                ss_y=self.static_shift
            )
            self.modem_data.mt_dict[self.station].Z = new_z_obj
            self.plot()
        else:
            pass

    def plot(self):
        """
        plot the data
        """

        if self.station is None:
            return

        z_obj = self.modem_data.mt_dict[self.station].Z
        t_obj = self.modem_data.mt_dict[self.station].Tipper
        period = self.modem_data.period_list

        # need to make sure that resistivity and phase is computed
        z_obj.compute_resistivity_phase()

        plt.rcParams["font.size"] = self.plot_settings.fs
        fontdict = {"size": self.plot_settings.fs + 2, "weight": "bold"}

        # --> make key word dictionaries for plotting
        kw_xx = {
            "color": self.plot_settings.cted,
            "marker": self.plot_settings.mted,
            "ms": self.plot_settings.ms,
            "ls": ":",
            "lw": self.plot_settings.lw,
            "e_capsize": self.plot_settings.e_capsize,
            "e_capthick": self.plot_settings.e_capthick,
            "picker": 3,
        }

        kw_yy = {
            "color": self.plot_settings.ctmd,
            "marker": self.plot_settings.mtmd,
            "ms": self.plot_settings.ms,
            "ls": ":",
            "lw": self.plot_settings.lw,
            "e_capsize": self.plot_settings.e_capsize,
            "e_capthick": self.plot_settings.e_capthick,
            "picker": 3,
        }

        # convert to apparent resistivity and phase
        if self.plot_z == True:
            scaling = np.zeros_like(z_obj.z)
            for ii in range(2):
                for jj in range(2):
                    scaling[:, ii, jj] = 1.0 / np.sqrt(z_obj.freq)
            plot_res = abs(z_obj.z.real * scaling)
            plot_res_err = abs(z_obj.z_err * scaling)
            plot_phase = abs(z_obj.z.imag * scaling)
            plot_phase_err = abs(z_obj.z_err * scaling)
            h_ratio = [1, 1, 0.5]

        elif self.plot_z == False:
            plot_res = z_obj.resistivity
            plot_res_err = z_obj.resistivity_err
            plot_phase = z_obj.phase
            plot_phase_err = z_obj.phase_err
            h_ratio = [1.5, 1, 0.5]

        # find locations where points have been masked
        nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
        nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
        nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
        nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
        ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
        nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]

        self.figure.clf()
        self.figure.suptitle(str(self.station), fontdict=fontdict)

        # set the grid of subplots
        if np.all(t_obj.tipper == 0.0) == True:
            self.plot_tipper = False
        else:
            self.plot_tipper = True

        gs = gridspec.GridSpec(3, 4, height_ratios=h_ratio)
        gs.update(
            wspace=self.plot_settings.subplot_wspace,
            left=self.plot_settings.subplot_left,
            top=self.plot_settings.subplot_top,
            bottom=self.plot_settings.subplot_bottom,
            right=self.plot_settings.subplot_right,
            hspace=self.plot_settings.subplot_hspace,
        )

        axrxx = self.figure.add_subplot(gs[0, 0])
        axrxy = self.figure.add_subplot(gs[0, 1], sharex=axrxx)
        axryx = self.figure.add_subplot(gs[0, 2], sharex=axrxx)
        axryy = self.figure.add_subplot(gs[0, 3], sharex=axrxx)

        axpxx = self.figure.add_subplot(gs[1, 0])
        axpxy = self.figure.add_subplot(gs[1, 1], sharex=axrxx)
        axpyx = self.figure.add_subplot(gs[1, 2], sharex=axrxx)
        axpyy = self.figure.add_subplot(gs[1, 3], sharex=axrxx)

        axtxr = self.figure.add_subplot(gs[2, 0], sharex=axrxx)
        axtxi = self.figure.add_subplot(gs[2, 1], sharex=axrxx)
        axtyr = self.figure.add_subplot(gs[2, 2], sharex=axrxx)
        axtyi = self.figure.add_subplot(gs[2, 3], sharex=axrxx)

        self.ax_list = [
            axrxx,
            axrxy,
            axryx,
            axryy,
            axpxx,
            axpxy,
            axpyx,
            axpyy,
            axtxr,
            axtxi,
            axtyr,
            axtyi,
        ]

        # plot data response
        erxx = mtplottools.plot_errorbar(
            axrxx, period[nzxx], plot_res[nzxx, 0, 0], plot_res_err[nzxx, 0, 0], **kw_xx
        )
        erxy = mtplottools.plot_errorbar(
            axrxy, period[nzxy], plot_res[nzxy, 0, 1], plot_res_err[nzxy, 0, 1], **kw_xx
        )
        eryx = mtplottools.plot_errorbar(
            axryx, period[nzyx], plot_res[nzyx, 1, 0], plot_res_err[nzyx, 1, 0], **kw_yy
        )
        eryy = mtplottools.plot_errorbar(
            axryy, period[nzyy], plot_res[nzyy, 1, 1], plot_res_err[nzyy, 1, 1], **kw_yy
        )
        # plot phase
        epxx = mtplottools.plot_errorbar(
            axpxx,
            period[nzxx],
            plot_phase[nzxx, 0, 0],
            plot_phase_err[nzxx, 0, 0],
            **kw_xx,
        )
        epxy = mtplottools.plot_errorbar(
            axpxy,
            period[nzxy],
            plot_phase[nzxy, 0, 1],
            plot_phase_err[nzxy, 0, 1],
            **kw_xx,
        )
        epyx = mtplottools.plot_errorbar(
            axpyx,
            period[nzyx],
            plot_phase[nzyx, 1, 0],
            plot_phase_err[nzyx, 1, 0],
            **kw_yy,
        )
        epyy = mtplottools.plot_errorbar(
            axpyy,
            period[nzyy],
            plot_phase[nzyy, 1, 1],
            plot_phase_err[nzyy, 1, 1],
            **kw_yy,
        )

        # plot tipper
        if self.plot_tipper == True:
            ertx = mtplottools.plot_errorbar(
                axtxr,
                period[ntx],
                t_obj.tipper[ntx, 0, 0].real,
                t_obj.tipper_err[ntx, 0, 0],
                **kw_xx,
            )
            erty = mtplottools.plot_errorbar(
                axtyr,
                period[nty],
                t_obj.tipper[nty, 0, 1].real,
                t_obj.tipper_err[nty, 0, 1],
                **kw_yy,
            )

            eptx = mtplottools.plot_errorbar(
                axtxi,
                period[ntx],
                t_obj.tipper[ntx, 0, 0].imag,
                t_obj.tipper_err[ntx, 0, 0],
                **kw_xx,
            )
            epty = mtplottools.plot_errorbar(
                axtyi,
                period[nty],
                t_obj.tipper[nty, 0, 1].imag,
                t_obj.tipper_err[nty, 0, 1],
                **kw_yy,
            )

        # ----------------------------------------------
        # get error bar list for editing later
        if self.plot_tipper == False:
            try:
                self._err_list = [
                    [erxx[1][0], erxx[1][1], erxx[2][0]],
                    [erxy[1][0], erxy[1][1], erxy[2][0]],
                    [eryx[1][0], eryx[1][1], eryx[2][0]],
                    [eryy[1][0], eryy[1][1], eryy[2][0]],
                    [epxx[1][0], epxx[1][1], epxx[2][0]],
                    [epxy[1][0], epxy[1][1], epxy[2][0]],
                    [epyx[1][0], epyx[1][1], epyx[2][0]],
                    [epyy[1][0], epyy[1][1], epyy[2][0]],
                ]
                line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
            except IndexError:
                print("Found no Z components for {0}".format(self.station))
                line_list = [[None], [None], [None], [None]]

                self._err_list = [
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                ]

        else:
            try:
                line_list = [
                    [erxx[0]],
                    [erxy[0]],
                    [eryx[0]],
                    [eryy[0]],
                    [ertx[0]],
                    [erty[0]],
                ]

                self._err_list = [
                    [erxx[1][0], erxx[1][1], erxx[2][0]],
                    [erxy[1][0], erxy[1][1], erxy[2][0]],
                    [eryx[1][0], eryx[1][1], eryx[2][0]],
                    [eryy[1][0], eryy[1][1], eryy[2][0]],
                    [epxx[1][0], epxx[1][1], epxx[2][0]],
                    [epxy[1][0], epxy[1][1], epxy[2][0]],
                    [epyx[1][0], epyx[1][1], epyx[2][0]],
                    [epyy[1][0], epyy[1][1], epyy[2][0]],
                    [ertx[1][0], ertx[1][1], ertx[2][0]],
                    [eptx[1][0], eptx[1][1], eptx[2][0]],
                    [erty[1][0], erty[1][1], erty[2][0]],
                    [epty[1][0], epty[1][1], epty[2][0]],
                ]
            except IndexError:
                print("Found no Z components for {0}".format(self.station))
                line_list = [[None], [None], [None], [None], [None], [None]]

                self._err_list = [
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                    [None, None, None],
                ]

        # ------------------------------------------
        # make things look nice
        # set titles of the Z components
        label_list = [["$Z_{xx}$"], ["$Z_{xy}$"], ["$Z_{yx}$"], ["$Z_{yy}$"]]
        for ax, label in zip(self.ax_list[0:4], label_list):
            ax.set_title(
                label[0], fontdict={"size": self.plot_settings.fs + 2, "weight": "bold"}
            )

        # set legends for tipper components
        # fake a line
        l1 = plt.Line2D([0], [0], linewidth=0, color="w", linestyle="None", marker=".")
        t_label_list = ["Re{$T_x$}", "Im{$T_x$}", "Re{$T_y$}", "Im{$T_y$}"]
        label_list += [["$T_{x}$"], ["$T_{y}$"]]
        for ax, label in zip(self.ax_list[-4:], t_label_list):
            ax.legend(
                [l1],
                [label],
                loc="upper left",
                markerscale=0.01,
                borderaxespad=0.05,
                labelspacing=0.01,
                handletextpad=0.05,
                borderpad=0.05,
                prop={"size": max([self.plot_settings.fs, 5])},
            )

        # --> set limits if input
        if self.plot_settings.res_xx_limits is not None:
            axrxx.set_ylim(self.plot_settings.res_xx_limits)
        if self.plot_settings.res_xy_limits is not None:
            axrxy.set_ylim(self.plot_settings.res_xy_limits)
        if self.plot_settings.res_yx_limits is not None:
            axryx.set_ylim(self.plot_settings.res_yx_limits)
        if self.plot_settings.res_yy_limits is not None:
            axryy.set_ylim(self.plot_settings.res_yy_limits)

        if self.plot_settings.phase_xx_limits is not None:
            axpxx.set_ylim(self.plot_settings.phase_xx_limits)
        if self.plot_settings.phase_xy_limits is not None:
            axpxy.set_ylim(self.plot_settings.phase_xy_limits)
        if self.plot_settings.phase_yx_limits is not None:
            axpyx.set_ylim(self.plot_settings.phase_yx_limits)
        if self.plot_settings.phase_yy_limits is not None:
            axpyy.set_ylim(self.plot_settings.phase_yy_limits)

        # set axis properties
        for aa, ax in enumerate(self.ax_list):
            ax.tick_params(axis="y", pad=self.plot_settings.ylabel_pad)
            ylabels = ax.get_yticks().tolist()
            if aa < 8:
                ylabels[-1] = ""
                ylabels[0] = ""
                ax.set_yticklabels(ylabels)
                plt.setp(ax.get_xticklabels(), visible=False)
                if self.plot_z == True:
                    ax.set_yscale("log", nonpositive="clip")

            else:
                ax.set_xlabel("Period (s)", fontdict=fontdict)

            if aa < 4 and self.plot_z is False:
                ax.set_yscale("log", nonpositive="clip")

            # set axes labels
            if aa == 0:
                if self.plot_z == False:
                    ax.set_ylabel(
                        "App. Res. ($\mathbf{\Omega \cdot m}$)", fontdict=fontdict
                    )
                elif self.plot_z == True:
                    ax.set_ylabel("Re[Z (mV/km nT)]", fontdict=fontdict)
            elif aa == 4:
                if self.plot_z == False:
                    ax.set_ylabel("Phase (deg)", fontdict=fontdict)
                elif self.plot_z == True:
                    ax.set_ylabel("Im[Z (mV/km nT)]", fontdict=fontdict)
            elif aa == 8:
                ax.set_ylabel("Tipper", fontdict=fontdict)

            if aa > 7:
                if self.plot_settings.tipper_limits is not None:
                    ax.set_ylim(self.plot_settings.tipper_limits)
                else:
                    pass

            ax.set_xscale("log", nonpositive="clip")
            ax.set_xlim(
                xmin=10 ** (np.floor(np.log10(period[0]))) * 1.01,
                xmax=10 ** (np.ceil(np.log10(period[-1]))) * 0.99,
            )
            ax.grid(True, alpha=0.25)

        # ----------------------------------------------
        # plot model response
        if self.modem_resp is not None:
            try:
                resp_z_obj = self.modem_resp.mt_dict[self.station].Z
                resp_z_err = np.nan_to_num((z_obj.z - resp_z_obj.z) / z_obj.z_err)
                resp_z_obj.compute_resistivity_phase()

                resp_t_obj = self.modem_resp.mt_dict[self.station].Tipper
                resp_t_err = np.nan_to_num(
                    (t_obj.tipper - resp_t_obj.tipper) / t_obj.tipper_err
                )
            except KeyError:
                print("Could not find {0} in .resp file".format(self.station))
                self.mpl_widget.draw()
                return

            try:
                resp_period = 1.0 / resp_z_obj.freq
            except TypeError:
                resp_period = 1.0 / resp_t_obj.freq

            # convert to apparent resistivity and phase
            if self.plot_z == True:
                scaling = np.zeros_like(resp_z_obj.z)
                for ii in range(2):
                    for jj in range(2):
                        scaling[:, ii, jj] = 1.0 / np.sqrt(resp_z_obj.freq)
                r_plot_res = abs(resp_z_obj.z.real * scaling)
                r_plot_phase = abs(resp_z_obj.z.imag * scaling)

            elif self.plot_z == False:
                r_plot_res = resp_z_obj.resistivity
                r_plot_phase = resp_z_obj.phase

            # find locations where points have been masked
            nzxx_r = np.nonzero(resp_z_obj.z[:, 0, 0])[0]
            nzxy_r = np.nonzero(resp_z_obj.z[:, 0, 1])[0]
            nzyx_r = np.nonzero(resp_z_obj.z[:, 1, 0])[0]
            nzyy_r = np.nonzero(resp_z_obj.z[:, 1, 1])[0]
            ntx_r = np.nonzero(resp_t_obj.tipper[:, 0, 0])[0]
            nty_r = np.nonzero(resp_t_obj.tipper[:, 0, 1])[0]

            rms_xx = resp_z_err[nzxx_r, 0, 0].std()
            rms_xy = resp_z_err[nzxy_r, 0, 1].std()
            rms_yx = resp_z_err[nzyx_r, 1, 0].std()
            rms_yy = resp_z_err[nzyy_r, 1, 1].std()

            # --> make key word dictionaries for plotting
            kw_xx = {
                "color": self.plot_settings.ctem,
                "marker": self.plot_settings.mtem,
                "ms": self.plot_settings.ms,
                "ls": ":",
                "lw": self.plot_settings.lw,
                "e_capsize": self.plot_settings.e_capsize,
                "e_capthick": self.plot_settings.e_capthick,
            }

            kw_yy = {
                "color": self.plot_settings.ctmm,
                "marker": self.plot_settings.mtmm,
                "ms": self.plot_settings.ms,
                "ls": ":",
                "lw": self.plot_settings.lw,
                "e_capsize": self.plot_settings.e_capsize,
                "e_capthick": self.plot_settings.e_capthick,
            }

            # plot data response
            rerxx = mtplottools.plot_errorbar(
                axrxx, resp_period[nzxx_r], r_plot_res[nzxx_r, 0, 0], None, **kw_xx
            )
            rerxy = mtplottools.plot_errorbar(
                axrxy, resp_period[nzxy_r], r_plot_res[nzxy_r, 0, 1], None, **kw_xx
            )
            reryx = mtplottools.plot_errorbar(
                axryx, resp_period[nzyx_r], r_plot_res[nzyx_r, 1, 0], None, **kw_yy
            )
            reryy = mtplottools.plot_errorbar(
                axryy, resp_period[nzyy_r], r_plot_res[nzyy_r, 1, 1], None, **kw_yy
            )
            # plot phase
            repxx = mtplottools.plot_errorbar(
                axpxx, resp_period[nzxx_r], r_plot_phase[nzxx_r, 0, 0], None, **kw_xx
            )
            repxy = mtplottools.plot_errorbar(
                axpxy, resp_period[nzxy_r], r_plot_phase[nzxy_r, 0, 1], None, **kw_xx
            )
            repyx = mtplottools.plot_errorbar(
                axpyx, resp_period[nzyx_r], r_plot_phase[nzyx_r, 1, 0], None, **kw_yy
            )
            repyy = mtplottools.plot_errorbar(
                axpyy, resp_period[nzyy_r], r_plot_phase[nzyy_r, 1, 1], None, **kw_yy
            )

            # plot tipper
            if self.plot_tipper == True:
                rertx = mtplottools.plot_errorbar(
                    axtxr,
                    resp_period[ntx_r],
                    resp_t_obj.tipper[ntx_r, 0, 0].real,
                    None,
                    **kw_xx,
                )
                rerty = mtplottools.plot_errorbar(
                    axtyr,
                    resp_period[nty_r],
                    resp_t_obj.tipper[nty_r, 0, 1].real,
                    None,
                    **kw_yy,
                )

                reptx = mtplottools.plot_errorbar(
                    axtxi,
                    resp_period[ntx_r],
                    resp_t_obj.tipper[ntx_r, 0, 0].imag,
                    None,
                    **kw_xx,
                )
                repty = mtplottools.plot_errorbar(
                    axtyi,
                    resp_period[nty_r],
                    resp_t_obj.tipper[nty_r, 0, 1].imag,
                    None,
                    **kw_yy,
                )

            if self.plot_tipper == False:
                line_list[0] += [rerxx[0]]
                line_list[1] += [rerxy[0]]
                line_list[2] += [reryx[0]]
                line_list[3] += [reryy[0]]
                label_list[0] += ["$Z^m_{xx}$ " + "rms={0:.2f}".format(rms_xx)]
                label_list[1] += ["$Z^m_{xy}$ " + "rms={0:.2f}".format(rms_xy)]
                label_list[2] += ["$Z^m_{yx}$ " + "rms={0:.2f}".format(rms_yx)]
                label_list[3] += ["$Z^m_{yy}$ " + "rms={0:.2f}".format(rms_yy)]
            else:
                line_list[0] += [rerxx[0]]
                line_list[1] += [rerxy[0]]
                line_list[2] += [reryx[0]]
                line_list[3] += [reryy[0]]
                line_list[4] += [rertx[0]]
                line_list[5] += [rerty[0]]
                label_list[0] += ["$Z^m_{xx}$ " + "rms={0:.2f}".format(rms_xx)]
                label_list[1] += ["$Z^m_{xy}$ " + "rms={0:.2f}".format(rms_xy)]
                label_list[2] += ["$Z^m_{yx}$ " + "rms={0:.2f}".format(rms_yx)]
                label_list[3] += ["$Z^m_{yy}$ " + "rms={0:.2f}".format(rms_yy)]
                label_list[4] += [
                    "$T^m_{x}$ " + "rms={0:.2f}".format(resp_t_err[ntx, 0, 0].std())
                ]
                label_list[5] += [
                    "$T^m_{y}$" + "rms={0:.2f}".format(resp_t_err[nty, 0, 1].std())
                ]

            legend_ax_list = self.ax_list[0:4]
            if self.plot_tipper == True:
                legend_ax_list += [self.ax_list[-4], self.ax_list[-2]]

            for aa, ax in enumerate(legend_ax_list):
                ax.legend(
                    line_list[aa],
                    label_list[aa],
                    loc=self.plot_settings.legend_loc,
                    bbox_to_anchor=self.plot_settings.legend_pos,
                    markerscale=self.plot_settings.legend_marker_scale,
                    borderaxespad=self.plot_settings.legend_border_axes_pad,
                    labelspacing=self.plot_settings.legend_label_spacing,
                    handletextpad=self.plot_settings.legend_handle_text_pad,
                    borderpad=self.plot_settings.legend_border_pad,
                    prop={"size": max([self.plot_settings.fs, 5])},
                )

        # make rectangular picker
        for ax, name in zip(
            self.ax_list,
            [
                "rxx",
                "rxy",
                "ryx",
                "ryy",
                "pxx",
                "pxy",
                "pyx",
                "pyy",
                "txr",
                "txi",
                "tyr",
                "tyi",
            ],
        ):
            setattr(
                self,
                f"rect_{name}",
                mplwidgets.RectangleSelector(
                    ax,
                    self.on_select_rect,
                    drawtype="box",
                    useblit=True,
                    interactive=True,
                    minspanx=5,
                    minspany=5,
                    spancoords="pixels",
                    button=[1],
                ),
            )

        self.mpl_widget.draw()

    def on_pick(self, event):
        """
        mask a data point when it is clicked on.
        """
        try:
            data_point = event.artist
            data_period = data_point.get_xdata()[event.ind]
            data_value = data_point.get_ydata()[event.ind]
        except AttributeError:
            print("Picked a bad point")
            return

        # get the indicies where the data point has been edited
        try:
            p_index = np.where(self.modem_data.period_list == data_period)[0][0]
            s_index = np.where(self.modem_data.data_array["station"] == self.station)[
                0
            ][0]
        except IndexError:
            return

        if self._key == "tip":
            data_value_2 = self.modem_data.mt_dict[self.station].Tipper.tipper[
                p_index, self._comp_index_x, self._comp_index_y
            ]
            if self._ax_index % 2 == 0:
                data_value_2 = data_value_2.imag
            else:
                data_value_2 = data_value_2.real

        elif self._key == "z":
            if self.plot_z == True:
                data_value_2 = self.modem_data.mt_dict[self.station].Z.z[
                    p_index, self._comp_index_x, self._comp_index_y
                ]

                if self._ax_index % 2 == 0:
                    data_value_2 = data_value_2.imag
                else:
                    data_value_2 = data_value_2.real
            elif self.plot_z == False and self._ax_index < 4:
                data_value_2 = self.modem_data.mt_dict[self.station].Z.phase[
                    p_index, self._comp_index_x, self._comp_index_y
                ]
            elif self.plot_z == False and self._ax_index >= 4:
                data_value_2 = self.modem_data.mt_dict[self.station].Z.resistivity[
                    p_index, self._comp_index_x, self._comp_index_y
                ]

        if event.mouseevent.button == 1:
            # mask the point in the data mt_dict

            self.modem_data.data_array[s_index][self._key][
                p_index, self._comp_index_x, self._comp_index_y
            ] = (0 + 0j)

            if self._key == "tip":
                self.modem_data.mt_dict[self.station].Tipper.tipper[
                    p_index, self._comp_index_x, self._comp_index_y
                ] = (0 + 0j)
            elif self._key == "z":
                self.modem_data.mt_dict[self.station].Z.z[
                    p_index, self._comp_index_x, self._comp_index_y
                ] = (0 + 0j)

            # plot the points as masked
            self._ax.plot(
                data_period,
                data_value,
                color=(0, 0, 0),
                marker="x",
                ms=self.plot_settings.ms * 2,
                mew=4,
            )

            self._ax2.plot(
                data_period,
                data_value_2,
                color=(0, 0, 0),
                marker="x",
                ms=self.plot_settings.ms * 2,
                mew=4,
            )

            self._ax.figure.canvas.draw()
            self._ax2.figure.canvas.draw()

        # Increase error bars
        if event.mouseevent.button == 3:
            # make sure just checking the top plots

            # put the new error into the error array
            if self._key == "tip":
                err = self.modem_data.mt_dict[self.station].Tipper.tipper_err[
                    p_index, self._comp_index_x, self._comp_index_y
                ]
                err = err + self.add_t_error
                self.modem_data.mt_dict[self.station].Tipper.tipper_err[
                    p_index, self._comp_index_x, self._comp_index_y
                ] = err

            if self._key == "z":
                err = self.modem_data.mt_dict[self.station].Z.z_err[
                    p_index, self._comp_index_x, self._comp_index_y
                ]
                err = err * self.add_z_error

                self.modem_data.mt_dict[self.station].Z.z_err[
                    p_index, self._comp_index_x, self._comp_index_y
                ] = err

            self.modem_data.data_array[s_index][self._key + "_err"][
                p_index, self._comp_index_x, self._comp_index_y
            ] = err

            # make error bar array
            try:
                e_index = event.ind[0]
                eb = self._err_list[self._ax_index][2].get_paths()[e_index].vertices
            except IndexError:
                return

            # make ecap array
            ecap_l = self._err_list[self._ax_index][0].get_data()[1][e_index]
            ecap_u = self._err_list[self._ax_index][1].get_data()[1][e_index]

            # change apparent resistivity error
            if self._key == "tip":
                neb_u = eb[0, 1] - self.add_t_error
                neb_l = eb[1, 1] + self.add_t_error
                ecap_l = ecap_l - self.add_t_error
                ecap_u = ecap_u + self.add_t_error
            elif self._key == "z":
                if self.plot_z:
                    neb_u = eb[0, 1] - self.add_z_error * abs(eb[0, 1]) / 2
                    neb_l = eb[1, 1] + self.add_z_error * abs(eb[1, 1]) / 2
                    ecap_l = ecap_l - self.add_z_error * abs(eb[0, 1]) / 2
                    ecap_u = ecap_u + self.add_z_error * abs(eb[1, 1]) / 2
                elif not self.plot_z:
                    if self._ax_index < 4:
                        neb_u = eb[0, 1] - self.add_z_error * np.sqrt(abs(eb[0, 1]))
                        neb_l = eb[1, 1] + self.add_z_error * np.sqrt(abs(eb[1, 1]))
                        ecap_l = ecap_l - self.add_z_error * np.sqrt(abs(eb[0, 1]))
                        ecap_u = ecap_u + self.add_z_error * np.sqrt(abs(eb[1, 1]))
                    else:
                        neb_u = eb[0, 1] - self.add_z_error / 100 * abs(eb[0, 1]) * 4
                        neb_l = eb[1, 1] + self.add_z_error / 100 * abs(eb[1, 1]) * 4
                        ecap_l = ecap_l - self.add_z_error / 100 * abs(eb[0, 1]) * 4
                        ecap_u = ecap_u + self.add_z_error / 100 * abs(eb[1, 1]) * 4

            # set the new error bar values
            eb[0, 1] = neb_u
            eb[1, 1] = neb_l

            # reset the error bars and caps
            ncap_l = self._err_list[self._ax_index][0].get_data()
            ncap_u = self._err_list[self._ax_index][1].get_data()
            ncap_l[1][e_index] = ecap_l
            ncap_u[1][e_index] = ecap_u

            # set the values
            self._err_list[self._ax_index][0].set_data(ncap_l)
            self._err_list[self._ax_index][1].set_data(ncap_u)
            self._err_list[self._ax_index][2].get_paths()[e_index].vertices = eb

            # need to redraw the figure
            self._ax.figure.canvas.draw()

    def in_axes(self, event):
        """
        figure out which axes you just chose the point from
        """

        ax_index_dict = {
            0: (0, 0),
            1: (0, 1),
            2: (1, 0),
            3: (1, 1),
            4: (0, 0),
            5: (0, 1),
            6: (1, 0),
            7: (1, 1),
            8: (0, 0),
            9: (0, 0),
            10: (0, 1),
            11: (0, 1),
        }

        ax_pairs = {
            0: 4,
            1: 5,
            2: 6,
            3: 7,
            4: 0,
            5: 1,
            6: 2,
            7: 3,
            8: 9,
            9: 8,
            10: 11,
            11: 10,
        }
        # make the axis an attribute
        self._ax = event.inaxes

        # find the component index so that it can be masked
        for ax_index, ax in enumerate(self.ax_list):
            if ax == event.inaxes:
                self._comp_index_x, self._comp_index_y = ax_index_dict[ax_index]
                self._ax_index = ax_index
                self._ax2 = self.ax_list[ax_pairs[ax_index]]
                if ax_index < 8:
                    self._key = "z"

                else:
                    self._key = "tip"

    def _get_frequency_range(self, period_01, period_02):

        fmin = min([period_01, period_02])
        fmax = max([period_01, period_02])
        prange = np.where(
            (self.modem_data.period_list >= fmin)
            & (self.modem_data.period_list <= fmax)
        )

        return prange

    def on_select_rect(self, eclick, erelease):
        x1 = eclick.xdata
        x2 = erelease.xdata

        f_idx = self._get_frequency_range(x1, x2)

        print(self._key, self._ax_index)
        for ff in f_idx:
            period = self.modem_data.period_list[ff]
            if self._key == "z":
                self._ax.plot(
                    period,
                    self.modem_data.mt_dict[self.station].Z.resistivity[
                        ff, self._comp_index_x, self._comp_index_y
                    ],
                    color=(0, 0, 0),
                    marker="x",
                    ms=self.plot_settings.ms * 2,
                    mew=4,
                )
                self._ax2.plot(
                    period,
                    self.modem_data.mt_dict[self.station].Z.phase[
                        ff, self._comp_index_x, self._comp_index_y
                    ],
                    color=(0, 0, 0),
                    marker="x",
                    ms=self.plot_settings.ms * 2,
                    mew=4,
                )

                self.modem_data.mt_dict[self.station].Z.z[
                    ff, self._comp_index_x, self._comp_index_y
                ] = (0.0 + 0.0 * 1j)
                self.modem_data.mt_dict[self.station].Z.z_err[
                    ff, self._comp_index_x, self._comp_index_y
                ] = 0.0
            elif self._key == "tip":
                self._ax.plot(
                    period,
                    self.modem_data.mt_dict[self.station]
                    .Tipper.tipper[ff, self._comp_index_x, self._comp_index_y]
                    .real,
                    color=(0, 0, 0),
                    marker="x",
                    ms=self.plot_settings.ms * 2,
                    mew=4,
                )
                self._ax2.plot(
                    period,
                    self.modem_data.mt_dict[self.station]
                    .Tipper.tipper[ff, self._comp_index_x, self._comp_index_y]
                    .imag,
                    color=(0, 0, 0),
                    marker="x",
                    ms=self.plot_settings.ms * 2,
                    mew=4,
                )

                self.modem_data.mt_dict[self.station].Tipper.tipper[
                    ff, self._comp_index_x, self._comp_index_y
                ] = (0.0 + 0.0 * 1j)
                self.modem_data.mt_dict[self.station].Tipper.tipper_err[
                    ff, self._comp_index_x, self._comp_index_y
                ] = 0.0
        self._ax.figure.canvas.draw()
        self._ax2.figure.canvas.draw()
