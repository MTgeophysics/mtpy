import numpy as np
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import pyqtSignal

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.gui.plot_parameter import COLORS, SIMPLE_COLORS
from mtpy.gui.SmartMT.ui_asset.groupbox_arrow import Ui_GroupBox_Arrow
from mtpy.gui.SmartMT.ui_asset.groupbox_ellipse import Ui_GroupBoxEllipse
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_index import Ui_GroupBox_Frequency_Period_Index
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_single import Ui_groupBoxFrequency_pereiod_single
from mtpy.gui.SmartMT.ui_asset.groupbox_linedir import Ui_GroupBox_Linedir
from mtpy.gui.SmartMT.ui_asset.groupbox_mesh_grid import Ui_GroupBox_mash_grid
from mtpy.gui.SmartMT.ui_asset.groupbox_padding import Ui_GroupBox_Padding
from mtpy.gui.SmartMT.ui_asset.groupbox_plot_control_mt_response import Ui_GroupBox_plot_control_mt_response
from mtpy.gui.SmartMT.ui_asset.groupbox_plot_control_resistivity_phase_pseudo_section import \
    Ui_GroupBox_plot_control_resistivity_phase_pseudo_section
from mtpy.gui.SmartMT.ui_asset.groupbox_rotation import Ui_GroupBox_Rotation
from mtpy.gui.SmartMT.ui_asset.groupbox_scale import Ui_GroupBox_Scale
from mtpy.gui.SmartMT.ui_asset.groupbox_station_select import Ui_GroupBox_Station_Select
from mtpy.gui.SmartMT.ui_asset.groupbox_stretch import Ui_GroupBox_Stretch
from mtpy.gui.SmartMT.ui_asset.groupbox_tolerance import Ui_GroupBoxTolerance
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_multiple import Ui_groupBoxZ_Component_Multiple
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_single import Ui_groupBoxZ_Component_Single
from mtpy.imaging.mtcolors import cmapdict


class ZComponentMultiple(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Multiple()
        self.ui.setupUi(self)
        # z-component checkbox logic
        self.ui.checkBox_zyx.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_zxy.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_det.stateChanged.connect(self._multiple_zcomponent_logic)

    def _multiple_zcomponent_logic(self, int):
        """
        set up at least one component selected
        :return:
        """
        counter = 0
        if self.ui.checkBox_det.isChecked() + self.ui.checkBox_zxy.isChecked() + self.ui.checkBox_zyx.isChecked() == 1:
            # only one checkbox is checked, lock the checked box
            if self.ui.checkBox_det.isChecked():
                self.ui.checkBox_det.setEnabled(False)
            elif self.ui.checkBox_zxy.isChecked():
                self.ui.checkBox_zxy.setEnabled(False)
            else:
                self.ui.checkBox_zyx.setEnabled(False)
        else:
            self.ui.checkBox_det.setEnabled(True)
            self.ui.checkBox_zxy.setEnabled(True)
            self.ui.checkBox_zyx.setEnabled(True)

    def get_selection(self):
        zcomponent = []
        if self.ui.checkBox_det.isChecked():
            zcomponent.append('det')
        if self.ui.checkBox_zxy.isChecked():
            zcomponent.append('zxy')
        if self.ui.checkBox_zyx.isChecked():
            zcomponent.append('zyx')
        return zcomponent


class ZComponentSingle(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Single()
        self.ui.setupUi(self)

    def get_selection(self):
        if self.ui.radioButton_det.isChecked():
            return 'det'
        elif self.ui.radioButton_zxy.isChecked():
            return 'zxy'
        elif self.ui.radioButton_zyx.isChecked():
            return 'zyx'


class FrequencySingle(QtGui.QGroupBox):
    """
    Frequency selection (single frequency)
    """

    _unit_period = 'second'
    _unit_frequency = 'Hz'
    _title_period = 'Period'
    _title_frequency = 'Frequency'

    def __init__(self, parent, use_period=False):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self.use_period = use_period
        self.ui = Ui_groupBoxFrequency_pereiod_single()
        self.ui.setupUi(self)
        self._histogram = FrequencySingle.FrequencyHistogram()
        self.set_use_period(self.use_period)
        # add matplotlib canvas
        self.ui.verticalLayoutFrequencyPeriod.addWidget(self._histogram)
        # connect components
        self.ui.comboBoxPeriod.currentIndexChanged.connect(self.update_histogram)
        self.ui.comboBoxPeriod.editTextChanged.connect(self.update_histogram)
        self._histogram.mpl_connect('button_release_event', self._mouse_pick)

    def toggle_time_scale(self, *args):
        self.use_period = not self.use_period
        self.set_use_period(self.use_period)

    def set_use_period(self, use_period=False):
        if use_period:
            self._histogram.set_unit(self._unit_period)
            self._histogram.set_title(self._title_period)
            title = '%s (%s)' % (self._title_period, self._unit_period)
        else:
            self._histogram.set_unit(self._unit_frequency)
            self._histogram.set_title(self._title_frequency)
            title = '%s (%s)' % (self._title_frequency, self._unit_frequency)
        self.setTitle(title)
        self._update_frequency()

    def _mouse_pick(self, event):
        if not event.inaxes:
            return
        x = event.xdata
        self.ui.comboBoxPeriod.setEditText("%.5f" % x)

    def get_frequency(self):
        return float(self.ui.comboBoxPeriod.currentText())

    def update_histogram(self):
        value = float(self.ui.comboBoxPeriod.currentText())
        self._histogram.set_current_frequency(value)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._update_frequency()

    def _update_frequency(self):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))

            if self.use_period:
                all_periods = 1.0 / np.array(all_freqs)
                # self._histogram.set_data(all_periods)
                all_unique = sorted(list(set(all_periods)))

            else:
                # self._histogram.set_data(all_freqs)
                all_unique = sorted(list(set(all_freqs)))
            self._histogram.set_data(all_unique)
            self._histogram.update_figure()
            # sort all frequencies in ascending order
            for period in all_unique:
                self.ui.comboBoxPeriod.addItem("%.5f" % period)
            self.ui.comboBoxPeriod.setCurrentIndex(0)
            self.update_histogram()

    class FrequencyHistogram(MPLCanvas):
        def __init__(self, parent=None, width=5, hight=2, dpi=100):
            self.artists = dict()
            self._frequency = None
            self._current_frequency = None
            self._title = None
            self._unit = None
            MPLCanvas.__init__(self, parent, width, hight, dpi)
            self._lx = None
            self.cursor = None

            # self.mpl_connect('motion_notify_event', self.cursor)
            self.mpl_connect('button_release_event', self.mouse_pick)
            self.setMinimumSize(200, 150)
            self.resize(self.sizeHint())

        # def mouse_move(self, event):
        #     if not event.inaxes:
        #         return
        #     x = event.xdata
        #     y = event.ydata
        #     if self._cursor_x is None:
        #         self._cursor_x = self._axes.axvline(linewidth=1, color="green")
        #     if self._cursor_text is None:
        #         self._cursor_text = self._axes.text(0.0, 0.0, '', fontsize=8)
        #     self._cursor_x.set_xdata(x)
        #     self._cursor_text.set_text('period=%.2f' % x)
        #     self._cursor_text.set_position((x, y))
        #     self.draw()

        def set_title(self, title):
            self._title = title

        def set_unit(self, unit):
            if unit != self._unit:
                self._unit = unit
                self.cursor = Cursor(self._axes, track_y=False, text_format="%f " + self._unit, useblit=True)

        def mouse_pick(self, event):
            if not event.inaxes:
                return
            x = event.xdata
            self.set_current_frequency(x)

        def compute_initial_figure(self):
            if self._frequency is not None:
                self._axes.tick_params(axis='both', which='major', labelsize=6)
                self._axes.tick_params(axis='both', which='minor', labelsize=4)
                self._axes.hist(self._frequency)  # , 50, normed=1)
                if self._title and self._unit:
                    self._axes.set_xlabel("%s (%s)" % (self._title, self._unit), fontsize=8)
                    self.figure.suptitle('%s Distribution in Selected Stations' % self._title, fontsize=8)

                self._fig.set_tight_layout(True)

        def set_data(self, frequency):
            self._frequency = frequency
            self._lx = None
            self._current_frequency = None

        def set_current_frequency(self, freq):
            self._current_frequency = freq
            if self._lx is None:
                self._lx = self._axes.axvline(linewidth=2, color="red")
            self._lx.set_xdata(self._current_frequency)
            # if self._fig.canvas.supports_blit:
            #     self._axes.draw_artist(self._lx)
            #     self._fig.canvas.blit(self._axes.bbox)
            # else:
            #     self._fig.canvas.draw_idle()
            self._fig.canvas.draw_idle()

        def update_figure(self):
            # clear figure
            self._axes.cla()
            self.compute_initial_figure()
            self.draw()


class Ellipse(QtGui.QGroupBox):
    """
    ellipse_dict defined for mtpy.imagining.phase_tensor_maps.PlogPhaseTensorMaps
    """
    _colorby = ['phimin', 'phimax', 'skew', 'skew_seg', 'normalized_skew', 'normalized_skew_seg', 'phidet',
                'ellipticity']
    _cmap = ['mt_yl2rd', 'mt_bl2yl2rd', 'mt_wh2bl', 'mt_rd2bl', 'mt_bl2wh2rd', 'mt_seg_bl2wh2rd', 'mt_rd2gr2bl']

    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxEllipse()
        self.ui.setupUi(self)
        # set tooltips for color by and cmap
        for i in range(self.ui.comboBoxColor_by.count()):
            self.ui.comboBoxColor_by.setItemData(i, self.ui.comboBoxColor_by.itemText(i), QtCore.Qt.ToolTipRole)
        for i in range(self.ui.comboBox_cmap.count()):
            self.ui.comboBox_cmap.setItemData(i, self.ui.comboBox_cmap.itemText(i), QtCore.Qt.ToolTipRole)

        self.ui.doubleSpinBox_min.editingFinished.connect(self._increase_max)
        self.ui.doubleSpinBox_max.editingFinished.connect(self._decrease_min)

    def _increase_max(self):
        value = self.ui.doubleSpinBox_min.value()
        if value > self.ui.doubleSpinBox_max.value():
            self.ui.doubleSpinBox_max.setValue(value)

    def _decrease_min(self):
        value = self.ui.doubleSpinBox_max.value()
        if value < self.ui.doubleSpinBox_min.value():
            self.ui.doubleSpinBox_min.setValue(value)

    def get_ellipse_dict(self):
        ellipse_dict = {
            'size': self.ui.doubleSpinBox_size.value(),
            'colorby': self._colorby[self.ui.comboBoxColor_by.currentIndex()],
            'range': (
                self.ui.doubleSpinBox_min.value(), self.ui.doubleSpinBox_max.value(),
                self.ui.doubleSpinBox_step.value()),
            'cmap': self._cmap[self.ui.comboBox_cmap.currentIndex()]
        }
        return ellipse_dict


class FrequencyTolerance(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxTolerance()
        self.ui.setupUi(self)

    def get_tolerance_in_float(self):
        return self.ui.doubleSpinBox.value() / 100.0


class Arrow(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Arrow()
        self.ui.setupUi(self)
        self._simple_color = simple_color
        if not self._simple_color:
            # use all colors available to matplot
            self.ui.comboBox_color_imaginary.clear()
            self.ui.comboBox_color_real.clear()
            cnames = [name for name, hex in COLORS]
            self.ui.comboBox_color_imaginary.addItems(cnames)
            self.ui.comboBox_color_real.addItems(cnames)

    _direction = [0, 1]

    def hide_size(self):
        self.ui.label_size.hide()
        self.ui.doubleSpinBox_size.hide()

    def hide_head_length(self):
        self.ui.label_head_length.hide()
        self.ui.doubleSpinBox_head_length.hide()

    def hide_head_width(self):
        self.ui.label_head_width.hide()
        self.ui.doubleSpinBox_head_width.hide()

    def hide_color_real(self):
        self.ui.label_color_real.hide()
        self.ui.comboBox_color_real.hide()

    def hide_color_imaginary(self):
        self.ui.label_color_imaginary.hide()
        self.ui.comboBox_color_imaginary.hide()

    def hide_threshold(self):
        self.ui.label_threshold.hide()
        self.ui.doubleSpinBox_threshold.hide()

    def hide_direction(self):
        self.ui.label_direction.hide()
        self.ui.comboBox_direction.hide()

    def get_arrow_dict(self):
        if self.ui.groupBox_advanced_options.isChecked():
            arrow_dict = {
                'size': self.ui.doubleSpinBox_size.value(),
                'head_length': self.ui.doubleSpinBox_head_length.value(),
                'head_width': self.ui.doubleSpinBox_head_width.value(),
                'lw': self.ui.doubleSpinBox_line_width.value(),
                'threshold': self.ui.doubleSpinBox_threshold.value(),
                'direction': self._direction[self.ui.comboBox_direction.currentIndex()]
            }
            if self._simple_color:
                arrow_dict['color'] = (SIMPLE_COLORS[self.ui.comboBox_color_real.currentIndex()],
                                       SIMPLE_COLORS[self.ui.comboBox_color_imaginary.currentIndex()])
            else:
                arrow_dict['color'] = (COLORS[self.ui.comboBox_color_real.currentIndex()][1],
                                       COLORS[self.ui.comboBox_color_imaginary.currentIndex()][1])
            return arrow_dict
        else:
            return None

    def get_plot_tipper(self):
        if self.ui.checkBox_real.isChecked() and self.ui.checkBox_imaginary.isChecked():
            return 'yri'
        elif self.ui.checkBox_real.isChecked():
            return 'yr'
        elif self.ui.checkBox_imaginary.isChecked():
            return 'ri'
        else:
            return 'n'


class Padding(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Padding()
        self.ui.setupUi(self)

    def get_x_pad(self):
        return self.ui.doubleSpinBox_x.value()

    def get_y_pad(self):
        return self.ui.doubleSpinBox_y.value()


class Scale(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Scale()
        self.ui.setupUi(self)

    _tscale = ['period', 'freq']
    _mapscale = ['deg', 'm', 'km']

    def get_tscale(self):
        return self._tscale[self.ui.comboBox_time.currentIndex()]

    def hide_mapscale(self):
        self.ui.label_map.hide()
        self.ui.comboBox_map.hide()

    def get_mapscale(self):
        return self._mapscale[self.ui.comboBox_map.currentIndex()]


class Stretch(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Stretch()
        self.ui.setupUi(self)
        self.ui.checkBox_x_range.stateChanged.connect(self._x_range_state_change)
        self.ui.checkBox_y_range.stateChanged.connect(self._y_range_state_change)

    def _x_range_state_change(self, p_int):
        if p_int == 0:
            self.ui.doubleSpinBox_x_min.setEnabled(False)
            self.ui.doubleSpinBox_x_max.setEnabled(False)
        else:
            self.ui.doubleSpinBox_x_min.setEnabled(True)
            self.ui.doubleSpinBox_x_max.setEnabled(True)

    def _y_range_state_change(self, p_int):
        if p_int == 0:
            self.ui.doubleSpinBox_y_min.setEnabled(False)
            self.ui.doubleSpinBox_y_max.setEnabled(False)
        else:
            self.ui.doubleSpinBox_y_min.setEnabled(True)
            self.ui.doubleSpinBox_y_max.setEnabled(True)

    def get_stretch(self):
        return self.ui.doubleSpinBox_x.value(), self.ui.doubleSpinBox_y.value()

    def get_x_limits(self):
        return self.ui.doubleSpinBox_x_min.value(), self.ui.doubleSpinBox_x_max.value()

    def get_y_limits(self):
        return self.ui.doubleSpinBox_y_min.value(), self.ui.doubleSpinBox_y_max.value()


class LineDir(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Linedir()
        self.ui.setupUi(self)

    def get_linedir(self):
        if self.ui.radioButton_ns.isChecked():
            return 'ns'
        elif self.ui.radioButton_ew.isChecked():
            return 'ew'
        else:
            return None


class FrequencyIndex(QtGui.QGroupBox):
    _unit_period = 'second'
    _unit_frequency = 'Hz'
    _title_period = 'Period'
    _title_frequency = 'Frequency'

    def __init__(self, parent, use_period=False):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Frequency_Period_Index()
        self.ui.setupUi(self)
        self._mt_objs = None
        self.use_period = use_period
        self.set_use_period(self.use_period)

    def set_use_period(self, use_period=False):
        self.use_period = use_period
        if self.use_period:
            title = '%s (%s)' % (self._title_period, self._unit_period)
        else:
            title = '%s (%s)' % (self._title_frequency, self._unit_frequency)
        self.setTitle(title)
        self._update_frequency()

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._update_frequency()

    def _update_frequency(self):
        if self._mt_objs:
            self.ui.listWidget_frequency_period.clear()

            all_freqs = self._mt_objs[0].Z.freq
            # print all_freqs

            if all([all_freqs.shape == mt_obj.Z.freq.shape and np.allclose(all_freqs, mt_obj.Z.freq) for mt_obj in
                    self._mt_objs[1:]]):

                if self.use_period:
                    all_freqs = 1.0 / np.array(all_freqs)

                self.ui.listWidget_frequency_period.addItems(
                    ["%.5f %s" % (value, self._unit_period if self.use_period else self._unit_frequency) for value in
                     all_freqs])
                self.ui.listWidget_frequency_period.setCurrentRow(0)  # select the first row by default
                self.ui.listWidget_frequency_period.setEnabled(True)
            else:
                self.ui.listWidget_frequency_period.addItem("ERROR: frequency lists from stations are not identical")
                self.ui.listWidget_frequency_period.setEnabled(False)

    def get_index_list(self):
        return sorted([index.row() for index in self.ui.listWidget_frequency_period.selectedIndexes()], reverse=False)


class UniqueFrequencies(FrequencyIndex):
    def __init__(self, parent, use_period=False):
        FrequencyIndex.__init__(self, parent, use_period)
        self.unique_freqs = None

    def _update_frequency(self):
        if self._mt_objs:
            self.ui.listWidget_frequency_period.clear()

            unique_freqs = set()
            for mt_obj in self._mt_objs:
                unique_freqs.update(mt_obj.Z.freq)
            self.unique_freqs = np.array(list(unique_freqs))
            if self.use_period:
                self.unique_freqs = 1.0 / self.unique_freqs
            self.unique_freqs.sort()

            self.ui.listWidget_frequency_period.addItems(
                [
                    "%.5f %s" % (value, self._unit_period if self.use_period else self._unit_frequency)
                    for value in self.unique_freqs
                ]
            )

    def get_index_list(self):
        """
        should net be used
        :return:
        """
        pass

    def get_frequency_list(self):
        return sorted([self.unique_freqs[index.row()] for index in self.ui.listWidget_frequency_period.selectedIndexes()], reverse=False)


class StationSelection(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Station_Select()
        self.ui.setupUi(self)
        self.mt_objs = None

        self.ui.comboBox_station.currentIndexChanged.connect(self._current_station_changed)

    def _current_station_changed(self):
        self.station_changed.emit()

    station_changed = pyqtSignal()

    def set_data(self, mt_objs):
        self.ui.comboBox_station.clear()
        self.mt_objs = []
        for mt_obj in mt_objs:
            self.mt_objs.append(mt_obj)
            self.ui.comboBox_station.addItem(mt_obj.station)

    def get_station(self):
        index = self.ui.comboBox_station.currentIndex()
        return self.mt_objs[index]


class Rotation(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Rotation()
        self.ui.setupUi(self)
        self.ui.dial_rotation.valueChanged.connect(self._dial_value_changed)
        self.ui.doubleSpinBox_rotation.valueChanged.connect(self._text_value_changed)

    def _dial_value_changed(self, p_int):
        degree = (p_int - 180) % 360
        self.ui.doubleSpinBox_rotation.setValue(degree)

    def _text_value_changed(self):
        degree = (int(self.ui.doubleSpinBox_rotation.value()) + 180) % 360
        if degree != self.ui.dial_rotation.value():
            self.ui.dial_rotation.setValue(degree)

    def get_rotation_in_degree(self):
        return self.ui.doubleSpinBox_rotation.value()


class PlotControlMTResponse(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_plot_control_mt_response()
        self.ui.setupUi(self)

    def get_plot_num(self):
        if self.ui.radioButton_1.isChecked():
            return 1
        elif self.ui.radioButton_2.isChecked():
            return 2
        elif self.ui.radioButton_3.isChecked():
            return 3
        else:
            return 0  # should never reach here

    def get_strike(self):
        strike = "y"
        if self.ui.checkBox_strike_t.isChecked():
            strike += 't'
        if self.ui.checkBox_strike_p.isChecked():
            strike += 'p'
        if self.ui.checkBox_strike_i.isChecked():
            strike += 'i'
        if len(strike) > 1:
            return strike
        else:
            return 'n'

    def get_skew(self):
        if self.ui.radioButton_skew_y.isChecked():
            return 'y'
        else:
            return 'n'

    def get_ellipses(self):
        if self.ui.radioButton_ellipses_y.isChecked():
            return 'y'
        else:
            return 'n'

    def get_style(self):
        if self.ui.radioButton_compare.isChecked():
            return "compare"
        else:
            return "all"

    def hide_plot_style(self):
        self.ui.groupBox_plot_style.hide()


class MeshGrid(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_mash_grid()
        self.ui.setupUi(self)

        # connect signal
        self.ui.radioButton_imshow.toggled.connect(self._imshow_toggled)

    _grid_types = ['imshow', 'pcolormesh']
    _interpolation_methods = ['none', 'nearest', 'bilinear', 'bicubic',
                              'spline16', 'spline36', 'hanning', 'hamming',
                              'hermite', 'kaiser', 'quadric', 'catrom',
                              'gaussian', 'bessel', 'mitchell', 'sinc',
                              'lanczos']

    def _imshow_toggled(self, checked):
        if checked:
            self.ui.groupBox_interpolation_method.setHidden(False)
        else:
            self.ui.groupBox_interpolation_method.setHidden(True)

    def get_grid_type(self):
        if self.ui.radioButton_imshow.isChecked():
            return self._grid_types[0]
        elif self.ui.radioButton_pcolormesh.isChecked():
            return self._grid_types[1]
        else:
            return None  # should never reach here

    def get_interpolation_method(self):
        if not self.ui.groupBox_interpolation_method.isHidden():
            return self._interpolation_methods[
                self.ui.comboBox_interpolation_method.currentIndex()
            ]
        else:
            return None


class PlotControlResistivityPhasePseudoSection(QtGui.QGroupBox):
    """
    plot settings for resistivity phase pseudo section that cannot be standardized
    """
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_plot_control_resistivity_phase_pseudo_section()
        # set up gui
        self.ui.setupUi(self)

        self.ui.comboBox_res_cmap.addItems(cmapdict.keys())
        self.ui.comboBox_res_cmap.setCurrentIndex(self.ui.comboBox_res_cmap.findText('mt_rd2gr2bl'))
        self.ui.comboBox_phase_cmap.addItems(cmapdict.keys())
        self.ui.comboBox_phase_cmap.setCurrentIndex(self.ui.comboBox_phase_cmap.findText('mt_bl2gr2rd'))

        self.ui.doubleSpinBox_res_max.setMaximum(np.inf)
        self.ui.doubleSpinBox_res_min.setMaximum(np.inf)
        self.ui.doubleSpinBox_period_min.setMaximum(np.inf)
        self.ui.doubleSpinBox_period_max.setMaximum(np.inf)

        # connect signals
        self.ui.doubleSpinBox_period_min.valueChanged.connect(self._period_min_changed)
        self.ui.doubleSpinBox_period_max.valueChanged.connect(self._period_max_changed)
        self.ui.doubleSpinBox_phase_min.valueChanged.connect(self._phase_min_changed)
        self.ui.doubleSpinBox_phase_max.valueChanged.connect(self._phase_max_changed)
        self.ui.doubleSpinBox_res_min.valueChanged.connect(self._res_min_changed)
        self.ui.doubleSpinBox_res_max.valueChanged.connect(self._res_max_changed)

        self.ui.checkBox_phase.stateChanged.connect(self._phase_state_changed)
        self.ui.checkBox_period.stateChanged.connect(self._period_state_changed)
        self.ui.checkBox_resistivity.stateChanged.connect(self._res_state_changed)

    def _res_min_changed(self, value):
        if value > self.ui.doubleSpinBox_res_max.value():
            self.ui.doubleSpinBox_res_max.blockSignals(True)
            self.ui.doubleSpinBox_res_max.setValue(value)
            self.ui.doubleSpinBox_res_max.blockSignals(False)

    def _res_max_changed(self, value):
        if value < self.ui.doubleSpinBox_res_min.value():
            self.ui.doubleSpinBox_res_min.blockSignals(True)
            self.ui.doubleSpinBox_res_min.setValue(value)
            self.ui.doubleSpinBox_res_min.blockSignals(False)

    def _res_state_changed(self, p_int):
        state = bool(p_int != 0)
        self.ui.doubleSpinBox_res_min.setEnabled(state)
        self.ui.doubleSpinBox_res_max.setEnabled(state)

    def _phase_state_changed(self, p_int):
        state = bool(p_int != 0)
        self.ui.doubleSpinBox_phase_min.setEnabled(state)
        self.ui.doubleSpinBox_phase_max.setEnabled(state)

    def _period_state_changed(self, p_int):
        state = bool(p_int != 0)
        self.ui.doubleSpinBox_period_min.setEnabled(state)
        self.ui.doubleSpinBox_period_max.setEnabled(state)

    def _period_min_changed(self, value):
        if value > self.ui.doubleSpinBox_period_max.value():
            self.ui.doubleSpinBox_period_max.blockSignals(True)
            self.ui.doubleSpinBox_period_max.setValue(value)
            self.ui.doubleSpinBox_period_max.blockSignals(False)

    def _period_max_changed(self, value):
        if value < self.ui.doubleSpinBox_period_min.value():
            self.ui.doubleSpinBox_period_min.blockSignals(True)
            self.ui.doubleSpinBox_period_min.setValue(value)
            self.ui.doubleSpinBox_period_min.blockSignals(False)

    def _phase_min_changed(self, value):
        if value > self.ui.doubleSpinBox_phase_max.value():
            self.ui.doubleSpinBox_phase_max.blockSignals(True)
            self.ui.doubleSpinBox_phase_max.setValue(value)
            self.ui.doubleSpinBox_phase_max.blockSignals(False)

    def _phase_max_changed(self, value):
        if value < self.ui.doubleSpinBox_phase_min.value():
            self.ui.doubleSpinBox_phase_min.blockSignals(True)
            self.ui.doubleSpinBox_phase_min.setValue(value)
            self.ui.doubleSpinBox_phase_min.blockSignals(False)

    def get_phase_limit(self):
        if self.ui.checkBox_phase.isChecked():
            return self.ui.doubleSpinBox_phase_min.value(), self.ui.doubleSpinBox_phase_max.value()
        else:
            return None

    def get_period_limit(self):
        if self.ui.checkBox_period.isChecked():
            return self.ui.doubleSpinBox_period_min, self.ui.doubleSpinBox_period_max
        else:
            return None

    def get_resistivity_limits(self):
        if self.ui.checkBox_resistivity.isChecked():
            return self.ui.doubleSpinBox_res_min, self.ui.doubleSpinBox_res_max
        else:
            return None

    def get_plot_xx(self):
        return 'y' if self.ui.checkBox_zxx.isChecked() else 'n'

    def get_plot_xy(self):
        return 'y' if self.ui.checkBox_zxy.isChecked() else 'n'

    def get_plot_yx(self):
        return 'y' if self.ui.checkBox_zyx.isChecked() else 'n'

    def get_plot_yy(self):
        return 'y' if self.ui.checkBox_zyy.isChecked() else 'n'

    def get_res_cmap(self):
        return cmapdict[str(self.ui.comboBox_res_cmap.currentText())]

    def get_phase_cmap(self):
        return cmapdict[str(self.ui.comboBox_phase_cmap.currentText())]

    def get_tickspace(self):
        return self.ui.spinBox_tickspace.value()