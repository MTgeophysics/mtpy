# -*- coding: utf-8 -*-
import numpy as np
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import pyqtSignal

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.gui.plot_parameter import COLORS, SIMPLE_COLORS
from mtpy.gui.SmartMT.ui_asset.groupbox_arrow import Ui_GroupBox_Arrow
from mtpy.gui.SmartMT.ui_asset.groupbox_ellipse import Ui_GroupBoxEllipse
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_index import Ui_GroupBox_Frequency_Period_Index
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_select import Ui_GroupBox_frequency_select
from mtpy.gui.SmartMT.ui_asset.groupbox_linedir import Ui_GroupBox_Linedir
from mtpy.gui.SmartMT.ui_asset.groupbox_mesh_grid import Ui_GroupBox_mash_grid
from mtpy.gui.SmartMT.ui_asset.groupbox_padding import Ui_GroupBox_Padding
from mtpy.gui.SmartMT.ui_asset.groupbox_rotation import Ui_GroupBox_Rotation
from mtpy.gui.SmartMT.ui_asset.groupbox_scale import Ui_GroupBox_Scale
from mtpy.gui.SmartMT.ui_asset.groupbox_station_select import Ui_GroupBox_Station_Select
from mtpy.gui.SmartMT.ui_asset.groupbox_stretch import Ui_GroupBox_Stretch
from mtpy.gui.SmartMT.ui_asset.groupbox_tolerance import Ui_GroupBoxTolerance
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_multiple import Ui_groupBoxZ_Component_Multiple
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_single import Ui_groupBoxZ_Component_Single
from mtpy.gui.SmartMT.ui_asset.groupbox_z_unit import Ui_GroupBox_z_unit


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


class FrequencySelect(QtGui.QGroupBox):
    """
    frequency selection
    """

    def __init__(self, parent, show_period=True, show_frequency=True, allow_range_select=True,
                 select_multiple=True):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self._unique_periods = None
        self._unique_frequencies = None
        self._periods = None
        self._frequencies = None
        self._allow_range = allow_range_select
        self._select_multiple = select_multiple
        self.ui = Ui_GroupBox_frequency_select()
        self.ui.setupUi(self)

        self.ui.label_place_holder.hide()
        self.model_selected = QtGui.QStandardItemModel()
        self.ui.listView_selected.setModel(self.model_selected)
        self.frequency_delegate = FrequencySelect.FrequencyDelegate(self.ui.listView_selected)
        self.ui.listView_selected.setItemDelegate(self.frequency_delegate)

        self.histogram = FrequencySelect.Histogram(self, allow_range_select=self._allow_range)
        self.histogram.set_unit(self._units[0])
        self.histogram.set_tol(self.ui.doubleSpinBox_tolerance.value())
        self.histogram.frequency_selected.connect(self._frequency_selected)
        self.histogram.frequency_range_selected.connect(self._frequency_selected)
        self.ui.widget_histgram.layout().addWidget(self.histogram)

        self.ui.radioButton_period.setChecked(show_period)
        self.ui.radioButton_frequency.setChecked(show_frequency)
        self.ui.doubleSpinBox_tolerance.setHidden(not self._allow_range)
        self.ui.checkBox_existing_only.setChecked(not self._allow_range)
        self.ui.checkBox_existing_only.setHidden(not self._allow_range)
        self.ui.label_tolerance.setHidden(not self._allow_range)
        self.ui.radioButton_period.setHidden(not (show_period and show_frequency))
        self.ui.radioButton_frequency.setHidden(not (show_period and show_frequency))
        if self.ui.radioButton_frequency.isHidden():
            self.setTitle(self._type[1])
        elif self.ui.radioButton_period.isHidden():
            self.setTitle(self._type[0])

        self.ui.radioButton_frequency.toggled.connect(self._frequency_toggled)
        self.ui.checkBox_existing_only.toggled.connect(self.histogram.select_existing)
        self.ui.checkBox_existing_only.toggled.connect(self.model_selected.clear)
        self.ui.checkBox_show_existing.toggled.connect(self.histogram.show_existing)
        self.ui.checkBox_x_log_scale.toggled.connect(self.histogram.set_x_log_scale)
        self.ui.checkBox_y_log_scale.toggled.connect(self.histogram.set_y_log_scale)
        self.ui.pushButton_clear.clicked.connect(self._clear_all)
        self.ui.pushButton_delete.clicked.connect(self._delete_selected)
        self.ui.doubleSpinBox_tolerance.valueChanged.connect(self.histogram.set_tol)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._unique_frequencies = None
        self._unique_periods = None
        self._update_frequency()

    def get_frequencies(self):
        frequencies = [self.model_selected.item(index).data(QtCore.Qt.DisplayRole).toPyObject()
                       for index in range(self.model_selected.rowCount())]
        if self._allow_range:
            frequencies = [(freq[0], freq[1]) if isinstance(freq, tuple) else freq for freq in frequencies]
        else:
            frequencies = [freq[3] if isinstance(freq, tuple) else freq for freq in frequencies if (isinstance(freq, tuple) and len(freq) == 5) or isinstance(freq, float)]
        # print frequencies
        if self._select_multiple:
            return frequencies
        else:
            return frequencies[0] if frequencies else None

    _units = ['Hz', 's']

    _type = ['Frequency', 'Period']

    def _clear_all(self):
        self.model_selected.clear()
        self.histogram.clear_all_drawing()

    def _delete_selected(self):
        for item in [self.model_selected.item(index.row()) for index in self.ui.listView_selected.selectedIndexes()]:
            x = item.data(QtCore.Qt.DisplayRole).toPyObject()
            self.model_selected.removeRow(self.model_selected.indexFromItem(item).row())
            self.histogram.remove_marker(x)

    def _frequency_selected(self, x):
        if not self._select_multiple:
            self.histogram.clear_all_drawing()
            self.model_selected.clear()
        for item in [self.model_selected.item(index) for index in range(self.model_selected.rowCount())]:
            value = item.data(QtCore.Qt.DisplayRole).toPyObject()
            if value == x:
                return
            elif isinstance(value, tuple) and isinstance(x, float) and value[0] <= x <= value[1]:
                return  # x already in interval
            elif isinstance(x, tuple) and isinstance(value, float) and x[0] <= value <= x[1]:
                # existing value in new interval
                self.model_selected.removeRow(self.model_selected.indexFromItem(item).row())
                self.histogram.remove_marker(value)
            elif isinstance(x, tuple) and isinstance(value, tuple):
                if min(x[1], value[1]) - max(x[0], value[0]) >= 0:
                    # there is intersection between intervals, so marge them
                    mi = min(x[0], value[0])
                    ma = max(x[1], value[1])
                    uniques = self._unique_frequencies if self.ui.radioButton_frequency.isChecked() else self._unique_periods
                    num = len([freq for freq in uniques if mi <= freq <= ma ])  # num of existing freqs in the new interval
                    x = (mi, ma, num)
                    # remove old interval
                    self.model_selected.removeRow(self.model_selected.indexFromItem(item).row())
                    self.histogram.remove_marker(value)
            else:
                prec = self.frequency_delegate.prec
                while np.all(np.isclose(value, x, pow(.1, prec))):
                    prec += 1
                self.frequency_delegate.prec = prec
        new_item = FrequencySelect.FrequencyItem()
        new_item.setData(x, QtCore.Qt.DisplayRole)
        # update graphic
        if isinstance(x, float):
            self.histogram.add_marker(x)
            # new_item.setData(x, QtCore.Qt.UserRole)
        elif isinstance(x, tuple):
            self.histogram.add_marker(x)
            # new_item.setData(x[0], QtCore.Qt.UserRole)
        # update model
        self.model_selected.appendRow(new_item)
        self.model_selected.sort(0)

    def show_period(self):
        self.ui.radioButton_period.setChecked(True)

    def show_frequency(self):
        self.ui.radioButton_frequency.setChecked(True)

    def _frequency_toggled(self, is_checked):
        self.histogram.set_unit(self._units[0] if is_checked else self._units[1])
        self._update_frequency()

    def _update_frequency(self):
        self.model_selected.clear()
        if self._mt_objs is not None:
            if self._unique_frequencies is None:
                self._frequencies = [freq for mt_obj in self._mt_objs for freq in list(mt_obj.Z.freq)]
                all_unique = set(self._frequencies)
                self._unique_frequencies = sorted(list(all_unique))
            if self.ui.radioButton_period.isChecked() and self._unique_periods is None:
                self._periods = 1. / np.array(self._frequencies)
                all_unique = set(self._periods)
                self._unique_periods = sorted(list(all_unique))
            self.histogram.set_data(
                self._periods if self.ui.radioButton_period.isChecked() else self._frequencies,
                self._unique_periods if self.ui.radioButton_period.isChecked() else self._unique_frequencies
            )
            self.frequency_delegate.freqs = self._unique_periods if self.ui.radioButton_period.isChecked() else self._unique_frequencies
            self.histogram.update_figure()

    class FrequencyItem(QtGui.QStandardItem):
        def __lt__(self, other):
            value = self.data(QtCore.Qt.DisplayRole).toPyObject()
            other_value = other.data(QtCore.Qt.DisplayRole).toPyObject()
            if isinstance(value, tuple):
                value = value[0]
            if isinstance(other_value, tuple):
                other_value = other_value[0]
            return value < other_value

    class FrequencyDelegate(QtGui.QStyledItemDelegate):
        _prec = 5  # decimal places

        def get_prec(self):
            return self._prec

        def set_prec(self, prec):
            self._prec = prec

        prec = property(get_prec, set_prec)

        def displayText(self, value, locale):
            py_obj = value.toPyObject()
            if isinstance(py_obj, float):
                return '{:.{prec}f}'.format(py_obj, prec=self._prec)
            elif isinstance(py_obj, tuple) and len(py_obj) == 3:  # (min, max, num)
                return '{}{}, {}{} ({num} selected)'.format(
                    '(' if py_obj[0] == -np.inf else '[',
                    '{:.{prec}f}'.format(py_obj[0], prec=self._prec),
                    '{:.{prec}f}'.format(py_obj[1], prec=self._prec),
                    ')' if py_obj[1] == np.inf else ']',
                    num=py_obj[2]
                )
            elif len(py_obj) == 5:  # (min, max, num, freq, tol)
                return u'{:.{prec}f} ±{tol}% ({num} selected)'.format(
                    py_obj[3], prec=self._prec, tol=py_obj[4], num=py_obj[2])
            # elif isinstance(py_obj, set):
            #     return '{{}}'.format(','.join(['{:.{prec}f}'.format(f, prec=self._prec) for f in py_obj if isinstance(f, float)]))
            return value

    class Histogram(MPLCanvas):
        def __init__(self, parent, y_log_scale=False, x_log_scale=False, allow_range_select=True):
            self._frequencies = None
            self._unique_frequencies = None
            self._title = None
            self._unit = None
            self._press = None
            self._tol = None
            MPLCanvas.__init__(self, parent, 5, 1.5)
            self._lx = {}
            self._cursor = None
            self._select_existing_only = False
            self._show_existing = False
            self._x_log_scale = x_log_scale
            self._y_log_scale = y_log_scale
            self._select_range = allow_range_select

            if self._select_range:
                self.mpl_connect('button_press_event', self.on_press)
            self.mpl_connect('button_release_event', self.on_release)

        def add_marker(self, x):
            if isinstance(x, float):
                lx = self._lx.setdefault(x, self._draw_v_line(x))
                # self._axes.draw_artist(lx)
                self.draw_idle()
            elif isinstance(x, tuple):
                if len(x) == 3:
                    lx = self._lx.setdefault(x, self._fill_v_area(x[0], x[1]))
                elif len(x) == 5:
                    lx = self._lx.setdefault(x, (
                        self._draw_v_line(x[3]),
                        self._fill_v_area(x[0], x[1])
                    ))
            else:
                raise NotImplemented
            self.draw_idle()

        def remove_marker(self, x):
            if x in self._lx:
                marker = self._lx[x]
                if isinstance(marker, tuple):
                    for m in marker:
                        m.remove()
                else:
                    marker.remove()
                self.draw_idle()
                del self._lx[x]

        def clear_all_drawing(self):
            for key in self._lx.keys():
                marker = self._lx[key]
                if isinstance(marker, tuple):
                    for m in marker:
                        m.remove()
                else:
                    marker.remove()
            self._lx.clear()
            self.draw_idle()

        def set_unit(self, unit):
            if unit != self._unit:
                self._unit = unit
                self._cursor = Cursor(self._axes,
                                      track_y=False,
                                      show_drag=self._select_range,
                                      text_format="%f" + self._unit,
                                      useblit=True)

        def select_existing(self, select_existing):
            self._select_existing_only = select_existing
            self.clear_all_drawing()

        def set_tol(self, tol):
            self._tol = tol

        def show_existing(self, show_existing):
            self._show_existing = show_existing
            self.update_figure()

        def set_data(self, frequencies, unique_frequencies=None):
            self._frequencies = frequencies
            if unique_frequencies is not None:
                self._unique_frequencies = unique_frequencies
            else:
                self._unique_frequencies = sorted(list(set(frequencies)))
            self._lx.clear()

        def set_y_log_scale(self, ischecked):
            self._y_log_scale = ischecked
            self.update_figure()

        def set_x_log_scale(self, isChecked):
            self._x_log_scale = isChecked
            self.update_figure()

        frequency_selected = pyqtSignal(float)
        frequency_range_selected = pyqtSignal(tuple)

        def _get_valid_cursor_loc(self, event):
            if not event.inaxes:
                pos = self._axes.get_position()
                if self.height() * pos.y0 < event.y < self.height() * pos.y1:
                    x = -np.inf if event.x < self.width() * pos.x0 else np.inf
                else:
                    x = None
            else:
                x = event.xdata
            return x

        def on_press(self, event):
            self._press = self._get_valid_cursor_loc(event)

        def on_release(self, event):
            x = self._get_valid_cursor_loc(event)
            if x:
                if self._press and self._press != x:  # emit (min, max, num)
                    if self._press < x:
                        self.frequency_range_selected.emit(
                            (
                                self._press,
                                x,
                                len([freq for freq in self._unique_frequencies if self._press <= freq <= x])
                            )
                        )
                    elif self._press > x:
                        self.frequency_range_selected.emit(
                            (
                                x,
                                self._press,
                                len([freq for freq in self._unique_frequencies if x <= freq <= self._press])
                            )
                        )
                elif not self._select_range or self._select_existing_only:
                    x = self._find_closest(x)
                    self.frequency_selected.emit(x)
                else:  # emit (min, max, num, freq, tol)
                    tol = x * self._tol/100.
                    min = x - tol
                    max = x + tol
                    self.frequency_range_selected.emit(
                        (
                            min,
                            max,
                            len([freq for freq in self._unique_frequencies if min <= freq <= max]),
                            x,
                            self._tol
                        )
                    )
            self._press = None

        def _find_closest(self, x):
            return min(self._frequencies, key=lambda freq: abs(freq - x))

        def compute_initial_figure(self):
            self._axes.tick_params(axis='both', which='major', labelsize=6)
            self._axes.tick_params(axis='both', which='minor', labelsize=4)
            if self._frequencies is not None:
                bins = np.array(self._unique_frequencies)  # get center of bins
                bins -= np.diff(np.r_[0, bins]).dot(.5)  # shift left
                self._axes.hist(self._frequencies, bins=bins)  # , 50, normed=1)
                if self._y_log_scale:
                    self._axes.set_yscale('log')
                if self._x_log_scale:
                    self._axes.set_xscale('log')
                if self._show_existing:
                    for freq in self._unique_frequencies:
                        self._axes.axvline(freq, linewidth=1, color='black', alpha=0.2)
            if self._title and self._unit:
                self._axes.set_xlabel("%s (%s)" % (self._title, self._unit), fontsize=8)
                self.figure.suptitle('%s Distribution in Selected Stations' % self._title, fontsize=8)

            self._fig.set_tight_layout(True)

        def update_figure(self):
            self._axes.cla()
            self.compute_initial_figure()
            for key in self._lx.keys():
                if isinstance(key, float):
                    self._lx[key] = self._draw_v_line(key)
                elif isinstance(key, tuple):
                    if len(key) == 3:
                        self._lx[key] = self._fill_v_area(key[0], key[1])
                    elif len(key) == 5:
                        self._lx[key] = (self._draw_v_line(key[3]), self._fill_v_area(key[0], key[1]))
            self.draw()

        def _draw_v_line(self, x):
            if x == -np.inf:
                x = self._axes.get_xlim()[0]
            if x == np.inf:
                x = self._axes.get_xlim()[1]
            return self._axes.axvline(x=x, linewidth=1, color="red")

        def _fill_v_area(self, x1, x2):
            if x1 == -np.inf:
                x1 = self._axes.get_xlim()[0]
            if x2 == np.inf:
                x2 = self._axes.get_xlim()[1]
            return self._axes.axvspan(x1, x2, alpha=0.5, color='red')


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
        return sorted(
            [self.unique_freqs[index.row()] for index in self.ui.listWidget_frequency_period.selectedIndexes()],
            reverse=False)


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


class ZUnit(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_z_unit()
        self.ui.setupUi(self)

    def get_unit(self):
        return 'km' if self.ui.radioButton_km.isChecked() else 'm'
