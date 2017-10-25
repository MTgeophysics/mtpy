# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
import numpy as np
from qtpy.QtCore import Signal
from qtpy.QtWidgets import QGroupBox, QStyledItemDelegate
from qtpy.QtGui import QStandardItemModel, QStandardItem
from qtpy import QtCore

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_select import Ui_GroupBox_frequency_select
from mtpy.gui.SmartMT.ui_asset.groupbox_select_periods_from_files import Ui_GroupBox_select_from_files
from mtpy.gui.SmartMT.utils.matplotlib_utils import gen_hist_bins


class FrequencySelection(QGroupBox):
    """
    frequency selection
    """

    def __init__(self, parent, show_period=True, show_frequency=True, allow_range_select=True,
                 select_multiple=True):
        QGroupBox.__init__(self, parent)
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
        self.model_selected = QStandardItemModel()
        self.ui.listView_selected.setModel(self.model_selected)
        self.frequency_delegate = FrequencySelection.FrequencyDelegate(self.ui.listView_selected)
        self.ui.listView_selected.setItemDelegate(self.frequency_delegate)

        self.histogram = FrequencySelection.Histogram(self, allow_range_select=self._allow_range)
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
        frequencies = [self.model_selected.item(index).data(QtCore.Qt.DisplayRole)
                       for index in range(self.model_selected.rowCount())]
        if self._allow_range:
            frequencies = [(freq[0], freq[1]) if isinstance(freq, tuple) else freq
                           for freq in frequencies]
        else:
            frequencies = [freq[3] if isinstance(freq, tuple) else freq
                           for freq in frequencies
                           if (isinstance(freq, tuple) and len(freq) == 5)
                           or isinstance(freq, float)]
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
        for item in [self.model_selected.item(index.row())
                     for index in self.ui.listView_selected.selectedIndexes()]:
            x = item.data(QtCore.Qt.DisplayRole)
            self.model_selected.removeRow(self.model_selected.indexFromItem(item).row())
            self.histogram.remove_marker(x)

    def _frequency_selected(self, x):
        if not self._select_multiple:
            self.histogram.clear_all_drawing()
            self.model_selected.clear()
        for item in [self.model_selected.item(index) for index in range(self.model_selected.rowCount())]:
            value = item.data(QtCore.Qt.DisplayRole)
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
                    uniques = self._unique_frequencies \
                        if self.ui.radioButton_frequency.isChecked() \
                        else self._unique_periods
                    num = len([freq for freq in uniques if mi <= freq <= ma])  # num of existing freqs in the new interval
                    x = (mi, ma, num)
                    # remove old interval
                    self.model_selected.removeRow(self.model_selected.indexFromItem(item).row())
                    self.histogram.remove_marker(value)
            else:
                prec = self.frequency_delegate.prec
                while np.all(np.isclose(value, x, pow(.1, prec))):
                    prec += 1
                self.frequency_delegate.prec = prec
        new_item = FrequencySelection.FrequencyItem()
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
                self._periods if self.ui.radioButton_period.isChecked()
                else self._frequencies,
                self._unique_periods if self.ui.radioButton_period.isChecked()
                else self._unique_frequencies
            )
            self.frequency_delegate.freqs = self._unique_periods \
                if self.ui.radioButton_period.isChecked() \
                else self._unique_frequencies
            self.histogram.update_figure()

    class FrequencyItem(QStandardItem):
        def __lt__(self, other):
            value = self.data(QtCore.Qt.DisplayRole)
            other_value = other.data(QtCore.Qt.DisplayRole)
            if isinstance(value, tuple):
                value = value[0]
            if isinstance(other_value, tuple):
                other_value = other_value[0]
            return value < other_value

    class FrequencyDelegate(QStyledItemDelegate):
        _prec = 5  # decimal places

        def get_prec(self):
            return self._prec

        def set_prec(self, prec):
            self._prec = prec

        prec = property(get_prec, set_prec)

        def displayText(self, value, locale):
            if isinstance(value, float):
                return '{:.{prec}f}'.format(value, prec=self._prec)
            elif isinstance(value, tuple) and len(value) == 3:  # (min, max, num)
                return '{}{}, {}{} ({num} selected)'.format(
                    '(' if value[0] == -np.inf else '[',
                    '{:.{prec}f}'.format(value[0], prec=self._prec),
                    '{:.{prec}f}'.format(value[1], prec=self._prec),
                    ')' if value[1] == np.inf else ']',
                    num=value[2]
                )
            elif len(value) == 5:  # (min, max, num, freq, tol)
                return u'{:.{prec}f} ±{tol}% ({num} selected)'.format(
                    value[3], prec=self._prec, tol=value[4], num=value[2])
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

        frequency_selected = Signal(float)
        frequency_range_selected = Signal(tuple)

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
                                len([freq for freq in self._unique_frequencies
                                     if self._press <= freq <= x])
                            )
                        )
                    elif self._press > x:
                        self.frequency_range_selected.emit(
                            (
                                x,
                                self._press,
                                len([freq for freq in self._unique_frequencies
                                     if x <= freq <= self._press])
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
                            len([freq for freq in self._unique_frequencies
                                 if min <= freq <= max]),
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
                bins = gen_hist_bins(self._unique_frequencies)
                self._axes.hist(self._frequencies, bins=bins)  # , 50, normed=1)
                if self._y_log_scale:
                    self._axes.set_yscale('log', nonposy='clip')
                if self._x_log_scale:
                    self._axes.set_xscale('log', nonposx='clip')
                if self._show_existing:
                    for freq in self._unique_frequencies:
                        self._axes.axvline(freq, linewidth=1, color='black', alpha=0.2)
            if self._title and self._unit:
                self._axes.set_xlabel("%s (%s)" % (self._title, self._unit), fontsize=8)
                self.figure.suptitle('%s Distribution in Selected Stations' %
                                     self._title, fontsize=8)

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


class FrequencySelectionFromFile(QGroupBox):
    """
    select frequencies/periods from the selected edi files
    """
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self.model_stations = QStandardItemModel()
        self.model_selected_frequencies = QStandardItemModel()

        # setup ui
        self.ui.setupUi(self)
        self.ui.listView_selected.setModel(self.model_selected_frequencies)
        self.ui.listWidget_stations.setModel(self.model_stations)

        # connect signals

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self.model_selected_frequencies.clear()

        self._update_stations()

        self.data_changed.emit()

    def _update_stations(self):
        if self._mt_objs is not None:
            self.model_stations.clear()
            for mt_obj in self._mt_objs:
                new_item = QStandardItem()
                new_item.setData(mt_obj.station, QtCore.Qt.DisplayRole)
                new_item.setData(mt_obj.fn, QtCore.Qt.ToolTipRole)
                self.model_stations.appendRow(new_item)




