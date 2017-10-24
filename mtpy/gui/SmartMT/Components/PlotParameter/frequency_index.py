# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
import numpy as np
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_index import Ui_GroupBox_Frequency_Period_Index


class FrequencyIndex(QGroupBox):
    _unit_period = 'second'
    _unit_frequency = 'Hz'
    _title_period = 'Period'
    _title_frequency = 'Frequency'

    def __init__(self, parent, use_period=False):
        QGroupBox.__init__(self, parent)
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
