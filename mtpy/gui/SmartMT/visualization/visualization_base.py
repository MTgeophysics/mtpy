# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""

import abc

from PyQt4 import QtGui, QtCore


class VisualizationBase(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self._groupBoxParameters = None

    def get_QGroupBox_parameters(self):
        return self._groupBoxParameters

    @staticmethod
    @abc.abstractmethod
    def get_plot_name():
        return VisualizationBase.__name__

    @staticmethod
    @abc.abstractmethod
    def get_plot_description():
        return VisualizationBase.__name__
