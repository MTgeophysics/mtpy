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

    @abc.abstractmethod
    def __init__(self, parent):
        self._parent = parent
        self._mt_objs = None

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self.update_ui()

    @abc.abstractmethod
    def update_ui(self):
        pass

    @abc.abstractproperty
    def parameter_ui(self):
        return "should not see this"

    @staticmethod
    @abc.abstractmethod
    def plot_name():
        return VisualizationBase.__name__

    @staticmethod
    @abc.abstractmethod
    def plot_description():
        return VisualizationBase.__name__
