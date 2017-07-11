# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""

import abc

from PyQt4 import QtGui, QtCore

from mtpy.utils.mtpylog import MtPyLog
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

class VisualizationBase(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, parent):
        self._parent = parent
        self._mt_objs = None
        self._fig = None
        self._logger = MtPyLog().get_mtpy_logger(__name__)

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

    @abc.abstractmethod
    def plot(self):
        pass

    def show_figure(self):
        self.plot()
        if self._fig:
            # self._fig.show()
            widget = QtGui.QWidget()
            layout = QtGui.QVBoxLayout()
            canvas = FigureCanvas(self._fig.get_figure())
            toolbar = NavigationToolbar(canvas, widget)
            layout.addWidget(toolbar)
            layout.addWidget(canvas)
            widget.setLayout(layout)
            self._parent._parent.create_subwindow(widget, self.plot_name(), overide=False)

