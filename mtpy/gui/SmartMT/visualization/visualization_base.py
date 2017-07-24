# -*- coding: utf-8 -*-
"""
    Description:
    base class used by plot/visualization plugins,
    where the standard methods are defined as abstract methods
    and standard utility methods are implemented

    Author: YingzhiGou
    Date: 20/06/2017
"""

import abc

import matplotlib.pyplot as plt
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from mtpy.gui.SmartMT.gui.plot_parameter import PlotParameter
from mtpy.gui.SmartMT.gui.progress_bar import ProgressBar
from mtpy.utils.mtpylog import MtPyLog


class VisualizationBase(object):
    """
    plugin base for data visualization
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, parent):
        """
        set up ui and attributes here
        this function sets up the empty gui area for adding gui components that defines information necessary
        for visualization
        :param parent:
        """
        self._parent = parent
        self._mt_objs = None
        self._fig = None
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self._parameter_ui = PlotParameter(self._parent)
        self._plotting_object = None
        # connect plot button
        QtCore.QObject.connect(self._parameter_ui.ui.pushButtonPlot, QtCore.SIGNAL("clicked()"), self.show_figure)

    def set_data(self, mt_objs):
        """
        set input data, default is a collection of mt_objs
        then call to update gui components that may requires mt_objs as input
        :param mt_objs:
        :return:
        """
        self._mt_objs = mt_objs
        self.update_ui()

    @property
    def parameter_ui(self):
        """
        returns GUI for selecting parameters of the plot
        :return:
        """
        return self._parameter_ui

    @abc.abstractmethod
    def update_ui(self):
        """
        all gui components of parameter selections that require data as input should be included in this
        function. this function will be called every time when the input data is changed
        :return:
        """
        pass

    @staticmethod
    @abc.abstractmethod
    def plot_name():
        """
        static methods that returns the name of the plot/visualization,
        which will be used to populate the plot option list in the GUI
        :return: name of the plot
        """
        return VisualizationBase.__name__

    @staticmethod
    @abc.abstractmethod
    def plot_description():
        """
        static method that returns the description and notes of the plot/visualization,
        which will be used to populate the plot description text area in the GUI
        HTML is preferred text format
        :return: description of the plot (ideally in html)
        """
        return VisualizationBase.__name__

        # _plottingFinished = QtCore.pyqtSignal()

        # def _create_plot(self):
        #     self.plot()
        # self._plottingFinished.emit()

    @abc.abstractmethod
    def plot(self):
        """
        run the routines to,
            1) collect parameters from the GUI
            2) generate plot (if any plotting object exists, assigned it to self._plotting_object)
            3) assign the matplotlib.figure to self._fig for display (MOST IMPORTANT)
        :return:
        """
        pass

    @abc.abstractmethod
    def get_parameter_str(self):
        """
        returns the string that describe the plotting parameter, which will be used in the
        tooltips of the popup image window containing the plot to help users identify plots
        :return:
        """
        pass

    def show_figure(self):
        """
        function that creates plot by calling self.plot(), then create subwindow for the created image,
        exception and error handling, progress indication are handled here
        :return:
        """
        # clear the figure if there is already one up
        plt.clf()
        # show progress bar
        progressbar = ProgressBar(title='Generating image...')
        progressbar.onStart()
        # self._plottingFinished.connect(progressbar.onFinished)
        # self._create_plot()
        try:
            self.plot()
            if self._fig:
                # self._fig.show()
                widget = MPLCanvasWidget(self._fig)

                progressbar.onFinished()

                self._parent._parent.create_subwindow(widget, "%s" % self.plot_name(), overide=False,
                                                      tooltip=self.get_parameter_str())
            else:
                progressbar.onFinished()
        except Exception as e:
            progressbar.onFinished()
            QtGui.QMessageBox.critical(self._parameter_ui, 'Plotting Error', e.message, QtGui.QMessageBox.Close)


class MPLCanvasWidget(QtGui.QWidget):
    def __init__(self, fig):
        QtGui.QWidget.__init__(self)
        self._fig = fig
        self._layout = QtGui.QVBoxLayout()
        self._canvas = FigureCanvas(self._fig)
        self._toolbar = NavigationToolbar(self._canvas, self)
        self._layout.addWidget(self._toolbar)
        self._layout.addWidget(self._canvas)
        self.setLayout(self._layout)

    def get_fig(self):
        return self._fig
