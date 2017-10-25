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
import inspect
import traceback

import matplotlib.pyplot as plt
from qtpy import QtCore, QT_VERSION
from qtpy.QtWidgets import QWidget, QVBoxLayout
from qtpy.QtCore import Signal
from matplotlib.figure import Figure
if QT_VERSION.startswith('4'):
    from matplotlib.backends.backend_qt4agg import FigureCanvas
    from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
else:
    from matplotlib.backends.backend_qt5agg import FigureCanvas
    from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from mtpy.gui.SmartMT.gui.plot_parameter import PlotParameter
from mtpy.utils.mtpylog import MtPyLog


class VisualizationBase(QtCore.QThread):
    """
    plugin base for data visualization
    """

    # __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, parent):
        """
        set up ui and attributes here
        this function sets up the empty gui area for adding gui components that defines information necessary
        for visualization
        :param parent:
        """
        QtCore.QThread.__init__(self, parent)
        self._parent = parent
        self._mt_objs = None
        self._fig = None
        self._plotting_object = None
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self._parameter_ui = PlotParameter(self._parent)

        # add plot common setting gui
        # self._common_ui = CommonSettings(self._parameter_ui)
        # apply default common settings
        self.default_common_settings()
        # self._parameter_ui.add_parameter_groubox(self._common_ui)

    def default_common_settings(self):
        """
        this function will be called to initialize the common_settings parameters of the plot.
        including titles (title text and its font and positions) and figure size.
        The default value should fit for the most purpose but override this function to change
        the default values when necessary.
        :return:
        """
        pass

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
    def get_plot_tooltip(self):
        """
        returns the string that describe the plotting parameter, which will be used in the
        tooltips of the popup image window containing the plot to help users identify plots
        :return:
        """
        pass

    # plotting_started = pyqtSignal()
    # plotting_finished = pyqtSignal()

    def run(self):
        # self.setTerminationEnabled(True)
        plt.clf()
        try:
            self.plot()
            # change size and title
            if self._parameter_ui.customized_figure_size():
                self._fig.set_size_inches(self._parameter_ui.get_size_inches_width(),
                                          self._parameter_ui.get_size_inches_height())
                self._fig.set_dpi(self._parameter_ui.get_dpi())
                self._fig.set_tight_layout(self._parameter_ui.get_layout())
            if self._parameter_ui.customized_figure_title():
                self._fig.suptitle(self._parameter_ui.get_title(),
                                   **self._parameter_ui.get_title_font_dict())

            self.plotting_completed.emit(self._fig)
        except Exception as e:
            frm = inspect.trace()[-1]
            mod = inspect.getmodule(frm[0])
            self.plotting_error.emit("{}: {}".format(mod.__name__, e.message), traceback.format_exc())

    plotting_error = Signal(str, str)
    plotting_completed = Signal(Figure)

    def get_fig(self):
        return self._fig


class MPLCanvasWidget(QWidget):
    def __init__(self, fig):
        QWidget.__init__(self)
        self._fig = fig
        self._layout = QVBoxLayout()
        self._canvas = FigureCanvas(self._fig)
        self._toolbar = NavigationToolbar(self._canvas, self)
        self._layout.addWidget(self._toolbar)
        self._layout.addWidget(self._canvas)
        self.setLayout(self._layout)

        # resize the figure
        # dpi = self._fig.get_dpi()
        # width = int(self._fig.get_figwidth() * dpi) + margins[0] + margins[2]
        # height = int(self._fig.get_figheight() * dpi) + margins[1] + margins[3]
        # self._canvas.resize(width, height)
        # self._canvas.resize(self._canvas.sizeHint())
        # this is a hack to keep the figure close to the right size (pixels)
        # everything above this code does not work
        size = self.sizeHint()
        margins = self._layout.getContentsMargins()
        width = size.width() + margins[0]
        height = size.height() + self._toolbar.sizeHint().height()
        # print width, height
        self.resize(width, height)

    def get_fig(self):
        return self._fig
