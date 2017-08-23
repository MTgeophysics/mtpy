# -*- coding: utf-8 -*-
"""
    Description:
        classes for embedding matplot figures
        based on https://matplotlib.org/examples/user_interfaces/embedding_in_qt4.html

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""

import matplotlib.pyplot as plt
from matplotlib.backends import qt_compat
from matplotlib.widgets import AxesWidget

use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE
if use_pyside:
    from PySide import QtGui, QtCore
else:
    from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MPLCanvas(FigureCanvas):
    """
    Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.).
    """

    def __init__(self, parent=None, width=5, hight=4, dpi=100):
        self._fig = Figure(figsize=(width, hight), dpi=dpi, facecolor='none')
        self._axes = self._fig.add_subplot(111)

        self.compute_initial_figure()

        FigureCanvas.__init__(self, self._fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass

    def update_figure(self):
        pass


class Cursor(AxesWidget):
    """
    inspaired by matplotlib.widgets.Cursor
    """
    def __init__(self, ax, track_x=True, track_y=True, show_coord=True, text_format="(%.2f, %.2f)", col='green', useblit=True, **lineprops):
        AxesWidget.__init__(self, ax)

        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('draw_event', self.clear)

        self.visible = True
        self.horizOn = track_y
        self.vertOn = track_x
        self.show_coord = show_coord
        self.text_format = text_format
        self.useblit = useblit and self.canvas.supports_blit

        if self.useblit:
            lineprops['animated'] = True
        if 'color' not in lineprops:
            lineprops['color'] = col
        self.lineh = ax.axhline(ax.get_ybound()[0], visible=False, **lineprops)
        self.linev = ax.axvline(ax.get_xbound()[0], visible=False, **lineprops)
        self.text = ax.text(ax.get_xbound()[0], ax.get_ybound()[0], '', fontsize=8, color=col, visible=False)

        self.background = None
        self.needclear = False

    def clear(self, event):
        """clear the cursor"""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.linev.set_visible(False)
        self.lineh.set_visible(False)
        self.text.set_visible(False)

    def onmove(self, event):
        """on mouse motion draw the cursor if visible"""
        if self.ignore(event):
            return
        if not self.canvas.widgetlock.available(self):
            return
        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)
            self.text.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        self.needclear = True
        if not self.visible:
            return
        self.linev.set_xdata((event.xdata, event.xdata))

        self.lineh.set_ydata((event.ydata, event.ydata))
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)

        if self.visible and self.show_coord:
            self.text.set_visible(True)
            if self.vertOn and self.horizOn:
                self.text.set_text(self.text_format % (event.xdata, event.ydata))
            elif self.vertOn:
                self.text.set_text(self.text_format % event.xdata)
            elif self.ly:
                self.text.set_text(self.text_format % event.ydata)
            self.text.set_position((event.xdata, event.ydata))
        else:
            self.text.set_visible(False)

        self._update()

    def _update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.ax.draw_artist(self.text)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()
        return False


class MathTextLabel(QtGui.QWidget):
    def __init__(self, mathText,  parent=None, **kwargs):
        QtGui.QWidget.__init__(self, parent, **kwargs)

        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        r, g, b, a = self.palette().base().color().getRgbF()

        self._figure = Figure(edgecolor=(r, g, b), facecolor=(r, g, b))
        self._canvas = FigureCanvas(self._figure)

        layout.addWidget(self._canvas)

        self._figure.clear()
        text = self._figure.suptitle(mathText,
                                     x=0.0,
                                     y=1.0,
                                     horizontalalignment='left',
                                     verticalalignment='top',
                                     size=QtGui.QApplication.font().pointSize()*2)
        self._canvas.draw()

        (x0, y0), (x1, y1) = text.get_window_extent().get_points()
        w = x1 - x0
        h = y1 - y0

        self._figure.set_size_inches(w/80, h/80)
        self.setFixedSize(w ,h)
