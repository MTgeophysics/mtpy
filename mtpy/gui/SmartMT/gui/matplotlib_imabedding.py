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
from matplotlib import patches
from matplotlib.backends import qt_compat
from matplotlib.widgets import AxesWidget

use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE
if use_pyside:
    from PySide import QtGui
else:
    from PyQt4 import QtGui
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

    def __init__(self, ax, track_x=True, track_y=True, show_coord=True, text_format="(%.2f, %.2f)", col='green',
                 useblit=True, show_drag=False, **lineprops):
        AxesWidget.__init__(self, ax)

        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('draw_event', self.clear)
        if show_drag:
            self.connect_event('button_press_event', self.on_press)
            self.connect_event('button_release_event', self.on_release)

        self.visible = True
        self.horizOn = track_y
        self.vertOn = track_x
        self.drag_on = show_drag
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
        self.area = patches.Rectangle((0, 0), 0, 0, fc=col, ec=None, alpha=0.2, visible=False)
        self.ax.add_patch(self.area)

        self._press = None
        self.background = None
        self.needclear = False

    def on_press(self, event):
        if not event.inaxes:
            return
        else:
            self._press = event.xdata, event.ydata
            x0, y0 = self._press
            self.area.set_visible(self.visible)
            if self.vertOn:
                self.area.set_x(x0)
            else:
                self.area.set_x(0)
            if self.horizOn:
                self.area.set_y(y0)
            else:
                self.area.set_y(0)

    def on_release(self, event):
        self._press = None
        self.area.set_visible(False)

    def clear(self, event):
        """clear the cursor"""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.linev.set_visible(False)
        self.lineh.set_visible(False)
        self.text.set_visible(False)
        self.area.set_visible(False)

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
            self.area.set_visible(False)

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
        if self.drag_on and self._press and self.visible:
            self.area.set_visible(True)
            x0, y0 = self._press
            if self.vertOn:
                dx = event.xdata - x0
                self.area.set_width(dx)
            else:
                self.area.set_width(self.ax.get_xlim()[1])
            if self.horizOn:
                dy = event.ydata - y0
                self.area.set_height(dy)
            else:
                self.area.set_height(self.ax.get_ylim()[1])

        if self.visible and self.show_coord:
            self.text.set_visible(True)
            if self.vertOn and self.horizOn:
                self.text.set_text(self.text_format % (event.xdata, event.ydata))
            elif self.vertOn:
                self.text.set_text(self.text_format % event.xdata)
            elif self.horizOn:
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
            self.ax.draw_artist(self.area)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()
        return False


class MathTextLabel(QtGui.QWidget):
    def __init__(self, parent=None, math_text=None, horizontalalignment='left', verticalalignment='top', **kwargs):
        QtGui.QWidget.__init__(self, parent, **kwargs)

        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        r, g, b, a = self.palette().base().color().getRgbF()

        self._figure = Figure(edgecolor=(r, g, b), facecolor=(r, g, b))
        self._canvas = FigureCanvas(self._figure)

        layout.addWidget(self._canvas)

        self._figure.clear()
        self._display_text = self._figure.suptitle(math_text,
                                                   x=0.0,
                                                   y=1.0,
                                                   horizontalalignment=horizontalalignment,
                                                   verticalalignment=verticalalignment,
                                                   size=QtGui.QApplication.font().pointSize() * 2)
        self._canvas.draw()

        (x0, y0), (x1, y1) = self._display_text.get_window_extent().get_points()
        w = x1 - x0
        h = y1 - y0

        self._figure.set_size_inches(w / 80, h / 80)
        self.setFixedSize(w, h)

    def set_math_text(self, math_text):
        self._display_text.set_text(math_text)
        self._canvas.draw_idle()
