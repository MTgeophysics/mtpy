import numpy as np
from PyQt4 import QtGui

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_single import (
    Ui_groupBoxFrequency_pereiod_single,
)


class FrequencySingle(QtGui.QGroupBox):
    """
    Frequency selection (single frequency)
    """

    _unit_period = "second"
    _unit_frequency = "Hz"
    _title_period = "Period"
    _title_frequency = "Frequency"

    def __init__(self, parent, use_period=False):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self.use_period = use_period
        self.ui = Ui_groupBoxFrequency_pereiod_single()
        self.ui.setupUi(self)
        self._histogram = FrequencySingle.FrequencyHistogram()
        self.set_use_period(self.use_period)
        # add matplotlib canvas
        self.ui.verticalLayoutFrequencyPeriod.addWidget(self._histogram)
        # connect components
        self.ui.comboBoxPeriod.currentIndexChanged.connect(self.update_histogram)
        self.ui.comboBoxPeriod.editTextChanged.connect(self.update_histogram)
        self._histogram.mpl_connect("button_release_event", self._mouse_pick)

    def toggle_time_scale(self, *args):
        self.use_period = not self.use_period
        self.set_use_period(self.use_period)

    def set_use_period(self, use_period=False):
        if use_period:
            self._histogram.set_unit(self._unit_period)
            self._histogram.set_title(self._title_period)
            title = "%s (%s)" % (self._title_period, self._unit_period)
        else:
            self._histogram.set_unit(self._unit_frequency)
            self._histogram.set_title(self._title_frequency)
            title = "%s (%s)" % (self._title_frequency, self._unit_frequency)
        self.setTitle(title)
        self._update_frequency()

    def _mouse_pick(self, event):
        if not event.inaxes:
            return
        x = event.xdata
        self.ui.comboBoxPeriod.setEditText("%.5f" % x)

    def get_frequency(self):
        return float(self.ui.comboBoxPeriod.currentText())

    def update_histogram(self):
        value = float(self.ui.comboBoxPeriod.currentText())
        self._histogram.set_current_frequency(value)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._update_frequency()

    def _update_frequency(self):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))

            if self.use_period:
                all_periods = 1.0 / np.array(all_freqs)
                # self._histogram.set_data(all_periods)
                all_unique = sorted(list(set(all_periods)))

            else:
                # self._histogram.set_data(all_freqs)
                all_unique = sorted(list(set(all_freqs)))
            self._histogram.set_data(all_unique)
            self._histogram.update_figure()
            # sort all frequencies in ascending order
            for period in all_unique:
                self.ui.comboBoxPeriod.addItem("%.5f" % period)
            self.ui.comboBoxPeriod.setCurrentIndex(0)
            self.update_histogram()

    class FrequencyHistogram(MPLCanvas):
        def __init__(self, parent=None, width=5, height=2, dpi=100):
            self.artists = dict()
            self._frequency = None
            self._current_frequency = None
            self._title = None
            self._unit = None
            MPLCanvas.__init__(self, parent, width, height, dpi)
            self._lx = None
            self.cursor = None

            # self.mpl_connect('motion_notify_event', self.cursor)
            self.mpl_connect("button_release_event", self.mouse_pick)
            self.setMinimumSize(200, 150)
            self.resize(self.sizeHint())

        # def mouse_move(self, event):
        #     if not event.inaxes:
        #         return
        #     x = event.xdata
        #     y = event.ydata
        #     if self._cursor_x is None:
        #         self._cursor_x = self._axes.axvline(linewidth=1, color="green")
        #     if self._cursor_text is None:
        #         self._cursor_text = self._axes.text(0.0, 0.0, '', fontsize=8)
        #     self._cursor_x.set_xdata(x)
        #     self._cursor_text.set_text('period=%.2f' % x)
        #     self._cursor_text.set_position((x, y))
        #     self.draw()

        def set_title(self, title):
            self._title = title

        def set_unit(self, unit):
            if unit != self._unit:
                self._unit = unit
                self.cursor = Cursor(
                    self._axes,
                    track_y=False,
                    text_format="%f " + self._unit,
                    useblit=True,
                )

        def mouse_pick(self, event):
            if not event.inaxes:
                return
            x = event.xdata
            self.set_current_frequency(x)

        def compute_initial_figure(self):
            if self._frequency is not None:
                self._axes.tick_params(axis="both", which="major", labelsize=6)
                self._axes.tick_params(axis="both", which="minor", labelsize=4)
                self._axes.hist(self._frequency)  # , 50, normed=1)
                if self._title and self._unit:
                    self._axes.set_xlabel(
                        "%s (%s)" % (self._title, self._unit), fontsize=8
                    )
                    self.figure.suptitle(
                        "%s Distribution in Selected Stations" % self._title, fontsize=8
                    )

                self._fig.set_tight_layout(True)

        def set_data(self, frequency):
            self._frequency = frequency
            self._lx = None
            self._current_frequency = None

        def set_current_frequency(self, freq):
            self._current_frequency = freq
            if self._lx is None:
                self._lx = self._axes.axvline(linewidth=2, color="red")
            self._lx.set_xdata(self._current_frequency)
            # if self._fig.canvas.supports_blit:
            #     self._axes.draw_artist(self._lx)
            #     self._fig.canvas.blit(self._axes.bbox)
            # else:
            #     self._fig.canvas.draw_idle()
            self._fig.canvas.draw_idle()

        def update_figure(self):
            # clear figure
            self._axes.cla()
            self.compute_initial_figure()
            self.draw()
