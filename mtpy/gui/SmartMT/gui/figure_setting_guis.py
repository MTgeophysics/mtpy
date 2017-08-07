from PyQt4 import QtGui

from mtpy.gui.SmartMT.gui.plot_parameter import COLORS, SIMPLE_COLORS
from mtpy.gui.SmartMT.ui_asset.groupbox_color_bar import Ui_GroupBox_ColorBar
from mtpy.gui.SmartMT.ui_asset.groupbox_common import Ui_GroupBox_common_settings
from mtpy.gui.SmartMT.ui_asset.groupbox_font import Ui_GroupBox_Font


class ColorBar(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_ColorBar()
        self.ui.setupUi(self)

        # connect event
        self.ui.horizontalSlider_x.valueChanged.connect(self._x_slider_value_changed)
        self.ui.horizontalSlider_y.valueChanged.connect(self._y_slider_value_changed)
        self.ui.horizontalSlider_width.valueChanged.connect(self._width_slider_value_changed)
        self.ui.horizontalSlider_height.valueChanged.connect(self._height_slider_value_changed)
        self.ui.comboBox_orientation.currentIndexChanged.connect(self._orientation_changed)
        self.ui.doubleSpinBox_x.editingFinished.connect(self._update_slider_x)
        self.ui.doubleSpinBox_y.editingFinished.connect(self._update_slider_y)
        self.ui.doubleSpinBox_width.editingFinished.connect(self._update_slider_width)
        self.ui.doubleSpinBox_height.editingFinished.connect(self._update_slider_height)

    def _width_slider_value_changed(self, value):
        self.ui.doubleSpinBox_width.setValue(value / 100.0)

    def _height_slider_value_changed(self, value):
        self.ui.doubleSpinBox_height.setValue(value / 100.0)

    def _x_slider_value_changed(self, value):
        self.ui.doubleSpinBox_x.setValue(value / 100.0)

    def _y_slider_value_changed(self, value):
        self.ui.doubleSpinBox_y.setValue(value / 100.0)

    def _orientation_changed(self, *args):
        x = self.ui.doubleSpinBox_x.value()
        y = self.ui.doubleSpinBox_y.value()
        width = self.ui.doubleSpinBox_width.value()
        height = self.ui.doubleSpinBox_height.value()
        self.ui.horizontalSlider_x.setValue(100 - x * 100)
        self.ui.horizontalSlider_y.setValue(100 - y * 100)
        self.ui.horizontalSlider_width.setValue(height * 100)
        self.ui.horizontalSlider_height.setValue(width * 100)

    def _update_slider_x(self):
        value = int(self.ui.doubleSpinBox_x.value() * 100)
        self.ui.horizontalSlider_x.setValue(value)

    def _update_slider_y(self):
        value = int(self.ui.doubleSpinBox_y.value() * 100)
        self.ui.horizontalSlider_y.setValue(value)

    def _update_slider_width(self):
        value = int(self.ui.doubleSpinBox_width.value() * 100)
        self.ui.horizontalSlider_width.setValue(value)

    def _update_slider_height(self):
        value = int(self.ui.doubleSpinBox_height.value() * 100)
        self.ui.horizontalSlider_height.setValue(value)

    _cb_orientation = ['vertical', 'horizontal']

    def get_colorbar_dict(self):
        if self.isChecked():
            cb_dict = {
                'orientataion': self._cb_orientation[self.ui.comboBox_orientation.currentIndex()],
                'position': (
                    self.ui.doubleSpinBox_x.value(),
                    self.ui.doubleSpinBox_y.value(),
                    self.ui.doubleSpinBox_width.value(),
                    self.ui.doubleSpinBox_height.value()
                )
            }
            return cb_dict
        else:
            return None


class Font(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Font()
        self.ui.setupUi(self)
        self.ui.checkBox_size.stateChanged.connect(self.size_state_changed)
        self.ui.checkBox_weight.stateChanged.connect(self.weight_state_changed)
        self.ui.checkBox_color.stateChanged.connect(self.color_state_changed)
        self._simple_color = simple_color
        if not self._simple_color:
            self.ui.comboBox_color.clear()
            cnames = [name for name, hex in COLORS]
            self.ui.comboBox_color.addItems(cnames)

    def size_state_changed(self, p_int):
        if p_int == 0:
            self.ui.spinBox_size.setEnabled(False)
        else:
            self.ui.spinBox_size.setEnabled(True)

    def weight_state_changed(self, p_int):
        if p_int == 0:
            self.ui.comboBox_weight.setEnabled(False)
        else:
            self.ui.comboBox_weight.setEnabled(True)

    def color_state_changed(self, p_int):
        if p_int == 0:
            self.ui.comboBox_color.setEnabled(False)
        else:
            self.ui.comboBox_color.setEnabled(True)

    def hide_size(self):
        self.ui.spinBox_size.hide()
        self.ui.checkBox_size.hide()

    def hide_weight(self):
        self.ui.comboBox_weight.hide()
        self.ui.checkBox_weight.hide()

    def hide_color(self):
        self.ui.comboBox_color.hide()
        self.ui.checkBox_color.hide()

    def get_size(self):
        return self.ui.spinBox_size.value()

    def get_weight(self):
        return str(self.ui.comboBox_weight.currentText())

    def get_color(self):
        if self._simple_color:
            return SIMPLE_COLORS[self.ui.comboBox_color.currentIndex()]
        else:
            return COLORS[self.ui.comboBox_color.currentIndex()][1]


class CommonSettings(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_common_settings()
        self.ui.setupUi(self)

        # dpi
        self.ui.spinBox_dpi.valueChanged.connect(self._dpi_changed)
        # inches
        self.ui.doubleSpinBox_width_inches.valueChanged.connect(self._width_inches_changed)
        self.ui.doubleSpinBox_height_inches.valueChanged.connect(self._height_inches_changed)
        # pixels
        self.ui.spinBox_width_pixels.valueChanged.connect(self._width_pixels_changed)
        self.ui.spinBox_height_pixels.valueChanged.connect(self._height_pixels_changed)
        # x
        self.ui.horizontalSlider_x.valueChanged.connect(self._x_slider_changed)
        self.ui.doubleSpinBox_x.valueChanged.connect(self._x_spinbox_changed)
        # y
        self.ui.horizontalSlider_y.valueChanged.connect(self._y_slider_changed)
        self.ui.doubleSpinBox_y.valueChanged.connect(self._y_spinbox_changed)

    def get_title(self):
        return str(self.ui.lineEdit_title.text())

    def set_title(self, title):
        self.ui.lineEdit_title.setText(title)

    def _x_slider_changed(self, value):
        self.ui.doubleSpinBox_x.blockSignals(True)
        self.ui.doubleSpinBox_x.setValue(value / 100.0)
        self.ui.doubleSpinBox_x.blockSignals(False)

    def _y_slider_changed(self, value):
        self.ui.doubleSpinBox_y.blockSignals(True)
        self.ui.doubleSpinBox_y.setValue(value / 100.0)
        self.ui.doubleSpinBox_y.blockSignals(False)

    def _x_spinbox_changed(self, value):
        self.ui.horizontalSlider_x.blockSignals(True)
        self.ui.horizontalSlider_x.setValue(value * 100)
        self.ui.horizontalSlider_x.blockSignals(False)

    def _y_spinbox_changed(self, value):
        self.ui.horizontalSlider_y.blockSignals(True)
        self.ui.horizontalSlider_y.setValue(value * 100)
        self.ui.horizontalSlider_y.blockSignals(False)

    def _dpi_changed(self, dpi):
        self.ui.doubleSpinBox_height_inches.blockSignals(True)
        self.ui.doubleSpinBox_width_inches.blockSignals(True)
        self.ui.spinBox_height_pixels.setValue(
            self.ui.doubleSpinBox_height_inches.value() * dpi
        )
        self.ui.spinBox_width_pixels.setValue(
            self.ui.doubleSpinBox_width_inches.value() * dpi
        )
        self.ui.doubleSpinBox_height_inches.blockSignals(False)
        self.ui.doubleSpinBox_width_inches.blockSignals(False)

    def _width_pixels_changed(self, width):
        self.ui.doubleSpinBox_width_inches.blockSignals(True)
        new_width_inches = width / float(self.ui.spinBox_dpi.value())
        self.ui.doubleSpinBox_width_inches.setValue(new_width_inches)
        self.ui.doubleSpinBox_width_inches.blockSignals(False)

    def _height_pixels_changed(self, height):
        self.ui.doubleSpinBox_height_inches.blockSignals(True)
        new_height_inches = height / float(self.ui.spinBox_dpi.value())
        self.ui.doubleSpinBox_height_inches.setValue(new_height_inches)
        self.ui.doubleSpinBox_height_inches.blockSignals(False)

    def _width_inches_changed(self, width):
        self.ui.spinBox_width_pixels.blockSignals(True)
        self.ui.spinBox_width_pixels.setValue(
            width * self.ui.spinBox_dpi.value()
        )
        self.ui.spinBox_width_pixels.blockSignals(False)

    def _height_inches_changed(self, height):
        self.ui.spinBox_height_pixels.blockSignals(True)
        self.ui.spinBox_height_pixels.setValue(
            height * self.ui.spinBox_dpi.value()
        )
        self.ui.spinBox_height_pixels.blockSignals(False)

    def customized_figure_title(self):
        return self.ui.groupBox_title.isChecked()

    def customized_figure_size(self):
        return self.ui.groupBox_figure_size.isChecked()

    def get_size_inches_width(self):
        return self.ui.doubleSpinBox_width_inches.value()

    def get_size_inches_height(self):
        return self.ui.doubleSpinBox_height_inches.value()

    def get_dpi(self):
        return self.ui.spinBox_dpi.value()

    def get_layout(self):
        if self.ui.checkBox_tight_layout.isChecked():
            return True
        else:
            return False

    def get_title_font_dict(self):
        font_properties = {
            'x': self.ui.doubleSpinBox_x.value(),
            'y': self.ui.doubleSpinBox_y.value(),
            'horizontalalignment': self._horizontalalignment[
                self.ui.comboBox_horizontal_alignment.currentIndex()],
            'verticalalignment': self._verticalalignment[
                self.ui.comboBox_vertical_alignment.currentIndex()
            ],
            'fontsize': self.ui.spinBox_fontsize.value()
        }
        return font_properties

    _horizontalalignment = ['right', 'center', 'left']
    _verticalalignment = ['top', 'center', 'bottom', 'baseline']
