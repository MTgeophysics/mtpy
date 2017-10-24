# -*- coding: utf-8 -*-
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_aspect import Ui_GroupBox_aspect
from mtpy.gui.SmartMT.ui_asset.groupbox_text_box import Ui_GroupBox_text_box


class AspectRatio(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_aspect()
        self.ui.setupUi(self)

        # connect signal
        self.ui.radioButton_aspect_float.toggled.connect(self._aspect_float_toggled)

    def _aspect_float_toggled(self, checked):
        self.ui.doubleSpinBox_aspect_float.setEnabled(checked)

    def get_aspect(self):
        if self.ui.radioButton_aspect_auto.isChecked():
            return 'auto'
        elif self.ui.radioButton_aspect_equal.isChecked():
            return 'equal'
        else:
            return self.ui.doubleSpinBox_aspect_float.value()


class TextBox(QGroupBox):
    def __init__(self, parent, point_size=True, key_size=False):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_text_box()
        self.ui.setupUi(self)

        self._point_size = point_size
        self._key_size = key_size
        self.ui.comboBox_size.model().item(len(self._size_keys)).setEnabled(self._point_size)

        # connect signal
        self.ui.checkBox_size.stateChanged.connect(self._size_state_changed)
        self.ui.checkBox_weight.stateChanged.connect(self._weight_state_changed)
        self.ui.comboBox_size.currentIndexChanged.connect(self._size_index_changed)

        self.ui.horizontalSlider_x.valueChanged.connect(self._x_slider_value_changed)
        self.ui.horizontalSlider_y.valueChanged.connect(self._y_slider_value_changed)
        self.ui.doubleSpinBox_x.editingFinished.connect(self._update_slider_x)
        self.ui.doubleSpinBox_y.editingFinished.connect(self._update_slider_y)

        self.ui.horizontalSlider_x_pad.valueChanged.connect(self._x_pad_slider_value_changed)
        self.ui.horizontalSlider_y_pad.valueChanged.connect(self._y_pad_slider_value_changed)
        self.ui.doubleSpinBox_x_pad.editingFinished.connect(self._update_slider_x_pad)
        self.ui.doubleSpinBox_y_pad.editingFinished.connect(self._update_slider_y_pad)

    _size_keys = [
        'xx-small',
        'x-small',
        'small',
        'medium',
        'large',
        'x-large',
        'xx-large'
    ]

    def _location_x_changed(self, p_int):
        self.ui.horizontalSlider_x.setEnabled(p_int != 0)
        self.ui.doubleSpinBox_x.setEnabled(p_int != 0)

    def _location_y_changed(self, p_int):
        self.ui.horizontalSlider_y.setEnabled(p_int != 0)
        self.ui.doubleSpinBox_y.setEnabled(p_int != 0)

    def _size_index_changed(self, p_int):
        self.ui.spinBox_size.setEnabled(p_int >= len(self._size_keys))

    def _size_state_changed(self, p_int):
        if self._key_size:
            self.ui.comboBox_size.setEnabled(p_int != 0)
        else:
            self.ui.comboBox_size.setCurrentIndex(len(self._size_keys))

    def _weight_state_changed(self, p_int):
        self.ui.comboBox_weight.setEnabled(p_int != 0)

    def _x_slider_value_changed(self, value):
        self.ui.doubleSpinBox_x.setValue(value / 100.0)

    def _y_slider_value_changed(self, value):
        self.ui.doubleSpinBox_y.setValue(value / 100.0)

    def _update_slider_x(self):
        value = int(self.ui.doubleSpinBox_x.value() * 100)
        self.ui.horizontalSlider_x.setValue(value)

    def _update_slider_y(self):
        value = int(self.ui.doubleSpinBox_y.value() * 100)
        self.ui.horizontalSlider_y.setValue(value)

    def _x_pad_slider_value_changed(self, value):
        self.ui.doubleSpinBox_x_pad.setValue(value / 100.0)

    def _y_pad_slider_value_changed(self, value):
        self.ui.doubleSpinBox_y_pad.setValue(value / 100.0)

    def _update_slider_x_pad(self):
        value = int(self.ui.doubleSpinBox_x_pad.value() * 100)
        self.ui.horizontalSlider_x_pad.setValue(value)

    def _update_slider_y_pad(self):
        value = int(self.ui.doubleSpinBox_y_pad.value() * 100)
        self.ui.horizontalSlider_y_pad.setValue(value)

    def get_size(self):
        if self.ui.checkBox_size.isChecked():
            return self._size_keys[self.ui.comboBox_size.currentIndex()] \
                if self.ui.comboBox_size.currentIndex() < len(self._size_keys) \
                else self.ui.spinBox_size.value()
        else:
            return None

    def get_weight(self):
        if self.ui.checkBox_weight.isChecked():
            return str(self.ui.comboBox_weight.currentText())
        else:
            return None

    def get_location(self):
        if self.ui.groupBox_location.isChecked():
            return self.ui.doubleSpinBox_x.value(), self.ui.doubleSpinBox_y.value()
        else:
            return None

    def get_xpad(self):
        if self.ui.groupBox_padding.isChecked():
            return self.ui.doubleSpinBox_x_pad.value()
        else:
            return None

    def get_ypad(self):
        if self.ui.groupBox_padding.isChecked():
            return self.ui.doubleSpinBox_y_pad.value()
        else:
            return None
