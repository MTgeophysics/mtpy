from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.Components import COLORS, SIMPLE_COLORS
from mtpy.gui.SmartMT.ui_asset.groupbox_font import Ui_GroupBox_Font


class Font(QGroupBox):
    def __init__(self, parent, simple_color=True, point_size=True, key_size=False):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Font()
        self.ui.setupUi(self)
        self.ui.checkBox_size.stateChanged.connect(self.size_state_changed)
        self.ui.checkBox_weight.stateChanged.connect(self.weight_state_changed)
        self.ui.checkBox_color.stateChanged.connect(self.color_state_changed)
        self.ui.comboBox_size.currentIndexChanged.connect(self.size_index_changed)
        self._simple_color = simple_color
        self._point_size = point_size
        self._key_size = key_size

        self.ui.comboBox_size.model().item(len(self._size_keys)).setEnabled(self._point_size)

        if not self._simple_color:
            self.ui.comboBox_color.clear()
            cnames = [name for name, hex in COLORS]
            self.ui.comboBox_color.addItems(cnames)

    _size_keys = [
        'xx-small',
        'x-small',
        'small',
        'medium',
        'large',
        'x-large',
        'xx-large'
    ]

    def size_index_changed(self, p_int):
        self.ui.spinBox_size.setEnabled(p_int >= len(self._size_keys))

    def size_state_changed(self, p_int):
        if self._key_size:
            self.ui.comboBox_size.setEnabled(p_int != 0)
        else:
            self.ui.comboBox_size.setCurrentIndex(len(self._size_keys))

    def weight_state_changed(self, p_int):
        self.ui.comboBox_weight.setEnabled(p_int != 0)

    def color_state_changed(self, p_int):
        self.ui.comboBox_color.setEnabled(p_int != 0)

    def hide_size(self):
        self.ui.spinBox_size.hide()
        self.ui.checkBox_size.hide()
        self.ui.comboBox_size.hide()

    def hide_weight(self):
        self.ui.comboBox_weight.hide()
        self.ui.checkBox_weight.hide()

    def hide_color(self):
        self.ui.comboBox_color.hide()
        self.ui.checkBox_color.hide()

    def get_size(self):
        if self.ui.checkBox_size.isChecked():
            return self._size_keys[self.ui.comboBox_size.currentIndex()] \
                if self.ui.comboBox_size.currentIndex() < len(self._size_keys) \
                else self.ui.spinBox_size.value()
        else:
            return None

    def get_weight(self):
        return str(self.ui.comboBox_weight.currentText())

    def get_color(self):
        if self._simple_color:
            return SIMPLE_COLORS[self.ui.comboBox_color.currentIndex()]
        else:
            return COLORS[self.ui.comboBox_color.currentIndex()][1]
