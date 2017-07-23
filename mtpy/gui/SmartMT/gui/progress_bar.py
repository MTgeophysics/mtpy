from PyQt4 import QtGui


class ProgressBar(QtGui.QWidget):
    def __init__(self, parent=None, title=None, minimum=0, maximum=0, value=0):
        super(ProgressBar, self).__init__(parent)
        layout = QtGui.QVBoxLayout(self)
        self.progressbar = QtGui.QProgressBar(self)
        self.progressbar.setMinimum(minimum)
        self.progressbar.setMaximum(maximum)
        self.progressbar.setValue(value)
        self.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.setFixedSize(256, 64)
        # self.progressbar.setValue(1)
        layout.addWidget(self.progressbar)
        if title:
            self.setWindowTitle(title)

    def setValue(self, value):
        self.progressbar.setValue(value)

    def setMaximumValue(self, value):
        self.progressbar.setMaximum(value)

    def incrementValue(self, increment=1):
        self.progressbar.setValue(self.progressbar.value()+increment)

    def onStart(self):
        self.show()

    def onFinished(self):
        self.hide()
