import numpy as np
from PyQt4 import QtGui, QtCore


class ProgressBar(QtGui.QWidget):
    style = """
    QProgressBar
    {
        border: 2px solid grey;
        border-radius: 5px;
        text-align: center;
    }
    """
    def __init__(self, parent=None, title=None, minimum=0, maximum=1, value=0):
        super(ProgressBar, self).__init__(parent)
        layout = QtGui.QGridLayout(self)

        self.progressbar = QtGui.QProgressBar(self)
        self.progressbar.setMinimum(minimum)
        self.progressbar.setMaximum(maximum)
        self.progressbar.setValue(value)
        self.progressbar.setStyleSheet(self.style)
        self.label = QtGui.QLabel("")
        self.label.setStyleSheet("Qlabel { font-size: 20px }")
        self.label.setAlignment(QtCore.Qt.AlignCenter)

        self.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.setFixedSize(256, 64)
        # self.progressbar.setValue(1)
        layout.addWidget(self.progressbar, 0, 0)
        layout.addWidget(self.label, 0, 0)
        self.setLayout(layout)

        if title:
            self.setWindowTitle(title)

    def setValue(self, value):
        self.progressbar.setValue(value)

    def setMaximumValue(self, value):
        self.progressbar.setMaximum(value)

    def incrementValue(self, increment=1):
        self.progressbar.setValue(self.progressbar.value() + increment)

    def onStart(self):
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.show()
        # self.setWindowState(self.windowState() | QtCore.Qt.WindowActive)

    def onFinished(self):
        self.hide()

    def updateIndicatorText(self, string):
        self.label.setText(string)


class BusyOverlay(QtGui.QWidget):
    """
    display an overlay to the window to indicate busy status
    this code is based on https://wiki.python.org/moin/PyQt/A%20full%20widget%20waiting%20indicator
    """

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        palette = QtGui.QPalette(self.palette())
        palette.setColor(palette.Background, QtCore.Qt.transparent)
        self.setPalette(palette)
        self.timer = None

    def paintEvent(self, event):
        painter = QtGui.QPainter()
        painter.begin(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.fillRect(event.rect(), QtGui.QBrush(QtGui.QColor(255, 255, 255, 127)))
        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))

        for i in range(6):
            if (self.counter / 5) % 6 == i:
                painter.setBrush(QtGui.QBrush(QtGui.QColor(127 + (self.counter % 5) * 32, 127, 127)))
            else:
                painter.setBrush(QtGui.QBrush(QtGui.QColor(127, 127, 127)))
            painter.drawEllipse(
                self.width() / 2 + 30 * np.cos(2 * np.pi * i / 6.0) - 10,
                self.height() / 2 + 30 * np.sin(2 * np.pi * i / 6.0) - 10,
                20, 20
            )

        painter.end()

    def showEvent(self, event):
        self.timer = self.startTimer(50)
        self.counter = 0

    def timerEvent(self, event):
        self.counter += 1
        self.update()

    def hideEvent(self, event):
        if self.timer is not None:
            self.killTimer(self.timer)
