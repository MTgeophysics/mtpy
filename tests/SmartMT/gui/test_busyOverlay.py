from PyQt4 import QtCore
from unittest import TestCase

import sys
from PyQt4.QtGui import QMainWindow, QWidget, QTextEdit, QGridLayout, QPushButton, QApplication
from PyQt4.QtTest import QTest

from mtpy.gui.SmartMT.gui.busy_indicators import BusyOverlay
from tests.SmartMT.gui import _click_area


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)

        widget = QWidget(self)
        self.editor = QTextEdit()
        self.editor.setPlainText("0123456789" * 100)
        layout = QGridLayout(widget)
        layout.addWidget(self.editor, 0, 0, 1, 3)
        self.button = QPushButton("Wait")
        layout.addWidget(self.button, 1, 1, 1, 1)

        self.setCentralWidget(widget)
        self.overlay = BusyOverlay(self.centralWidget())
        self.overlay.hide()
        self.button.clicked.connect(self.overlay.show)

    def resizeEvent(self, event):
        self.overlay.resize(event.size())
        event.accept()


class TestBusyOverlay(TestCase):
    def test(self):
        app = QApplication(sys.argv)
        window = MainWindow()
        window.show()
        QTest.qWaitForWindowShown(window)
        _click_area(window.button)
        # QTest.mouseClick(window.button, QtCore.Qt.LeftButton)
        QTest.qWait(5000)
        window.overlay.hide()
