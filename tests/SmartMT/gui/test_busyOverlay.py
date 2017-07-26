from unittest import TestCase

import sys
from PyQt4.QtGui import QMainWindow, QWidget, QTextEdit, QGridLayout, QPushButton, QApplication
from PyQt4.QtTest import QTest

from mtpy.gui.SmartMT.gui.busy_indicators import BusyOverlay


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)

        widget = QWidget(self)
        self.editor = QTextEdit()
        self.editor.setPlainText("0123456789" * 100)
        layout = QGridLayout(widget)
        layout.addWidget(self.editor, 0, 0, 1, 3)
        button = QPushButton("Wait")
        layout.addWidget(button, 1, 1, 1, 1)

        self.setCentralWidget(widget)
        self.overlay = BusyOverlay(self.centralWidget())
        self.overlay.hide()
        button.clicked.connect(self.overlay.show)

    def resizeEvent(self, event):
        self.overlay.resize(event.size())
        event.accept()


class TestBusyOverlay(TestCase):
    def test(self):
        app = QApplication(sys.argv)
        window = MainWindow()
        window.show()
        QTest.qWaitForWindowShown(window)
        # QTest.qWait(10000)
        app.exec_()
