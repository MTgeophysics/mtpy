import random

from PyQt4 import QtCore
from PyQt4.QtTest import QTest


def _click_area(qobj, pos=None, offset=None):
    geom = qobj.geometry()
    if pos is None:
        pos = QtCore.QPoint(geom.width() * random.random(), geom.height() * random.random())
    if offset is not None:
        pos += offset
    QTest.mouseClick(qobj, QtCore.Qt.LeftButton, pos=pos)

