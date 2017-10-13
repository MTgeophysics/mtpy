import random

from PyQt4 import QtCore
from PyQt4.QtTest import QTest


def _click_area(qobj, pos=None, offset=None):
    geom = qobj.geometry()
    if pos is None:
        x = int(geom.width() * random.uniform(0.1, 0.9))  # avid to click on the edge of widgets
        y = int(geom.height() * random.uniform(0.1, 0.9))
        pos = QtCore.QPoint(x, y)
    if offset is not None:
        pos += offset
    QTest.mouseClick(qobj, QtCore.Qt.LeftButton, pos=pos)
    print(pos.x(), pos.y())
