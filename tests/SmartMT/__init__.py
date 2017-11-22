from __future__ import print_function

import os
import pytest

qtpy = pytest.importorskip("qtpy")
pytestmark = pytest.mark.skipif(os.name == "posix" and 'DISPLAY' not in os.environ)

import random
import sys

import matplotlib
import sip
from qtpy import QtCore
from qtpy.QtTest import QTest
from qtpy.QtWidgets import QApplication

from mtpy.utils.mtpylog import MtPyLog

sip.setdestroyonexit(False)

app = QApplication(sys.argv)

import matplotlib.pyplot as plt
plt.ion()

MtPyLog.get_mtpy_logger(__name__).info("Testing using matplotlib backend {}".format(matplotlib.rcParams['backend']))


def _click_area(qobj, pos=None, offset=None, modifier=QtCore.Qt.NoModifier):
    geom = qobj.geometry()
    if pos is None:
        x = int(geom.width() * random.uniform(0.2, 0.8))  # avid to click on the edge of widgets
        y = int(geom.height() * random.uniform(0.2, 0.8))
        pos = QtCore.QPoint(x, y)
    if offset is not None:
        pos += offset
    QTest.mouseClick(qobj, QtCore.Qt.LeftButton, modifier=modifier, pos=pos)
    # print(pos.x(), pos.y())
