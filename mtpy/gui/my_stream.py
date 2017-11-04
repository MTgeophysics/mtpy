# -*- coding: utf-8 -*-
"""
Created on Fri Sep 02 12:33:27 2016

@author: jpeacock
"""

try:
    from PyQt4 import QtCore
except ImportError:
    try:
        from PyQt5 import QtCore
    except ImportError:
        raise ImportError('Cannot find PyQt4 or PyQt5')


class MyStream(QtCore.QObject):
    """
    this class will emit a signal
    """
    message = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(MyStream, self).__init__(parent)

    def write(self, message):
        self.message.emit(str(message))