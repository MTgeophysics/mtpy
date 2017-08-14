import os
import re

from PyQt4 import QtGui


class FileValidator(QtGui.QValidator):
    def validate(self, QString, p_int):
        QString = str(QString)
        basename = os.path.basename(QString)
        dir = os.path.dirname(QString)
        if os.path.isfile(QString):
            return QtGui.QValidator.Acceptable, p_int
        elif dir == "" or (os.path.isdir(dir) and any([item.startswith(basename) for item in os.listdir(dir)])):
            return QtGui.QValidator.Intermediate, p_int
        else:
            return QtGui.QValidator.Invalid, p_int


class DirectoryValidator(QtGui.QValidator):
    def validate(self, QString, p_int):
        QString = str(QString)
        basename = os.path.basename(QString)
        dir = os.path.dirname(QString)
        if os.path.isdir(QString):
            return QtGui.QValidator.Acceptable, p_int
        elif dir == "" or (os.path.isdir(dir) and any([item.startswith(basename)
                                                       for item in os.listdir(dir)
                                                       if os.path.isdir(os.path.join(dir, item))])):
            return QtGui.QValidator.Intermediate, p_int
        else:
            return QtGui.QValidator.Invalid, p_int


class FloatValidator(QtGui.QValidator):
    def validate(self, string, position):
        # print "validating"
        string = str(string)
        if valid_float_string(string):
            state = QtGui.QValidator.Acceptable
        elif string == "" or string[position - 1] in 'e.-+':
            state = QtGui.QValidator.Intermediate
        else:
            state = QtGui.QValidator.Invalid
        return state, position

    def fixup(self, text):
        text = str(text)
        match = _float_re.search(text)
        # print match.groups()[0]
        return match.groups()[0] if match else ""


_float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')


def valid_float_string(string):
    match = _float_re.search(string)
    return match.groups()[0] == string if match else False
