import os
import re

from qtpy.QtGui import QValidator


class FileValidator(QValidator):
    def validate(self, QString, p_int):
        QString = str(QString)
        basename = os.path.basename(QString)
        dir = os.path.dirname(QString)
        if os.path.isfile(QString):
            return QValidator.Acceptable, p_int
        elif dir == "" or (os.path.isdir(dir) and any([item.startswith(basename) for item in os.listdir(dir)])):
            return QValidator.Intermediate, p_int
        else:
            return QValidator.Invalid, p_int


class DirectoryValidator(QValidator):
    def validate(self, QString, p_int):
        QString = str(QString)
        basename = os.path.basename(QString)
        dir = os.path.dirname(QString)
        if os.path.isdir(QString):
            return QValidator.Acceptable, p_int
        elif dir == "" or (os.path.isdir(dir) and any([item.startswith(basename)
                                                       for item in os.listdir(dir)
                                                       if os.path.isdir(os.path.join(dir, item))])):
            return QValidator.Intermediate, p_int
        else:
            return QValidator.Invalid, p_int


class FloatValidator(QValidator):
    def validate(self, string, position):
        # print "validating"
        string = str(string)
        if valid_float_string(string):
            state = QValidator.Acceptable
        elif string == "" or string[position - 1] in 'e.-+':
            state = QValidator.Intermediate
        else:
            state = QValidator.Invalid
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
