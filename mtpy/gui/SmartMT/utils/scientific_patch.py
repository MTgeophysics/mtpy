# -*- coding: utf-8 -*-
"""
    Description:
        monkey patches for QSpinbox instances
    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
import re
import types

import numpy as np

# Regular expression to find floats. Match groups are the whole string, the
# whole coefficient, the decimal part of the coefficient, and the exponent
# part.
from mtpy.gui.SmartMT.utils.validator import FloatValidator


def patch_QDoubleSpinBox(spinbox):
    """
    this function patches the q double spain box (instances) to display numbers in scientific notation
    :param spinbox:
    :return:
    """
    spinbox.setMinimum(-np.inf)
    spinbox.setMaximum(np.inf)
    spinbox.validator = FloatValidator()
    spinbox.lineEdit().setValidator(spinbox.validator)
    spinbox.setDecimals(1000)

    def validate(self, text, position):
        return self.validator.validate(text, position)

    spinbox.validate = types.MethodType(validate, spinbox)

    def fixup(self, text):
        return self.validator.fixup(text)

    spinbox.fixup = types.MethodType(fixup, spinbox)

    def valueFromText(self, text):
        text = str(text)
        return float(text)

    spinbox.valueFromText = types.MethodType(valueFromText, spinbox)

    def textFromValue(self, value):
        return format_float(value)

    spinbox.textFromValue = types.MethodType(textFromValue, spinbox)

    def stepBy(self, steps):
        text = self.cleanText()
        groups = _float_re.search(text).groups()
        decimal = float(groups[1])
        decimal += steps
        new_string = '{:g}'.format(decimal) + (groups[3] if groups[3] else "")
        self.lineEdit().setText(new_string)

    spinbox.stepBy = types.MethodType(stepBy, spinbox)

    return spinbox


def format_float(value):
    """modified form of 'g' format specifier."""
    string = "{:g}".format(value).replace("e+", "e")
    string = re.sub("e(-?)0*(\d+)", r"e\1\2", string)
    print string
    return string
