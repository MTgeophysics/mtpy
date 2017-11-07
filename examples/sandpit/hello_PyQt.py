#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Here we provide the necessary imports.
# The basic GUI widgets are located in QtGui module. 
import sys
from PyQt4.QtGui import QApplication, QWidget

# Every PyQt4 application must create an application object.
# The application object is located in the QtGui module.
a = QApplication(sys.argv)

# The QWidget widget is the base class of all user interface objects in PyQt4.
# We provide the default constructor for QWidget. The default constructor has no parent.
# A widget with no parent is called a window. 
w = QWidget()

w.resize(320, 240)  # The resize() method resizes the widget.
w.setWindowTitle("Hello, World!")  # Here we set the title for our window.
w.show()  # The show() method displays the widget on the screen.

sys.exit(a.exec_())  # Finally, we enter the mainloop of the application.
