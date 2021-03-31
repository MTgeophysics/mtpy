# import os
# os.environ['QT_API'] = 'pyqt'  # use pyqt4 api

from qtpy import QT_VERSION
import matplotlib

if QT_VERSION.startswith("4"):
    matplotlib.use("Qt4Agg")
elif QT_VERSION.startswith("5"):
    matplotlib.use("Qt5Agg")
