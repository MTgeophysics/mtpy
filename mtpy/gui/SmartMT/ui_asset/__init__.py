"""
This file is responsible for dynamic generation of ui_asset/*.py and ui_asset/*.pyc
when the start.py is first run. It will pick pyqt4 or pyqt5 whichever is installed.
"""
import os
import sys
import threading

from mtpy import MtPyLog

_thread_lock = threading.Lock()
_thread_lock.acquire()

_GUI_RESOURCE_BUILT = False

if not _GUI_RESOURCE_BUILT:
    logger = MtPyLog.get_mtpy_logger(__name__)

    from qtpy.uic import compileUiDir
    from qtpy import QT_VERSION

    if getattr(sys, 'frozen', False):
        # running in a pyinstaller bundle
        _dir = os.path.join(sys._MEIPASS, "mtpy/gui/SmartMT/ui_asset")
        sys.path.append(sys._MEIPASS)
    else:
        _dir = os.path.dirname(__file__)

    if QT_VERSION.startswith("4"):
        _resource_suffix = "_py{python_version}_qt4_rc".format(python_version=sys.version_info[0])
    else:
        _resource_suffix = "_qt5_rc"

    logger.debug("Setting PyQt resource suffix to {}".format(_resource_suffix))
    logger.debug("compiling Qt .ui designer files ...")
    compileUiDir(_dir, recurse=True, resource_suffix=_resource_suffix)
    del _dir, _resource_suffix
    _GUI_RESOURCE_BUILT = True

_thread_lock.release()
del _thread_lock
