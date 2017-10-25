import os
import sys
from qtpy.uic import compileUiDir
from qtpy import QT_VERSION
_dir = os.path.dirname(__file__)

if QT_VERSION.startswith("4"):
    _resource_suffix = "_py{python_version}_qt4_rc".format(python_version=sys.version_info[0])
else:
    _resource_suffix = "_qt5_rc"
compileUiDir(_dir, recurse=True, resource_suffix=_resource_suffix)
del _dir, _resource_suffix
