import os
import sys
from qtpy.uic import compileUiDir

_dir = os.path.dirname(__file__)
_resource_suffix = "_py{python_version}_rc".format(python_version=sys.version_info[0])
compileUiDir(_dir, recurse=True, resource_suffix=_resource_suffix)
del _dir, _resource_suffix
