import functools
import inspect
import os
from mtpy.utils.mtpylog import MtPyLog


class deprecated(object):
    """
        Description:
            used to mark functions, methods and classes deprecated, and prints warning message when it called
            decorators based on https://stackoverflow.com/a/40301488

        Usage:
            todo: write usage

        Author: YingzhiGou
        Date: 20/06/2017
    """
    def __init__(self, reason):  # pragma: no cover
        if inspect.isclass(reason) or inspect.isfunction(reason):
            raise TypeError("Reason for deprecation must be supplied")
        self.reason = reason

    def __call__(self, cls_or_func):  # pragma: no cover
        if inspect.isfunction(cls_or_func):
            if hasattr(cls_or_func, 'func_code'):
                _code = cls_or_func.__code__
            else:
                _code = cls_or_func.__code__
            fmt = "Call to deprecated function or method {name} ({reason})."
            filename = _code.co_filename
            lineno = _code.co_firstlineno + 1

        elif inspect.isclass(cls_or_func):
            fmt = "Call to deprecated class {name} ({reason})."
            filename = cls_or_func.__module__
            lineno = 1

        else:
            raise TypeError(type(cls_or_func))

        msg = fmt.format(name=cls_or_func.__name__, reason=self.reason)

        @functools.wraps(cls_or_func)
        def new_func(*args, **kwargs):  # pragma: no cover
            import warnings
            warnings.simplefilter('always', DeprecationWarning)  # turn off filter
            warnings.warn_explicit(msg, category=DeprecationWarning, filename=filename, lineno=lineno)
            warnings.simplefilter('default', DeprecationWarning)  # reset filter
            return cls_or_func(*args, **kwargs)

        return new_func


class gdal_data_check(object):
    _has_checked = False
    _gdal_data_found = False
    _logger = MtPyLog.get_mtpy_logger(__name__)

    def __init__(self, func, raise_error=False):
        """
        this decorator should only be used for the function that requres gdal and gdal-data
        to function correctly.

        the decorator will check if the GDAL_DATA is set and the path in GDAL_DATA is exist.
        If GDAL_DATA is not set, then try to use external program "gdal-config --datadir" to
        findout where the data files are installed.

        If failed to find the data file, then ImportError will be raised.

        :param func: function to be decorated
        """
        self._func = func
        if not self._has_checked:
            self._gdal_data_found = self._check_gdal_data()
            self._has_checked = True
        if not self._gdal_data_found:
            if(raise_error):
                raise ImportError("GDAL  is NOT installed correctly")
            else:
                print ("Ignore GDAL as it is not working. Will use pyproj")

    def __call__(self, *args, **kwargs):  # pragma: no cover
        return self._func(*args, **kwargs)

    def _check_gdal_data(self):
        if 'GDAL_DATA' not in os.environ:
            # gdal data not defined, try to define
            from subprocess import Popen, PIPE
            self._logger.warning("GDAL_DATA environment variable is not set  Please see https://trac.osgeo.org/gdal/wiki/FAQInstallationAndBuilding#HowtosetGDAL_DATAvariable ")
            try:
                # try to find out gdal_data path using gdal-config
                self._logger.info("Trying to find gdal-data path ...")
                process = Popen(['gdal-config', '--datadir'], stdout=PIPE)
                (output, err) = process.communicate()
                exit_code = process.wait()
                output = output.strip()
                if exit_code == 0 and os.path.exists(output):
                    os.environ['GDAL_DATA'] = output
                    self._logger.info("Found gdal-data path: {}".format(output))
                    return True
                else:
                    self._logger.error(
                        "\tCannot find gdal-data path. Please find the gdal-data path of your installation and set it to "
                        "\"GDAL_DATA\" environment variable. Please see "
                        "https://trac.osgeo.org/gdal/wiki/FAQInstallationAndBuilding#HowtosetGDAL_DATAvariable for "
                        "more information.")
                    return False
            except Exception:
                return False
        else:
            if os.path.exists(os.environ['GDAL_DATA']):
                self._logger.info("GDAL_DATA is set to: {}".format(os.environ['GDAL_DATA']))

                try:
                    from osgeo import osr
                    from osgeo.ogr import OGRERR_NONE
                except:
                    self._logger.error("Failed to load module osgeo; looks like GDAL is NOT working")
                    # print ("Failed to load module osgeo !!! ")

                    return False
                # end try

                return True
            else:
                self._logger.error("GDAL_DATA is set to: {}, but the path does not exist.".format(os.environ['GDAL_DATA']))
                return False
