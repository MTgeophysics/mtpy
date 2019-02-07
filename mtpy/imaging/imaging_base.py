"""
    this file contains the base class of all imaging classes, where the few common routines are defined as well as
    some abstract member functions such as plot() and set_data()
    Author: YingzhiGou
    Date: 20/06/2017
"""

import abc,six
import mtpy.core.mt as mt
# get a logger object for this module, using the utility class MtPyLog to
# config the logger
from mtpy.utils.mtpylog import MtPyLog

@six.add_metaclass(abc.ABCMeta)
class ImagingBase():
    """
    Description:
        This is the base class for all the imaging classes, with standardized API (as abstract methods)
        Also some common functionality should be implemented here.
    """

    def __init__(self):
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)
        self._data = None
        self._fig = None

    @abc.abstractmethod
    def plot(self, **kwargs):
        pass

    def show(self, block=True):
        """
        display the image
        :return:
        """
        if self._fig is None:
            self.plot()
        if block:
            import matplotlib.pyplot as plt
            plt.show()
        else:
            self._fig.show()

    # @abc.abstractmethod
    # def set_param(self, **kwargs):
    #     pass

    def export_image(self, fn, **kwargs):
        if self._fig is None:
            self.plot()
        self._fig.savefig(fn, **kwargs)

    # @abc.abstractmethod
    def set_data(self, data):
        pass

    def _set_edis(self, edis):
        """
        set edi (mt.MT) objects for plotting
        :param edis:
        :return:
        """
        if isinstance(edis, list):
            if all(isinstance(edi, mt.MT) for edi in edis):
                self._data = edis
                self._reset_fig()
            else:
                # todo: raise an exception
                pass
        elif isinstance(edis, mt.MT):
            # if its a single mt object, put it in an list
            self._data = [edis]
            self._reset_fig()
        else:
            # todo raise an exception
            pass

    def _set_edi(self, edi):
        """
        set a single data source as the source for plotting
        :param edi:
        :return:
        """
        if isinstance(edi, mt.MT):
            self._data = edi
            self._reset_fig()
        else:
            # todo raise an exception
            pass

    def close(self):
        """
        close the figure
        :return:
        """
        self._reset_fig()

    def _reset_fig(self):
        if self._fig is not None:
            self._fig.clf()
            self._fig = None

    def get_figure(self):
        if self._fig is None:
            self.plot()
        return self._fig

    def get_data(self):
        return self._data

    # ========================================
    # set properties
    # ========================================
    data = property(set_data, get_data, doc="the data (mt objects) that are to be plotted")
    fig = property(None, get_figure, doc="matplotlib fig object")


class ImagingError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class ParameterError(ImagingError):
    def __init__(self, *args, **kwargs):
        ImagingError.__init__(self, *args, **kwargs)
