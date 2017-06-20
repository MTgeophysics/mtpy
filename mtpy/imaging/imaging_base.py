"""
    Description:
        todo: write description

    Usage:
        todo: write usage

    Author: YingzhiGou
    Date: 20/06/2017
"""

import abc
import mtpy.core.mt as mt


class ImagingBase:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def plot(self, **kwargs):
        pass

    def show(self, block=True):
        """
        display the image
        :return:
        """
        if self._fig == None:
            self.plot()
        if block:
            import matplotlib.pyplot as plt
            plt.show()
        else:
            self._fig.show()


    def export_image(self, fn, dpi=800):
        if self._fig == None:
            self.plot()
        self._fig.savefile(fn, dpi)

    @abc.abstractmethod
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
            else:
                # todo: raise an exception
                pass
        elif isinstance(edis, mt.MT):
            # if its a single mt object, put it in an list
            self._data = [edis]
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
            # if its a single mt object, put it in an list
            self._data = edi
        else:
            # todo raise an exception
            pass


    def clear(self):
        """
        clear the figure
        :return:
        """
        # todo decide if this is necessory
        raise NotImplemented
