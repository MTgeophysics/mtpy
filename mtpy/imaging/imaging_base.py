"""
    Description:
        todo: write description

    Usage:
        todo: write usage

    Author: YingzhiGou
    Date: 20/06/2017
"""

from abc import ABCMeta
import mtpy.core.mt as mt


class ImagingBase:
    __metaclass__ = ABCMeta

    @abstractmethod
    def plot(self, **kwargs):
        pass

    def show(self):
        """
        display the image
        :return:
        """
        if self._fig == None:
            self.plot()
        self._fig.show()


    @abstractmethod
    def export_image(self, fn, dpi=800):
        if self.fig == None:
            self.plot()
        self._fig.savefile(fn, dpi)

    @abstractmethod
    def set_edi(self, edis):
        if isinstance(edis, list):
            if all(isinstance(edi, mt.MT) for edi in edis):
                self._edis = edis
            else:
                # todo: raise an exception
                pass
        elif isinstance(edis, mt.MT):
            # if its a single mt object, put it in an list
            self._edis = [edis]
        else:
            # todo raise an exception
            pass

