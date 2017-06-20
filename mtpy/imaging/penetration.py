"""
    Description:
        todo: write description

    Usage:
        todo: write usage

    Author: YingzhiGou
    Date: 20/06/2017
"""

from imaging_base import ImagingBase
import numpy as np
import matplotlib.pyplot as plt
from utils.decorator import deprecated

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
from mtpy.utils.mtpylog import MtPyLog
logger = MtPyLog().get_mtpy_logger(__name__)


class Depth1D(ImagingBase):
    def __init__(self, edis = None):
        self._edis = None
        self._fig = None
        self.set_edi(edis)

    def export_image(self):
        raise NotImplemented

    def plot(self, rholist=set('zxy', 'zyx', 'det')):
        if not isinstance(rholist, set):
            rholist = set(rholist)

        self._fig = plt.figure()
        self._fig.grid(True)

        for edi in self._edis:
            zeta = edi.Z            # the attribute Z represent the impedance tensor 2X2 matrix
            freqs = zeta.freq       # frequencies
            scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

            # The periods array
            periods = 1.0 / freqs
            legendh = []

            if 'zxy' in rholist:
                # One of the 4-components: XY
                penetration_depth = scale_param * \
                                    np.sqrt(zeta.resistivity[:, 0, 1] * periods)

                # pen_zxy, = plt.semilogx(periods, -penetration_depth, '-*',label='Zxy')
                pen_zxy, = self._fig.semilogx(
                    periods, -penetration_depth, color='#000000', marker='*', label='Zxy')
                # See
                # http://matplotlib.org/1.3.1/examples/pylab_examples/line_styles.html

                legendh.append(pen_zxy)

            if 'zyx' in rholist:
                penetration_depth = scale_param * \
                                    np.sqrt(zeta.resistivity[:, 1, 0] * periods)

                pen_zyx, = self._fig.semilogx(
                    periods, -penetration_depth, color='g', marker='o', label='Zyx')
                legendh.append(pen_zyx)

            if 'det' in rholist:
                # determinant
                det2 = np.abs(zeta.det[0])
                det_penetration_depth = scale_param * \
                                        np.sqrt(0.2 * periods * det2 * periods)

                # pen_det, = plt.semilogx(periods, -det_penetration_depth, '-^', label='Determinant')
                pen_det, = self._fig.semilogx(
                    periods, -det_penetration_depth, color='b', marker='^', label='Determinant')
                legendh.append(pen_det)

            self._fig.legend(
                handles=legendh,
                bbox_to_anchor=(
                    0.1,
                    0.5),
                loc=3,
                ncol=1,
                borderaxespad=0.)

            self._fig.title("Penetration Depth for file %s" % edi._get_fn())
            self._fig.xlabel("Log Period (seconds)", fontsize=16)
            self._fig.ylabel("Penetration Depth (meters)", fontsize=16)


### ======================== old APIs ==============================

def plot_edi_dir(edi_path, rholist=['zxy', 'zyx', 'det']):
    """ plot edi files from the input directory edi_path
    """
    raise NotImplemented

def plot_edi_file(edifile, rholist=['zxy', 'zyx', 'det'], savefile=None):
    """
    Plot the input edi_file
    Args:
        edi_file: path2edifile
        rholist: a list of the rho to be used.
        savefile: path2savefig, not save if None
    Returns:
    """
    raise NotImplemented

### ======================== main ============================

###############################################################################
# plot one-by-one edi files in a given dirpath
# How to Run:
#    export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH
#    python  mtpy/imaging/penetration_depth1d.py data/edi_files/
# python  mtpy/imaging/penetration_depth1d.py
# tests/data/edifiles/15125A.edi

if __name__ = '__main__':
    import sys, os

    if len(sys.argv) < 2:
        print (
            "\n please provide path to edi files\n USAGE:  %s path2edifile" %
            sys.argv[0])
        sys.exit(1)
    else:
        edi_path = sys.argv[1]

        if os.path.isfile(edi_path):
            plot_edi_file(edi_path, savefile='C:/temp/pen_depth.jpg')
            # rholist can be any of ['zxy','zyx','det'], default all of them
            elif os.path.isdir(edi_path):  # choose a suitable function below at run
            # plot_edi_dir(edi_path )
            plot_edi_dir(edi_path, rholist=['det'])
        else:
            logger.error("Usage %s %s", sys.argv[0], "path2edi")