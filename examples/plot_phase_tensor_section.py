"""
Plot phase tensor ellipses as a pseudo section (distance along profile vs period)

"""

import os
import sys
import glob

from mtpy.imaging.plotptpseudosection import PlotPhaseTensorPseudoSection


def main(edi_path, save_file=None):
    """
    Plot Phase Tensor Pseudo Section
    :param edi_path: path2edi dir
    :param save_file: save file of the plot
    :return:
    """

    edi_files = glob.glob(os.path.join(edi_path, "*.edi"))

    ptpObj = PlotPhaseTensorPseudoSection(fn_list=edi_files,
                                     tscale='period',
                                     #ylim=(1e-1, 1e3),  # orig period range to plot
                                     ylim=(0, 10000),  # period range to plot
                                     # xlim = (0,10000),
                                     stretch=(2000, 40),  # determines (x,y) aspect ratio of plot
                                     station_id=(0, 10),  # indices for showing station names
                                     ellipse_dict={'size': 6},
                                     plot_tipper='yri',
                                     arrow_dict={'size': 5, 'head_length': 0.2,
                                                 'head_width': 0.1, 'lw': 0.5},
                                     # arrow parameters, adjust as necessary. lw = linewidth
                                     font_size=4,
                                     dpi=300)



    ptpObj.plot()

    return


###################################################################################################
# How to Run:
# cd path2/mtpy2
# export PYTHONPATH=/path2/mtpy2   # the full path to your repo dir: mtpy2
# python examples/plot_phase_tensor_section.py ./examples/data/edi_files/georgina
# python examples/plot_phase_tensor_section.py ./examples/data/edi_files
# python examples/plot_phase_tensor_section.py ./tests/data/edifiles
##################################################################################################
if __name__ == '__main__':
    edi_path = sys.argv[1]
    main(edi_path)
