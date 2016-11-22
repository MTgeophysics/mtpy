# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots phase tensor ellipses as a pseudo section (distance along profile vs period)
"""
import sys
import os
import glob

import mtpy.imaging.plotptpseudosection as ptp

def main(edi_path):
    """Plot Phase Tensor Pseudo Section
    Args:
        edi_path: path to edi files

    Returns:
    """

    elst = glob.glob(os.path.join(edi_path, "*.edi"))

    ptp.PlotPhaseTensorPseudoSection(fn_list=elst,
                                     tscale='period',
                                     ylim=(1e-1, 1e3),  # period range to plot
                                     # xlim = (0,10000),
                                     stretch=(50, 16),  # determines (x,y) aspect ratio of plot
                                     station_id=(0, 10),  # indices for showing station names
                                     ellipse_dict={'size': 3},
                                     plot_tipper='yr',
                                     arrow_dict={'size': 3, 'head_length': 0.1,
                                                 'head_width': 0.1, 'lw': 0.5},
                                     # arrow parameters, adjust as necessary. lw = linewidth
                                     font_size=4,
                                     dpi=300)

    return


#########################################################
# How to Run:
# export PYTHONPATH=/path2/mtpy2   # the full path to your repo dir: mtpy2
# python examples/plot_phase_tensor_section.py /g/data/ha3/fxz547/mtpy2/examples/data/edi_files/georgin
#########################################################
if __name__ == '__main__':
    edi_path = sys.argv[1]
    main(edi_path)
