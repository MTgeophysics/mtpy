# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 07:29:58 2013

@author: Alison Kirkby

plots phase tensor ellipses as a map for a given frequency
"""
import sys
import os
import glob
import matplotlib.pyplot as plt
import mtpy.imaging.plotptmaps as pptmaps

def main(edi_path, save_path=None):
    """Plot Phase Tensor Map
    Args:
        edi_path: path to edi files
        save_path: None  not to save the figure to file; or /tmp/georgina
    Returns:
    """
    # gets edi file names as a list
    elst = glob.glob(os.path.join(edi_path, "*.edi"))

    # frequency to plot
    #plot_freq = 9.4
    plot_freq = 100.0  # check the freq range in your input edi files

    # parameters describing ellipses
    ellipse_dict = {'size': 0.01, 'colorby': 'phimin', 'range': (0, 90, 1), 'cmap': 'mt_bl2gr2rd'}

    # parameters describing the induction vector arrows
    arrow_dict = {'size': 0.2,
                  'lw': 0.01,
                  'head_width': 0.002,
                  'head_length': 0.002,
                  'threshold': 0.8,
                  'direction': 0}

    # parameters describing the arrow legend (should be self explanatory)
    arrow_legend_dict = {'position': 'upper right',
                         'fontpad': 0.0025,
                         'xborderpad': 0.07,
                         'yborderpad': 0.015}

    m = pptmaps.PlotPhaseTensorMaps(fn_list=elst,
                                    plot_freq=plot_freq,
                                    arrow_legend_dict=arrow_legend_dict,
                                    ftol=0.2,
                                    xpad=0.02,
                                    plot_tipper='yr',
                                    arrow_dict=arrow_dict,
                                    #ellipse_dict=ellipse_dict,
                                    fig_size=(4, 4),
                                    mapscale='deg',  #deg or km
                                    save_fn=save_path)

    m.redraw_plot()

    # if save_path is not None:
    #     plt.savefig(save_path, dpi=300)

    return

###################################################################################################
# How to Run:
# cd /path2/mtpy2
# export PYTHONPATH=/path2/mtpy2
# python examples/plot_phase_tensor_map.py ./examples/data/edi_files/georgina ./localdir/mtpy_map
# python examples/plot_phase_tensor_map.py ./examples/data/edi_files ./localdir/mtpy_map
# python examples/plot_phase_tensor_map.py tests/data/edifiles/ ./localdir/mtpy_map
###################################################################################################
if __name__ == '__main__':

    edi_path = sys.argv[1]

    if len(sys.argv) > 2:
        save_file = sys.argv[2]
    else:
        save_file = None

    main(edi_path, save_path=save_file)
