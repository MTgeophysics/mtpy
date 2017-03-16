"""
plots phase tensor ellipses as a map for a given frequency.
write to an image file in [jpg, png]-formats if OutputDir/File is specified

"""

import os
import sys
import glob

from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps


def plot_pt(edi_file_list, freq, save_path=None):
    """Plot Phase Tensor Map in lat-lon coordinate
    Args:
        edi_file_list: a list of edi files
        save_path: None not to save.
    Returns:
    """

    # parameters describing ellipses, differ for different map scales: deg, m, km
    ellipse_dict = {'size': 0.03, 'colorby': 'phimin', 'range': (0, 90, 1), 'cmap': 'mt_bl2gr2rd'}

    # parameters describing the induction vector arrows
    arrow_dict = {'size': 0.05,
                  'lw': 0.005,
                  'head_width': 0.002,
                  'head_length': 0.002,
                  'threshold': 0.8,
                  'direction': 0}

    # parameters describing the arrow legend (should be self explanatory)
    arrow_legend_dict = {'position': 'upper right',
                         'fontpad': 0.0025,
                         'xborderpad': 0.07,
                         'yborderpad': 0.015}

    ptObj = PlotPhaseTensorMaps(fn_list=edi_file_list,
                                    plot_freq=freq,
                                    # arrow_legend_dict=arrow_legend_dict,
                                    ftol=0.2,
                                    xpad=0.06,  # plot margin; change according to lat-lon of your dataset
                                    ypad=0.06,
                                    plot_tipper='yri',
                                    arrow_dict=arrow_dict,
                                    ellipse_dict=ellipse_dict,
                                    #fig_size=(6, 5),
                                    mapscale='deg',  # deg or m, or km
                                    save_fn=save_path)  #fig_dpi=300 OK)

    # do the plot and save figure is save_path provided
    ptObj.plot(save_path=save_path)

    # if save_path is not None:
    #     ptObj.save_figure(save_path, fig_dpi=300)

    return

def plot_pt_utm(edi_file_list, freq, save_path=None):
    """
    UnboundLocalError: local variable 'zone1' referenced before assignment for GA_UA latest survey

    :param edi_file_list:
    :param freq:
    :param save_path:
    :return:
    """

    ellipse_dict = {'size': 20, 'colorby': 'phimin', 'range': (0, 90, 1), 'cmap': 'mt_bl2gr2rd'}

    pt2 = PlotPhaseTensorMaps(
        fn_list=edi_file_list,
        plot_freq=freq,
        save_fn=save_path,
        xpad=20,  # plot margin; change according to your dataset and mapscale km/m
        ypad=20,
        mapscale='km',   #mapscale='m' can cause big numbers in ticks labels
        ellipse_dict=ellipse_dict,
        # plot_tipper='yri',
        #plot_title='??Customisable Title?? '
    )

    pt2.plot()

    if save_path is not None:
        pt2.save_figure(save_path, fig_dpi=300, file_format='jpg')  #pdf)


###################################################################################################
# How to Run:
# cd /path2/mtpy2
# export PYTHONPATH=/path2/mtpy2
# python examples/modem_phase_tensor_maps.py ./examples/data/edi_files/georgina 10 /e/MTPY2_Outputs/
# python examples/modem_phase_tensor_maps.py ./examples/data/edi_files 10 /e/MTPY2_Outputs/
# python examples/modem_phase_tensor_maps.py tests/data/edifiles/ 10 /e/MTPY2_Outputs/
# python examples/modem_phase_tensor_maps.py E:/Datasets/MT_Datasets/GA_UA_edited_10s-10000s 0.0625 E:/MTPY2_Outputs
###################################################################################################
if __name__ == '__main__':

    # the MT edi dir
    edi_path = sys.argv[1]
    # get edi file names as a list
    edifiles = glob.glob(os.path.join(edi_path, "*.edi"))
    if edifiles is None or (edifiles)==0:
        raise Exception("No EDI files found!!!  Please check the input dir")
        sys.exit(1)


    # the MT frequency
    # check the freq range in your input edi files: 10 for georgina tests/data/edifiles
    freq=0.0625
    if len(sys.argv)>2:
        freq=float(sys.argv[2])

    if len(sys.argv) > 3:
        save_file = sys.argv[3]
    else:
        save_file = None


    #plot_pt(edifiles, freq, save_path=save_file)

    # KM map scale may have bugs which make the MT stations too close: E:/Datasets/MT_Datasets/GA_UA_edited_10s-10000s 0.0625
    plot_pt_utm(edifiles, freq, save_path=save_file)

