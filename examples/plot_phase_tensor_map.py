"""
plots phase tensor ellipses as a map for a given frequency.
write to an image file in [jpg, png]-formats if OutputDir/File is specified

"""

import glob
import os
import sys

from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps


def plot_pt(edi_file_list, freq, save_path=None):
    """
    Plot Phase Tensor Map in Lat-Long (unprojected) Coordinate system
    :param edi_file_list: a list of edi files
    :param freq: the frequency value (within a certain tolerance ftol)
    :param save_path: path to file/dir to save the plot as jpg or png file
    :return:
    """

    # 1) Define plots params
    # parameters describing ellipses, differ for different map scales: deg, m, km
    # Try different size to find a suitable value for your case. as a
    # guidance: 1 degree=100KM
    ellipse_dict = {
        'size': 0.2,
        'colorby': 'phimin',
        'range': (
            0,
            90,
            1),
        'cmap': 'mt_bl2gr2rd'}

    # adjust to suitable size: parameters describing the induction vector arrows
    arrow_dict = {'size': 0.5,
                  'lw': 0.2,
                  'head_width': 0.04,
                  'head_length': 0.04,
                  'threshold': 0.8,
                  'direction': 0}

    # parameters describing the arrow legend (not necessarily used)
    # arrow_legend_dict = {'position': 'upper right',
    #                      'fontpad': 0.0025,
    #                      'xborderpad': 0.07,
    #                      'yborderpad': 0.015}

    # 2) Construct plotting object
    pt_obj = PlotPhaseTensorMaps(fn_list=edi_file_list,
                                 plot_freq=freq,
                                 ftol=0.2,  # fre tolerance considered as equal
                                 mapscale='deg',  # deg or m, or km
                                 xpad=0.4,  # plot margin; change according to lat-lon in edifiles
                                 ypad=0.4,  # ~ 2* ellipse size
                                 ellipse_dict=ellipse_dict,
                                 plot_tipper='yr',
                                 arrow_dict=arrow_dict,
                                 # arrow_legend_dict=arrow_legend_dict,
                                 # fig_spython examples/plot_phase_tensor_map.py tests/data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.pngize=(6, 5),
                                 # fig_dpi=300, the default is OK. Higher dpi
                                 # may distort figure
                                 save_fn=save_path)

    # 3) do the plot and save figure - if the param save_path provided
    path2figure = pt_obj.plot(save_path=save_path)

    pt_obj.export_pt_params_to_file(save_path=save_path)

    print ("Please check your output figure: %s" % path2figure)

    return


def plot_pt_utm(edi_file_list, freq, save_path=None):
    """
    UnboundLocalError: local variable 'zone1' referenced before assignment for GA_UA latest survey

    :param edi_file_list:
    :param freq:
    :param save_path:
    :return:
    """

    ellipse_dict = {
        'size': 20,
        'colorby': 'phimin',
        'range': (
            0,
            90,
            1),
        'cmap': 'mt_bl2gr2rd'}

    pt2 = PlotPhaseTensorMaps(
        fn_list=edi_file_list,
        plot_freq=freq,
        save_fn=save_path,
        xpad=20,  # plot margin; change according to your dataset and mapscale km/m
        ypad=20,
        mapscale='km',  # mapscale='m' can cause big numbers in ticks labels
        ellipse_dict=ellipse_dict,
        # plot_tipper='yri',
        # plot_title='??Customisable Title?? '
    )

    pt2.plot()

    if save_path is not None:
        pt2.save_figure(save_path, fig_dpi=300, file_format='jpg')  # pdf)


##########################################################################
# How to Run:
# cd /path2/mtpy2;  export PYTHONPATH=/path2/mtpy2
# python examples/plot_phase_tensor_map.py  examples/data/edi_files/georgina 10 /e/MTPY2_Outputs/
# python examples/plot_phase_tensor_map.py  examples/data/edi_files 10 /e/MTPY2_Outputs/
# python examples/plot_phase_tensor_map.py  tests/data/edifiles/ 10 /e/MTPY2_Outputs/
# python examples/plot_phase_tensor_map.py  E:/Datasets/MT_Datasets/GA_UA_edited_10s-10000s 0.0625 E:/MTPY2_Outputs
##########################################################################
if __name__ == '__main__':

    # the MT edi dir
    edi_path = sys.argv[1]
    # get edi file names as a list
    edifiles = glob.glob(os.path.join(edi_path, "*.edi"))
    if edifiles is None or (edifiles) == 0:
        print("No EDI files found!!!  Please check the input dir")
        sys.exit(1)

    # the MT frequency
    # check the freq range in your input edi files: 10 for georgina tests/data/edifiles
    # freq = 0.0625
    if len(sys.argv) > 2:
        input_freq = float(sys.argv[2])

    if len(sys.argv) > 3:
        save_file = sys.argv[3]
    else:
        save_file = None

    plot_pt(edifiles, input_freq, save_path=save_file)

    # KM map scale may have bugs which make the MT stations too close:
    # E:/Datasets/MT_Datasets/GA_UA_edited_10s-10000s 0.0625
    # plot_pt_utm(edifiles, freq, save_path=save_file)
