"""
plot phase tensor map and save to file
modified python mtpy/imaging/plot_phase_tensor_maps_with_tipper.py
"""

import os, sys
import glob

#from mtpy.tests.common import MTPY_DEVELOPMENT_TEST_DATA
#from mtpy.imaging.plot_phase_tensor_maps_with_tipper import  PlotPhaseTensorMaps
import mtpy.imaging.plotptmaps as pptmaps
###############################################################
# How2Run:

#   export PYTHONPATH=/Softlab/Githubz/mtpy2
#    python plot1_phase_tensor_map.py data/edi2/

#   python plot1_phase_tensor_map.py /Softlab/Githubz/mtpy1/examples/data/edi_files/georgina
#   python plot1_phase_tensor_map.py /Softlab/Githubz/mtpy1/examples/data/edi_files
# compare to  
# /Softlab/Githubz/mtpy2$  python examples/plot_phase_tensor_map.py examples/data/edi_files/georgina
# /Softlab/Githubz/mtpy2$  python examples/plot_phase_tensor_map.py examples/data/edi_files/

#  python plot1_phase_tensor_map.py /Softlab/Data/MT_datasets/GA_UA_edited_10s-10000s/ 0.1

#--------------------------------------------------------------
if __name__ == "__main__":


    # edi_file_list = glob.glob(os.path.join(MTPY_DEVELOPMENT_TEST_DATA,'*.edi'))

    edi_path = sys.argv[1]

    edi_file_list = glob.glob(os.path.join(edi_path, '*.edi'))

    print(edi_file_list)

    freq=0.0625  #10.0
    if len(sys.argv)>2:
        freq=float(sys.argv[2])

    pt1 = pptmaps.PlotPhaseTensorMaps(fn_list=edi_file_list,
                                    plot_freq=freq,) #, plot_freq=250, plot_tipper='yri')
    
    
    # pt1.arrow.arrow_size = .2
    # pt1.arrow.arrow_head_width = 0.00125
    # pt1.arrow.arrow_head_length = 0.0025
    # pt1.arrow.arrow_lw = .01
    # pt1.redraw_plot()


# map scale km
    pt1 = pptmaps.PlotPhaseTensorMaps(
        fn_list=edi_file_list,
        plot_freq=freq,
        plot_tipper='yri',
        mapscale='km',
        plot_title='??Customisable Title?? '
    )

    
    
    # pt1.arrow.arrow_size = .2
    # pt1.arrow.arrow_head_width = 0.00125
    # pt1.arrow.arrow_head_length = 0.0025
    # pt1.arrow.arrow_lw = .01
    # pt1.redraw_plot()


    pt1.save_figure(
        save_file_name_or_dir=edi_path,
        file_name_if_dir_provided='pt_map_with_tipper',
        file_format='png')  #pdf
