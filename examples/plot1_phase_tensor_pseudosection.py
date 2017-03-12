"""
extension from the mtpy1/mtpy/imaging/phase_tensor_pseudo_section_plot.py

"""
import os, sys
import glob

import mtpy.imaging.mtplot as mtplot
from mtpy.imaging.mtplottools import MTArrows, MTEllipse
from mtpy.tests.common import MTPY_DEVELOPMENT_TEST_DATA

from mtpy.imaging.phase_tensor_pseudo_section_plot import PlotPhaseTensorPseudoSection
import matplotlib.pyplot as plt


# Define the ellipse and arrow properties
ellipse_dict = {'range': (20, 70), 'cmap': 'mt_bl2gr2rd',
                'colorby': 'phimin', 'size': 10}
ellipse = MTEllipse(ellipse_dict=ellipse_dict)
arrow = MTArrows({'size': 60, 'head_length': 4})

def test_plot1(edi_file_list):
    """ from  mtpy1/mtpy/imaging/phase_tensor_pseudo_section_plot.py
    """

    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=1,  # plot an ellipse at every frequency value
        plot_tipper='yri',
        stretch=(1500, 35),
        scale_arrow=False
        )

    plt.rcdefaults()

# Why the plot below becomes smaller?
    pt1.plot_tipper = 'yri'
    pt1.ellipse_freq = 2  # plot ellipse at every second frequency
    pt1.arrow.arrow_size = 200
    pt1.redraw_plot()

    # plt.show()

    pt1.arrow.arrow_size = 100
    edi_file_list = glob.glob(os.path.join(MTPY_DEVELOPMENT_TEST_DATA, '*.edi'))
    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=3,  # plot every 3rd ellipse
        plot_tipper='yri',
        stretch=(1500, 35),
        scale_arrow=False
        )

    return pt1

def test_plot2(edi_file_list):
    """Adapted from demo_scripts/imaging/phase_tensor_pseudo_section_plot.ipynb
    """

    # Define the ellipse and arrow properties
    # ellipse_dict = {'range': (20, 70), 'cmap': 'mt_bl2gr2rd',
    #                 'colorby': 'phimin', 'size': 10}
    # ellipse = MTEllipse(ellipse_dict=ellipse_dict)
    # arrow = MTArrows({'size': 60, 'head_length': 4})

    # Plot the phase tensor pseudo section
    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=1,  # plot an ellipse at every frequency value
        plot_tipper='yri',
        stretch=(1500, 35),
        scale_arrow=False,
        fig_size=[10, 12]
        )


    # Change some properties and replot.  But why the figure window become smaller?
    pt1.ellipse_freq = 2  # plot ellipse at every second frequency
    pt1.arrow.arrow_size = 50  # change arrow size to 50
    pt1.fig_size=[10, 12]
    pt1.font_size=14
    pt1.redraw_plot()

    return pt1

def test_plot3(edi_file_list):
    """Adapted from demo_scripts/imaging/phase_tensor_pseudo_section_plot.ipynb
    """
    arrow.arrow_size = 100
    arrow.arrow_head_width = 3
    arrow.arrow_head_length = 4

    # plot every 3rd ellipse
    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=5,  #=3 plot every 3rd ellipse
        plot_tipper='yri', # plot real and imaginary tipper arrows
        stretch=(1500, 35),
        scale_arrow=False,
        fig_size=[10, 12],
        #font_size=14,
        )

    return pt1

def test_plot4(edi_file_list):
    """ colorby different
    """

    # Colorby 'skew'
    ellipse.ellipse_colorby = 'skew'
    arrow.arrow_size = 50
    # plot every 3rd ellipse
    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=3,  # plot every 3rd ellipse
        plot_tipper='yri', # plot real and imaginary tipper arrows
        stretch=(1500, 35),
        scale_arrow=False,
        fig_size=[10, 12],
        #font_size=14   # too big?
       )


    # Colorby 'normalized_skew'
    ellipse.ellipse_colorby = 'normalized_skew'
    # change arrow size
    arrow.arrow_size = 40
    # plot every 4th ellipse
    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=4,  # plot every 4th ellipse
        plot_tipper='yri', # plot real and imaginary tipper arrows
        stretch=(1500, 35),
        scale_arrow=False,
        fig_size=[10, 12],
        #font_size=14
        )


    # Colorby 'ellipticity'
    ellipse.ellipse_colorby = 'ellipticity'

    # plot every 4th ellipse
    pt1 = PlotPhaseTensorPseudoSection(
        edi_file_list,
        data_type='z',
        ellipse=ellipse,
        arrow=arrow,
        fn_list=edi_file_list,
        tscale='frequency',
        ellipse_freq=4,  # plot every 4th ellipse
        plot_tipper='yri', # plot real and imaginary tipper arrows
        stretch=(1500, 35),
        scale_arrow=False,
        fig_size=[10, 12],
        #font_size=14
        )

    return pt1


###################################################
#1) set env variables PYTHONPATH and MTPYPATH before starting jupyther notebook,
# OR alternatively
#2) do the following two statements. Then it will all works fine.
# sys.path.insert(0,'/Softlab/Githubz/mtpy1')
# os.environ['MTPYPATH']='/Softlab/Githubz/mtpy1'
# 
# python testz/test_plot_phase_tensor_pseudosection.py MTPy_development/test_data/
# python testz/test_plot_phase_tensor_pseudosection.py examples/data/edi_files/
# python testz/test_plot_phase_tensor_pseudosection.py examples/data/edi_files/georgina
#
# compare to Alison script examples/plot_phase_tensor_section.py examples/data/edi_files/georgina
# They all have changing figure size.
#----------------------------------------------------
if __name__ == '__main__':
    """the script commandline run entry point.
    """

    #edi_file_list = glob.glob(os.path.join(MTPY_DEVELOPMENT_TEST_DATA, '*.edi'))
    
    edi_path = sys.argv[1]

    edi_file_list = glob.glob(os.path.join(edi_path, '*.edi'))

    print(edi_file_list)

    print("test_plot1.................")
    test_plot1(edi_file_list)

    print("test_plot2.................")
    test_plot2(edi_file_list)

    print("test_plot3.................")
    test_plot3(edi_file_list)

    print("test_plot4.................")
    test_plot4(edi_file_list)    