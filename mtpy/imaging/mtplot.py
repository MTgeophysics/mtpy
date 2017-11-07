"""
========
mtplot
========

**Provides**

    1. Different plotting options to represent the MT response.
    2. Ability to create text files of the plots for further analysis
    3. Class object that contains all the important information for an MT
       station.
============================= =================================================
Functions                      Description
============================= =================================================
plot_mt_response              plots resistivity and phase for a single station
                              Options include tipper, strike and skew.
plot_multiple_mt_responses    plots multiple stations at once with options
                              of plotting in single figure, all in one
                              figure as subplots or all in one plot for
                              direct comparison.
plot_pt                       plots the phase tensor ellipses and parameters
                              in one plot including strike angle, minimum
                              and maximum phase, skew angle and ellipticity
plot_pt_pseudosection         plots a pseudo section of phase
                              tensor ellipses assuming the
                              stations are along a profile line.
                              Options to plot induction arrows.
plot_mt_map                   plots phase tensor ellipses in map view for
                              a single frequency.  Options to plot
                              induction arrows.
plot_strike                   plots strike angle estimated from the
                              invariants of the impedance tensor defined
                              by Weaver et al. [2000,2003], strike angle
                              from the phase tensor and option to plot
                              strike estimated from the induction arrows.
plot_residual_pt_maps         plots the residual phase tensor between two
                              surveys in map view.
plot_residual_pt_ps           plots the residual phase tensor between two
                              surveys as a pseudo section.
============================= =================================================

All plot function return plot classes where the important properties are made
attributes which can be manipulated by the user.  All classes have been
written with the basic input being edi files.  This was assumed to be the
standard MT response file, but turns out to be not as widely used as thought.
So the inputs can be other arrays and class objects (see MTplot doc string for
details).  If you have a data file format you can create a class using the
objects in mtpy.core to create an input, otherwise contact us and we can try
to build something.

A typical use might be loading in all the .edi files in and plotting them in
different modes, like apparent resistivity and phase, phase tensor pseudo
section and strike angle.

:Example: ::

    >>> import mtpy.imaging.mtplot as mtplot
    >>> import os
    >>> import matplotlib.pyplot as plt
    >>> edipath = r"/home/MT/EDIfiles"
    >>> #--> create a list of full paths to the edi files
    >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
    >>> ...        if edi.find('.edi')>0]
    >>> #--> plot apparent resisitivity, phase and induction arrows
    >>> rpm = mtplot.plot_multiple_mt_responses(fn_lst=edilst, plot_style='1',
    >>> ...                                     plot_tipper='yr')
    >>> #--> close all the plots after done looking at them
    >>> plt.close('all')
    >>> #--> plot phase tensor pseudo section with induction arrows
    >>> pts = mtplot.plot_pt_pseudosection(fn_lst=edilst,
    >>> ...                                plot_tipper='yr')
    >>> #--> write out the phase tensor parameter values to files
    >>> pts.export_pt_params_to_file()
    >>> #--> change coloring scheme to color by skew and a segmented colormap
    >>> pts.ellipse_colorby = 'skew_seg'
    >>> pts.ellipse_cmap = 'mt_seg_bl2wh2rd'
    >>> pts.ellipse_range = (-9, 9, 3)
    >>> pts.redraw_plot()

:Authors:
    Lars Krieger,
    Jared Peacock, and
    Kent Invariarty


:Version: 0.0.1 of 2013

"""
# ==============================================================================
import mtpy.imaging.plotnresponses as plotnresponses
import mtpy.imaging.plotpseudosection as plotrpps
import mtpy.imaging.plotpt as plotpt
import mtpy.imaging.phase_tensor_maps as plotptmaps
import mtpy.imaging.phase_tensor_pseudosection as plotptps
import mtpy.imaging.plotresidualptmaps as plotresidualptmaps
import mtpy.imaging.plotresidualptps as plotresidualptps
import mtpy.imaging.plot_mt_response as plotresponse
import mtpy.imaging.plotstations as plotstations
import mtpy.imaging.plotstrike as plotstrike

# reload(plotpt)
# reload(plotptmaps)

#==============================================================================



def plot_mt_response(**kwargs):
    """
    plots the MT response for a single station.

    """

    return plotresponse.PlotMTResponse(**kwargs)

def plot_multiple_mt_responses(**kwargs):
    """
    plot multiple MT responses

    """

    return plotnresponses.PlotMultipleResponses(**kwargs)

def plot_pt(**kwargs):
    """
    plots the phase tensor ellipses along with the strike, minimum phase,
    maximum phase, skew and ellipticity.

    """

    return plotpt.PlotPhaseTensor(**kwargs)

def plot_pt_pseudosection(**kwargs):
    """
    plots the phase tensor ellipses as a pseudo section.

    """

    return plotptps.PlotPhaseTensorPseudoSection(**kwargs)

def plot_pt_map(**kwargs):

    """
    plots a map of phase tensor ellipses for a given frequency.

    """

    return plotptmaps.PlotPhaseTensorMaps(**kwargs)

def plot_strike(**kwargs):
    """
    plots the strike angle.

    """

    return plotstrike.PlotStrike(**kwargs)

def plot_resphase_pseudosection(**kwargs):
    """
    plots resistivity and phase as a pseudo section

    """

    return plotrpps.PlotResPhasePseudoSection(**kwargs)

def plot_station_locations(**kwargs):
    """
    Plot station locations in map view.

    """

    return plotstations.PlotStations(**kwargs)

def plot_residual_pt_maps(fn_list1, fn_list2, **kwargs):
    """
    plot residual pt between two measurements in map view

    """

    return plotresidualptmaps.PlotResidualPTMaps(fn_list1, fn_list2, **kwargs)

def plot_residual_pt_ps(fn_list1, fn_list2, **kwargs):
    """
    plot residual ps between two measurements as a pseudo section

    """

    return plotresidualptps.PlotResidualPTps(fn_list1, fn_list2, **kwargs)


# reset the doc strings of these helper functions to that of the class
# there is probably a more elegant way to do this, but for now, this
# works
plot_mt_response.__doc__ = plotresponse.PlotMTResponse.__doc__
plot_multiple_mt_responses.__doc__ = \
                                plotnresponses.PlotMultipleResponses.__doc__
plot_pt.__doc__ = plotpt.PlotPhaseTensor.__doc__
plot_pt_pseudosection.__doc__ = plotptps.PlotPhaseTensorPseudoSection.__doc__
plot_pt_map.__doc__ = plotptmaps.PlotPhaseTensorMaps.__doc__
plot_strike.__doc__ = plotstrike.PlotStrike.__doc__
plot_resphase_pseudosection.__doc__ = \
                                    plotrpps.PlotResPhasePseudoSection.__doc__
plot_station_locations.__doc__ = plotstations.PlotStations.__doc__
plot_residual_pt_maps.__doc__ = plotresidualptmaps.PlotResidualPTMaps.__doc__
plot_residual_pt_ps.__doc__ = plotresidualptps.PlotResidualPTps.__doc__
