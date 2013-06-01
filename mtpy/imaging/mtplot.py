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
    >>> #create a list of full paths to the edi files
    >>> edilst = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
    >>> ...        if edi.find('.edi')>0]
    >>> #plot apparent resisitivity, phase and induction arrows as individual
    >>> #figures
    >>> rpm = mtplot.plot_multiple_mt_responses(fn_lst=edilst, plot_style='1',
    >>> ...                                     plot_tipper='yr')
    >>> #close all the plots after done looking at them
    >>> plt.close('all')
    >>> #plot phase tensor pseudo section with induction arrows
    >>> pts = mtplot.plot_pt_pseudosection(fn_lst=edilst, 
    >>> ...                                plot_tipper='yr')
    >>> # write out the phase tensor parameter values to files
    >>> pts.writeTextFiles()
    >>> #change coloring scheme to color by skew and a segmented colormap
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
#==============================================================================

from mtpy.imaging.plotnreponses import PlotMultipleResponses as plotnresponses
from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection as plotrpps
from mtpy.imaging.plotpt import PlotPhaseTensor as plotpt
from mtpy.imaging.plotptpseudosection import PlotPhaseTensorPseudoSection as plotptps
from mtpy.imaging.plotptmaps import PlotPhaseTensorMaps as plotptmaps
from mtpy.imaging.plotresponse import PlotResponse as plotresponse
from mtpy.imaging.plotstrike import PlotStrike as plotstrike

#==============================================================================


def plot_mt_response(filename=None, z_array=None, z_err_array=None, 
                     period=None, fignum=1, plotnum=1, title=None, dpi=300, 
                     rot_z=0, plot_yn='y', plot_tipper='n', plot_strike='n',
                     plot_skew='n', tipper_array=None, tipper_err_array=None, 
                     tipper_object=None, res_array=None, res_err_array=None,
                     phase_array=None, phase_err_array=None, plot_pt='n',
                     res_phase_object=None, z_object=None, mt_object=None):
                         
    """
    plots the MT response for a single station.  
    
    """
                         
    kwargs = {'filename' : filename,
              'z_array' : z_array,
              'z_err_array' : z_err_array,
              'period' : period,
              'fignum' : fignum,
              'plotnum' : plotnum,
              'title' : title,
              'dpi' : dpi,
              'rot_z' : rot_z,
              'plot_yn' : plot_yn,
              'plot_tipper' : plot_tipper,
              'plot_strike' : plot_strike,
              'plot_skew' : plot_skew,
              'plot_pt' : plot_pt,
              'tipper_array' : tipper_array,
              'tipper_err_array' : tipper_err_array,
              'tipper_object' : tipper_object,
              'res_array' : res_array,
              'res_err_array' : res_err_array,
              'phase_array' : phase_array,
              'phase_err_array' : phase_err_array,
              'res_phase_object' : res_phase_object,
              'z_object' : z_object,
              'mt_object' : mt_object}
              
    
    return plotresponse(**kwargs)

def plot_multiple_mt_responses(fn_lst=None, res_object_lst=None,
                               z_object_lst=None, tipper_object_lst=None,
                               mt_object_lst=None, fignum=1, dpi=300, rot_z=0,
                               plot_num=1, plot_style='1', plot_yn='y', 
                               plot_tipper='n',plot_strike='n', plot_skew='n',
                               title=None):
    """
    plot multiple MT responses    
    
    """
    
    kwargs = {'fn_lst' : fn_lst,
              'res_object_lst' : res_object_lst,
              'z_object_lst' : z_object_lst,
              'tipper_object_lst' : tipper_object_lst,
              'mt_object_lst' : mt_object_lst,
              'fignum' : fignum,
              'dpi' : dpi,
              'rot_z' : rot_z,
              'plot_num' : plot_num,
              'plot_style' : plot_style,
              'plot_yn' : plot_yn,
              'plot_tipper' : plot_tipper,
              'plot_strike' : plot_strike,
              'plot_skew' : plot_skew,
              'title' : title}
              
    return plotnresponses(**kwargs)

def plot_pt(filename=None, z_object=None, mt_object=None, 
            pt_object=None, fignum=1, dpi=300, rot_z=0, plot_yn='y',
            ellipse_dict=None):
                
    """
    plots the phase tensor ellipses along with the strike, minimum phase,
    maximum phase, skew and ellipticity.    
    
    """
                
    kwargs = {'filename' : filename,
              'z_object' : z_object,
              'mt_object' : mt_object,
              'pt_object' : pt_object,
              'fignum' : fignum,
              'dpi' : dpi,
              'rot_z' : rot_z,
              'plot_yn' : plot_yn,
              'ellipse_dict' : ellipse_dict}
              
    return plotpt(**kwargs)

def plot_pt_pseudosection(fn_lst=None, res_object_lst=None,
                          z_object_lst=None, tipper_object_lst=None, 
                          mt_object_lst=None, ellipse_dict={}, 
                          stretch=(50,25), stationid=(0,4), title=None,
                          cb_dict={}, linedir='ns', fignum=1, rot_z=0, 
                          fig_size=[6,6], dpi=300, plot_tipper='n', 
                          arrow_dict={}, tscale='period', 
                          font_size=7, plot_yn='y', xlim=None, ylim=None):
    """
    plots a pseudo section of phase tensor ellipses for a given profile.    
    
    """
    
    kwargs = {'fn_lst' : fn_lst,
              'res_object_lst' : res_object_lst,
              'z_object_lst' : z_object_lst,
              'tipper_object_lst' : tipper_object_lst,
              'mt_object_lst' : mt_object_lst,
              'ellipse_dict' : ellipse_dict,
              'stretch' : stretch,
              'stationid' : stationid,
              'cb_dict' : cb_dict,
              'linedir' : linedir,
              'arrow_dict' : arrow_dict,
              'rot_z' : rot_z,
              'fig_size' : fig_size,
              'tscale' : tscale,
              'fignum' : fignum,
              'plot_yn' : plot_yn,
              'font_size' : font_size,
              'dpi' : dpi,
              'title' : title,
              'plot_tipper' : plot_tipper,
              'xlim' : xlim,
              'ylim' : ylim}
              
    return plotptps(**kwargs)

def plot_pt_map(fn_lst=None, res_object_lst=None,
                z_object_lst=None, tipper_object_lst=None, mt_object_lst=None,
                plot_freq=1, ellipse_dict={}, cb_dict={},
                arrow_dict={}, xpad=.2, ypad=.2, rot_z=0,
                fig_size=[8,8], station_dict=None, tscale='period',
                mapscale='latlon', fignum=1, image_dict=None, plot_yn='y',
                arrow_legend_dict={}, font_size=7, dpi=300, title=None,
                reference_point=(0,0), plot_tipper='n', ftol=.1):
                    
    """
    plots a map of phase tensor ellipses for a given frequency.
    
    """
    
    kwargs = {'fn_lst' : fn_lst,
              'res_object_lst' : res_object_lst,
              'tipper_object_lst' : tipper_object_lst,
              'mt_object_lst' : mt_object_lst,
              'plot_freq' : plot_freq,
              'ellipse_dict' : ellipse_dict,
              'cb_dict' : cb_dict,
              'arrow_dict' : arrow_dict,
              'xpad' : xpad,
              'ypad' : ypad,
              'rot_z' : rot_z,
              'fig_size' : fig_size,
              'station_dict' : station_dict,
              'tscale' : tscale,
              'mapscale' : mapscale,
              'fignum' : fignum,
              'image_dict' : image_dict,
              'plot_yn' : plot_yn,
              'arrow_legend_dict' : arrow_legend_dict,
              'font_size' : font_size,
              'dpi' : dpi,
              'title' : title,
              'reference_point' : reference_point,
              'plot_tipper' : plot_tipper,
              'ftol' : ftol}
    
    return plotptmaps(**kwargs) 

def plot_strike(fn_lst=None, z_object_lst=None, tipper_object_lst=None,
                mt_object_lst=None, fignum=1, font_size=10, dpi=300, rot_z=0,
                period_tolerance=.05, text_dict={}, plot_range='data',
                plot_type=1, plot_tipper='n', pt_error_floor=None,
                plot_yn='y', fold=True, bin_width=5):
                    
    kwargs = {'fn_lst' : fn_lst,
              'z_object_lst' : z_object_lst,
              'tipper_object_lst' : tipper_object_lst,
              'mt_object_lst' : mt_object_lst,
              'fignum' : fignum,
              'font_size' : font_size,
              'dpi' : dpi,
              'rot_z' : rot_z,
              'period_tolerance' : period_tolerance,
              'text_dict' : text_dict,
              'plot_range' : plot_range,
              'plot_type' : plot_type,
              'plot_tipper' : plot_tipper,
              'pt_error_floor' : pt_error_floor,
              'plot_yn' : plot_yn,
              'fold' : fold,
              'bin_width' : bin_width}
              
    return plotstrike(**kwargs)

def plot_resphase_pseudosection():
    pass
