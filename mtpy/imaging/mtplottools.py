# -*- coding: utf-8 -*-
"""
===========
mtplottools
===========

Contains helper functions and classes for plotting



@author: jpeacock-pr
"""
# ==============================================================================

import numpy as np
# import mtpy.core.mt
from mtpy.core import mt
import mtpy.core.z as mtz
import mtpy.utils.exceptions as mtex
import mtpy.utils.gis_tools as gis_tools
import matplotlib.mlab as mlab

# ==============================================================================


# define text formating for plotting

ckdict = {'phiminang': r'$\Phi_{min}$ (deg)',
          'phimin': r'$\Phi_{min}$ (deg)',
          'phimaxang': r'$\Phi_{max}$ (deg)',
          'phimax': r'$\Phi_{max}$ (deg)',
          'phidet': r'Det{$\Phi$} (deg)',
          'skew': r'Skew (deg)',
          'normalized_skew': r'Normalized Skew (deg)',
          'ellipticity': r'Ellipticity',
          'skew_seg': r'Skew (deg)',
          'normalized_skew_seg': r'Normalized Skew (deg)',
          'geometric_mean': r'$\sqrt{\Phi_{min} \cdot \Phi_{max}}$',
          'strike': r'Azimuth (deg)',
          'azimuth': r'Azimuth (deg)'}

labeldict = dict([(ii, '$10^{' + str(ii) + '}$') for ii in range(-20, 21)])


# ==============================================================================
# Arrows properties for induction vectors               
# ==============================================================================
class MTArrows(object):
    """
    Helper class to read a dictionary of arrow properties
    
    Arguments:
    -----------
        **arrow_dict** : dictionary for arrow properties
                        * 'size' : float
                                  multiplier to scale the arrow. *default* is 5
                        * 'head_length' : float
                                         length of the arrow head *default* is 
                                         1.5
                        * 'head_width' : float
                                        width of the arrow head *default* is 
                                        1.5
                        * 'lw' : float
                                line width of the arrow *default* is .5
                                
                        * 'color' : tuple (real, imaginary)
                                   color of the arrows for real and imaginary
                                   
                        * 'threshold': float
                                      threshold of which any arrow larger than
                                      this number will not be plotted, helps 
                                      clean up if the data is not good. 
                                      *default* is 1, note this is before 
                                      scaling by 'size'
                                      
                        * 'direction : [ 0 | 1 ]
                                     - 0 for arrows to point toward a conductor
                                     - 1 for arrow to point away from conductor
    
    Attributes:
    -----------
    
        -arrow_color_imag     color of imaginary induction arrow
        -arrow_color_real     color of real induction arrow
        -arrow_direction      convention of arrows pointing to or away from 
                              conductors, see above.
        -arrow_head_length    length of arrow head in relative points
        -arrow_head_width     width of arrow head in relative points
        -arrow_lw             line width of arrows
        -arrow_size           scaling factor to multiple arrows by to be visible
        -arrow_threshold      threshold for plotting arrows, anything above 
                              this number will not be plotted.
                              
    """

    def __init__(self, **kwargs):
        super(MTArrows, self).__init__(**kwargs)

        self.arrow_size = 2.5
        self.arrow_head_length = .15 * self.arrow_size
        self.arrow_head_width = .1 * self.arrow_size
        self.arrow_lw = .5 * self.arrow_size
        self.arrow_threshold = 2
        self.arrow_color_imag = 'b'
        self.arrow_color_real = 'k'
        self.arrow_direction = 0

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if(v in list(kwargs.keys())):
                setattr(self, v, kwargs.pop(v, None))

    def _read_arrow_dict(self, arrow_dict):

        for key in list(arrow_dict.keys()):
            setattr(self, key, arrow_dict[key])


# ==============================================================================
#  ellipse properties           
# ==============================================================================
class MTEllipse(object):
    """
    helper class for getting ellipse properties from an input dictionary
    
    Arguments:
    -------------
        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor 
                          ellipses with keys:
                              
                          * 'size' -> size of ellipse in points 
                                     *default* is .25
                          
                          * 'colorby' : [ 'phimin' | 'phimax' | 'beta' | 
                                    'skew_seg' | 'phidet' | 'ellipticity' ]
                                    
                                    - 'phimin' -> colors by minimum phase
                                    - 'phimax' -> colors by maximum phase
                                    - 'skew' -> colors by skew
                                    - 'skew_seg' -> colors by skew in 
                                                   discrete segments 
                                                   defined by the range
                                    - 'normalized_skew' -> colors by 
                                                    normalized_skew
                                                    see Booker, 2014
                                    - 'normalized_skew_seg' -> colors by 
                                                   normalized_skew
                                                   discrete segments 
                                                   defined by the range
                                    - 'phidet' -> colors by determinant of
                                                 the phase tensor
                                    - 'ellipticity' -> colors by ellipticity
                                    *default* is 'phimin'
                            
                          * 'range' : tuple (min, max, step)
                                     Need to input at least the min and max
                                     and if using 'skew_seg' to plot
                                     discrete values input step as well
                                     *default* depends on 'colorby'
                                     
                          * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 
                                       'mt_wh2bl' | 'mt_rd2bl' | 
                                       'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' | 
                                       'mt_rd2gr2bl']
                                      
                                   - 'mt_yl2rd'       --> yellow to red
                                   - 'mt_bl2yl2rd'    --> blue to yellow to red
                                   - 'mt_wh2bl'       --> white to blue
                                   - 'mt_rd2bl'       --> red to blue
                                   - 'mt_bl2wh2rd'    --> blue to white to red
                                   - 'mt_bl2gr2rd'    --> blue to green to red
                                   - 'mt_rd2gr2bl'    --> red to green to blue
                                   - 'mt_seg_bl2wh2rd' --> discrete blue to 
                                                           white to red
    
    Attributes:
    ------------
    
        -ellipse_cmap         ellipse color map, see above for options
        -ellipse_colorby      parameter to color ellipse by
        -ellipse_range        (min, max, step) values to color ellipses
        -ellipse_size         scaling factor to make ellipses visible
                                                           
                                                           
    """

    def __init__(self, **kwargs):
        super(MTEllipse, self).__init__()
        self.ellipse_size = 2
        self.ellipse_colorby = 'phimin'
        self.ellipse_range = (0, 90, 10)
        self.ellipse_cmap = 'mt_bl2gr2rd'

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if(v in list(kwargs.keys())):
                setattr(self, v, kwargs.pop(v, None))


    def _read_ellipse_dict(self, ellipse_dict):
        """
        read in dictionary and set default values if no entry given
        """
        # check all values are populated:
        default_dict = {'size': 2,
                        'ellipse_range':[0,0],
                        'ellipse_colorby':'skew',
                        'ellipse_cmap':'mt_bl2gr2rd'}
        for key in list(default_dict.keys()):
            if key not in list(ellipse_dict.keys()):
                ellipse_dict[key] = default_dict[key]
        
        # --> set the ellipse properties
        for key in list(ellipse_dict.keys()):
            setattr(self, key, ellipse_dict[key])

        try:
            self.ellipse_range[2]
        except IndexError:
            self.ellipse_range = (self.ellipse_range[0],
                                  self.ellipse_range[1],
                                  1)

        # set color ranges
        if(self.ellipse_range[0] == self.ellipse_range[1]): # override default-dict values
            if self.ellipse_colorby == 'skew' or \
                            self.ellipse_colorby == 'skew_seg' or \
                            self.ellipse_colorby == 'normalized_skew' or \
                            self.ellipse_colorby == 'normalized_skew_seg':

                self.ellipse_range = (-9, 9, 3)
            elif self.ellipse_colorby == 'ellipticity':
                self.ellipse_range = (0, 1, .1)
            else:
                self.ellipse_range = (0, 90, 5)
        # end if
        # only one colormap valid for skew_seg at this point in time
        if self.ellipse_colorby == 'skew_seg' or \
                        self.ellipse_colorby == 'normalized_skew_seg':
            print("Updating colormap to mt_seg_bl2wh2rd as this is the only available segmented colormap at this time")
            self.ellipse_cmap = 'mt_seg_bl2wh2rd'


        # set colormap to yellow to red
        '''
        if self.ellipse_colorby == 'skew' or \
                        self.ellipse_colorby == 'normalized_skew':
            self.ellipse_cmap = 'mt_bl2wh2rd'

        elif self.ellipse_colorby == 'skew_seg' or \
                        self.ellipse_colorby == 'normalized_skew_seg':
            self.ellipse_cmap = 'mt_seg_bl2wh2rd'

        else:
            self.ellipse_cmap = 'mt_bl2gr2rd'
        '''


# ==============================================================================
# Plot settings
# ==============================================================================
class PlotSettings(MTArrows, MTEllipse):
    """
    Hold all the plot settings that one might need
    """

    def __init__(self, **kwargs):
        super(PlotSettings, self).__init__(**kwargs)

        # figure properties:
        self.fig_num = 1
        self.fig_dpi = 150
        self.fig_size = None

        self.font_size = 7
        self.marker_size = 4
        self.marker_lw = .75
        self.lw = 1
        self.plot_title = None

        # line styles:
        self.xy_ls = ':'
        self.yx_ls = ':'
        self.det_ls = ':'
        self.skew_ls = ':'
        self.strike_ls = ':'

        # marker styles:
        self.xy_marker = 's'
        self.yx_marker = 'o'
        self.det_marker = 'v'
        self.skew_marker = 'd'
        self.strike_inv_marker = 'v'
        self.strike_pt_marker = '^'
        self.strike_tip_marker = '>'

        # marker color styles:
        self.xy_color = (0, 0, .75)
        self.yx_color = (.75, 0, 0)
        self.det_color = (0, .75, 0)
        self.skew_color = (.85, .35, 0)
        self.strike_inv_color = (.2, .2, .7)
        self.strike_pt_color = (.7, .2, .2)
        self.strike_tip_color = (.2, .7, .2)

        # marker face color styles:
        self.xy_mfc = (0, 0, .75)
        self.yx_mfc = (.75, 0, 0)
        self.det_mfc = (0, .75, 0)
        self.skew_mfc = (.85, .35, 0)
        self.strike_inv_mfc = (.2, .2, .7)
        self.strike_pt_mfc = (.7, .2, .2)
        self.strike_tip_mfc = (.2, .7, .2)

        # plot limits
        self.x_limits = None
        self.res_limits = None
        self.phase_limits = None
        self.tipper_limits = None
        self.strike_limits = None
        self.skew_limits = None
        self.pt_limits = None

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if(v in list(kwargs.keys())):
                setattr(self, v, kwargs.pop(v, None))


# ==============================================================================
# object for computing resistivity and phase
# ==============================================================================
class ResPhase(object):
    """
    Helper class to create a data type for just the resistivity and phase
    
    Arguments:
    ----------
        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None
                      
        **res_array** : np.ndarray((nf,2,2))
                        Array of resistivity values on a linear scale with 
                        the shape being [[rho_xx, rho_xy],[rho_yx, rho_yy]]
                        *default* is None
                        
        **res_err_array** : np.ndarray((nf,2,2))
                        Array of resistivity error values on a linear scale,  
                        the shape being [[rho_xx, rho_xy],[rho_yx, rho_yy]]
                        *default* is None
                        
        **phase_array** : np.ndarray((nf,2,2))
                        Array of resistivity values on a linear scale with 
                        the shape being [[phi_xx, phi_xy],[phi_yx, phi_yy]]
                        *default* is None
                        
        **phase_err_array** : np.ndarray((nf,2,2))
                        Array of resistivity error values on a linear scale,  
                        the shape being [[phi_xx, phi_xy],[phi_yx, phi_yy]]
                        *default* is None
                        
        **rot_z** : float
                    angle to rotate data by assuming North is 0 and measured
                    clockwise positive
                    
        **period** : np.ndarray(nf)
                     period array same length as res_array.  Needs to be input 
                     if inputing res and phase.
                        
    
    Attributes:
    ------------
    
        -period       array of periods corresponding to elements in res/phase
        
        -phase         phase array of shape (nf, 2, 2)
        -phase_err     phase error array of shape (nf, 2, 2)
        -phasedet      determinant of phase with shape (nf)
        -phasedet_err  determinant error of phase with shape (nf)
        -phasexx       xx component of phase with shape (nf)
        -phasexx_err   xx component of phase error with shape (nf)
        -phasexy       xy component of phase with shape (nf)
        -phasexy_err   xy component of phase error with shape (nf)
        -phaseyx       yx component of phase with shape (nf)
        -phaseyx_err   yx component of phase error with shape (nf)
        -phaseyy       yy component of phase with shape (nf)
        -phaseyy_err   yy component of phase error with shape (nf)
        
        -phase_quadrant [ 1 | 3 ] 
                        * 1 for both phases to be in 0 to 90, 
                        * 3 for xy to be in 0-90 and yx to be in -180 to 270            
        
        -res           res array of shape (nf, 2, 2)
        -res_err       res error array of shape (nf, 2, 2)
        -resdet        determinant of res with shape (nf)
        -resdet_err    determinant error of res with shape (nf)
        -resxx         xx component of res with shape (nf)
        -resxx_err     xx component of res error with shape (nf)
        -resxy         xy component of res with shape (nf)
        -resxy_err     xy component of res error with shape (nf)
        -resyx         yx component of res with shape (nf)
        -resyx_err     yx component of res error with shape (nf)
        -resyy         yy component of res with shape (nf)
        -resyy_err     yy component of res error with shape (nf)
        
    Methods:
    ---------
    
        -compute_res_phase   computes the above attributes on call
        -rotate              rotates the data
     
    """

    def __init__(self, z_object=None, res_array=None, res_err_array=None,
                 phase_array=None, phase_err_array=None, rot_z=0, period=None,
                 phase_quadrant=1):

        self._Z = z_object
        self.res = res_array
        self.res_err = res_err_array
        self.phase = phase_array
        self.phase_err = phase_err_array
        self.period = period
        self.phase_quadrant = phase_quadrant

        # check to make sure they are the same size
        if self.res is not None or self.phase is not None:
            if self.res.shape != self.phase.shape:
                raise mtex.MTpyError_Z('res_array and phase_array ' + \
                                       'are not the same shape')

            if self._Z is None:
                self._Z = mtz.Z()
                self._Z.set_res_phase(res_array, phase_array,
                                      res_err_array=res_err_array,
                                      phase_err_array=phase_err_array)
                if self.period is None:
                    raise mtex.MTpyError_Z('Need to input period to ' + \
                                           'compute z.')
                self._Z.freq = 1. / period

        try:
            if len(z_object.freq) == 0:
                raise mtex.MTpyError_Z('Need to set Z.freq to an' + \
                                       ' array that cooresponds to Z.z')
        except TypeError:
            raise mtex.MTpyError_Z('Need to set Z.freq to an' + \
                                   ' array that cooresponds to Z.z')
        if period is not None:
            self._Z.freq = 1. / period
        else:
            self.period = 1. / self._Z.freq

        if rot_z != 0:
            self.rotate(rot_z)

        # compute the resistivity and phase components
        self.compute_res_phase()

    def compute_res_phase(self):
        """
        computes the resistivity and phase and sets each component as an 
        attribute
        """

        if self._Z is not None:
            self._Z.compute_resistivity_phase()
            self.res = self._Z.resistivity
            self.phase = self._Z.phase
            self.res_err = self._Z.resistivity_err
            self.phase_err = self._Z.phase_err

        # check to see if a res_err_array was input if not set to zeros
        if self.res_err is None:
            self.res_err = np.zeros_like(self.res)

        else:
            # check to see if res and res_err are the same shape
            if self.res.shape != self.res_err.shape:
                mtex.MTpyError_inputarguments
                raise mtex.MTpyError_inputarguments('res_array and res_err_array ' + \
                                                    'are not the same shape')

        # check to see if a phase_err_array was input, if not set to zeros
        if self.phase_err is None:
            self.phase_err = np.zeros_like(self.phase)

        else:
            # check to see if res and res_err are the same shape
            if self.phase.shape != self.phase_err.shape:
                raise mtex.MTpyError_inputarguments('phase_array and ' + \
                                                    'phase_err_array are not the same shape')

        # --> set the attributes of the class to the components of each
        self.resxx = self.res[:, 0, 0]
        self.resxy = self.res[:, 0, 1]
        self.resyx = self.res[:, 1, 0]
        self.resyy = self.res[:, 1, 1]

        self.resxx_err = self.res_err[:, 0, 0]
        self.resxy_err = self.res_err[:, 0, 1]
        self.resyx_err = self.res_err[:, 1, 0]
        self.resyy_err = self.res_err[:, 1, 1]

        self.phasexx = self.phase[:, 0, 0]
        self.phasexy = self.phase[:, 0, 1]
        self.phaseyx = self.phase[:, 1, 0]
        self.phaseyy = self.phase[:, 1, 1]

        self.phasexx_err = self.phase_err[:, 0, 0]
        self.phasexy_err = self.phase_err[:, 0, 1]
        self.phaseyx_err = self.phase_err[:, 1, 0]
        self.phaseyy_err = self.phase_err[:, 1, 1]

        if self.phase_quadrant == 1:
            if self.phaseyx.mean() > 180:
                self.phaseyx -= 180
            else:
                self.phaseyx += 180

        # calculate determinant values
        zdet = np.array([np.linalg.det(zz) ** .5 for zz in self._Z.z])
        if self._Z.z_err is not None:
            zdetvar = np.array([np.linalg.det(zzv) ** .5 for zzv in self._Z.z_err])
        else:
            zdetvar = np.ones_like(zdet)

        # apparent resistivity
        self.resdet = 0.2 * (1. / self._Z.freq) * abs(zdet) ** 2
        self.resdet_err = 0.2 * (1. / self._Z.freq) * \
                          np.abs(zdet + zdetvar) ** 2 - self.resdet

        # phase
        self.phasedet = np.arctan2(zdet.imag, zdet.real) * (180 / np.pi)
        self.phasedet_err = np.arcsin(zdetvar / abs(zdet)) * (180 / np.pi)

    def rotate(self, rot_z):
        """
        rotate the impedance tensor by the angle rot_z (deg) assuming 
        0 is North and angle is positve clockwise.
        """

        self.rot_z = rot_z

        # need to fill the rotation angles in Z for this to work
        try:
            len(self.rot_z)
        except TypeError:
            self._Z.rotation_angle = np.array([rot_z
                                               for rr in range(self.res.shape[0])])
        if self._Z is not None:
            self._Z.rotate(self.rot_z)

        else:
            raise mtex.MTpyError_Z('Cannot rotate just resistivity and ' + \
                                   'phase data, need to input z')

        self.compute_res_phase()


# ==============================================================================
#  Define a tipper object
# ==============================================================================
class Tipper(object):
    """
    Helper class to put tipper data into a usable format
    
    Arguments:
    ----------
        **tipper_object** : class mtpy.core.z.Tipper
    
        **tipper_array** : np.ndarray((nf, 1, 2))
                           array of the complex tipper as [tx, ty]
        **tipper_err_array** : np.ndarray((nf, 1, 2))
                               array of the tipper error as [tx, ty]
   
   Attributes:
    -----------
    
        -mag_real   magnitude of real induction arrow
        -mag_imag   magnitude of imaginary induction arrow
        -ang_real   angle of real induction arrow assuming 0 is North positive 
                    clockwise
        _ang_imag   angle of imaginary induction arrow assuming 0 is North 
                    positive clockwise
        
    Methods:
    --------
    
        -compute_components   computes above attributes on call
        -rotate     rotates the data assuming 0 is North positive 
                    clockwise  
    """

    def __init__(self, tipper_object=None, tipper_array=None,
                 tipper_err_array=None, rot_t=0, freq=None):

        if tipper_object is not None:
            self._Tipper = tipper_object
        else:
            self._Tipper = mtz.Tipper(tipper_array=tipper_array,
                                      tipper_err_array=tipper_err_array,
                                      freq=freq)

        self.freq = freq

        if rot_t != 0:
            self.rotate(rot_t)

        else:
            self.compute_components()

    def compute_components(self):

        if self._Tipper.tipper is None:
            self.mag_imag = np.zeros_like(self.freq)
            self.mag_real = np.zeros_like(self.freq)
            self.ang_real = np.zeros_like(self.freq)
            self.ang_imag = np.zeros_like(self.freq)
        else:
            self.mag_real = self._Tipper.mag_real
            self.ang_real = self._Tipper.angle_real
            self.mag_imag = self._Tipper.mag_imag
            self.ang_imag = self._Tipper.angle_imag

    def rotate(self, rot_t):

        self.rot_t = rot_t
        self._Tipper.rotate(self.rot_t)

        self.compute_components()


# ==============================================================================
#  make an MT object that has all the important information and methods               
# ==============================================================================

class MTplot(mt.MT):
    """
    This class will be able to read in the imporant information from either
    an .edi file or information input by hand and can be extended to other 
    file types in the future.  This is a helper class to get all the 
    information needed in one place.
    
    The normal usage is to input an edi file from which all the information
    is read in.  However, it would seem that not everyone uses the .edi format
    so an option is to write you're own class for that particular file type, or
    give it to us to deal with.  Or write enough code to get it into one of 
    the suppported forms as arrays of z or resistivity and phase, or put those
    into a z_object or tipper_object.  
    
    
    Arguments:
    ----------
    
        **fn** : string
                       full path to file to be read in.  At the moment only
                       .edi type files are supported. *default* is None
        
        **z** : np.array((nf, 2, 2), dtype='complex')
                impedance tensor with length of nf -> the number of freq
                *default* is None
                
        **z_err** : np.array((nf, 2, 2), dtype='real')
                    impedance tensor error estimates, same shape as z.
                    *default* is None
                    
        **res_array** : np.array((nf, 2, 2))
                        array of resistivity values in linear scale.
                        *default* is None
                        
        **res_err_array** : np.array((nf, 2, 2))
                            array of resistivity error estimates, same shape 
                            as res_array. *default* is None
                            
        **phase_array** : np.array((nf, 2, 2))
                          array of phase values in degrees, same shape as 
                          res_array. *default* is None
                          
        **phase_err_array** : np.array((nf, 2, 2))
                              array of phase error estimates, same shape as 
                              phase_array. *default* is None
                              
        **tipper_array** : np.array((nf, 1, 2), dtype='complex')
                           array of tipper values for tx, ty. *default* is None
                           
        **tipper_err_array** : np.array((nf, 1, 2))
                               array of tipper error estimates, same shape as
                               tipper_array. *default* is None
        
        **station** : string
                      name of the station to be plotted. *default* is None
                      
        **period** : np.array(nf)
                     array of periods that coorespond to the components of z
                     *default* is None
                      
        **lat** : float
                 latitude of the station to be plotted in decimal degrees.
                 *default* is None
                 
        **lon** : float
                 longitude of the station to be plotted in decimal degrees.
                 *default* is None
                 
        **elev** : float
                   elevation of the station to be plotted in meters.
                   *default* is None
                              
        **rot_z** : float or np.array(nf)
                    angle (deg) to rotate the data assuming North is 0 and 
                    angle is positive clockwise.  Can be input as an array to
                    rotate different periods by different angles. 
                    *default* is 0
                    
        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None
                      
        **tipper_object** : class mtpy.core.z.Tipper
                            object of mtpy.core.z. If this is input be sure the
                            attribute z.freq is filled.  
                            *default* is None 
    Attributes:
    -----------
        -z             impedance tensor as an np.array((nf, 2, 2))
            
        -z_err         estimates of impedance tensor error same shape as z
        
        -tipper        tipper elements in an np.array((nf, 1, 2))
        
        -tipper_err    estimates of tipper error, same shape as tipper
        
        -station       station name
        
        -period        period as an np.array(nf)
        
        -lat           latitude in decimal degrees
        
        -lon           longitude in decimal degrees
        
        -elev          elevation in meters
        
        -fn            filename read from
        
     
     These can be get/set by simple dot syntax.  
        
    :Example: ::
        
        >>> mt1 = mtplot.MTplot(fn=r'/home/mt/edifiles/mt01.edi')
        >>> mt1.station
        >>> 'mt01'
        >>> mt1.station = 'pb075'
        >>> mt1.station
        >>> 'pb075'
        
    **Note:** that the format of z and tipper are:
        
        ::
            
          z[ii, :, :] = np.array([[z_xx, Z_xy], [Z_yx, Z_yy]])
          z[ii, 0, 0] = z_xx
          z[ii, 0, 1] = z_xy
          z[ii, 1, 0] = z_yx
          z[ii, 1, 1] = z_yy
          tipper[ii, :, :] = np.array([[tx],[ty]])
          tipper[ii, 0, 0] = tx
          tipper[ii, 0, 1] = ty
        
    Methods:
    --------
        -get_ResPhase        returns a ResPhase object
        -get_PhaseTensor     returns a PhaseTensor object
        -get_Tipper          returns a Tipper object
        -get_Zinvariants     returns a Zinvariants object
        
    :Example: ::
        
        >>> import mtpy.imaging.mtplot as mtplot
        >>> mt1 = mtplot.MTplot(fn=r'/home/mt/edifiles/mt01.edi')
        >>> # if you don't have an .edi file but res and phase
        >>> mt1 = mtplot.MTplot(res_array=res, phase_array=phase, 
        >>> ...                 period=period, station='mt01')
        
        
        
    """

    def __init__(self, fn=None, z_object=None, tipper_object=None, **kwargs):

        super(MTplot, self).__init__()

        self._fn = fn
        if self.fn is not None:
            self.read_mt_file(self.fn)

        if z_object is not None:
            self.Z = z_object

        if tipper_object is not None:
            self.Tipper = tipper_object

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if(v in list(kwargs.keys())):
                setattr(self, v, kwargs.pop(v, None))


    @property
    def period(self):
        if np.all(self.Z.z == 0 + 0j):
            if np.all(self.Tipper.tipper == 0 + 0j):
                return None
            else:
                return 1. / self.Tipper.freq
        else:
            return 1. / self.Z.freq


# ==============================================================================
# define a class object that contains list of mt_objects
# ==============================================================================
class MTplot_list(object):
    """
    manipulates a list of MTplot objects
    
    """

    def __init__(self, fn_list=None, res_object_list=None, z_object_list=None,
                 tipper_object_list=None, mt_object_list=None):

        self._fn_list = fn_list
        self._res_object_list = res_object_list
        self._z_object_list = z_object_list
        self._tipper_object_list = tipper_object_list
        self.mt_list = mt_object_list

        if self.mt_list is None:
            self.mt_list = get_mtlist(fn_list=self._fn_list,
                                      res_object_list=self._res_object_list,
                                      z_object_list=self._z_object_list,
                                      tipper_object_list=self._tipper_object_list)

    def sort_by_offsets_profile(self, line_direction='ew'):
        """
        get list of offsets to sort the mt list
        
        """

        mm = sort_by_offsets(self.mt_list, line_direction=line_direction)

        self.mt_list_sort = mm[0]
        self.station_list = mm[1]
        self.offset_list = mm[2]

    def get_station_locations(self, map_scale='latlon', ref_point=(0, 0)):
        """
        creates a dictionary where the keys are the stations and the values
        are the index in the plot_mesh grid for the station location.
        
        *Note*: the handling of zone changes in UTM coordinates is rough and 
        needs to be changed.  If there are zone changes in your survey, stick
        to latlon.
        
        Arguments:
        ----------
            **map_scale**: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
            
            **ref_point**: (map_scale_x, map_scale_y)
                           reference point to center the map on, needs to be in
                           map_scale coordinates.
                           
        Returns:
        --------
            **map_dict**: dictionary
                          where keys are station names, and the values are the
                          index values (x, y) for the plot_meshgrid
                          
            **map_xarr**: np.ndarray(num_stations)
                          east-west values to of station to plot. 
                          
            **map_yarr**: np.ndarray(num_stations)
                          north-south values to of station to plot.  
        """

        mm = get_station_locations(self.mt_list, map_scale=map_scale,
                                   ref_point=ref_point)

        self.map_dict = mm[0]
        self.map_xarr = mm[1]
        self.map_yarr = mm[2]

    def get_rp_arrays(self, plot_period, sort_by='line', line_direction='ew',
                      map_scale='latlon', ref_point=(0, 0), ftol=.1):
        """
        get resistivity and phase values in the correct order according to 
        offsets and periods for either map view or pseudosection.
    
        Attributes:
        -----------
            **sort_by**: [ 'line' | 'map' ]
                         * 'line' --> sort the station distances into a line 
                                      according to line_direction
                         * 'map' --> sort the station distances into map 
                                     coordinates
            
            **line_direction**: [ 'ew' | 'ns' ]
            
            **map_scale**: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
            
            **ref_point**: (x, y)
                           reference point to center the plot on, this point 
                           needs to be in map coordinates
                           
            **ftol**: float
                      tolerance to match periods in mt_list with plot_period
                      
            **plot_period**: np.ndarray(nt)
                             array of periods in seconds to get data for.
                             
        Returns:
        --------
            Returns the individual components of resisitivity (in log scale) 
            and phase (deg) with a shape according to sort_by.  
            
            * If sort_by == 'line', the returned shape is (num_periods, 
                                                           num_stations)
            * If sort_by == 'map', the returned shape is (num_periods, 
                                                          num_stations, 
                                                          num_stations)
                                                           
        
            **resxx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       apparent resistivity (log 10 scale) for xx component
            **resxy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       apparent resistivity (log 10 scale) for xy component
            **resyx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       apparent resistivity (log 10 scale) for yx component
            **resyy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       apparent resistivity (log 10 scale) for yy component
                       
            **phasexx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       phase (deg) for xx component
            **phasexy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       phase (deg) for xy component
            **phaseyx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       phase (deg) for yx component
            **phaseyy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       phase (deg) for yy component
        
        """

        mm = get_rp_arrays(self.mt_list, sort_by=sort_by,
                           line_direction=line_direction,
                           map_scale=map_scale,
                           ref_point=ref_point,
                           ftol=ftol,
                           plot_period=plot_period)

        if sort_by == 'line':
            self.resxx_ps = mm[0]
            self.resxy_ps = mm[1]
            self.resyx_ps = mm[2]
            self.resyy_ps = mm[3]

            self.phasexx_ps = mm[4]
            self.phasexy_ps = mm[5]
            self.phaseyx_ps = mm[6]
            self.phaseyy_ps = mm[7]

            self.station_list = mm[8]
            self.offset_list = mm[9]

        if sort_by == 'map':
            self.resxx_map = mm[0]
            self.resxy_map = mm[1]
            self.resyx_map = mm[2]
            self.resyy_map = mm[3]

            self.phasexx_map = mm[4]
            self.phasexy_map = mm[5]
            self.phaseyx_map = mm[6]
            self.phaseyy_map = mm[7]

            self.map_x = mm[8]
            self.map_y = mm[9]

    def get_pt_arrays(self, plot_period, sort_by='line', line_direction='ew',
                      map_scale='latlon', ref_point=(0, 0), ftol=.1):
        """
        get resistivity and phase values in the correct order according to 
        offsets and periods for either map view or pseudosection.
    
        Attributes:
        -----------
            **sort_by**: [ 'line' | 'map' ]
                         * 'line' --> sort the station distances into a line 
                                      according to line_direction
                         * 'map' --> sort the station distances into map 
                                     coordinates
            
            **line_direction**: [ 'ew' | 'ns' ]
            
            **map_scale**: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
            
            **ref_point**: (x, y)
                           reference point to center the plot on, this point needs
                           to be in map coordinates
                           
            **ftol**: float
                      tolerance to match periods in mt_list with plot_period
                      
            **plot_period**: np.ndarray(nt)
                             array of periods in seconds to get data for.
                             
        Returns:
        --------
            Returns the individual parameters of the phase tensor (deg)
            
            * If sort_by == 'line', the returned shape is (num_periods, 
                                                           num_stations)
            * If sort_by == 'map', the returned shape is (num_periods, 
                                                          num_stations, 
                                                          num_stations)
                                                           
        
            **phimin**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                       minimum phase or 2nd principal component of phase tensor
                       
            **phimax**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                        maximum phase or 1st principal component of phase tensor
                       
            **skew**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                      skew angle of phase tensor
                       
            **azimuth**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                         regional strike direction estimated from phase tensor,
                         with a 90 degree ambiguity
                       
            **ellipticity**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                             ratio of phimin to phimax suggesting dimensionality.
        
        """

        mm = get_pt_arrays(self.mt_list,
                           sort_by=sort_by,
                           line_direction=line_direction,
                           map_scale=map_scale,
                           ref_point=ref_point,
                           ftol=ftol,
                           plot_period=plot_period)

        if sort_by == 'line':
            self.phimin_ps = mm[0]
            self.phimax_ps = mm[1]
            self.skew_ps = mm[2]
            self.azimuth_ps = mm[3]
            self.ellipticity_ps = mm[4]

            self.station_list = mm[5]
            self.offset_list = mm[6]

        elif sort_by == 'map':
            self.phimin_map = mm[0]
            self.phimax_map = mm[1]
            self.skew_map = mm[2]
            self.azimuth_map = mm[3]
            self.ellipticity_map = mm[4]

            self.map_xarr = mm[5]
            self.map_yarr = mm[6]


# ==============================================================================
# get list of mt objects     
# ==============================================================================
def get_mtlist(fn_list=None, res_object_list=None, z_object_list=None,
               tipper_object_list=None, mt_object_list=None):
    """
    gets a list of mt objects from the inputs  

    Arguments:     
    -----------
        **fn_list** : list of strings
                          full paths to .edi files to plot
                          
        **res_object_list** : list of mtplot.ResPhase objects
                             *default* is none
                          
        **z_object_list** : list of class mtpy.core.z.Z
                           object of mtpy.core.z.  If this is input be sure the
                           attribute z.freq is filled.  *default* is None
                      
        **mt_object_list** : list of class mtpy.imaging.mtplot.MTplot
                            object of mtpy.imaging.mtplot.MTplot
                            *default* is None
                            
    Returns:
    ---------
    
        **mt_list** : list of MTplot instances
    """

    # first need to find something to loop over

    if fn_list is not None:
        ns = len(fn_list)
        mt_list = [MTplot(fn=fn) for fn in fn_list]
        print('Reading {0} stations'.format(ns))
        return mt_list

    elif mt_object_list is not None:
        return mt_object_list

    elif z_object_list is not None:
        ns = len(z_object_list)
        mt_list = [MTplot(z_object=z_obj) for z_obj in z_object_list]
        try:
            nt = len(tipper_object_list)
            if nt != ns:
                raise mtex.MTpyError_inputarguments('length ' + \
                                                    ' of z_list is not equal to tip_list' + \
                                                    '; nz={0}, nt={1}'.format(ns, nt))
            for mt, tip_obj in zip(mt_list, tipper_object_list):
                mt._Tipper = tip_obj
        except TypeError:
            pass
        print('Reading {0} stations'.format(ns))
        return mt_list


# elif tipper_object_list is not None:
#        
#    elif type(fn_list[0]) is MTplot:
#        return mt_list
#        
#    else:
#        try:
#            ns = len(fn_list)
#            mt_list = [MTplot(fn=fn) for fn in fn_list]
#            print 'Reading {0} stations'.format(ns)
#            return mt_list
#        except TypeError:
#            try:
#                ns = len(res_object_list)
#                mt_list = [MTplot(res_phase_object=res_obj) 
#                            for res_obj in res_object_list]
#                try:
#                    nt = len(tipper_object_list)
#                    if nt != ns:
#                        raise mtex.MTpyError_inputarguments('length '+\
#                              ' of z_list is not equal to tip_list'+\
#                              '; nz={0}, nt={1}'.format(ns, nt))
#                    for mt,tip_obj in zip(mt_list,tipper_object_list):
#                        mt._Tipper = tip_obj 
#                except TypeError:
#                    pass
#                print 'Reading {0} stations'.format(ns)
#                return mt_list
#            except TypeError:
#                try: 
#                    ns = len(z_object_list)
#                    mt_list = [MTplot(z_object=z_obj) for z_obj in z_object_list]
#                    try:
#                        nt = len(tipper_object_list)
#                        if nt != ns:
#                            raise mtex.MTpyError_inputarguments('length '+\
#                                  ' of z_list is not equal to tip_list'+\
#                                  '; nz={0}, nt={1}'.format(ns, nt))
#                        for mt,tip_obj in zip(mt_list,tipper_object_list):
#                            mt._Tipper = tip_obj 
#                    except TypeError:
#                        pass
#                    print 'Reading {0} stations'.format(ns)
#                    return mt_list
#                    
#                except TypeError:
#                    try:
#                        ns = len(mt_object_list)
#                        print 'Reading {0} stations'.format(ns)
#                        return mt_list
#                    except TypeError:
#                        raise IOError('Need to input an iteratable list')

# ==============================================================================
# sort an mt_list by offset values in a particular direction                  
# ==============================================================================
def sort_by_offsets(mt_list, line_direction='ew'):
    """
    get list of offsets for the given line_direction.
    
    Arguments:
    ----------
        **mt_list**: list
                    list of MTplot objects
        
        **line_direction**: [ 'ew' | 'ns' ]
        
    Returns:
    --------
        **sort_mt_list**: list
                         list of MTplot objects sorted by offset in the
                         line_direction
                         
        **station_list**: list
                         list of stations sorted by offset
        
        **offset_list**: np.ndarray(num_stations)
                        array of sorted offset values corresponding to station
                        in station_list
    """

    dtype = [('station', 'S10'), ('offset', float), ('spot', int)]
    slist = []
    # get offsets
    for ii, mt in enumerate(mt_list):
        # get offsets between stations
        if ii == 0:
            east0 = mt.lon
            north0 = mt.lat
            offset = 0.0
        else:
            east = mt.lon
            north = mt.lat
            # if line is predominantly e-w
            if line_direction == 'ew':
                if east0 < east:
                    offset = np.sqrt((east0 - east) ** 2 + (north0 - north) ** 2)
                elif east0 > east:
                    offset = -1 * np.sqrt((east0 - east) ** 2 + (north0 - north) ** 2)
                else:
                    offset = 0
            # if line is predominantly n-s
            elif line_direction == 'ns':
                if north0 < north:
                    offset = np.sqrt((east0 - east) ** 2 + (north0 - north) ** 2)
                elif north0 > north:
                    offset = -1 * np.sqrt((east0 - east) ** 2 + (north0 - north) ** 2)
                else:
                    offset = 0
        # append values to list for sorting
        slist.append((mt.station, offset, ii))

    # create a structured array according to the data type and values
    v_array = np.array(slist, dtype=dtype)

    # sort the structured array by offsets
    sorted_array = np.sort(v_array, order=['offset'])

    # create an offset list as an attribute
    offset_list = np.array([ss[1] for ss in sorted_array])

    # create a station list as an attribute
    station_list = np.array([ss[0] for ss in sorted_array])

    # create an index list of the sorted index values
    index_list = [ss[2] for ss in sorted_array]

    # create a new mt_list according to the offsets from the new index_list
    sort_mt_list = [mt_list[ii] for ii in index_list]

    return sort_mt_list, station_list, offset_list


# ==============================================================================
# get map values
# ==============================================================================
def get_station_locations(mt_list, map_scale='latlon', ref_point=(0, 0)):
    """
    creates a dictionary where the keys are the stations and the values
    are the index in the plot_mesh grid for the station location.
    
    *Note*: the handling of zone changes in UTM coordinates is rough and 
    needs to be changed.  If there are zone changes in your survey, stick
    to latlon.
    
    Arguments:
    ----------
        **map_scale**: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
        
        **ref_point**: (map_scale_x, map_scale_y)
                       reference point to center the map on, needs to be in
                       map_scale coordinates.
 
                       
    Returns:
    --------
        **map_dict**: dictionary
                      where keys are station names, and the values are the
                      index values (x, y) for the plot_meshgrid
                      
        **plot_meshgrid**: np.ndarray(num_stations, num_stations)
                           a meshgrid (x, y) for plotting map view in 
                           map_scale coordinates.  
    """

    # make some empty arrays
    lat_list = np.zeros(len(mt_list))
    lon_list = np.zeros(len(mt_list))
    elev_list = np.zeros(len(mt_list))
    x_arr = np.zeros(len(mt_list))
    y_arr = np.zeros(len(mt_list))

    if map_scale == 'eastnorth':
        dscale = 1.
    elif map_scale == 'eastnorthkm':
        dscale = 1000.

    map_station_dict = {}
    # need to sort by station
    for ii, mt in enumerate(mt_list):
        lat_list[ii] = mt.lat
        lon_list[ii] = mt.lon
        elev_list[ii] = mt.elev

        # if map scale is lat lon set parameters
        if map_scale == 'latlon':
            x = mt.lon - ref_point[0]
            y = mt.lat - ref_point[1]

        # if map scale is in meters easting and northing
        elif map_scale == 'eastnorth' or map_scale == 'eastnorthkm':
            east, north, zone = gis_tools.project_point_ll2utm(mt.lat, mt.lon)

            east /= dscale
            north /= dscale

            # set the first point read in as a refernce other points
            if ii == 0:
                zone1 = zone
                x = east - ref_point[0]
                y = north - ref_point[1]

            # read in all the other point
            else:
                # check to make sure the zone is the same this needs
                # to be more rigorously done
                if zone1 != zone:
                    print('Zone change at station ' + mt.station)
                    if zone1[0:2] == zone[0:2]:
                        pass
                    elif int(zone1[0:2]) < int(zone[0:2]):
                        east += 500000
                    else:
                        east -= -500000

                    x = east - ref_point[0]
                    y = north - ref_point[1]
                else:
                    x = east - ref_point[0]
                    y = north - ref_point[1]
        else:
            raise NameError('mapscale not recognized')

        # put the location of each ellipse into an array in x and y
        x_arr[ii] = x
        y_arr[ii] = y

        map_station_dict[mt.station] = (x, y, mt.elev)

    return map_station_dict, x_arr, y_arr


# ==============================================================================
# grid data onto a map view
# ==============================================================================
def grid_data(data_array, x, y, nx=None, ny=None):
    """
    Project data onto a regular grid for plotting.
    
    
    Arguments:
    -----------
        **data_array**: np.ndarray (len(x), len(y))
                        array of data values to be gridded
                        
        **x**: np.ndarray(len(x))
               array of values that coorespond  
    
        **nx**: int
                number of cells in the x-direction.  If none, 2 times the 
                number of x components
                
        **ny**: int
                number of cells in the x-direction.  If none, 2 times the 
                number of y components
                
    Returns:
    ---------
        **grid_array**: np.ndarray(nx, ny)
                        array of data set on a regular grid
        
        **xg**: np.ndarray(nx, ny)
                array of x-grid values
                
        **yg**: np.ndarray(nx, ny)
                array of y-grid values
                
        
    """

    if nx is None:
        nx = 2 * len(x)

    if ny is None:
        ny = 2 * len(y)

    # create evenly spaced intervals to grid over
    xi = np.linspace(x.min(), x.max(), num=nx, endpoint=True)
    yi = np.linspace(y.min(), y.max(), num=ny, endpoint=True)

    xg, yg = np.meshgrid(xi, yi)

    grid_array = mlab.griddata(x, y, data_array, xg, yg)

    return grid_array, xg, yg


# ==============================================================================
# get resistivity and phase arrays for plotting
# ==============================================================================
def get_rp_arrays(mt_list, plot_period, sort_by='line', line_direction='ew',
                  map_scale='latlon', ref_point=(0, 0), ftol=.1):
    """
    get resistivity and phase values in the correct order according to 
    offsets and periods for either map view or pseudosection.

    Attributes:
    -----------
        **sort_by**: [ 'line' | 'map' ]
                     * 'line' --> sort the station distances into a line 
                                  according to line_direction
                     * 'map' --> sort the station distances into map 
                                 coordinates
        
        **line_direction**: [ 'ew' | 'ns' ]
        
        **map_scale**: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
        
        **ref_point**: (x, y)
                       reference point to center the plot on, this point needs
                       to be in map coordinates
                       
        **ftol**: float
                  tolerance to match periods in mt_list with plot_period
                  
        **plot_period**: np.ndarray(nt)
                         array of periods in seconds to get data for.
                         
    Returns:
    --------
        Returns the individual components of resisitivity (in log scale) and 
        phase (deg) with a shape according to sort_by.  
        
        * If sort_by == 'line', the returned shape is (num_periods, 
                                                       num_stations)
        * If sort_by == 'map', the returned shape is (num_periods, 
                                                      num_stations, 
                                                      num_stations)
                                                       
    
        **resxx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   apparent resistivity (log 10 scale) for xx component
        **resxy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   apparent resistivity (log 10 scale) for xy component
        **resyx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   apparent resistivity (log 10 scale) for yx component
        **resyy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   apparent resistivity (log 10 scale) for yy component
                   
        **phasexx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   phase (deg) for xx component
        **phasexy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   phase (deg) for xy component
        **phaseyx**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   phase (deg) for yx component
        **phaseyy**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   phase (deg) for yy component
    
    """
    if plot_period is None:
        raise mtex.MTpyError_inputarguments('Need to input an array of ' + \
                                            'periods')

    ns = len(mt_list)
    nt = len(plot_period)

    # make a dictionary of the periods to plot for a reference
    period_dict = dict([(key, vv)
                        for vv, key in enumerate(plot_period)])

    # get arrays in pseudosection format
    if sort_by == 'line':
        # sort the data by offset
        mt_list_sort, station_list, offset_list = sort_by_offsets(mt_list,
                                                                  line_direction=line_direction)

        # create empty arrays to put data into
        resxx = np.zeros((nt, ns))
        resxy = np.zeros((nt, ns))
        resyx = np.zeros((nt, ns))
        resyy = np.zeros((nt, ns))

        phasexx = np.zeros((nt, ns))
        phasexy = np.zeros((nt, ns))
        phaseyx = np.zeros((nt, ns))
        phaseyy = np.zeros((nt, ns))

        for ii, mt in enumerate(mt_list_sort):
            # get resisitivity and phase in a dictionary and append to a list
            rp = mt.Z

            for rr, rper in enumerate(plot_period):
                jj = None
                for kk, iper in enumerate(mt.period):
                    if iper == rper:
                        jj = period_dict[rper]
                        resxx[jj, ii] = np.log10(rp.res_xx[kk])
                        resxy[jj, ii] = np.log10(rp.res_xy[kk])
                        resyx[jj, ii] = np.log10(rp.res_yx[kk])
                        resyy[jj, ii] = np.log10(rp.res_yy[kk])

                        phasexx[jj, ii] = rp.phase_xx[kk]
                        phasexy[jj, ii] = rp.phase_xy[kk]
                        phaseyx[jj, ii] = rp.phase_yx[kk]
                        phaseyy[jj, ii] = rp.phase_yy[kk]

                        break

                    elif rper * (1 - ftol) <= iper and \
                                    iper <= rper * (1 + ftol):
                        jj = period_dict[rper]
                        resxx[jj, ii] = np.log10(rp.res_xx[kk])
                        resxy[jj, ii] = np.log10(rp.res_xy[kk])
                        resyx[jj, ii] = np.log10(rp.res_yx[kk])
                        resyy[jj, ii] = np.log10(rp.res_yy[kk])

                        phasexx[jj, ii] = rp.phase_xx[kk]
                        phasexy[jj, ii] = rp.phase_xy[kk]
                        phaseyx[jj, ii] = rp.phase_yx[kk]
                        phaseyy[jj, ii] = rp.phase_yy[kk]

                        break
                    else:
                        pass

                if jj is None:
                    print('did not find period {0:.6g} (s) for {1}'.format(
                        rper, mt.station))
        return resxx, resxy, resyx, resyy, phasexx, phasexy, phaseyx, phaseyy, \
               station_list, offset_list

    elif sort_by == 'map':
        map_dict, x, y = get_station_locations(mt_list,
                                               map_scale=map_scale,
                                               ref_point=ref_point)

        resxx = np.zeros((nt, ns))
        resxy = np.zeros((nt, ns))
        resyx = np.zeros((nt, ns, ns))
        resyy = np.zeros((nt, ns, ns))

        phasexx = np.zeros((nt, ns))
        phasexy = np.zeros((nt, ns))
        phaseyx = np.zeros((nt, ns))
        phaseyy = np.zeros((nt, ns))

        for ii, mt in enumerate(mt_list):
            # get resisitivity and phase in a dictionary and append to a list
            rp = mt.get_ResPhase()
            #            rp = mt.Z

            for rr, rper in enumerate(plot_period):
                jj = None
                for kk, iper in enumerate(mt.period):
                    if iper == rper:
                        jj = period_dict[rper]

                        resxx[jj, ii] = np.log10(rp.res_xx[kk])
                        resxy[jj, ii] = np.log10(rp.res_xy[kk])
                        resyx[jj, ii] = np.log10(rp.res_yx[kk])
                        resyy[jj, ii] = np.log10(rp.res_yy[kk])

                        phasexx[jj, ii] = rp.phase_xx[kk]
                        phasexy[jj, ii] = rp.phase_xy[kk]
                        phaseyx[jj, ii] = rp.phase_yx[kk]
                        phaseyy[jj, ii] = rp.phase_yy[kk]

                        break

                    elif rper * (1 - ftol) <= iper and \
                                    iper <= rper * (1 + ftol):
                        jj = period_dict[rper]

                        resxx[jj, ii] = np.log10(rp.res_xx[kk])
                        resxy[jj, ii] = np.log10(rp.res_xy[kk])
                        resyx[jj, ii] = np.log10(rp.res_yx[kk])
                        resyy[jj, ii] = np.log10(rp.res_yy[kk])

                        phasexx[jj, ii] = rp.phase_xx[kk]
                        phasexy[jj, ii] = rp.phase_xy[kk]
                        phaseyx[jj, ii] = rp.phase_yx[kk]
                        phaseyy[jj, ii] = rp.phase_yy[kk]

                        break
                    else:
                        pass

                if jj is None:
                    print('did not find period {0:.6g} (s) for {1}'.format(
                        rper, mt.station))
        return resxx, resxy, resyx, resyy, + \
            phasexx, phasexy, phaseyx, phaseyy, x, y, map_dict


# ==============================================================================
# get phase tensor arrays for plotting
# ==============================================================================
def get_pt_arrays(mt_list, plot_period, sort_by='line', line_direction='ew',
                  map_scale='latlon', ref_point=(0, 0), ftol=.1):
    """
    get resistivity and phase values in the correct order according to 
    offsets and periods for either map view or pseudosection.

    Attributes:
    -----------
        **sort_by**: [ 'line' | 'map' ]
                     * 'line' --> sort the station distances into a line 
                                  according to line_direction
                     * 'map' --> sort the station distances into map 
                                 coordinates
        
        **line_direction**: [ 'ew' | 'ns' ]
        
        **map_scale**: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
        
        **ref_point**: (x, y)
                       reference point to center the plot on, this point needs
                       to be in map coordinates
                       
        **ftol**: float
                  tolerance to match periods in mt_list with plot_period
                  
        **plot_period**: np.ndarray(nt)
                         array of periods in seconds to get data for.
                         
    Returns:
    --------
        Returns the individual parameters of the phase tensor (deg)
        
        * If sort_by == 'line', the returned shape is (num_periods, 
                                                       num_stations)
        * If sort_by == 'map', the returned shape is (num_periods, 
                                                      num_stations, 
                                                      num_stations)
                                                       
    
        **phimin**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                   minimum phase or 2nd principal component of phase tensor
                   
        **phimax**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                    maximum phase or 1st principal component of phase tensor
                   
        **skew**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                  skew angle of phase tensor
                   
        **azimuth**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                     regional strike direction estimated from phase tensor,
                     with a 90 degree ambiguity
                   
        **ellipticity**: np.ndarray(nt, ns) or np.ndarray(nt, ns, ns)
                         ratio of phimin to phimax suggesting dimensionality.
    
    """
    if plot_period is None:
        raise mtex.MTpyError_inputarguments('Need to input an array of ' + \
                                            'periods')

    ns = len(mt_list)
    nt = len(plot_period)

    # make a dictionary of the periods to plot for a reference
    period_dict = dict([(key, vv)
                        for vv, key in enumerate(plot_period)])

    # get arrays in pseudosection format
    if sort_by == 'line':

        mt_list_sort, slist, olist = sort_by_offsets(mt_list,
                                                     line_direction=line_direction)

        # create empty arrays to put data into need to reset to zero in case
        # something has changed


        phimin = np.zeros((nt, ns))
        phimax = np.zeros((nt, ns))
        skew = np.zeros((nt, ns))
        azimuth = np.zeros((nt, ns))
        ellipticity = np.zeros((nt, ns))

        for ii, mt in enumerate(mt_list_sort):
            # get resisitivity and phase in a dictionary and append to a list
            pt = mt.pt

            for rr, rper in enumerate(plot_period):
                jj = None
                for kk, iper in enumerate(mt.period):
                    if iper == rper:
                        jj = period_dict[rper]
                        phimin[jj, ii] = pt.phimin[kk]
                        phimax[jj, ii] = pt.phimax[kk]
                        skew[jj, ii] = pt.beta[kk]
                        azimuth[jj, ii] = pt.azimuth[kk]
                        ellipticity[jj, ii] = pt.ellipticity[kk]

                        break

                    elif rper * (1 - ftol) <= iper and \
                                    iper <= rper * (1 + ftol):
                        jj = period_dict[rper]
                        phimin[jj, ii] = pt.phimin[kk]
                        phimax[jj, ii] = pt.phimax[kk]
                        skew[jj, ii] = pt.beta[kk]
                        azimuth[jj, ii] = pt.azimuth[kk]
                        ellipticity[jj, ii] = pt.ellipticity[kk]

                        break
                    else:
                        pass

                if jj is None:
                    print('did not find period {0:.6g} (s) for {1}'.format(
                        rper, mt.station))
        return phimin, phimax, skew, azimuth, ellipticity, slist, olist

    elif sort_by == 'map':
        map_dict, x, y = get_station_locations(mt_list,
                                               map_scale=map_scale,
                                               ref_point=ref_point)

        phimin = np.zeros((nt, ns))
        phimax = np.zeros((nt, ns))
        skew = np.zeros((nt, ns))
        azimuth = np.zeros((nt, ns))
        ellipticity = np.zeros((nt, ns))

        for ii, mt in enumerate(mt_list):
            # get resisitivity and phase in a dictionary and append to a list
            pt = mt.pt

            for rr, rper in enumerate(plot_period):
                jj = None
                for kk, iper in enumerate(mt.period):
                    if iper == rper:
                        jj = period_dict[rper]

                        phimin[jj, ii] = pt.phimin[0][kk]
                        phimax[jj, ii] = pt.phimax[0][kk]
                        skew[jj, ii] = pt.beta[0][kk]
                        azimuth[jj, ii] = pt.azimuth[0][kk]
                        ellipticity[jj, ii] = pt.ellipticity[0][kk]

                        break

                    elif rper * (1 - ftol) <= iper and \
                                    iper <= rper * (1 + ftol):
                        jj = period_dict[rper]
                        phimin[jj, ii] = pt.phimin[0][kk]
                        phimax[jj, ii] = pt.phimax[0][kk]
                        skew[jj, ii] = pt.beta[0][kk]
                        azimuth[jj, ii] = pt.azimuth[0][kk]
                        ellipticity[jj, ii] = pt.ellipticity[0][kk]

                        break
                    else:
                        pass

                if jj is None:
                    print('did not find period {0:.6g} (s) for {1}'.format(
                        rper, mt.station))
        return phimin, phimax, skew, azimuth, ellipticity, x, y, map_dict


# ==============================================================================
# function for writing values to file
# ==============================================================================
def make_value_str(value, value_list=None, spacing='{0:^8}',
                   value_format='{0: .2f}', append=False, add=False):
    """
    helper function for writing values to a file, takes in a value and either
    appends or adds value to value_list according to the spacing and format of 
    the string.
    
    Arguments:
    ----------
        **value** : float
        
        **value_list** : list of values converted to strings
        
        **spacing** : spacing of the string that the value will be converted
                      to.
                      
        **value_format** : format of the string that the value is being 
                            coverted to.
        
        **append** : [ True | False]
                     if True then appends the value to value list
        
        **add** : [ True | False ]
                  if True adds value string to the other value strings in
                  value_list
    
    Returns:
    --------
        **value_list** : the input value_list with the new value either 
                        added or appended.
        or
        
        **value_str** : value string if add and append are false
    """

    value_str = spacing.format(value_format.format(value))

    if append is True:
        value_list.append(value_str)
        return value_list
    if add is True:
        value_list += value_str
        return value_list

    if append == False and add == False:
        return value_str

    return value_list


# ==============================================================================
# function for error bar plots 
# ==============================================================================
def plot_errorbar(ax, x_array, y_array, y_error=None, x_error=None,
                  color='k', marker='x', ms=2, ls=':', lw=1, e_capsize=2,
                  e_capthick=.5, picker=None):
    """
    convinience function to make an error bar instance
    
    Arguments:
    ------------
        **ax** : matplotlib.axes instance 
                 axes to put error bar plot on
    
        **x_array** : np.ndarray(nx)
                      array of x values to plot
                      
        **y_array** : np.ndarray(nx)
                      array of y values to plot
                      
        **y_error** : np.ndarray(nx)
                      array of errors in y-direction to plot
        
        **x_error** : np.ndarray(ns)
                      array of error in x-direction to plot
                      
        **color** : string or (r, g, b)
                    color of marker, line and error bar
                    
        **marker** : string
                     marker type to plot data as
                     
        **ms** : float
                 size of marker
                 
        **ls** : string
                 line style between markers
                 
        **lw** : float
                 width of line between markers
        
        **e_capsize** : float
                        size of error bar cap
        
        **e_capthick** : float
                         thickness of error bar cap
        
        **picker** : float
                     radius in points to be able to pick a point. 
        
        
    Returns:
    ---------
        **errorbar_object** : matplotlib.Axes.errorbar 
                              error bar object containing line data, 
                              errorbars, etc.
    """
    # this is to make sure error bars plot in full and not just a dashed line
    if x_error is not None:
        #        x_err_high = np.array(x_error)
        #        x_err_low = np.array(x_err_high)
        #        x_err_low[x_err_high>=x_array] = x_array[x_err_high>=x_array]*.9999
        #        x_err = [x_err_low, x_err_high]
        x_err = x_error
    else:
        x_err = None

    if y_error is not None:
        #        y_err_high = np.array(y_error)
        #        y_err_low = np.array(y_err_high)
        #        y_err_low[y_err_high>=y_array] = y_array[y_err_high>=y_array]*.9999
        #        y_err = [y_err_low, y_err_high]
        y_err = y_error
    else:
        y_err = None

    errorbar_object = ax.errorbar(x_array,
                                  y_array,
                                  marker=marker,
                                  ms=ms,
                                  mfc='None',
                                  mew=lw,
                                  mec=color,
                                  ls=ls,
                                  xerr=x_err,
                                  yerr=y_err,
                                  ecolor=color,
                                  color=color,
                                  picker=picker,
                                  lw=lw,
                                  elinewidth=lw,
                                  capsize=e_capsize,
                                  #                                  capthick=e_capthick
                                  )
    return errorbar_object
