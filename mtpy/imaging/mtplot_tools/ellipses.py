# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:19:16 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

# =============================================================================


class MTEllipse:
    """
    helper class for getting ellipse properties from an input dictionary

    Arguments:
    -------------

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

    """

    def __init__(self, **kwargs):
        self.ellipse_size = 2
        self.ellipse_colorby = "phimin"
        self.ellipse_range = (0, 90, 10)
        self.ellipse_cmap = "mt_bl2gr2rd"
        self.ellipse_spacing = 1

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if v in list(kwargs.keys()):
                setattr(self, v, kwargs.pop(v, None))
        self.get_range()
        self.get_color_map()

    def get_color_map(self):
        """
        get a color map
        """
        if self.ellipse_colorby in ["skew_seg", "normalized_skew_seg"]:
            self.ellipse_cmap = "mt_seg_bl2wh2rd"

    def get_range(self):
        """
        get an appropriate range for the colorby
        """
        # set color ranges
        if self.ellipse_range[0] == self.ellipse_range[1]:
            if self.ellipse_colorby in [
                "skew",
                "skew_seg",
                "normalized_skew",
                "normalized_skew_seg",
            ]:

                self.ellipse_range = (-9, 9, 3)
            elif self.ellipse_colorby == "ellipticity":
                self.ellipse_range = (0, 1, 0.1)
            else:
                self.ellipse_range = (0, 90, 5)

    @property
    def ellipse_cmap_n_segments(self):
        return float(
            (self.ellipse_range[1] - self.ellipse_range[1])
            / (2 * self.ellipse_range[2])
        )

    @property
    def ellipse_cmap_bounds(self):
        try:
            return np.arange(
                self.ellipse_range[0],
                self.ellipse_range[1] + self.ellipse_range[2],
                self.ellipse_range[2],
            )
        except IndexError:
            return None

    def get_pt_color_array(self, pt_object):
        """
        Get the appropriat color by array
        """

        # get the properties to color the ellipses by
        if (
            self.ellipse_colorby == "phiminang"
            or self.ellipse_colorby == "phimin"
        ):
            color_array = pt_object.phimin
        elif (
            self.ellipse_colorby == "phimaxang"
            or self.ellipse_colorby == "phimax"
        ):
            color_array = pt_object.phimax
        elif self.ellipse_colorby == "phidet":
            color_array = np.sqrt(abs(pt_object.det)) * (180 / np.pi)
        elif (
            self.ellipse_colorby == "skew"
            or self.ellipse_colorby == "skew_seg"
        ):
            color_array = pt_object.beta
        elif self.ellipse_colorby == "ellipticity":
            color_array = pt_object.ellipticity
        elif self.ellipse_colorby in ["strike", "azimuth"]:
            color_array = pt_object.azimuth % 180
            color_array[np.where(color_array > 90)] -= 180
        else:
            raise NameError(self.ellipse_colorby + " is not supported")
        return color_array

    @property
    def ellipse_properties(self):
        return {
            "size": self.ellipse_size,
            "range": self.ellipse_range,
            "cmap": self.ellipse_cmap,
            "colorby": self.ellipse_colorby,
            "spacing": self.ellipse_spacing,
        }
