# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:16:31 2022

@author: jpeacock
"""

# ==============================================================================
# Arrows properties for induction vectors
# ==============================================================================
class MTArrows:
    """
    Helper class to read a dictionary of arrow properties

    Arguments:
    -----------
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

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.arrow_size = 2.5
        self.arrow_head_length = 0.15 * self.arrow_size
        self.arrow_head_width = 0.1 * self.arrow_size
        self.arrow_lw = 0.5 * self.arrow_size
        self.arrow_threshold = 2
        self.arrow_color_imag = "c"
        self.arrow_color_real = "k"
        self.arrow_direction = 0

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if v in list(kwargs.keys()):
                setattr(self, v, kwargs.pop(v, None))
