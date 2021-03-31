# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:19:35 2017

@author: u64125
"""

from mtpy.core.mt import MT
import mtpy.analysis.geometry as mtg
import numpy as np

mtObj = MT(r"C:\Git\mtpy\examples\data\edi_files\pb42c.edi")

dimensionality_result = np.array(
    [
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
    ]
)

dimensionality = mtg.dimensionality(z_object=mtObj.Z)

assert np.all(np.abs(dimensionality - dimensionality < 1e-8))
