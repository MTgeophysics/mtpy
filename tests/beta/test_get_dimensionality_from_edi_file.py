# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:19:35 2017

@author: u64125
"""

import numpy as np

import mtpy.analysis.geometry as mtg
from mtpy.core.mt import MT


def test_fun():
    """
    test function
    :return: T/F
    """
    # mtObj = MT(r'C:\Git\mtpy\examples\data\edi_files\pb42c.edi')
    mtObj = MT(r'E:/Githubz/mtpy/examples/data/edi_files/pb42c.edi')
    dimensionality_result = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                                      2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

    dimensionality = mtg.dimensionality(z_object=mtObj.Z)
    differ = np.abs(dimensionality - dimensionality_result)

    print differ

    assert np.all(np.abs(dimensionality - dimensionality_result) < 1e-8)

if __name__ == "__main__":
    test_fun()
