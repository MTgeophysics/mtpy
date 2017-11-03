# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:35:39 2013

@author: Alison Kirkby

plots resistivity and phase as a coloured pseudo section (distance along profile vs period) 
"""
import os

from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection
import os.path as op
import matplotlib.pyplot as plt
plt.ion() # make figure disappear automatically:
#plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc

def test_func():
    """
    test fun
    :return:
    """

    # path to edis
    epath = r'E:\Githubz\mtpy\examples\data\edi_files'

    save_path = r'E:\Githubz\mtpy\temp\resphase.png'

    elst = [op.join(epath,edi) for edi in os.listdir(epath) if edi.endswith('.edi')][::4]

    print (elst)
    resphase = PlotResPhasePseudoSection(fn_list=elst)

    resphase.save_plot(save_path)

    assert (os.path.exists(save_path))

if __name__ == "__main__":
    test_func()
