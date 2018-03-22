# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: u64125
"""

# import mtpy.modeling.modem_new as mtmn
import os.path as op

from mtpy.imaging.modem_phase_tensor_maps import PlotPTMaps

workdir = r'/Softlab/Githubz/mtpy2/examples/data/ModEM_files'
workdir = r'/g/data/ha3/fxz547/Githubz/mtpy2/examples/data/ModEM_files'
workdir = r'C:\mtpywin\mtpy\examples\data\ModEM_files'


modeldir = op.join(workdir, 'VicSynthetic05')

ptmap = PlotPTMaps(data_fn=op.join(modeldir, 'ModEM_Data_noise10inv.dat'),
                   ellipse_size=40,lw=0.5,color='k',
                   resp_fn=op.join(modeldir, 'Modular_MPI_NLCG_NLCG_015.dat'))

ptmap.plot()
