# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: u64125
"""

import os.path as op

from mtpy.imaging.modem_phase_tensor_maps import PlotPTMaps

workdir = r'C:\mtpywin\mtpy\examples\data\ModEM_files\VicSynthetic05'


ptmap = PlotPTMaps(data_fn=op.join(workdir, 'ModEM_Data_noise10inv.dat'),
                   resp_fn=op.join(workdir, 'Modular_MPI_NLCG_NLCG_015.dat'), # comment out to plot data only
                   ellipse_size=40 # scaling factor for ellipses
                   )

ptmap.plot(period=20, # index of period to plot
           edgecolor='k' # colour for edge of ellipses
           )
