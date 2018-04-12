# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: u64125
"""

import os.path as op

from mtpy.modeling.modem.phase_tensor_maps import PlotPTMaps

scripts = op.dirname(__file__)
examples = op.dirname(scripts)
data = op.join(examples, 'data')

modemFiles = op.join(data, 'ModEM_files')
workdir = op.join(modemFiles, 'VicSynthetic05')
savepath='/tmp'

ptmap = PlotPTMaps(data_fn=op.join(workdir, 'ModEM_Data_noise10inv.dat'),
                   resp_fn=op.join(workdir, 'Modular_MPI_NLCG_NLCG_015.dat'), # comment out to plot data only
                   #cb_pt_pad=0.1,
                   cb_tick_step=9,
                   ellipse_dict =  {'size': 40,
                                    'ellipse_range':[-9, 9],
                                    'ellipse_colorby':'skew',
                                    'ellipse_cmap':'mt_bl2gr2rd'},
                   residual_cmap='mt_wh2or'
                   )

ptmap.plot(period=20, # index of period to plot
           edgecolor='none', # colour for edge of ellipses
           )

ptmap.save_figure(op.join(savepath,'PTMap'),
                  file_format='png',
                  fig_dpi=400) # change to your desired resolution