# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 13:13:29 2016

@author: Alison Kirkby

plot phase tensor ellipses (data, model, residual from outputs from ModEM)

"""
import os
import os.path as op

os.chdir(r'C:\mtpywin\mtpy') # change to path to your mtpy installation

from mtpy.modeling.modem.phase_tensor_maps import PlotPTMaps

scripts = op.dirname(__file__)
examples = op.dirname(scripts)
data = op.join(examples, 'data')

modemFiles = op.join(data, 'ModEM_files')
workdir = op.join(modemFiles, 'VicSynthetic05')
savepath=r'C:\tmp'

ptmap = PlotPTMaps(data_fn=op.join(workdir, 'ModEM_Data_noise10inv.dat'),
                   resp_fn=op.join(workdir, 'Modular_MPI_NLCG_NLCG_015.dat'), # comment out to plot data only
                   #cb_pt_pad=0.1,
                   ellipse_dict =  {'size': 40,
                                    'ellipse_range':[0, 90],
                                    'ellipse_colorby':'phimin',
                                    'ellipse_cmap':'mt_bl2gr2rd'},
                   residual_cmap='mt_wh2or'
                   )

ptmap.plot(period=20, # index of period to plot
           edgecolor='k', # colour for edge of ellipses
           lw = 0.5 # linewidth of edge of ellipses
           )

ptmap.save_figure(save_path=savepath,
                  file_format='png',
                  fig_dpi=400) # change to your desired resolution