# -*- coding: utf-8 -*-
"""
plots edi files as resistivity and phase vs period
Only works for Windows installation due to hard-coded paths

CreatedOn:      Wed Sep 18 15:35:39 2013
CreatedBy:      Alison Kirkby

LastUpdated:    2017-01-24
UpdatedBy:      fei.zhang@ga.gov.au

LastUpdated:    2018-03-22  AK updated to use as script rather than command line

"""
import os
import os.path as op

os.chdir(r'C:\mtpywin\mtpy') # change to your path to your mtpy installation

from mtpy.core import mt


edi_path = r'C:\mtpywin\mtpy\examples\data\edi_files'

edi_list = [op.join(edi_path,ff) for ff in os.listdir(edi_path) if ff.endswith('.edi')]

savepath = r'C:\tmp'

for edi_file in edi_list:
    mt_obj = mt.MT(edi_file)
    pt_obj = mt_obj.plot_mt_response(plot_yn='n')
    pt_obj.plot()
    pt_obj.save_plot(op.join(savepath,
                             op.basename(edi_file)[:-4]+'.png'),
                     fig_dpi=400) # change to your preferred file resolution


