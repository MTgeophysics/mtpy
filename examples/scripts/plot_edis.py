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
import tempfile

from mtpy.core import mt

try:
    # PACK_ROOT = os.environ['PACK_ROOT']
    # mtpy_path = os.path.join(PACK_ROOT, 'mtpy')
    mtpy_path = os.environ['MTPY_ROOT']
except:
    print("Define environment variable MTPY_ROOT to be the mtpy source code (clone) directory.")
    raise Exception("MTPY_ROOT var not defined")

os.chdir(mtpy_path) # change to your path to your mtpy installation



edi_path = os.path.join(mtpy_path, 'examples', 'data', 'edi_files')

edi_list = [os.path.join(edi_path, ff) for ff in os.listdir(edi_path) if ff.endswith('.edi')]

temp_dir = tempfile.gettempdir()
print('Using temporary directory ' + temp_dir)
savepath = temp_dir

for edi_file in edi_list:
    mt_obj = mt.MT(edi_file)
    pt_obj = mt_obj.plot_mt_response(plot_yn='n')
    pt_obj.plot()
    pt_obj.save_plot(os.path.join(savepath,
                             os.path.basename(edi_file)[:-4]+'.png'),
                     fig_dpi=400) # change to your preferred file resolution
