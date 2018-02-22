# -*- coding: utf-8 -*-
"""
plots phase tensor ellipses as a pseudo section (distance along profile vs period)

CreatedOn:      Wed Sep 18 15:35:39 2013
CreatedBy:      Alison Kirkby


"""


import os
import os.path as op
from mtpy.core import mt


edi_path = r'C:\Git\mtpy\data\edifiles2'
edi_list = [op.join(edi_path,ff) for ff in os.listdir(edi_path) if ff.endswith('.edi')]

for edi_file in edi_list:
    mt_obj = mt.MT(edi_file)
    pt_obj = mt_obj.plot_mt_response(plot_yn='n')
    pt_obj.plot()

