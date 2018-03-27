# -*- coding: utf-8 -*-
"""
plots edi files as resistivity and phase vs period

CreatedOn:      Wed Sep 18 15:35:39 2013
CreatedBy:      Alison Kirkby

LastUpdated:    2017-01-24
UpdatedBy:      fei.zhang@ga.gov.au

LastUpdated:    2018-03-22  AK updated to use as script rather than command line

"""

from mtpy.core import mt


edi_file = r'C:\mtpywin\mtpy\examples\data\edi_files_2\Synth06.edi'

mt_obj = mt.MT(edi_file)
pt_obj = mt_obj.plot_mt_response(plot_yn='n')
pt_obj.plot()



