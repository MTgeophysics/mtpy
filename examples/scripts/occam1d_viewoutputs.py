# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 08:29:24 2016

@author: Alison Kirkby

View the data, response and model output from Occam1D

"""

import mtpy.modeling.occam1d as mtoc1d
import os.path as op

# working directory
wd = r'C:\mtpywin\mtpy\examples\model_files\Occam1d'

# model and data file names
modelfn=op.join(wd,'Model1D')
datafn=op.join(wd,'Occam1d_DataFile_DET.dat')

iterfn  = op.join(wd,'ITER_97.iter')
respfn  = op.join(wd,'ITER_97.resp')

# read in the model, don't need these lines to view the model but useful if you want to analyse the data
oc1m = mtoc1d.Model(model_fn=modelfn)
oc1m.read_iter_file(iterfn)


# read in the data file
oc1d = mtoc1d.Data(data_fn=datafn)
oc1d.read_data_file(data_fn=datafn)
oc1d.read_resp_file(resp_fn=respfn,data_fn=datafn)

# plot the model and response
pr = mtoc1d.Plot1DResponse(data_te_fn = datafn,
                           data_tm_fn = datafn,
                           model_fn = modelfn,
                           resp_te_fn = respfn,
                           iter_te_fn = iterfn,
                           resp_tm_fn = respfn,
                           iter_tm_fn = iterfn,
                           depth_limits = (0,10)
                           )