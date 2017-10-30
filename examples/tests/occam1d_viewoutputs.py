"""
Script to visualise the 1D model output.

"""
import os
#os.chdir(r'C:\Users\alktrt\Documents\Git\mtpy')
import mtpy.modeling.occam1d as mtoc1d
import os.path as op

import numpy as np
import matplotlib.pyplot as plt

# working directory
wd = r'M:\Presentation\Cover thickness workshop_MT\Occam1d\occam1d_files\model01'

# inversion directory
idir = wd

# model and data file names
modelfn=op.join(idir,'Model1D')
datafn=op.join(idir,'Occam1d_DataFile_DET.dat')

# list of all the 'smooth' files, exclude the log file
fn_list = np.array([ff for ff in os.listdir(idir)])

#go to model results directory and find the latest iteration file
#iterfile = 'ITER_135.iter'
#respfile = 'ITER_135.resp'
iterfile = 'Smooth_156.iter'
respfile = 'Smooth_156.resp'


# get maximum iteration
iterfn  = op.join(idir,iterfile)
respfn  = op.join(idir,respfile)

# read in the model, don't need these lines to view the output but useful if you want to analyse the data
oc1m = mtoc1d.Model(model_fn=modelfn)
oc1m.read_iter_file(iterfn)

# read in the data file
oc1d = mtoc1d.Data(data_fn=datafn)
oc1d.read_data_file(data_fn=datafn)
oc1d.read_resp_file(resp_fn=respfn,data_fn=datafn)

# plot the model output
pr = mtoc1d.Plot1DResponse(data_te_fn = datafn,
                           data_tm_fn = datafn,
                           model_fn = modelfn,
                           resp_te_fn = respfn,
                           iter_te_fn = iterfn,
                           resp_tm_fn = respfn,
                           iter_tm_fn = iterfn,
                           depth_limits = (0,1),

                           )
pr.axm.set_xlim(1e-1,1e3)
plt.savefig(op.join(idir,'occam1d_07-E1_minRMS.png'),dpi=300)