"""
Script to visualise the 1D model output.

"""
import mtpy.modeling.occam1d as mtoc1d
import os.path as op
import numpy as np
import glob


# directory containing inputs and outputs
idir = r'/home/workshop/MT/Occam1d/07E1_input_output'

# Iteration file name. If not defined, the code will select the last iteration
iterfile = None     #  User can specify the iteration file to plot, e.g.'output_170.iter'   

# plot file name to save to
plot_file_name = '07E1_smooth.png'

# model and data file names
modelfn=op.join(idir,'Model1D')
datafn=op.join(idir,'Occam1d_DataFile_DET.dat')
if iterfile is not None:
	respfile = iterfile.replace('.iter','.resp')

#if not already defined, go to model results directory and find the latest iteration file
if iterfile is None:
    respFiles = glob.glob('%s/*.resp'%(idir))
    iterFiles = glob.glob('%s/*.iter'%(idir))
    
    respDict = {}
    iterDict = {}
    for rf, iterf in zip(respFiles, iterFiles):
    	
    	iterNo = int(rf.split('_')[-1].split('.')[0])
    	respDict[iterNo] = rf
    
    	iterNo = int(iterf.split('_')[-1].split('.')[0])
    	iterDict[iterNo] = iterf
    
    iterfile = iterDict[np.amax(iterDict.keys())]
    respfile = respDict[np.amax(respDict.keys())]
    
# get iteration file to plot
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
                           fig_size=(10,5),
                           depth_units = 'm',
                           override_legend_subscript = 'DET',
                           depth_limits = (0,500),

                           )
pr.axm.set_xlim(1e-1,1e5)
pr.save_figure(op.join(idir,plot_file_name))
