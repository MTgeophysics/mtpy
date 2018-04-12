
import os
import os.path as op
os.chdir(r'C:/mtpywin/mtpy') # change to your path to your mtpy installation to ensure you are using the correct version.

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

from mtpy.modeling.modem import PlotSlices
wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM'

savepath = r'C:\tmp' # change to your desired save path

model_fn = op.join(wd,'Modular_MPI_NLCG_004.rho')
data_fn = op.join(wd,'ModEM_Data.dat')
fs = 8 # fontsize on plot

ps = PlotSlices(model_fn=model_fn, data_fn=data_fn,
                save_path=wd,
                plot_yn='n')

fig = plt.figure(figsize=[6,3])
ax = plt.subplot(111)

gd, gz, gv = ps.get_slice("STA", nsteps=1000)

ci =ax.pcolormesh(gd, gz, gv,
              vmin=1,vmax=1e4,
              norm=colors.LogNorm(),
               cmap='bwr_r', 
               rasterized=True)

ax.invert_yaxis()
ax.set_aspect(1)
ax.set_ylim(10)
plt.setp(ax.get_xticklabels(),fontsize=fs)
plt.setp(ax.get_yticklabels(),fontsize=fs)
plt.xlabel('Distance, km',fontsize=fs)
plt.ylabel('Depth, km',fontsize=fs)

plt.savefig(op.join(savepath,'DepthSlice.png'),
            dpi=400) # change to your desired figure resolution