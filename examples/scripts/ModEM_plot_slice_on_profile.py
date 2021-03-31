import matplotlib.pyplot as plt
from matplotlib import colors

from mtpy.modeling.modem import PlotSlices
from tests import MODEM_DIR, TEST_TEMP_DIR

save_path = TEST_TEMP_DIR.joinpath("ModEM")
if not save_path.exists():
    save_path.mkdir()

model_fn = MODEM_DIR.joinpath("Modular_MPI_NLCG_004.rho")
data_fn = MODEM_DIR.joinpath("ModEM_Data.dat")
fs = 8  # fontsize on plot

ps = PlotSlices(model_fn=model_fn, data_fn=data_fn, save_path=save_path, plot_yn="n")

fig = plt.figure(figsize=[6, 3])
ax = plt.subplot(111)

gd, gz, gv = ps.get_slice("STA", nsteps=1000)

ci = ax.pcolormesh(
    gd, gz, gv, vmin=1, vmax=1e4, norm=colors.LogNorm(), cmap="bwr_r", rasterized=True
)

ax.invert_yaxis()
ax.set_aspect(1)
ax.set_ylim(10)
plt.setp(ax.get_xticklabels(), fontsize=fs)
plt.setp(ax.get_yticklabels(), fontsize=fs)
plt.xlabel("Distance, km", fontsize=fs)
plt.ylabel("Depth, km", fontsize=fs)

# plt.savefig(op.join(savepath,'DepthSlice.png'),
#             dpi=400) # change to your desired figure resolution
