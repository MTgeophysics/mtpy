import os
import os.path as op

os.chdir(r"C:\mtpywin\mtpy")  # change to path to your mtpy installation

from mtpy.modeling.modem import PlotSlices

wd = r"C:\mtpywin\mtpy\examples\model_files\ModEM_2"

savepath = r"C:\tmp"

model_fn = op.join(wd, "Modular_MPI_NLCG_004.rho")
data_fn = op.join(wd, "ModEM_Data.dat")


pObj = PlotSlices(
    model_fn=model_fn,
    data_fn=data_fn,
    save_path=wd,
    cmap="seismic_r",
    plot_stations=True,
    station_id=[5, 8],  # indices (start,finish) of station label to plot as labels
    ew_limits=[
        -220,
        220,
    ],  # option to specify limits, if not provided will auto calculate from data
    ns_limits=[
        -170,
        170,
    ],  # option to specify limits, if not provided will auto calculate from data
    font_size=6,
    fig_size=(6, 3),
    plot_yn="y",
    fig_dpi=400,  # change to your preferred file resolution
)

# pObj.save_path = savepath
# figs,fpaths = pObj.export_slices(plane='N-E', # options are 'N-Z', 'E-Z', and 'N-E'
#                                 indexlist=[32], # depth (or east/west) index to plot
#                                 station_buffer=20e3,
#                                 save=True,
#                                 )
