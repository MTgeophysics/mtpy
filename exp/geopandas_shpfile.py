import os,sys
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'

mpl.rcParams['figure.figsize']=[30,10]

import geopandas as gpd

def plot_shapefile(shpfile):
    shpf = gpd.GeoDataFrame.from_file(shpfile)

    print("columns and shape", shpf.columns, shpf.shape)

    shpf.plot(figsize=[10,8],linewidth=2)
    plt.show()


if __name__ == "__main__":


    if len(sys.argv)<2:
        sys.exit(1)
    elif os.path.isfile(sys.argv[1]):
        path2shp = sys.argv[1]
        plot_shapefile(path2shp)

    elif os.path.isdir(sys.argv[1]):
        shpfiles= glob.glob(sys.argv[1]+"/PT*.shp")
        for afile in shpfiles:
            plot_shapefile(afile)
    else:
        print("invalid input", sys.argv[1])
