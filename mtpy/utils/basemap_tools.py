
import numpy as np
import matplotlib.pyplot as plt

def add_basemap_frame(basemap,tick_interval=2, coastline_kwargs={},states_kwargs={},
                      mlabels=[False,False,False,True],plabels=[True,False,False,False]):
    """
    add a standard map frame (lat/lon labels and tick marks, coastline and states) to basemap
    
    :param tick_interval: tick interval in degrees
    :param coastline_kwargs: dictionary containing arguments to pass into the drawcoastlines function
    :param states_kwargs: dictionary containing arguments to pass into the drawstates function
    :param mlabels: where to place meridian (longitude) labels on plot (list containing True/False for [left,right,top,bottom])
    :param plabels: where to place parallels (latitudes) labels on plot (list containing True/False for [left,right,top,bottom])
    """
    
    basemap.drawmeridians(np.arange(np.floor(basemap.lonmin),np.ceil(basemap.lonmax),2),labels=mlabels)#
    basemap.drawparallels(np.arange(np.floor(basemap.latmin),np.ceil(basemap.latmax),2),labels=plabels)#

    basemap.drawcoastlines(**coastline_kwargs)
    basemap.drawstates(**states_kwargs)


def plot_data(x,y,values,basemap=None,cbar=False,**param_dict):
    """
    plot array data, either 1d or 2d
    
    :param x: x position of points
    :param y: y position of points
    :param values: values to plot, if 1D, a scatter plot will be made, if 2D, a pcolormesh plot will be made
    :param basemap: supply a basemap, if None, data will be plotted on current axes
    :param cbar: True/False, whether or not to show a colorbar
    
    """
    

    if len(np.shape(values)) == 1:
        if basemap is None:
            # plot a scatter plot with values coloured
            plt.scatter(x,y,c=values,**param_dict)
        else:
            x,y = basemap(x,y)
            basemap.scatter(x,y,c=values,**param_dict)
    elif len(np.shape(values)) == 2:
        if basemap is None:
            # plot a pcolormesh plot
            plt.pcolormesh(x,y,values,**param_dict)
        else:
            x,y = basemap(x,y)
            basemap.pcolormesh(x,y,values,**param_dict)

    plt.gca().set_aspect(1)
    if cbar:
        plt.colorbar(shrink=0.5)