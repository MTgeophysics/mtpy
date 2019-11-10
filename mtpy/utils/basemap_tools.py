
import numpy as np
import matplotlib.pyplot as plt
from mtpy.modeling.modem import Data
from mtpy.utils.calculator import nearest_index

        

def get_latlon_extents_from_modem_data(stations_obj):

    return stations_obj.lon.min(),stations_obj.lon.max(),stations_obj.lat.min(),stations_obj.lat.max()
        

def compute_tick_interval_from_map_extent(lonMin,lonMax,latMin,latMax):
    """
    estimate an even tick interval based on map extent based on some sensible options
    
    """
    tick_options = np.array([1.,2.,5.])
    tick_options = np.hstack([tick_options*0.01,
                              tick_options*0.1,
                              tick_options,
                              tick_options*10.])
    
    # return nearest interval from list of options, to the average of latitude
    # extent and longitude extent
    return tick_options[nearest_index((lonMax-lonMin+latMax-latMin)/10.,tick_options)]
    


def compute_map_extent_from_modem_data(stations_obj,buffer=None,buffer_factor=0.1):
    """
    compute extent for a plot from data extent from ModEM data file
    
    :param data_fn: full path to modem data file
    :param buffer: optional argument; buffer in latitude/longitude (if not provided, 
    this is assumed to be a fraction of the maximum of the north-south or east-west extent)
    :param buffer_factor: fraction of north-south or east-west extent for buffer (if buffer not provided)
    
    """

    lonMin, lonMax, latMin, latMax = get_latlon_extents_from_modem_data(stations_obj)
    
    # compute buffer
    if buffer is None:
        buffer = max([(lonMax-lonMin)*buffer_factor,(latMax-latMin)*buffer_factor])
        
    return lonMin - buffer, lonMax + buffer, latMin - buffer, latMax + buffer
    

def compute_lonlat0_from_modem_data(stations_obj):
    """
    compute lat0 and lon0 for creating a basemap, using data centre point in modem data file
    """
    
    return stations_obj.center_point['lon'], stations_obj.center_point['lat']


def initialise_basemap(stations_obj,buffer=None,**basemap_kwargs):
    """
    create a new basemap instance
    
    """
    
    from mpl_toolkits.basemap import Basemap

    lonMin, lonMax, latMin, latMax = compute_map_extent_from_modem_data(stations_obj,buffer=buffer)
    lon_0, lat_0 = compute_lonlat0_from_modem_data(stations_obj)
        
    # update basemap arguments with defaults if not provided
    basemap_kwargs['llcrnrlon'] = basemap_kwargs.pop('llcrnrlon',lonMin)
    basemap_kwargs['urcrnrlon'] = basemap_kwargs.pop('urcrnrlon',lonMax)
    basemap_kwargs['llcrnrlat'] = basemap_kwargs.pop('llcrnrlat',latMin)
    basemap_kwargs['urcrnrlat'] = basemap_kwargs.pop('urcrnrlat',latMax)
    basemap_kwargs['lat_0'] = basemap_kwargs.pop('lat_0',lat_0)
    basemap_kwargs['lon_0'] = basemap_kwargs.pop('lon_0',lon_0)
    basemap_kwargs['resolution'] = basemap_kwargs.pop('resolution','l')
    basemap_kwargs['projection'] = basemap_kwargs.pop('projection','cyl')

    return Basemap(**basemap_kwargs)



def add_basemap_frame(basemap,tick_interval=None, coastline_kwargs={},states_kwargs={},
                      mlabels=[False,False,False,True],plabels=[True,False,False,False]):
    """
    add a standard map frame (lat/lon labels and tick marks, coastline and states) to basemap
    
    :param tick_interval: tick interval in degrees
    :param coastline_kwargs: dictionary containing arguments to pass into the drawcoastlines function
    :param states_kwargs: dictionary containing arguments to pass into the drawstates function
    :param mlabels: where to place meridian (longitude) labels on plot (list containing True/False for [left,right,top,bottom])
    :param plabels: where to place parallels (latitudes) labels on plot (list containing True/False for [left,right,top,bottom])
    """
    if tick_interval is None:
        tick_interval = compute_tick_interval_from_map_extent(basemap.lonmin,basemap.lonmax,basemap.latmin,basemap.latmax)
        print("tick_interval",tick_interval)
    
    basemap.drawmeridians(np.arange(np.floor(basemap.lonmin),np.ceil(basemap.lonmax),tick_interval),labels=mlabels)#
    basemap.drawparallels(np.arange(np.floor(basemap.latmin),np.ceil(basemap.latmax),tick_interval),labels=plabels)#

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