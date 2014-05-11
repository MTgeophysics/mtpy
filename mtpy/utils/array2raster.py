# -*- coding: utf-8 -*-
"""
Created on Sun May 11 12:15:37 2014

@author: jrpeacock
"""
try:
    from osgeo import ogr, gdal, osr
except ImportError:
    raise ImportError('Did not find GDAL, be sure it is installed correctly and '
          'all the paths are correct')
import numpy as np
import mtpy.modeling.modem_new as modem
import os

ogr.UseExceptions()

class ModEM2Raster(object):
    """
    create a raster image of a model slice from a ModEM model
    
    """

    def __init__(self, **kwargs):
        self.model_fn = kwargs.pop('model_fn', None)
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.projection = kwargs.pop('projection', 'WGS84')
        self.origin = kwargs.pop('origin', None)
        
        self.pad_east = None
        self.pad_north = None
        self.res_array = None
        self.cell_size_east = None
        self.cell_size_north = None
        
    def _get_model(self):
        """
        get model to put into array
        """

        model_obj = modem.Model()
        model_obj.model_fn = self.model_fn
        model_obj.read_model_file()
        
        self.cell_size_east = np.median(model_obj.nodes_east)
        self.cell_size_north = np.median(model_obj.nodes_north)
        
        self.pad_east = np.where(model_obj.nodes_east[0:10] > 
                                                    self.cell_size_east*1.1)
        self.pad_north = np.where(model_obj.nodes_north[0:10] > 
                                                    self.cell_size_north*1.1)
        self.grid_z = model_obj.grid_z.copy()                                            
        self.res_array = model_obj.res_model
        
    def write_raster_files(self, save_path=None):
        """
        write a raster file for each layer
        
        """
        if self.origin is None:
            raise ValueError('Need to input an origin as (lon, lat) of the'
                             'southeast corner of the station grid.')
        if save_path is not None:
            self.save_path = save_path
            
        self._get_model()
        
        for ii in range(self.res_array.shape[2]):
            d = self.grid_z[ii]
            raster_fn = os.path.join(save_path, 'Depth_{0:.2f}.tif'.format(d))
            array2raster(raster_fn, self.origin, self.cell_size_east, 
                         self.cell_size_north, self.res_array[:,:,ii],
                         self.projection)            
        
        
#==============================================================================
# create a raster from an array     
#==============================================================================
    
def array2raster(raster_fn, origin, cell_width, cell_height, res_array,
                 projection='WGS84'):
    """
    converts an array into a raster file that can be read into a GIS program.
    
    Arguments:
    --------------
        **raster_fn** : string
                        full path to save raster file to
                        
        **origin** : (lon, lat)
                     longitude and latitude of southwest corner of the array
        
        **cell_width** : float (in meters)
                         size of model cells in east-west direction
                         
        **cell_height** : float (in meters)
                         size of model cells in north-south direction
                         
        **res_array** : np.ndarray(east, north)
                        resistivity array in linear scale.
        
        **projection** : string
                        name of the projection datum
                        
    Output:
    ----------
        * creates a geotiff file projected into projection in UTM.  The 
          values are in log scale for coloring purposes.
    
    """
    res_array = np.log10(res_array[::-1])

    ncols = res_array.shape[1]
    nrows = res_array.shape[0]
    
    utm_cs, utm_point = transform_ll_to_utm(origin[0], origin[1], projection)
    
    origin_east = utm_point[0]
    origin_north = utm_point[1]
    
    # set drive to make a geo tiff
    driver = gdal.GetDriverByName('GTiff')

    # make a raster with the shape of the array to be written    
    out_raster = driver.Create(raster_fn, ncols, nrows, 1, gdal.GDT_Float32)
    out_raster.SetGeoTransform((origin_east, cell_width, 0, 
                                origin_north, 0, cell_height))

    # create a band for the raster data to be put in
    outband = out_raster.GetRasterBand(1)
    outband.WriteArray(res_array)
    
    # geo reference the raster
    #out_raster_georef = osr.SpatialReference()
    #out_raster_georef.ImportFromEPSG(4326)
    out_raster.SetProjection(utm_cs.ExportToWkt())
    
    # be sure to flush the data  
    outband.FlushCache()
 
#==============================================================================
#  transform coordinate systems
#==============================================================================
def transform_ll_to_utm(lon, lat, reference_ellipsoid='WGS84'):    
    def get_utm_zone(longitude):
        return (int(1+(longitude+180.0)/6.0))
    
    def is_northern(latitude):
        """
        Determines if given latitude is a northern for UTM
        """
        if (latitude < 0.0):
            return 0
        else:
            return 1
            
    utm_coordinate_system = osr.SpatialReference()
    # Set geographic coordinate system to handle lat/lon  
    utm_coordinate_system.SetWellKnownGeogCS(reference_ellipsoid) 
    utm_coordinate_system.SetUTM(get_utm_zone(lon), is_northern(lat))
    
    # Clone ONLY the geographic coordinate system 
    ll_coordinate_system = utm_coordinate_system.CloneGeogCS() 
    # create transform component
    ll_to_utm_geo_transform = osr.CoordinateTransformation(ll_coordinate_system, 
                                                          utm_coordinate_system)
                                                                                  
    
    utm_point = ll_to_utm_geo_transform.TransformPoint(lon, lat, 0)
        
    # returns easting, northing, altitude  
    return utm_coordinate_system, utm_point
#==============================================================================
# example test   
#==============================================================================
#rasterOrigin = (-119.000, 37.80)
#pixelWidth = 500
#pixelHeight = 500
#newRasterfn = r'c:\Users\jrpeacock\Documents\test.tif'
#
#
#    
#array = np.array([[ 10, .001, .01, .1, 1, 10, 100, 1000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
#                  [ 1, 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
#                  [ 1, 0.1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 1, 1, 1],
#                  [ 1, 0.5, 1, 1000, 1, 1, 1, 0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1, 1, 1],
#                  [ .1, 0.5, .1, 0.5, 0.5, .1, .1, 0.5, .1, 0.5, .1, 0.5, 0.5, 0.5, .1, 0.5, .1, .1, .1],
#                  [ 1, 0.5, 1, 1, 0.5, 1, 1, 0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1, 1, 1],
#                  [ 1, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 1, 0.5, 1, 0.5, 0.5, 0.5, 1],
#                  [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 20, 22, 24, 26, 27, 30, 40, 10, 1, 1],
#                  [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
#                  [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
#                  
#array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array)
#==============================================================================
# modem test
#==============================================================================
mfn = r"c:\Users\jrpeacock\Documents\ModEM_Mesh"

m_obj = ModEM2Raster()
m_obj.model_fn = mfn
m_obj.origin = (-119.1, 37.85)
m_obj.write_raster_files(save_path=r"c:\Users\jrpeacock\Documents\RasterTest")
