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
import scipy.interpolate as interpolate

ogr.UseExceptions()

class ModEM2Raster(object):
    """
    create a raster image of a model slice from a ModEM model
    
    :Example: ::
        >>> import mtpy.utils.array2raster as a2r
        >>> mfn = r"/home/ModEM/Inv1/Modular_NLCG_110.rho"
        >>> m_obj = a2r.ModEM2Raster()
        >>> m_obj.model_fn = mfn
        >>> m_obj.lower_left_corner = (-119.11, 37.80)
        >>> m_obj.write_raster_files(save_path=r"/home/ModEM/Inv1/GIS_depth_slices")

    
    """

    def __init__(self, **kwargs):
        self.model_fn = kwargs.pop('model_fn', None)
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.projection = kwargs.pop('projection', 'WGS84')
        self.lower_left_corner = kwargs.pop('lower_left_corner', None)
        
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
                                    self.cell_size_east*1.1)[0][-1]
        self.pad_north = np.where(model_obj.nodes_north[0:10] > 
                                    self.cell_size_north*1.1)[0][-1]
        self.grid_z = model_obj.grid_z.copy()                                            
        self.res_array = model_obj.res_model[self.pad_north:-self.pad_north,
                                             self.pad_east:-self.pad_east,
                                             :]
                                             
    def interpolate_grid(self, pad_east=None, pad_north=None, cell_size=None):
        """
        interpolate the irregular model grid onto a regular grid.
        
        """
        
        model_obj = modem.Model()
        model_obj.model_fn = self.model_fn
        model_obj.read_model_file()

        self.grid_z = model_obj.grid_z.copy()                                            

        if cell_size is not None:
            self.cell_size_east = cell_size
            self.cell_size_north = cell_size
        else:
            self.cell_size_east = np.median(model_obj.nodes_east)
            self.cell_size_north = np.median(model_obj.nodes_north)
        
        if pad_east is not None:
            self.pad_east = pad_east
        else:
            self.pad_east = np.where(model_obj.nodes_east[0:10] > 
                                     self.cell_size_east*1.1)[0][-1]
        if pad_north is not None:
            self.pad_north = pad_north
        else:
            self.pad_north = np.where(model_obj.nodes_north[0:10] > 
                                    self.cell_size_north*1.1)[0][-1]
        
            
        new_east = np.arange(model_obj.grid_east[self.pad_east],
                             model_obj.grid_east[-self.pad_east-1],
                             self.cell_size_east)
        new_north = np.arange(model_obj.grid_north[self.pad_north],
                             model_obj.grid_north[-self.pad_north-1],
                             self.cell_size_north)
            
        model_n, model_e = np.broadcast_arrays(model_obj.grid_north[:, None], 
                                                      model_obj.grid_east[None, :])

                                             
        new_res_arr = np.zeros((new_north.shape[0],
                                new_east.shape[0],
                                model_obj.grid_z.shape[0]))
                                
        for z_index in range(model_obj.grid_z.shape[0]):
            res = model_obj.res_model[:, :, z_index]
            new_res_arr[:, :, z_index] = interpolate.griddata(
                                         (model_n.ravel(), model_e.ravel()),
                                         res.ravel(), 
                                         (new_north[:, None], new_east[None, :]))
            
#        #1) first need to make x, y, z have dimensions (nx, ny, nz), similar to res
#        north, east, vert = np.broadcast_arrays(model_obj.grid_north[:, None, None], 
#                                                model_obj.grid_east[None, :, None], 
#                                                model_obj.grid_z[None, None, :])
#        
#        #2) next interpolate ont the new mesh
#        new_res = interpolate.griddata((north.ravel(), 
#                                        east.ravel(), 
#                                        vert.ravel()),
#                                        model_obj.res_model.ravel(),
#                                        (new_north[:, None, None], 
#                                         new_east[None, :, None], 
#                                         model_obj.grid_z[None, None, :]),
#                                         method='linear')
        self.res_array = new_res_arr
        
    def write_raster_files(self, save_path=None, pad_east=None, 
                           pad_north=None, cell_size=None):
        """
        write a raster file for each layer
        
        """
        if self.lower_left_corner is None:
            raise ValueError('Need to input an lower_left_corner as (lon, lat)')
        if save_path is not None:
            self.save_path = save_path
            
        self.interpolate_grid(pad_east=pad_east, pad_north=pad_north, 
                              cell_size=cell_size)
        
        for ii in range(self.res_array.shape[2]):
            d = self.grid_z[ii]
            raster_fn = os.path.join(save_path, 'Depth_{0:.2f}_{1}.tif'.format(d, 
                                     self.projection))
            array2raster(raster_fn, 
                         self.lower_left_corner, 
                         self.cell_size_east, 
                         self.cell_size_north, 
                         np.log10(self.res_array[:,:,ii]),
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
    res_array = np.flipud(res_array[::-1])

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
#mfn = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\Modular_NLCG_110.rho"
#
#m_obj = ModEM2Raster()
#m_obj.model_fn = mfn
#m_obj.origin = (-119.11, 37.80)
#m_obj.write_raster_files(save_path=r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\GIS_depth_slices")
