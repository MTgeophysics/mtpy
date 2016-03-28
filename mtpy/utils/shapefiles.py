# -*- coding: utf-8 -*-
"""
Create shape files for phase tensor ellipses.

Created on Sun Apr 13 12:32:16 2014

@author: jrpeacock
"""
try:
    from osgeo import ogr, gdal, osr
except ImportError:
    raise ImportError('Did not find GDAL, be sure it is installed correctly and '
          'all the paths are correct')
import numpy as np
import os
import mtpy.core.mt as mt
import mtpy.modeling.modem_new as modem
import mtpy.analysis.pt as mtpt

ogr.UseExceptions()

class PTShapeFile(object):
    """
    write shape file for GIS plotting programs
    
    ======================== ==================================================
    key words/attributes      Description
    ======================== ==================================================
    edi_list                 list of edi files, full paths
    ellipse_size             size of normalized ellipse in map scale
                             *default* is .01
    mt_obj_list              list of mt.MT objects
                             *default* is None, filled if edi_list is given
    plot_period              list or value of period to convert to shape file
                             *default* is None, which will write a file for
                             every period in the edi files 
    ptol                     tolerance to look for given periods
                             *default* is .05
    pt_dict                  dictionary with keys of plot_period.  Each
                             dictionary key is a structured array containing
                             the important information for the phase tensor.
    projection               projection of coordinates see EPSG for all options
                             *default* is WSG84 in lat and lon
    save_path                path to save files to
                             *default* is current working directory.
    ======================== ==================================================
    
    
    ======================== ==================================================
    Methods                   Description
    ======================== ==================================================
    _get_plot_period         get a list of all frequencies possible from
                             input files
    _get_pt_array            get phase tensors from input files and put the
                             information into a structured array 
    write_shape_files        write shape files based on attributes of class
    ======================== ==================================================

    * This will project the data into UTM WSG84
    :Example: ::
        >>> edipath = r"/home/edi_files_rotated_to_geographic_north"
        >>> edilist = [os.path.join(edipath, edi) \
                      for edi in os.listdir(edipath)\
                      if edi.find('.edi')>0]
        >>> pts = PTShapeFile(edilist, save_path=r"/home/gis")
        >>> pts.write_shape_files()
        
    * To project into another datum, set the projection attribute
    :Example: ::
        >>> pts = PTShapeFile(edilist, save_path=r"/home/gis")
        >>> pts.projection = 'NAD27'
        >>> pts.write_shape_files()
        
  
    """
    
    def __init__(self, edi_list=None, **kwargs):
        self.edi_list = edi_list
        self.projection = 'WGS84'
        self.plot_period = None
        self.save_path = os.getcwd()
        self.ellipse_size = 500.0
        self._theta = np.arange(0, 2*np.pi, np.pi/180.)
        self.ptol = .05
        
        self.mt_obj_list = None
        self.pt_dict = None
        
        if self.edi_list is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edi_list]
            
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        
        if self.mt_obj_list is not None:
            self._get_plot_period()
            self._get_pt_array()
        
        self._proj_dict = {'WGS84':4326, 'NAD27':4267}
        
        self.utm_cs = None
        self._rotation_angle = 0.0
        
    def _set_rotation_angle(self, rotation_angle):
        """
        rotate all mt_objs to rotation angle
        """
        
        self._rotation_angle = float(rotation_angle)
        
        for mt_obj in self.mt_obj_list:
            mt_obj.rotation_angle = float(self._rotation_angle)
            
    def _get_rotation_angle(self):
        return self._rotation_angle
        
    rotation_angle = property(_get_rotation_angle, _set_rotation_angle,
                              doc="rotation angle of Z and Tipper")
    
    def _get_plot_period(self):
        """
        from the list of edi's get a frequency list from all possible 
        frequencies.
        
        """

        if self.plot_period is None:
            #get all frequencies from all edi files
            all_freqs = []
            for mt_obj in self.mt_obj_list:
                all_freqs.extend(list(mt_obj.Z.freq))
            
            #sort all frequencies so that they are in descending order,
            #use set to remove repeats and make an array
            self.plot_period = 1./np.array(sorted(list(set(all_freqs)),
                                        reverse=True))
        else:
            if type(self.plot_period) is list:
                pass
            if type(self.plot_period) is int or type(self.plot_period) is float:
                self.plot_period = [self.plot_period]
                
    def _get_pt_array(self):
        """
        get the phase tensor information into a form that is more structured
        to manipulate easier later.
        
        make a dictionary with keys being the plot period values and each
        key has a structured array that contains all the important information
        collected from each station.
        """
        self.pt_dict = {}  
        if self.plot_period is None:
            self._get_plot_period()
            
        for plot_per in self.plot_period:
            self.pt_dict[plot_per] = []
            for mt_obj in self.mt_obj_list:
                try:
                    p_index = [ff for ff, f2 in enumerate(1./mt_obj.Z.freq) 
                               if (f2 > plot_per*(1-self.ptol)) and
                                  (f2 < plot_per*(1+self.ptol))][0]
                    if self.projection is None:
                        east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                        self.utm_cs = osr.SpatialReference()
                        # Set geographic coordinate system to handle lat/lon  
                        self.utm_cs.SetWellKnownGeogCS(self.projection)
                    else:
                        self.utm_cs, utm_point = transform_ll_to_utm(mt_obj.lon, 
                                                                mt_obj.lat,
                                                                self.projection)
                        east, north, elev = utm_point

                    pt_tuple = (mt_obj.station, east, north,
                                mt_obj.pt.phimin[0][p_index],
                                mt_obj.pt.phimax[0][p_index],
                                mt_obj.pt.azimuth[0][p_index],
                                mt_obj.pt.beta[0][p_index],
                                2*mt_obj.pt.beta[0][p_index])           
                    self.pt_dict[plot_per].append(pt_tuple)
                except IndexError:
                    pass
                
            self.pt_dict[plot_per] = np.array(self.pt_dict[plot_per],
                                              dtype=[('station', '|S15'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('azimuth', np.float),
                                                     ('skew', np.float),
                                                     ('n_skew', np.float)])
            
                
    def write_shape_files(self, ):
        """
        write shape file from given attributes
        """
        
        self._get_pt_array()
            
        for plot_per in self.plot_period:
            #shape file path
            shape_fn = os.path.join(self.save_path, 
                                    'PT_{0:.5g}s_{1}.shp'.format(plot_per,
                                    self.projection))
            
            # remove the shape file if it already exists, has trouble over writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)
            
            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')
            
            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)
                
            #create shape file
            data_source = driver.CreateDataSource(shape_fn)
            
#            ##if you read from a raster get the georeference point otherwise create one
#            spatial_ref = osr.SpatialReference()
#            #this puts it in the wsg84 reference frame.
#            spatial_ref.ImportFromEPSG(self._proj_dict[self.projection]) 
            
            ##create a layer to put the ellipses onto
            layer = data_source.CreateLayer('PT', self.utm_cs, ogr.wkbPolygon)
            
            #make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)
            
            field_phimin = ogr.FieldDefn('phi_min', ogr.OFTReal)
            layer.CreateField(field_phimin)
            
            field_phimax = ogr.FieldDefn('phi_max', ogr.OFTReal)
            layer.CreateField(field_phimax)
            
            field_skew = ogr.FieldDefn('skew', ogr.OFTReal)
            layer.CreateField(field_skew)
            
            field_normalized_skew = ogr.FieldDefn('n_skew', ogr.OFTReal)
            layer.CreateField(field_normalized_skew)
            
            poly_list = []
            phimax = self.pt_dict[plot_per]['phimax'].max()
            for pt_array in self.pt_dict[plot_per]:
                
                #need to make an ellipse first using the parametric equation
                azimuth = -np.deg2rad(pt_array['azimuth'])
                width = self.ellipse_size*(pt_array['phimax']/phimax)
                height = self.ellipse_size*(pt_array['phimin']/phimax) 
                x0 = pt_array['east']
                y0 = pt_array['north']
                
                x = x0+height*np.cos(self._theta)*np.cos(azimuth)-\
                       width*np.sin(self._theta)*np.sin(azimuth)  
                y = y0+height*np.cos(self._theta)*np.sin(azimuth)+\
                       width*np.sin(self._theta)*np.cos(azimuth) 
                       
                #1) make a geometry shape of the ellipse
                ellipse = ogr.Geometry(ogr.wkbLinearRing)
                
                for ii, jj in zip(x, y):
                    ellipse.AddPoint(np.round(ii, 6), np.round(jj, 6))
                    
                ellipse.CloseRings()
                    
                #2) make a polygon
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ellipse)
                
                poly_list.append(poly)
            
            
                ##4) this part is confusing but we need to create a feature that has the
                ##   same definition as the layer that we created.
                # get the layer definition
                feature_def = layer.GetLayerDefn()
                
                #create a new feature
                new_feature = ogr.Feature(feature_def)
                #set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                #create the feature in the layer. 
                layer.CreateFeature(new_feature)
                
                #
                ###5) create a field to color by
                new_feature.SetField("Name", pt_array['station'])
                new_feature.SetField("phi_min", pt_array['phimin'])
                new_feature.SetField("phi_max", pt_array['phimax'])
                new_feature.SetField("skew", pt_array['skew'])
                new_feature.SetField("n_skew", pt_array['n_skew'])
            
                #add the new feature to the layer.
                layer.SetFeature(new_feature)
                
                #apparently need to destroy the feature
                new_feature.Destroy()
                
                
            # Need to be sure that all the new info is saved to 
            data_source.SyncToDisk()
            
            #write a projection file
#            spatial_ref.MorphToESRI()
            self.utm_cs.MorphToESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
#            prj_file.write(spatial_ref.ExportToWkt())
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()
            
            data_source.Destroy()
            
            print 'Wrote shape file to {0}'.format(shape_fn)
            
    def write_data_pt_shape_files_modem(self, modem_data_fn, 
                                        rotation_angle=0.0):
        """
        write pt files from a modem data file.
        
        """
            
        modem_obj = modem.Data()
        modem_obj.read_data_file(modem_data_fn)
        
        self.plot_period = modem_obj.period_list.copy()
        self.mt_obj_list = [modem_obj.mt_dict[key] 
                            for key in modem_obj.mt_dict.keys()]
                                
        self._set_rotation_angle(rotation_angle)
        
        
        self.write_shape_files()
        
    def write_resp_pt_shape_files_modem(self, modem_data_fn, modem_resp_fn,
                                        rotation_angle=0.0):
        """
        write pt files from a modem response file where ellipses are normalized
        by the data file.
        
        """

        #first get the data and response and place them in array for later use
        modem_data_obj = modem.Data()
        modem_data_obj.read_data_file(modem_data_fn)
        
        self.plot_period = modem_data_obj.period_list.copy()
        self.mt_obj_list = [modem_data_obj.mt_dict[key] 
                            for key in modem_data_obj.mt_dict.keys()]
        self._get_pt_array()
        
        self._set_rotation_angle(rotation_angle)
            
        modem_resp_obj = modem.Data()
        modem_resp_obj.read_data_file(modem_resp_fn)
        
        #rotate model response
        for r_key in modem_resp_obj.mt_dict.keys():
            modem_resp_obj.mt_dict[r_key].rotation_angle = float(rotation_angle)
            
        resp_pt_dict = {}        
        for p_index, plot_per in enumerate(self.plot_period):
            resp_pt_dict[plot_per] = []
            for key in modem_data_obj.mt_dict.keys():
                mt_obj = modem_data_obj.mt_dict[key]
                if self.projection is None:
                    east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                    self.utm_cs = osr.SpatialReference()
                    # Set geographic coordinate system to handle lat/lon  
                    self.utm_cs.SetWellKnownGeogCS(self.projection)
                else:
                    self.utm_cs, utm_point = transform_ll_to_utm(mt_obj.lon, 
                                                            mt_obj.lat,
                                                            self.projection)
                    east, north, elev = utm_point
                
                #get pt objects from data and model response
                try:
                    mpt = modem_resp_obj.mt_dict[key].pt
                
                    pt_tuple = (mt_obj.station, east, north,
                                mpt.phimin[0][p_index],
                                mpt.phimax[0][p_index],
                                mpt.azimuth[0][p_index],
                                mpt.beta[0][p_index],
                                2*mpt.beta[0][p_index])
                except KeyError:
                    pt_tuple = (mt_obj.station, east, north,
                                0,
                                0,
                                0,
                                0,
                                0)
                resp_pt_dict[plot_per].append(pt_tuple)

            #now make each period an array for writing to file    
            resp_pt_dict[plot_per] = np.array(resp_pt_dict[plot_per],
                                              dtype=[('station', '|S15'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('azimuth', np.float),
                                                     ('skew', np.float),
                                                     ('n_skew', np.float)])
                
        #write files
        for plot_per in self.plot_period:
            #shape file path
            shape_fn = os.path.join(self.save_path, 
                                    'Resp_PT_{0:.5g}s_{1}.shp'.format(plot_per,
                                    self.projection))
            
            # remove the shape file if it already exists, has trouble over writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)
            
            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')
            
            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)
                
            #create shape file
            data_source = driver.CreateDataSource(shape_fn)
            
            ##create a layer to put the ellipses onto
            layer = data_source.CreateLayer('RPT', self.utm_cs, ogr.wkbPolygon)
            
            #make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)
            
            field_phimin = ogr.FieldDefn('phi_min', ogr.OFTReal)
            layer.CreateField(field_phimin)
            
            field_phimax = ogr.FieldDefn('phi_max', ogr.OFTReal)
            layer.CreateField(field_phimax)
            
            field_skew = ogr.FieldDefn('skew', ogr.OFTReal)
            layer.CreateField(field_skew)
            
            field_normalized_skew = ogr.FieldDefn('n_skew', ogr.OFTReal)
            layer.CreateField(field_normalized_skew)
            
            poly_list = []
            phimax = self.pt_dict[plot_per]['phimax'].max()
            for pt_array in resp_pt_dict[plot_per]:
                
                #need to make an ellipse first using the parametric equation
                azimuth = -np.deg2rad(pt_array['azimuth'])
                width = self.ellipse_size*(pt_array['phimax']/phimax)
                height = self.ellipse_size*(pt_array['phimin']/phimax) 
                x0 = pt_array['east']
                y0 = pt_array['north']
                
                x = x0+height*np.cos(self._theta)*np.cos(azimuth)-\
                       width*np.sin(self._theta)*np.sin(azimuth)  
                y = y0+height*np.cos(self._theta)*np.sin(azimuth)+\
                       width*np.sin(self._theta)*np.cos(azimuth) 
                       
                #1) make a geometry shape of the ellipse
                ellipse = ogr.Geometry(ogr.wkbLinearRing)
                
                for ii, jj in zip(x, y):
                    ellipse.AddPoint(np.round(ii, 6), np.round(jj, 6))
                    
                ellipse.CloseRings()
                    
                #2) make a polygon
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ellipse)
                
                poly_list.append(poly)
            
            
                ##4) this part is confusing but we need to create a feature that has the
                ##   same definition as the layer that we created.
                # get the layer definition
                feature_def = layer.GetLayerDefn()
                
                #create a new feature
                new_feature = ogr.Feature(feature_def)
                #set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                #create the feature in the layer. 
                layer.CreateFeature(new_feature)
                
                #
                ###5) create a field to color by
                new_feature.SetField("Name", pt_array['station'])
                new_feature.SetField("phi_min", pt_array['phimin'])
                new_feature.SetField("phi_max", pt_array['phimax'])
                new_feature.SetField("skew", pt_array['skew'])
                new_feature.SetField("n_skew", pt_array['n_skew'])
            
                #add the new feature to the layer.
                layer.SetFeature(new_feature)
                
                #apparently need to destroy the feature
                new_feature.Destroy()
                
                
            # Need to be sure that all the new info is saved to 
            data_source.SyncToDisk()
            
            #write a projection file
            self.utm_cs.MorphToESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()
            
            data_source.Destroy()
            
            print 'Wrote shape file to {0}'.format(shape_fn)
        
    def write_residual_pt_shape_files_modem(self, modem_data_fn, modem_resp_fn,
                                            rotation_angle=0.0, normalize='1'):
        """
        write residual pt shape files from ModEM output
        
        normalize [ '1' | 'all' ]
                   * '1' to normalize the ellipse by itself, all ellipses are 
                         normalized to phimax, thus one axis is of length 
                         1*ellipse_size
                   * 'all' to normalize each period by the largest phimax
        
        """
        
        #first get the data and response and place them in array for later use
        modem_data_obj = modem.Data()
        modem_data_obj.read_data_file(modem_data_fn)
        
        self.plot_period = modem_data_obj.period_list.copy()
        self.mt_obj_list = [modem_data_obj.mt_dict[key] 
                            for key in modem_data_obj.mt_dict.keys()]
        self._get_pt_array()
        
        self._set_rotation_angle(rotation_angle)
            
        modem_resp_obj = modem.Data()
        modem_resp_obj.read_data_file(modem_resp_fn)
        
        #rotate model response
        for r_key in modem_resp_obj.mt_dict.keys():
            modem_resp_obj.mt_dict[r_key].rotation_angle = float(rotation_angle)
            
        residual_pt_dict = {}        
        for p_index, plot_per in enumerate(self.plot_period):
            residual_pt_dict[plot_per] = []
            for key in modem_data_obj.mt_dict.keys():
                mt_obj = modem_data_obj.mt_dict[key]
                if self.projection is None:
                    east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                    self.utm_cs = osr.SpatialReference()
                    # Set geographic coordinate system to handle lat/lon  
                    self.utm_cs.SetWellKnownGeogCS(self.projection)
                else:
                    self.utm_cs, utm_point = transform_ll_to_utm(mt_obj.lon, 
                                                            mt_obj.lat,
                                                            self.projection)
                    east, north, elev = utm_point
                
                #get pt objects from data and model response
                dpt = modem_data_obj.mt_dict[key].pt
                mpt = modem_resp_obj.mt_dict[key].pt

                #calculate the residual pt
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt, 
                                                   pt_object2=mpt)
                    rpt = rpt.residual_pt
                    rpt_mean = .25*np.linalg.norm(rpt.pt[p_index], ord='fro')
#                    rpt_mean = .25*np.sqrt(abs(rpt.pt[p_index, 0, 0])**2+
#                                          abs(rpt.pt[p_index, 0, 1])**2+
#                                          abs(rpt.pt[p_index, 1, 0])**2+
#                                          abs(rpt.pt[p_index, 1, 1])**2)
                    pt_tuple = (mt_obj.station, east, north,
                                rpt.phimin[0][p_index],
                                rpt.phimax[0][p_index],
                                rpt.azimuth[0][p_index],
                                rpt.beta[0][p_index],
                                rpt_mean)
#                                np.sqrt(abs(rpt.phimin[0][p_index]*
#                                            rpt.phimax[0][p_index])))           
                    residual_pt_dict[plot_per].append(pt_tuple)
                except mtpt.MTex.MTpyError_PT:
                    print key, dpt.pt.shape, mpt.pt.shape
                    pt_tuple = (mt_obj.station, east, north,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0)
                    residual_pt_dict[plot_per].append(pt_tuple)
            #now make each period an array for writing to file    
            residual_pt_dict[plot_per] = np.array(residual_pt_dict[plot_per],
                                              dtype=[('station', '|S15'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('azimuth', np.float),
                                                     ('skew', np.float),
                                                     ('geometric_mean', np.float)])
        #return residual_pt_dict
        #write files
        for plot_per in self.plot_period:
            #shape file path
            shape_fn = os.path.join(self.save_path, 
                                    'ResidualPT_{0:.5g}s_{1}.shp'.format(plot_per,
                                    self.projection))
            
            # remove the shape file if it already exists, has trouble over writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)
            
            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')
            
            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)
                
            #create shape file
            data_source = driver.CreateDataSource(shape_fn)
            
            ##create a layer to put the ellipses onto
            layer = data_source.CreateLayer('RPT', self.utm_cs, ogr.wkbPolygon)
            
            #make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)
            
            field_phimin = ogr.FieldDefn('phi_min', ogr.OFTReal)
            layer.CreateField(field_phimin)
            
            field_phimax = ogr.FieldDefn('phi_max', ogr.OFTReal)
            layer.CreateField(field_phimax)
            
            field_skew = ogr.FieldDefn('skew', ogr.OFTReal)
            layer.CreateField(field_skew)
            
            field_geometric_mean = ogr.FieldDefn('mean', ogr.OFTReal)
            layer.CreateField(field_geometric_mean)
            
            poly_list = []
            phimax = self.pt_dict[plot_per]['phimax'].max()
            for pt_array in residual_pt_dict[plot_per]:
                #need to make an ellipse first using the parametric equation
                azimuth = -np.deg2rad(pt_array['azimuth'])
                if normalize == '1':
                    width = self.ellipse_size*(pt_array['phimax']/pt_array['phimax'])
                    height = self.ellipse_size*(pt_array['phimin']/pt_array['phimax']) 
                elif normalize == 'all':
                    width = self.ellipse_size*(pt_array['phimax']/phimax)
                    height = self.ellipse_size*(pt_array['phimin']/phimax) 
                x0 = pt_array['east']
                y0 = pt_array['north']
                
                x = x0+height*np.cos(self._theta)*np.cos(azimuth)-\
                       width*np.sin(self._theta)*np.sin(azimuth)  
                y = y0+height*np.cos(self._theta)*np.sin(azimuth)+\
                       width*np.sin(self._theta)*np.cos(azimuth) 
                       
                #1) make a geometry shape of the ellipse
                ellipse = ogr.Geometry(ogr.wkbLinearRing)
                
                for ii, jj in zip(x, y):
                    ellipse.AddPoint(np.round(ii, 6), np.round(jj, 6))
                    
                ellipse.CloseRings()
                    
                #2) make a polygon
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ellipse)
                
                poly_list.append(poly)
            
            
                ##4) this part is confusing but we need to create a feature that has the
                ##   same definition as the layer that we created.
                # get the layer definition
                feature_def = layer.GetLayerDefn()
                
                #create a new feature
                new_feature = ogr.Feature(feature_def)
                #set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                #create the feature in the layer. 
                layer.CreateFeature(new_feature)
                
                #
                ###5) create a field to color by
                new_feature.SetField('Name', pt_array['station'])
                new_feature.SetField('phi_min', pt_array['phimin'])
                new_feature.SetField('phi_max', pt_array['phimax'])
                new_feature.SetField('skew', pt_array['skew'])
                new_feature.SetField('mean', pt_array['geometric_mean'])
            
                #add the new feature to the layer.
                layer.SetFeature(new_feature)
                
                #apparently need to destroy the feature
                new_feature.Destroy()
                
                
            # Need to be sure that all the new info is saved to 
            data_source.SyncToDisk()
            
            #write a projection file
            self.utm_cs.MorphToESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()
            
            data_source.Destroy()
            
            print 'Wrote shape file to {0}'.format(shape_fn)
        
#==============================================================================
# Tipper arrows
#==============================================================================
class TipperShapeFile(object):
    """
    write shape file for GIS plotting programs.
    
    currently only writes the real induction vectors.
    
    ======================== ==================================================
    key words/attributes      Description
    ======================== ==================================================
    arrow_direction          [ 1 | -1 ] 1 for Weise convention --> point 
                             toward conductors. *default* is 1 
                             (-1 is not supported yet)
    arrow_head_height        height of arrow head in map units
                             *default* is .002
    arrow_head_width         width of arrow head in map units
                             *default* is .001
    arrow_lw                 width of arrow in map units
                             *default* is .0005
                        
    arrow_size               size of normalized arrow length in map units
                             *default* is .01
                             
    edi_list                 list of edi files, full paths
    mt_obj_list              list of mt.MT objects
                             *default* is None, filled if edi_list is given
    plot_period              list or value of period to convert to shape file
                             *default* is None, which will write a file for
                             every period in the edi files 
    ptol                     tolerance to look for given periods
                             *default* is .05
    pt_dict                  dictionary with keys of plot_period.  Each
                             dictionary key is a structured array containing
                             the important information for the phase tensor.
    projection               projection of coordinates see EPSG for all options
                             *default* is WSG84 
    save_path                path to save files to
                             *default* is current working directory.
    ======================== ==================================================
    
    ======================== ==================================================
    Methods                   Description
    ======================== ==================================================
    _get_plot_period         get a list of all possible frequencies from data
    _get_tip_array           get Tipper information from data and put into 
                             a structured array for easy manipulation
    write_real_shape_files   write real induction arrow shape files
    write_imag_shape_files   write imaginary induction arrow shape files
    ======================== ==================================================
                
    :Example: ::
        >>> edipath = r"/home/edi_files_rotated_to_geographic_north"
        >>> edilist = [os.path.join(edipath, edi) \
                      for edi in os.listdir(edipath)\
                      if edi.find('.edi')>0]
        >>> tps = TipperShapeFile(edilist, save_path=r"/home/gis")
        >>> tps.arrow_head_height = .005
        >>> tps.arrow_lw = .0001
        >>> tps.arrow_size = .05
        >>> tps.write_shape_files()
        
    """
    
    def __init__(self, edi_list=None, **kwargs):
        self.edi_list = edi_list
        self.projection = 'WGS84'
        self.plot_period = None
        self.save_path = os.getcwd()
        self.arrow_size = 1000
        self.arrow_direction = 1
        self.ptol = .05
        self.arrow_head_width = 50
        self.arrow_head_height = 100
        self.arrow_lw = 20
        self.utm_cs = None
        
        self.mt_obj_list = None
        self.tip_dict = None
        
        if self.edi_list is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edi_list]
            
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
        if self.mt_obj_list is not None:
            self._get_plot_period()
            self._get_tip_array()
        
        self._proj_dict = {'WGS84':4326, 'NAD27':4267}
        self._rotation_angle = 0.0
        
    def _set_rotation_angle(self, rotation_angle):
        """
        rotate all mt_objs to rotation angle
        """
        
        self._rotation_angle = float(rotation_angle)
        
        for mt_obj in self.mt_obj_list:
            mt_obj.rotation_angle = float(self._rotation_angle)
            
    def _get_rotation_angle(self):
        return self._rotation_angle
        
    rotation_angle = property(_get_rotation_angle, _set_rotation_angle,
                              doc="rotation angle of Z and Tipper")
    
    def _get_plot_period(self):
        """
        from the list of edi's get a frequency list to invert for.
        
        """

        if self.plot_period is None:
            #get all frequencies from all edi files
            all_freqs = []
            for mt_obj in self.mt_obj_list:
                all_freqs.extend(list(mt_obj.Z.freq))
            
            #sort all frequencies so that they are in descending order,
            #use set to remove repeats and make an array
            self.plot_period = 1./np.array(sorted(list(set(all_freqs)),
                                        reverse=True))
        else:
            if type(self.plot_period) is list:
                pass
            if type(self.plot_period) is int or type(self.plot_period) is float:
                self.plot_period = [self.plot_period]
                
    def _get_tip_array(self):
        """
        get the phase tensor information into a form that is more structured
        to manipulate easier later.
        
        make a dictionary with keys being the plot period values and each
        key has a structured array that contains all the important information
        collected from each station.
        """
        self.tip_dict = {}        
        for plot_per in self.plot_period:
            self.tip_dict[plot_per] = []
            for mt_obj in self.mt_obj_list:
                mt_obj.Tipper._compute_mag_direction()
                try:
                    p_index = [ff for ff, f2 in enumerate(1./mt_obj.Z.freq) 
                               if (f2 > plot_per*(1-self.ptol)) and
                                  (f2 < plot_per*(1+self.ptol))][0]
                    if self.projection is None:
                        east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                        self.utm_cs = osr.SpatialReference()
                        # Set geographic coordinate system to handle lat/lon  
                        self.utm_cs.SetWellKnownGeogCS(self.projection)
                    else:
                        self.utm_cs, utm_point = transform_ll_to_utm(mt_obj.lon, 
                                                                mt_obj.lat,
                                                                self.projection)
                        east, north, elev = utm_point             
                        
                    if mt_obj.Tipper.tipper is not None: 
                        
                        if mt_obj.Tipper.tipper[p_index].all() != 0.0:
                            tp_tuple = (mt_obj.station, 
                                        east,
                                        north,
                                        mt_obj.Tipper.mag_real[p_index],
                                        mt_obj.Tipper.mag_imag[p_index],
                                        mt_obj.Tipper.angle_real[p_index],
                                        mt_obj.Tipper.angle_imag[p_index])           
                            self.tip_dict[plot_per].append(tp_tuple)
                        else:
                            tp_tuple = (mt_obj.station, 
                                        east,
                                        north,
                                        0,
                                        0,
                                        0,
                                        0)
                            self.tip_dict[plot_per].append(tp_tuple)
                            
                except IndexError:
                    pass
                
            self.tip_dict[plot_per] = np.array(self.tip_dict[plot_per],
                                              dtype=[('station', '|S15'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('mag_real', np.float),
                                                     ('mag_imag', np.float),
                                                     ('ang_real', np.float),
                                                     ('ang_imag', np.float)])
            
                
    def write_real_shape_files(self):
        """
        write shape file from given attributes
        """
        
        self._get_tip_array()
            
        for plot_per in self.plot_period:
            #shape file path
            shape_fn = os.path.join(self.save_path, 
                                    'Tip_{0:.5g}s_{1}_real.shp'.format(plot_per,
                                    self.projection))
            
            # remove the shape file if it already exists, has trouble over writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)
            
            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')
            
            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)
                
            #create shape file
            data_source = driver.CreateDataSource(shape_fn)
            
            ##if you read from a raster get the georeference point otherwise create one
            #spatial_ref = osr.SpatialReference()
            #this puts it in the wsg84 reference frame.
            #spatial_ref.ImportFromEPSG(self._proj_dict[self.projection]) 
            
            ##create a layer to put the ellipses onto
            layer = data_source.CreateLayer('TIPPER', self.utm_cs,
                                            ogr.wkbPolygon)
            
            #make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)
            
            field_mag_real = ogr.FieldDefn('mag_real', ogr.OFTReal)
            layer.CreateField(field_mag_real)
            
            field_ang_real = ogr.FieldDefn('ang_real', ogr.OFTReal)
            layer.CreateField(field_ang_real)
            
            for tp_arr in self.tip_dict[plot_per]:
                cos_t = np.cos(-np.deg2rad(tp_arr['ang_real']))
                sin_t = np.sin(-np.deg2rad(tp_arr['ang_real']))
                # calculate the points to make the line
                txr = 0
                tyr = tp_arr['mag_real']*self.arrow_size
                
                # make an arrow by drawing an outline.  have the arrow point 
                # north to start and then rotate later with the rotation 
                # matrix to properly orient it. 
                x0 = 0
                y0 = 0
                
                x1 = x0+self.arrow_lw
                y1 = y0
                
                x2 = x0+self.arrow_lw
                y2 = y0+tyr-self.arrow_head_height
                
                x3 = x0+self.arrow_lw+self.arrow_head_width                 
                y3 = y2
                
                x4 = x0+txr
                y4 = y0+tyr
                
                x7 = x0-self.arrow_lw
                y7 = y0
                
                x6 = x0-self.arrow_lw                  
                y6 = y0+tyr-self.arrow_head_height
                
                x5 = x0-self.arrow_lw-self.arrow_head_width              
                y5 = y6
                
                x = np.array([x0, x1, x2, x3, x4, x5, x6, x7]) 
                y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                
                rot_matrix = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
                
                # rotate the arrow to be properly oriented
                xy = np.array([x, y])
                rot_xy = np.dot(rot_matrix, xy)
                
                # shift the arrow to be centered on the station.
                x = tp_arr['east']+rot_xy[0]
                y = tp_arr['north']+rot_xy[1]
                       
                #1) make a geometry shape line
                arrow = ogr.Geometry(ogr.wkbLinearRing)
                for ii, jj in zip(x, y):
                    arrow.AddPoint(np.round(ii, 6), np.round(jj, 6))

                arrow.CloseRings()
                
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(arrow)
                ##4) this part is confusing but we need to create a 
                ##   feature that has the
                ##   same definition as the layer that we created.
                #    get the layer definition
                feature_def = layer.GetLayerDefn()
                
                #create a new feature
                new_feature = ogr.Feature(feature_def)
                #set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                #create the feature in the layer. 
                layer.CreateFeature(new_feature)
                
                #
                ###5) create a field to color by
                new_feature.SetField("Name", tp_arr['station'])
                new_feature.SetField("mag_real", tp_arr['mag_real'])
                new_feature.SetField("ang_real", tp_arr['ang_real'])
            
                #add the new feature to the layer.
                layer.SetFeature(new_feature)
                
                #apparently need to destroy the feature
                new_feature.Destroy()
                
                
            # Need to be sure that all the new info is saved to 
            data_source.SyncToDisk()
            
            #write a projection file
            #spatial_ref.MorphFromESRI() 
            self.utm_cs.MorphFromESRI() 
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()
            
            data_source.Destroy()
            
            print 'Wrote shape file to {0}'.format(shape_fn)
            
    def write_imag_shape_files(self):
        """
        write shape file from given attributes
        """
        
        self._get_tip_array()
            
        for plot_per in self.plot_period:
            #shape file path
            shape_fn = os.path.join(self.save_path, 
                                    'Tip_{0:.5g}s_{1}_imag.shp'.format(plot_per,
                                    self.projection))
            
            # remove the shape file if it already exists, has trouble over writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)
            
            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')
            
            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)
                
            #create shape file
            data_source = driver.CreateDataSource(shape_fn)
            
            ##if you read from a raster get the georeference point otherwise create one
            #spatial_ref = osr.SpatialReference()
            #this puts it in the wsg84 reference frame.
            #spatial_ref.ImportFromEPSG(self._proj_dict[self.projection])
            
            ##create a layer to put the ellipses onto
            layer = data_source.CreateLayer('TIPPER', self.utm_cs,
                                            ogr.wkbPolygon)
            
            #make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)
            
            field_mag_imag = ogr.FieldDefn('mag_imag', ogr.OFTReal)
            layer.CreateField(field_mag_imag)
            
            field_ang_imag = ogr.FieldDefn('ang_imag', ogr.OFTReal)
            layer.CreateField(field_ang_imag)
            
            for tp_arr in self.tip_dict[plot_per]:
                cos_t = np.cos(-np.deg2rad(tp_arr['ang_imag']))
                sin_t = np.sin(-np.deg2rad(tp_arr['ang_imag']))
                # calculate the points to make the line
                txr = 0
                tyr = tp_arr['mag_imag']*self.arrow_size
                
                # make an arrow by drawing an outline.  have the arrow point 
                # north to start and then rotate later with the rotation 
                # matrix to properly orient it.
                x0 = 0
                y0 = 0
                
                x1 = x0+self.arrow_lw
                y1 = y0
                
                x2 = x0+self.arrow_lw
                y2 = y0+tyr-self.arrow_head_height
                
                x3 = x0+self.arrow_lw+self.arrow_head_width                 
                y3 = y2
                
                x4 = x0+txr
                y4 = y0+tyr
                
                x7 = x0-self.arrow_lw
                y7 = y0
                
                x6 = x0-self.arrow_lw                  
                y6 = y0+tyr-self.arrow_head_height
                
                x5 = x0-self.arrow_lw-self.arrow_head_width              
                y5 = y6
                
                x = np.array([x0, x1, x2, x3, x4, x5, x6, x7]) 
                y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                
                rot_matrix = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
                
                # rotate the arrow to be properly oriented
                xy = np.array([x, y])
                rot_xy = np.dot(rot_matrix, xy)
                
                # shift the arrow to be centered on the station
                x = tp_arr['east']+rot_xy[0]
                y = tp_arr['north']+rot_xy[1]
                       
                #1) make a geometry shape line
                arrow = ogr.Geometry(ogr.wkbLinearRing)
                for ii, jj in zip(x, y):
                    arrow.AddPoint(np.round(ii, 6), np.round(jj, 6))

                arrow.CloseRings()
                
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(arrow)
                ##4) this part is confusing but we need to create a 
                ##   feature that has the
                ##   same definition as the layer that we created.
                #    get the layer definition
                feature_def = layer.GetLayerDefn()
                
                #create a new feature
                new_feature = ogr.Feature(feature_def)
                #set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                #create the feature in the layer. 
                layer.CreateFeature(new_feature)
                
                #
                ###5) create a field to color by
                new_feature.SetField("Name", tp_arr['station'])
                new_feature.SetField("mag_imag", tp_arr['mag_imag'])
                new_feature.SetField("ang_imag", tp_arr['ang_imag'])
            
                #add the new feature to the layer.
                layer.SetFeature(new_feature)
                
                #apparently need to destroy the feature
                new_feature.Destroy()
                
                
            # Need to be sure that all the new info is saved to 
            data_source.SyncToDisk()
            
            #write a projection file
            #spatial_ref.MorphFromESRI() 
            self.utm_cs.MorphFromESRI() 
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()
            
            data_source.Destroy()
            
            print 'Wrote shape file to {0}'.format(shape_fn)
            
    def write_tip_shape_files_modem(self, modem_data_fn, rotation_angle=0.0):
        """
        write tip files from a modem data file.
        
        """

        modem_obj = modem.Data()
        modem_obj.read_data_file(modem_data_fn)
        
        self.plot_period = modem_obj.period_list.copy()
        self.mt_obj_list = [modem_obj.mt_dict[key] 
                            for key in modem_obj.mt_dict.keys()]
                                
        self._set_rotation_angle(rotation_angle)
        
        self.write_imag_shape_files()
        self.write_real_shape_files()
        
    def write_tip_shape_files_modem_residual(self, modem_data_fn, 
                                             modem_resp_fn,
                                             rotation_angle):
        """
        write residual tipper files for modem
        
        """
        modem_data_obj = modem.Data()
        modem_data_obj.read_data_file(modem_data_fn)
        
        modem_resp_obj = modem.Data()
        modem_resp_obj.read_data_file(modem_resp_fn)
        

        
        self.plot_period = modem_data_obj.period_list.copy()
        mt_keys = sorted(modem_data_obj.mt_dict.keys())
        self.mt_obj_list = [modem_data_obj.mt_dict[key] 
                            for key in mt_keys]
                                
        self._set_rotation_angle(rotation_angle)
        
        for mt_obj, key in zip(self.mt_obj_list, mt_keys):
            resp_tipper = modem_resp_obj.mt_dict[key].Tipper.tipper
            mt_obj.Tipper.tipper[:, :, :] -= resp_tipper[:, :, :]
            
        self.write_imag_shape_files()
        self.write_real_shape_files()
        

#==============================================================================
# reproject a layer DOESNT WORK YET
#==============================================================================
def reproject_layer(in_shape_file, out_shape_file=None, out_proj='WGS84'):
    """
    reproject coordinates into a different coordinate system
    
    """    
    
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # input SpatialReference
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(2927)
    
    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(4326)
    
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    # get the input layer
    inDataSet = driver.Open(r'c:\data\spatial\basemap.shp')
    inLayer = inDataSet.GetLayer()
    
    # create the output layer
    outputShapefile = r'c:\data\spatial\basemap_4326.shp'
    if os.path.exists(outputShapefile):
        driver.DeleteDataSource(outputShapefile)
    outDataSet = driver.CreateDataSource(outputShapefile)
    outLayer = outDataSet.CreateLayer("basemap_4326", geom_type=ogr.wkbMultiPolygon)
    
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = inLayer.GetNextFeature()
    
    # close the shapefiles
    inDataSet.Destroy()
    outDataSet.Destroy() 

#==============================================================================
# create a raster from an array     
#==============================================================================
    
def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

#==============================================================================
#     
#==============================================================================
def transform_utm_to_ll(easting, northing, zone, 
                           reference_ellipsoid='WGS84'):
    
    utm_coordinate_system = osr.SpatialReference() 
    # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetWellKnownGeogCS(reference_ellipsoid)
    is_northern = northing > 0    
    utm_coordinate_system.SetUTM(zone, is_northern)
    
    # Clone ONLY the geographic coordinate system 
    ll_coordinate_system = utm_coordinate_system.CloneGeogCS() 
    
    # create transform component
    utm_to_ll_geo_transform = osr.CoordinateTransformation(utm_coordinate_system, 
                                                           ll_coordinate_system)
     # returns lon, lat, altitude
    return utm_to_ll_geo_transform.TransformPoint(easting, northing, 0)
    

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
# test
#==============================================================================
##edipath = r"c:\Users\jrpeacock\Documents\Mendenhall\MonoBasin\EDI_Files\GeographicNorth"
##edilst = [os.path.join(edipath, edi) for edi in os.listdir(edipath)
##          if edi.find('.edi') > 0]
##edilst.remove(os.path.join(edipath, 'mb035.edi'))
##
##pts = PTShapeFile(edilst, save_path=r"c:\Users\jrpeacock")
##pts.projection = 'NAD27'
##pts.ellipse_size = 1200
##pts.write_shape_files()
#
#
##tps = TipperShapeFile(edilst, save_path=r"c:\Users\jrpeacock")
##tps.projection = 'NAD27'
##tps.arrow_lw = 30
##tps.arrow_head_height = 100
##tps.arrow_head_width = 70
##tps.write_real_shape_files()
##tps.write_imag_shape_files()
#    
#mfn = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\Modular_NLCG_110.dat"
#sv_path = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\GIS_Tip_Response"
##sv_path = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\GIS_PT_Response"
##pts = PTShapeFile(save_path=sv_path)
##pts.projection = 'NAD27'
##pts.ellipse_size = 1200
##pts.write_pt_shape_files_modem(mfn)
#
#tps = TipperShapeFile(save_path=sv_path)
#tps.projection = 'NAD27'
#tps.arrow_lw = 30
#tps.arrow_head_height = 100
#tps.arrow_head_width = 70
#tps.write_tip_shape_files_modem(mfn)