# -*- coding: utf-8 -*-
"""
.. module:: mt_xml
   :synopsis: Deal with XML MT files 

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>

.. note:: This module is written to align with the tools written by
          Anna Kelbert <akelbert@usgs.gov>
"""

#==============================================================================
# Imports
#==============================================================================
import os
import datetime
import copy

import numpy as np
import xml.etree.cElementTree as ET
from xml.dom import minidom

import mtpy.core.z as mtz

dt_fmt = '%Y-%m-%d %H:%M:%S'

#==============================================================================
# Generic object to hold information
#==============================================================================
class XML_element(object):
    """
    Basically an ET element.  The key components are
        * 'name'  --> name of the element
        * 'attr'  --> attribute information of the element
        * 'value' --> value of the element
        
    Used the property function here to be sure that these 3 cannot be set 
    through the common k.value = 10, just in case there are similar names in 
    the xml file.  This seemed to be the safest to avoid those cases.
    """
    def __init__(self, name, attr, value, **kwargs):
        self._name = name
        self._attr = attr
        self._value = value
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, value):
        try:
            self._value = str(value)
        except ValueError:
            print('Cannot set {0} to string, set to None'.format(value))
            self._value = None
            
    @property
    def attr(self):
        return self._attr
    
    @attr.setter
    def attr(self, attr):
        if type(attr) is not dict:
            raise ValueError('attr needs to be a dictionary, not {0}'.format(type(attr)))
        else:
            self._attr = attr
            
    @property
    def name(self):
        return self._name
    
#==============================================================================
# General information
#==============================================================================
conditions_of_use = ''.join(['All data and metadata for this survey are ',
                             'available free of charge and may be copied ',
                             'freely, duplicated and further distributed ',
                             'provided this data set is cited as the ',
                             'reference. While the author(s) strive to ',
                             'provide data and metadata of best possible ',
                             'quality, neither the author(s) of this data ',
                             'set, not IRIS make any claims, promises, or ', 
                             'guarantees about the accuracy, completeness, ',
                             'or adequacy of this information, and expressly ',
                             'disclaim liability for errors and omissions in ',
                             'the contents of this file. Guidelines about ',
                             'the quality or limitations of the data and ',
                             'metadata, as obtained from the author(s), are ',
                             'included for informational purposes only.'])
                     
estimates = [XML_element('Estimate', {'type':'real', 'name':'VAR'}, None,
                   **{'Description':XML_element('Description', None, 'Variance'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Error Estimate'),
                      'Tag':XML_element('Tag', None, 'Variance')}),
            
             XML_element('Estimate', {'type':'complex', 'name':'COV'}, None,  
                   **{'Description':XML_element('Description', None, 
                                          'Full covariance between each two TF components'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Error Estimate'),
                      'Tag':XML_element('Tag', None, 'Covariance')}),
             
             XML_element('Estimate', {'type':'complex', 'name':'INVSIGCOV'}, None,  
                   **{'Description':XML_element('Description', None, 
                                          'Inverse Coherent Signal Power Matrix S'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Signal Power Estimate'),
                      'Tag':XML_element('Tag', None, 'inverse_signal_covariance')}),
             
             XML_element('Estimate',{'type':'complex', 'name':'RESIDCOV'}, None,  
                   **{'Description':XML_element('Description',None, 
                                          'Residual Covariance N'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Error Estimate'),
                      'Tag':XML_element('Tag', None, 'Coherence')}),
             
             XML_element('Estimate', {'type':'complex', 'name':'COH'}, None,  
                   **{'Description':XML_element('Description', None, 'Coherence'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Signal Coherence'),
                      'Tag':XML_element('Tag', None, 'Coherence')}),
             
             XML_element('Estimate', {'type':'complex', 'name':'PREDCOH'}, None,  
                   **{'Description':XML_element('Description', None, 'Multiple Coherence'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Signal Coherence'),
                      'Tag':XML_element('Tag', None, 'Multiple_Coherence')}),
             
             XML_element('Estimate', {'type':'complex', 'name':'SIGAMP'}, None,  
                   **{'Description':XML_element('Description', None, 'Signal Amplitude'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Signal Power Estimates'),
                      'Tag':XML_element('Tag', None, 'Signal_Amplitude')}),
             
             XML_element('Estimate', {'type':'complex', 'name':'SIGNOISE'}, None,  
                   **{'Description':XML_element('Description', None, 'Signal Noise'),
                      'ExternalUrl':XML_element('ExternalUrl', None, None),
                      'Intention':XML_element('Intention', None, 'Error Estimates'),
                      'Tag':XML_element('Tag', None, 'Signal_Noise')})]
                       
data_types = [XML_element('DataType',
                    {'units':'[mV/km]/[nT]', 
                     'name':'Z', 
                     'input':'H',
                     'output':'E'},
                     None, 
                     **{'Description':XML_element('Description', None, 'MT impedance'),
                       'ExternalUrl':XML_element('ExternalUrl', None, None),
                       'Intention':XML_element('Intention', None, 'primary data type'),
                       'Tag':XML_element('Tag', None, 'impedance')}),
              
              XML_element('DataType',
                    {'units':'[]', 
                     'name':'T', 
                     'input':'H',
                     'output':'H'},
                     None, 
                     **{'Description':XML_element('Description', None, 
                                            'Tipper-Vertical Field Transfer Function'),
                       'ExternalUrl':XML_element('ExternalUrl', None, None),
                       'Intention':XML_element('Intention', None, 'primary data type'),
                       'Tag':XML_element('Tag', None, 'tipper')})]
#==============================================================================
# Useful Functions
#==============================================================================                
class XML_Config(object):
    """
    Class to deal with configuration files for xml.
    
    Includes all the important information for the station and how data was 
    processed.  
    
    Key Information includes:
        
    ======================== ==================================================
    Name                     Purpose
    ======================== ==================================================
    ProductID                Station name
    ExternalUrl              External URL to link to data
    Notes                    Any important information on station, 
                             data collection.
    TimeSeriesArchived       Information on Archiving time series including 
                             URL.   
    Image                    A location to an image of the station or
                             the MT response.
    ======================== ==================================================    
    * ProductID --> station name
        * ExternalUrl --> external url to link to data
        * Notes --> any 
    
    
    """    
    
    def __init__(self, **kwargs):
        self.cfg_fn = None
        
        # Initialize the default attributes and values
        self.Description = XML_element('Description', 
                                 None, 
                                 'Magnetotelluric Transfer Functions')
                                 
        self.ProductId = XML_element('ProductID', None, None)
        
        self.Project = XML_element('Project', None, None)
        
        self.Survey = XML_element('Survey', None, None)
        
        self.Country = XML_element('Country', None, None)
       
        self.SubType = XML_element('SubType', None, 'MT_FT')
       
        self.Notes = XML_element('Notes', None, None)
        
        self.Tags = XML_element('Tags', None, 'impedance, tipper')
        
        
        self.Image = XML_element('Image', None, None, 
                           **{'PrimaryData':XML_element('PrimaryData', None, None),
                              'Filename':XML_element('Filename', None, None)})
                              
        
        self.Original = XML_element('Original', None, None,
                              **{'Attachment':XML_element('Attachment', None, None),
                                 'Filename':XML_element('Filename', None, None)})
        
       
        self.TimeSeriesArchived = XML_element('TimeSeriesArchived', None, None, 
                                        **{'Value':XML_element('Value', None, 0), 
                                           'URL':XML_element('URL', None, None)})
        
        self.ExternalUrl = XML_element('ExternalUrl', None, None, 
                                 **{'Description':XML_element('Description', None, None),
                                    'Url':XML_element('Url', None, None)})
        
        self.PrimaryData = XML_element('PrimaryData', None, None, 
                                 **{'Filename':XML_element('Filename', None, None)})
                                    
        
        self.Attachment = XML_element('Attachment', None, None, 
                                **{'Filename':XML_element('Filename', None, None),
                                   'Description':XML_element('Description', None, 
                                                       'Original file use to produce XML')})
                                   
        
        self.Provenance = XML_element('Provenance', None, None,
                                **{'CreationTime':XML_element('CreationTime', None, 
                                                     datetime.datetime.strftime(
                                                     datetime.datetime.utcnow(), 
                                                     dt_fmt)),
                                    'CreatingApplication':XML_element('CreatingApplication', 
                                                                None, 
                                                                'MTpy.core.mtxml'),
                                    'Submitter':XML_element('Submitter', None, None,
                                                      **{'Name':XML_element('Name', None, None),
                                                         'Email':XML_element('Email', None, None),
                                                         'Org':XML_element('Org', None, None),
                                                         'OrgURL':XML_element('OrgURL', None, None)}),
                                    'Creator':XML_element('Creator', None, None,
                                                    **{'Name':XML_element('Name', None, None),
                                                       'Email':XML_element('Email', None, None),
                                                       'Org':XML_element('Org', None, None),
                                                       'OrgURL':XML_element('OrgURL', None, None)})})
                                                       
        self.Copyright = XML_element('Copyright', None, None,
                               **{'Citation':XML_element('Citation', None, None,
                                                   **{'Title':XML_element('Title', None, None),
                                                      'Authors':XML_element('Authors', None, None),
                                                      'Year':XML_element('Year', None, None),
                                                      'Journal':XML_element('Journal', None, None),
                                                      'Volume':XML_element('Volume', None, None),
                                                      'DOI':XML_element('DOI', None, None)}),
                                  'ReleaseStatus':XML_element('ReleaseStatus', None, 'Closed'),
                                  'ConditionsOfUse':XML_element('ConditionsOfUse', None, conditions_of_use),
                                  'AdditionalInfo':XML_element('AdditionalInfo', None, None)})

        self.Site = XML_element('Site', None, None,
                          **{'Project':XML_element('Project', None, None),
                             'Survey':XML_element('Survey', None, None),
                             'YearCollected':XML_element('YearCollected', None, None),
                             'Id':XML_element('Id', None, None),
                             'Location':XML_element('Location', None, None,
                                              **{'Latitude':XML_element('Latitude', None, None),
                                                 'Longitude':XML_element('Longitude', None, None),
                                                 'Elevation':XML_element('Elevation', {'units':'meters'}, None),
                                                 'Declination':XML_element('Declination', {'epoch':'1995'}, None)}),
                            'Orientation':XML_element('Orientation', {'angle_to_geographic_north':'0.0'}, None), 
                            'AcquiredBy':XML_element('AcquiredBy', None, None),
                             'Start':XML_element('Start', None, None),
                             'End':XML_element('End', None, None),
                             'RunList':XML_element('RunList', None, None)})
                             
        self.FieldNotes = XML_element('FieldNotes', None, None,
                                **{'Instrument':XML_element('Instrument', None, None,
                                                      **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None),
                                                         'Settings':XML_element('Settings', None, None)}),
                                    'Dipole':XML_element('Dipole', {'name':'EX'}, None,
                                                      **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None),
                                                         'Length':XML_element('Length', {'units':'meters'}, None),
                                                         'Azimuth':XML_element('Azimuth', {'units':'degrees'}, None),
                                                         'Channel':XML_element('Channel', None, None)}),
                                    'Dipole_00':XML_element('Dipole', {'name':'EY'}, None,
                                                      **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None),
                                                         'Length':XML_element('Length', {'units':'meters'}, None),
                                                         'Azimuth':XML_element('Azimuth', {'units':'degrees'}, None),
                                                         'Channel':XML_element('Channel', None, None)}),
                                    'Magnetometer':XML_element('Magnetometer', {'name':'HX'}, None,
                                                        **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None),
                                                         'Azimuth':XML_element('Azimuth', {'units':'degrees'}, None),
                                                         'Channel':XML_element('Channel', None, None)}),
                                    'Magnetometer_00':XML_element('Magnetometer', {'name':'HY'}, None,
                                                        **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None),
                                                         'Azimuth':XML_element('Azimuth', {'units':'degrees'}, None),
                                                         'Channel':XML_element('Channel', None, None)}),
                                    'Magnetometer_01':XML_element('Magnetometer', {'name':'HZ'}, None,
                                                        **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None),
                                                         'Azimuth':XML_element('Azimuth', {'units':'degrees'}, None),
                                                         'Channel':XML_element('Channel', None, None)}),
                                    'DataQualityNotes':XML_element('DataQualityNotes', None, None,
                                                             **{'Rating':XML_element('Rating', None, None),
                                                                 'GoodFromPeriod':XML_element('GoodFromPeriod', None, None),
                                                                 'GoodToPeriod':XML_element('GoodToPeriod', None, None),
                                                                 'Comments':XML_element('Comments', None, None)}),
                                    'DataQualityWarnings':XML_element('DataQualityWarnings', None, None,
                                                                **{'Flag':XML_element('Flag', None, 0),
                                                                    'Comments':XML_element('Comments', None, None)})})
                                  
        self.ProcessingInfo = XML_element('ProcessingInfo', None, None,
                                    **{'ProcessedBy':XML_element('ProcessedBy', None, None),
                                       'ProcessingSoftware':XML_element('ProcessingSoftware', None, None,
                                                                  **{'Name':XML_element('Name', None, None),
                                                                     'LastMod':XML_element('LastMod', None, None),
                                                                     'Version':XML_element('Version', None, None),
                                                                     'Author':XML_element('Author', None, None)}),
                                        'SignConvention':XML_element('SignConvention', None, r'exp(+i\omega t)'),
                                        'RemoteRef':XML_element('RemoteRef', {'type':'Robust Remote Processing'}, None),
                                        'RemoteInfo':XML_element('RemoteInfo', None, None, 
                                                           **{'Project':XML_element('Project', None, None),
                                                              'Survey':XML_element('Survey', None, None),
                                                              'ID':XML_element('ID', None, None),
                                                              'AcquiredBy':XML_element('AcquiredBy', None, None),
                                                              'Name':XML_element('Name', None, None),
                                                              'YearCollected':XML_element('YearCollected', None, None),
                                                              'Location':XML_element('Location', {'datum':'WGS84'}, None,
                                                                               **{'Latitude':XML_element('Latitude', None, None),
                                                                                  'Longitude':XML_element('Longitude', None, None),
                                                                                  'Elevation':XML_element('Elevation', {'units':'meters'}, None)})
                                                               })})
        self.SiteLayout = XML_element('SiteLayout', None, None,
                                      **{'InputChannels': XML_element('InputChannels', 
                                                                      {'ref':'site', 'units':'m'}, 
                                                                      None,
                                              **{'Magnetic_hx':XML_element('Magnetic',
                                                                        {'name':"Hx",
                                                                         'orientation': '0.0', 
                                                                         'x':'0.0',
                                                                         'y':'0.0',
                                                                         'z':'0.0'},
                                                                         None),
                                                 'Magnetic_hy':XML_element('Magnetic',
                                                                        {'name':"Hy",
                                                                         'orientation': '0.0', 
                                                                         'x':'0.0',
                                                                         'y':'0.0',
                                                                         'z':'0.0'},
                                                                         None)}),
                                        'OutputChannels': XML_element('OutputChannels', 
                                                                      {'ref':'site', 'units':'m'}, 
                                                                      None,
                                              **{'Magnetic_hz':XML_element('Magnetic',
                                                                        {'name':"Hz",
                                                                         'orientation': '0.0', 
                                                                         'x':'0.0',
                                                                         'y':'0.0',
                                                                         'z':'0.0'},
                                                                         None),
                                                 'Electric_ex':XML_element('Electric',
                                                                        {'name':"Ex",
                                                                         'orientation': '0.0', 
                                                                         'x':'0.0',
                                                                         'y':'0.0',
                                                                         'z':'0.0',
                                                                         'x2':'0.0',
                                                                         'y2':'0.0',
                                                                         'z2':'0.0'},
                                                                         None),
                                                'Electric_ey':XML_element('Electric',
                                                                        {'name':"Ey",
                                                                         'orientation': '0.0', 
                                                                         'x':'0.0',
                                                                         'y':'0.0',
                                                                         'z':'0.0',
                                                                         'x2':'0.0',
                                                                         'y2':'0.0',
                                                                         'z2':'0.0'},
                                                                         None)})})
                        
#        self.InputChannels = XML_element('InputChannels', {'ref':'site', 'units':'m'}, None)
#        self.OutputChannels = XML_element('OutputChannels', {'ref':'site', 'units':'m'}, None)
        self.Data = XML_element('Data', {'count':0}, None)
        self.PeriodRange = XML_element('PeriodRange', None, None)
                                        
        self.Datum = XML_element('Datum', None, 'WGS84')
        self.Declination = XML_element('Declination', None, None)
                 
        self.StatisticalEstimates = XML_element('StatisticalEstimates', None, None)                           
        for ii, estimate in enumerate(estimates):
            setattr(self.StatisticalEstimates, 
                    'Estimate_{0:02}'.format(ii),
                    estimate)
        
        self.DataTypes = XML_element('DataTypes', None, None)
        for ii, d_type in enumerate(data_types):
            setattr(self.DataTypes, 'DataType_{0:02}'.format(ii), d_type)
                                   
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
    def read_cfg_file(self, cfg_fn=None):
        """
        Read in a cfg file making all key = value pairs attribures of 
        XML_Config.  Being sure all new attributes are XML_element objects.
        
        The assumed structure of the xml.cfg file is similar to:
            ``# XML Configuration File MTpy
 
            Attachement.Description = Original file use to produce XML
            Attachment.Filename = None
             
            Copyright.Citation.Authors = None
            Copyright.Citation.DOI = None
            Copyright.Citation.Journal = None
            Copyright.Citation.Title = None
            Copyright.Citation.Volume = None
            Copyright.Citation.Year = None
            
            PeriodRange(max=0)(min=0) = None``
            
        where the heirarchy of information is separated by a . and if the
        information has attribures they are in the name with (key=value) 
        syntax.
        
        Arguments
        -------------
            **cfg_fn** : string
                         full path to cfg file to read in
                         
        Example
        ---------
        :Read in xml.cfg file: ::
        
            >>> import mtpy.core.mtxml as mtxml
            >>> cfg_obj = mtxml.XML_Config()
            >>> cfg_obj.read_cfg_file(r"/home/MT/xml.cfg")
            
        """        
        
        if cfg_fn is not None:        
            self.cfg_fn = cfg_fn
            
        if not os.path.isfile(self.cfg_fn):
            raise NameError('Could not find {0}'.format(self.cfg_fn))
            
        with open(self.cfg_fn, 'r') as fid:
            lines = fid.readlines()
  
        for line in lines:
            # skip comments
            if line[0] == '#':
                pass
            
            # skip blank lines
            elif len(line.strip()) < 3:
                pass
            # else assume the line is metadata separated by = 
            else:
                self._read_cfg_line(line)
                    

    def _read_cfg_line(self, line):
        """
        read a configuration file line to make the appropriate attribute
        have the correct values and attributes.
        
        porbably should think of a better name for XML_element objects that are
        attributes of self.
        """
        
        # split the line by the last =
        line_list = self._split_cfg_line(line)
        # split the keys by . to get the different attributes
        key_list = line_list[0].strip().split('.')
        value = line_list[1].strip()
        if value in ['none', 'None']:
            value = None
        
        # loop over the keys to set them appropriately 
        for ii, key in enumerate(key_list):
            # get the name of the key and any attributes it might have
            name, attr = self._read_cfg_key(key)
            
            # if its the first key, see if its been made an attribute yet
            if ii == 0: 
                if not hasattr(self, name):
                    setattr(self, name, XML_element(name, None, None))
                # for looping purposes we need to get the current XML_element object
                cfg_attr = getattr(self, name)
                # be sure to set any attributes, need to do this here because
                # the test for hasattr will only make a new one if there
                # isn't one already, but since most things in the cfg file
                # are already attributes of self, they already exist.
                cfg_attr._attr = attr
            else:
                if not hasattr(cfg_attr, name):
                    setattr(cfg_attr, name, XML_element(name, None, None))
                cfg_attr = getattr(cfg_attr, name) 
                cfg_attr._attr = attr
        
        # set the value of the current XML_element object
        cfg_attr._value = value

    def _split_cfg_line(self, line):
        """
        split a cfg line by the last =, otherwise you get strings that are
        split in the wrong place.  For instance k.l(a=b) = None will be split
        at ['k.l(a', 'b)', None] but we want [k.l(a=b), None]
        """
        
        equal_find = -(line[::-1].find('='))
        line_list = [line[0:equal_find-1],
                     line[equal_find+1:]]
                     
        return line_list
        
            
    def _read_cfg_key(self, key):
        """
        read a key from a cfg file and check to see if has any attributes
        in the form of:
            parent.name(attribute=0)(attribute=2) = value
        """
        attr = {}
        if '(' and ')' in key:
            key_list = key.split('(')
            key = key_list[0]
            for key_attr in key_list[1:]:
                k_list = key_attr.replace(')', '').split('=')
                attr[k_list[0].strip()] = k_list[1].strip()
                
        return key, attr
        
    def write_cfg_file(self, cfg_fn=None):
        """
        Write out configuration file in the style of:
        parent.attribute = value
        
        Arguments
        --------------
            **cfg_fn** : string
                         full path to write the configuration file
                         
        """
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn
            
        # a list of lines that will be joined to write the file
        line_list = ['# XML Configuration File MTpy']
        
        # loop over attribute names 
        for attr_00_name in sorted(self.__dict__.keys()):
            # skip the data attribute cause we don't need that now
            if attr_00_name in ['Data', 'DataTypes', 'StatisticalEstimates']:
                continue
            
            # get the given attribute
            attr_00 = getattr(self, attr_00_name)
            
            # make sure it is of XML_element instance
            if isinstance(attr_00, XML_element):
                # be sure to add a new line for each parent attribute
                line_list.append(' ')
                # get the attributes associated with the parent
                attr_00_keys = self._get_attr_keys(attr_00)
            
                # if there are no attributes write the line
                if len(attr_00_keys) == 0:
                    line_list.append(self._write_cfg_line(attr_00))
                
                # otherwise loop through each attribute checking if there
                # are more attributes or no
                else:
                    for attr_01_name in attr_00_keys:
                        attr_01 = getattr(attr_00, attr_01_name)
                        attr_01_keys = self._get_attr_keys(attr_01)
                        
                        if len(attr_01_keys) == 0:
                            line_list.append(self._write_cfg_line(attr_01, attr_00))
                        else:
                            for attr_02_name in attr_01_keys:
                                attr_02 = getattr(attr_01, attr_02_name)
                                attr_02_keys = self._get_attr_keys(attr_02)
                                    
                                if len(attr_02_keys) == 0:
                                    line_list.append(self._write_cfg_line(attr_02, 
                                                                          [attr_00, attr_01]))
                                else:
                                    for attr_03_name in attr_02_keys:
                                        attr_03 = getattr(attr_02,
                                                          attr_03_name)
                                        attr_03_keys = self._get_attr_keys(attr_03)
                                        
                                        if len(attr_03_keys) == 0:
                                            line_list.append(self._write_cfg_line(attr_03, 
                                                                                  [attr_00, attr_01, attr_02]))
                                          
                                        else:
                                            for attr_04_name in attr_03_keys:
                                                attr_04 = getattr(attr_03, attr_04_name)
                                                line_list.append(self._write_cfg_line(attr_04,
                                                                                      [attr_00, attr_01, attr_02, attr_03]))
            else:
                print('Not including: {0}'.format(attr_00_name))                     

        # write the file
        with open(self.cfg_fn, 'w') as fid:
            fid.write('\n'.join(line_list))
            
        # show the user something happened
        print('-'*50)
        print('    Wrote xml configuration file to {0}'.format(self.cfg_fn))
        print('-'*50)
        
        

    def _write_cfg_line(self, XML_element_obj, parent=None):
        """
        write a configuration file line in the format of:
        parent.attribute = value
        """        
        
        if parent is None:
            parent_str = ''
            
        elif type(parent) is list:
            parent_str = '.'.join([p._name for p in parent]+[''])
            
        elif isinstance(parent, XML_element):
            parent_str = '{0}.'.format(parent._name)
        
        if XML_element_obj._attr is not None:
            attr_str = ''.join(['({0}={1})'.format(a_key, XML_element_obj._attr[a_key]) 
                                 for a_key in list(XML_element_obj._attr.keys())])
        else:
            attr_str = ''
        return '{0}{1}{2} = {3}'.format(parent_str,
                                        XML_element_obj._name, 
                                        attr_str, 
                                        XML_element_obj._value)
                                     
    def _get_attr_keys(self, attribute):
        return [a_key for a_key in sorted(attribute.__dict__.keys()) 
                if a_key not in ['_name', '_value', '_attr']]
#==============================================================================
#  EDI to XML
#==============================================================================
class MT_XML(XML_Config):
    """
    Class to read and write MT information from XML format.  This tries to 
    follow the format put forward by Anna Kelbert for archiving MT response 
    data.
    
    A configuration file can be read in that might make it easier to write
    multiple files for the same survey.  
    
    .. seealso:: mtpy.core.mt_xml.XML_Config
   
    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    Z               object of type mtpy.core.z.Z 
    Tipper          object of type mtpy.core.z.Tipper
    =============== ===========================================================
   
    .. note:: All other attributes are of the same name and of type XML_element,
              where attributes are name, value and attr.  Attr contains any 
              tag information.  This is left this way so that mtpy.core.mt.MT
              can read in the information.  **Use mtpy.core.mt.MT for 
              conversion between data formats.**
   
    =============== ===========================================================
    Methods         Description
    =============== =========================================================== 
    read_cfg_file   Read a configuration file in the format of XML_Config
    read_xml_file   Read an xml file
    write_xml_file  Write an xml file
    =============== ===========================================================
    
    :Example: ::
        >>> import mtpy.core.mt_xml as mtxml
        >>> x = mtxml.read_xml_file(r"/home/mt_data/mt01.xml")
        >>> x.read_cfg_file(r"/home/mt_data/survey_xml.cfg")
        >>> x.write_xml_file(r"/home/mt_data/xml/mt01.xml")
        
    """
    
    def __init__(self, **kwargs):
        
        XML_Config.__init__(self, **kwargs )
        self.edi_fn = None
        self.xml_fn = None
        self.cfg_fn = None
        
        self.parent_element = None

        self._Z = mtz.Z()
        self._Tipper = mtz.Tipper()
        
        self._order_list = ['Description',
                            'ProductId',
                            'SubType',
                            'Notes',
                            'Tags',
                            'ExternalUrl',
                            'PrimaryData',
                            #'TimeSeriesArchived',
                            #'Image',
                            #'Original',
                            'Attachment',
                            'Provenance',
                            'Copyright',
                            'Site',
                            'FieldNotes',
                            'ProcessingInfo',
                            'StatisticalEstimates',
                            'DataTypes',
                            'SiteLayout',
                            'Data',
                            'PeriodRange']
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
        
    def _get_name(self, name):
        if name.find('(') > 0:
            l = name.split('(')
            name = l[0]
            meta_dict = {}
            for ll in l[1:]:
                ll = ll.split('=')
                meta_dict[ll[0]] = ll[1].replace(')', '')
        else:
            meta_dict = {}
        
        return name, meta_dict
        
    def _format_data(self):
        """
        format the Z and tipper data apporpriately
        """
        # --> useful variables
        comp_dict_z = {(0, 0):('Zxx', 'Hx', 'Ex'), 
                       (0, 1):('Zxy', 'Hy', 'Ex'), 
                       (1, 0):('Zyx', 'Hx', 'Ey'),
                       (1, 1):('Zyy', 'Hy', 'Ey')}
                       
        header_dict = {}
        header_dict['Z'] = XML_element('Z',{'units':'[mV/km]/[nT]', 
                                      'type':'complex', 
                                      'size':'2 2'}, None)
        header_dict['Z.VAR'] = XML_element('Z.VAR', {'type':'real', 'size':'2 2'},
                                     None)
        header_dict['Z.INVSIGCOV'] = XML_element('Z.INVSIGCOV',{'type':'complex', 
                                                                'size':'2 2'}, 
                                                 None)
        header_dict['Z.RESIDCOV'] = XML_element('Z.RESIDCOV',{'type':'complex', 
                                                                'size':'2 2'}, 
                                                 None)
        header_dict['T'] = XML_element('T', {'units':'[]', 
                                       'type':'complex',
                                       'size':'1 2'}, None)
        header_dict['T.VAR'] = XML_element('T.VAR', {'type':'real', 
                                               'size':'1 2'}, None)
        header_dict['T.INVSIGCOV'] = XML_element('T.INVSIGCOV',{'type':'complex', 
                                                                'size':'2 2'}, 
                                                 None)
        header_dict['T.RESIDCOV'] = XML_element('T.RESIDCOV',{'type':'complex', 
                                                                'size':'1 1'}, 
                                                 None)
        nf = self.Z.freq.size 
        
        # determine whether to write tipper data or not
        if self.Tipper.tipper is not None:
            nz_tipper = np.any(self.Tipper.tipper) == 0
            if nz_tipper == True:
                write_tipper = False
            else:
                write_tipper = True
        else: 
            write_tipper = False
            
        # set the estimates to write out
        if write_tipper == True:
            estimates = ['Z', 'Z.VAR', 'Z.INVSIGCOV', 'Z.RESIDCOV',
                         'T', 'T.VAR', 'T.INVSIGCOV', 'T.RESIDCOV']
        else:
            estimates = ['Z', 'Z.VAR', 'Z.INVSIGCOV', 'Z.RESIDCOV']

        # make the data element
        
        self.Data = XML_element('Data', {'count':str(nf)}, None)  
        
        # loop through each period and add appropriate information
        for f_index, freq in enumerate(self.Z.freq):
            f_name = 'Period_{0:02}'.format(f_index)
            # set attribute period name with the index value
            # we are setting _name to have the necessary information so
            # we can name the attribute whatever we want.
            setattr(self.Data, f_name,
                    XML_element('Period', {'value':'{0:.6g}'.format(1./freq),
                                     'units':'seconds'}, None))
            d_attr = getattr(self.Data, f_name)
            # Get information from data
            for estimate in estimates:
                attr_name = estimate.replace('.', '_').replace('VAR', 'err').lower()
                estimate_name = estimate.replace('.', '')
                # need to make sure the attribute value is a copy otherwise it
                # will continue to rewrite itself.
                setattr(d_attr, estimate_name,
                        copy.deepcopy(header_dict[estimate]))
                c_attr = getattr(d_attr, estimate_name)
                if 'z' in attr_name:
                    count = 0
                    try:
                        z_arr = getattr(self.Z, attr_name)
                    except AttributeError:
                        #print 'No {0} information'.format(attr_name)
                        continue
                    for e_index in range(2):
                        for h_index in range(2):
                            c = comp_dict_z[(e_index, h_index)]
                            c_dict = {'name':c[0], 'input':c[1], 'output':c[2]}
                            z_value = z_arr[f_index, e_index, h_index]
                            if attr_name == 'z_err':
                                c_value = '{0:<+.6e}'.format(z_value)
                            else:
                                c_value = '{0:<+.6e} {1:<+.6e}'.format(z_value.real, 
                                                                       z_value.imag)
                            
                            setattr(c_attr, 
                                    'value_{0:02}'.format(count),
                                    XML_element('value', c_dict, c_value)) 

                            count += 1
                            
                if 't' in attr_name and write_tipper == True:
                    attr_name = attr_name.replace('t', 'tipper')
                    count = 0
                    if attr_name.lower() in ['tipper', 'tipper_err']:
                        tx = 1
                        ty = 2
                        comp_dict_t = {(0, 0):('Tx', 'Hx', 'Hz'),
                                       (0, 1):('Ty', 'Hy', 'Hz')}
                    elif attr_name.lower() == 'tipper_invsigcov':
                        tx = 2
                        ty = 2
                        comp_dict_t = {(0, 0):('', 'Hx', 'Hx'),
                                       (0, 1):('', 'Hx', 'Hy'),
                                       (1, 0):('', 'Hy', 'Hx'),
                                       (1, 1):('', 'Hy', 'Hy')}
                    elif attr_name.lower() == 'tipper_residcov':
                        tx = 1
                        ty = 1
                        comp_dict_t = {(0, 0):('', 'Hz', 'Hz')}
                        
                    for e_index in range(tx):
                        for h_index in range(ty):
                            c = comp_dict_t[(e_index, h_index)]
                            c_dict = {'name':c[0], 'input':c[1], 'output':c[2]}
                            try:
                                t_arr = getattr(self.Tipper, attr_name)
                            except AttributeError:
                                continue
                            t_value = t_arr[f_index, e_index, h_index]
                            if attr_name ==  'tipper_err': 
                                c_value = '{0:<+.6e}'.format(t_value)
                            else:
                                c_value = '{0:<+.6e} {1:<+.6e}'.format(t_value.real, 
                                                                       t_value.imag)
                             
                            setattr(c_attr, 
                                    'value_{0:02}'.format(count),
                                    XML_element('value', c_dict, c_value)) 
                            count += 1
                            
        self.PeriodRange._attr = {'min':'{0:.6e}'.format(1./self.Z.freq.min()),
                                  'max':'{0:.6e}'.format(1./self.Z.freq.max())}
        
    def _write_element(self, parent_et, XML_element_obj):
        """
        make a new element 
        """
        if XML_element_obj._attr is None:
            XML_element_obj._attr = {}
        else:
            for key in list(XML_element_obj._attr.keys()):
                XML_element_obj._attr[key] = str(XML_element_obj._attr[key])
#        if XML_element_obj._name is None:
#            XML_element_obj._name = 'None'
#        if XML_element_obj._value is None:
#            XML_element_obj.value = 'None'
            
        new_element = ET.SubElement(parent_et, 
                                    XML_element_obj._name,
                                    XML_element_obj._attr)
        new_element.text = XML_element_obj._value
        #new_element.tail = '\n'
        return new_element

    def write_xml_file(self, xml_fn, cfg_fn=None):
        """
        write xml from edi
        """
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn
            self.read_cfg_file()
        
        self.xml_fn = xml_fn

        # get data inot xml format
        self._format_data()

        # make the top of the tree element        
        emtf = ET.Element('EM_TF')
        
        # loop over the important information sections
        for element in self._order_list:
            # get the information for the given element
            d_00_obj = getattr(self, element)
            
            element_00 = self._write_element(emtf, d_00_obj)

            key_00_list = self._get_attr_keys(d_00_obj)
            if len(key_00_list) !=0:
                for key_00 in key_00_list:
                    d_01_obj = getattr(d_00_obj, key_00)
                    element_01 = self._write_element(element_00, d_01_obj)
                    
                    key_01_list = self._get_attr_keys(d_01_obj)
                    if len(key_01_list) !=0:
                        for key_01 in key_01_list:
                            d_02_obj = getattr(d_01_obj, key_01)
                            element_02 = self._write_element(element_01, d_02_obj)
                            
                            key_02_list = self._get_attr_keys(d_02_obj)
                            if len(key_02_list) !=0:
                                for key_02 in key_02_list:
                                    d_03_obj = getattr(d_02_obj, key_02)
                                    element_03 = self._write_element(element_02, d_03_obj)
                                    
                                    key_03_list = self._get_attr_keys(d_03_obj)
                                    if len(key_03_list) !=0:
                                        for key_03 in key_03_list:
                                            d_04_obj = getattr(d_03_obj, key_03)
                                            element_04 = self._write_element(element_03, d_04_obj)
                                    
                                            key_04_list = self._get_attr_keys(d_04_obj)
                                            if len(key_04_list) !=0:
                                                for key_04 in key_04_list:
                                                    d_05_obj = getattr(d_04_obj, key_04)
                                                    element_05 = self._write_element(element_04, d_05_obj)         
  
        #--> write xml file
        xmlstr = minidom.parseString(ET.tostring(emtf)).toprettyxml(indent="   ")
        with open(self.xml_fn, 'w') as fid:
            fid.write(xmlstr)
        
        print('-'*72)
        print('    Wrote xml file to: {0}'.format(self.xml_fn))
        print('-'*72)
        
        return self.xml_fn
       
    def read_xml_file(self, xml_fn):
        """
        read in an xml file and set attributes appropriately.
        
        Arguments
        --------------
            **xml_fn** : string
                         full path of xml file to read in
        """
        
        self.xml_fn = xml_fn
        
        et_xml = ET.parse(xml_fn)

        root = et_xml.getroot()
        for element_00 in root.getchildren():
            setattr(self, element_00.tag, self._read_element(element_00))
            
        try:
            setattr(self.FieldNotes,
                    'DataQualityNotes',
                    self.Site.DataQualityNotes)
            delattr(self.Site, 'DataQualityNotes')
        except AttributeError:
            pass
        
        try:
            setattr(self.FieldNotes,
                    'DataQualityWarnings',
                    self.Site.DataQualityWarnings)
            delattr(self.Site, 'DataQualityWarnings')
        except AttributeError:
            pass
            
        # set Z and Tipper
        nf = int(self.Data.attr['count'])
        period = np.zeros(nf, dtype=np.float)
        z = np.zeros((nf, 2, 2), dtype=np.complex)
        z_err = np.zeros((nf, 2, 2), dtype=np.float)
        
        self._Z = mtz.Z(z_array=z, z_err_array=z_err, freq=period)
        self._Z.z_invsigcov =  np.zeros((nf, 2, 2), dtype=np.complex)
        self._Z.z_residcov = np.zeros((nf, 2, 2), dtype=np.complex)
        
        input_dict = {'hx': 0, 'hy': 1, 'ex': 0, 'ey': 1}
        output_dict = input_dict
        
        p_count = 0
        for per_attr in dir(self.Data):
            if 'period' in per_attr.lower():
                p_obj = getattr(self.Data, per_attr)
                p1 = 1./float(p_obj.attr['value'])
                self._Z.freq[p_count] = 1./p1
                if hasattr(p_obj, 'Z'):
                    z_block = getattr(p_obj, 'Z')
                    for z_attr in dir(z_block):
                        if 'value_' in z_attr:
                            z_comp = getattr(z_block, z_attr)
                            ii = output_dict[z_comp.attr['output'].lower()]
                            jj = input_dict[z_comp.attr['input'].lower()]
                            z_value = [float(zz) for zz in z_comp.value.strip().split()]
                            z_value = complex(z_value[0], z_value[1])
                            self._Z.z[p_count, ii, jj] = z_value
                
                if hasattr(p_obj, 'ZVAR'):
                    z_block = getattr(p_obj, 'ZVAR')
                    for z_attr in dir(z_block):
                        if 'value_' in z_attr:
                            z_comp = getattr(z_block, z_attr)
                            ii = output_dict[z_comp.attr['output'].lower()]
                            jj = input_dict[z_comp.attr['input'].lower()]
                            z_value = float(z_comp.value)
                            self._Z.z_err[p_count, ii, jj] = z_value
                
                if hasattr(p_obj, 'ZINVSIGCOV'):
                    z_block = getattr(p_obj, 'ZINVSIGCOV')
                    for z_attr in dir(z_block):
                        if 'value_' in z_attr:
                            z_comp = getattr(z_block, z_attr)
                            ii = output_dict[z_comp.attr['output'].lower()]
                            jj = input_dict[z_comp.attr['input'].lower()]
                            z_value = [float(zz) for zz in z_comp.value.strip().split()]
                            z_value = complex(z_value[0], z_value[1])
                            self._Z.z_invsigcov[p_count, ii, jj] = z_value
                            
                if hasattr(p_obj, 'ZRESIDCOV'):
                    z_block = getattr(p_obj, 'ZRESIDCOV')
                    for z_attr in dir(z_block):
                        if 'value_' in z_attr:
                            z_comp = getattr(z_block, z_attr)
                            ii = output_dict[z_comp.attr['output'].lower()]
                            jj = input_dict[z_comp.attr['input'].lower()]
                            z_value = [float(zz) for zz in z_comp.value.strip().split()]
                            z_value = complex(z_value[0], z_value[1])
                            self._Z.z_residcov[p_count, ii, jj] = z_value
                p_count += 1
                
        # Fill Tipper
        t = np.zeros((nf, 1, 2), dtype=np.complex)
        t_err = np.zeros((nf, 1, 2), dtype=np.float) 
        
        self._Tipper = mtz.Tipper(tipper_array=t, tipper_err_array=t_err, freq=period)
        self._Tipper.tipper_invsigcov = np.zeros((nf, 2, 2), dtype=np.complex)
        self._Tipper.tipper_residcov = np.zeros((nf, 1, 1), dtype=np.complex)
        
        input_dict = {'hz':0, 'hx': 0, 'hy': 1}
        output_dict = {'hz':0, 'hx':0, 'hy':1}
        
        p_count = 0
        for per_attr in dir(self.Data):
            if 'period' in per_attr.lower():
                p_obj = getattr(self.Data, per_attr)
                p1 = float(p_obj.attr['value'])
                self._Tipper.freq[p_count] = 1./p1
                if hasattr(p_obj, 'T'):
                    t_block = getattr(p_obj, 'T')
                    for t_attr in dir(t_block):
                        if 'value_' in t_attr:
                            t_comp = getattr(t_block, t_attr)
                            ii = output_dict[t_comp.attr['output'].lower()]
                            jj = input_dict[t_comp.attr['input'].lower()]
                            t_value = [float(tt) for tt in t_comp.value.strip().split()]
                            t_value = complex(t_value[0], t_value[1])
                            self._Tipper.tipper[p_count, ii, jj] = t_value
                
                if hasattr(p_obj, 'TVAR'):
                    t_block = getattr(p_obj, 'TVAR')
                    for t_attr in dir(t_block):
                        if 'value_' in t_attr:
                            t_comp = getattr(t_block, t_attr)
                            ii = output_dict[t_comp.attr['output'].lower()]
                            jj = input_dict[t_comp.attr['input'].lower()]
                            t_value = float(t_comp.value)
                            self._Tipper.tipper_err[p_count, ii, jj] = t_value
                
                if hasattr(p_obj, 'TINVSIGCOV'):
                    t_block = getattr(p_obj, 'TINVSIGCOV')
                    for t_attr in dir(t_block):
                        if 'value_' in t_attr:
                            t_comp = getattr(t_block, t_attr)
                            ii = output_dict[t_comp.attr['output'].lower()]
                            jj = input_dict[t_comp.attr['input'].lower()]
                            t_value = [float(tt) for tt in t_comp.value.strip().split()]
                            t_value = complex(t_value[0], t_value[1])
                            self._Tipper.tipper_invsigcov[p_count, ii, jj] = t_value

                if hasattr(p_obj, 'TRESIDCOV'):
                    t_block = getattr(p_obj, 'TRESIDCOV')
                    for t_attr in dir(t_block):
                        if 'value_' in t_attr:
                            t_comp = getattr(t_block, t_attr)
                            ii = output_dict[t_comp.attr['output'].lower()]
                            jj = input_dict[t_comp.attr['input'].lower()]
                            t_value = [float(tt) for tt in t_comp.value.strip().split()]
                            t_value = complex(t_value[0], t_value[1])
                            self._Tipper.tipper_residcov[p_count, ii, jj] = t_value

                p_count += 1
            
    def _get_info_from_element(self, element):
        """
        Get information from an element, including name, attr, value
        
        Arguments
        ------------
            **element** : ET.Element
            
        Returns
        ---------
            **XML_element** XML_element Object
        """
        
        return_obj = XML_element(None, None, None)
        return_obj._name = element.tag
        try:
            return_obj._value = element.text.strip()
        except AttributeError:
            return_obj._value = None
        return_obj._attr = element.attrib
        
        return return_obj
        
    def _get_attr_name(self, parent, attr_name):
        """
        make attribute name, if one already exists, then add a number to it
        
        ex. attribute_01
        """
        attr_name = attr_name.replace('.', '')
        if attr_name == 'Period':
            if not hasattr(parent, 'Period_00'):
                return 'Period_00'
            else:
                for ii in range(0, 100):
                    new_attr_name = '{0}_{1:02}'.format(attr_name, ii)
                    if not hasattr(parent, new_attr_name):
                        break
                return new_attr_name
            
        if hasattr(parent, attr_name):
            for ii in range(0, 100):
                new_attr_name = '{0}_{1:02}'.format(attr_name, ii)
                if not hasattr(parent, new_attr_name):
                    break
    
        else:
            new_attr_name = attr_name
        
        return new_attr_name
        
    def _read_element(self, element):
        """
        read a given element and return something useful
        """

        child = self._get_info_from_element(element)
        
        children = element.getchildren()
        if len(children) > 0:
            for child_00 in children:
                attr_name = self._get_attr_name(child, child_00.tag)
                setattr(child, attr_name, self._get_info_from_element(child_00))
                
                children_01 = child_00.getchildren()
                if len(children_01) > 0:
                    parent_01 = getattr(child, attr_name)
                    for child_01 in children_01:
                        attr_01_name = self._get_attr_name(parent_01,
                                                     child_01.tag)
                                                     
                        setattr(parent_01, 
                                attr_01_name,
                                self._get_info_from_element(child_01))
                                
                        children_02 = child_01.getchildren()
                        if len(children_02) > 0:
                            parent_02 = getattr(parent_01, attr_01_name)
                            for child_02 in children_02:
                                attr_02_name = self._get_attr_name(parent_02, 
                                                             child_02.tag)
                                                             
                                setattr(parent_02,
                                        attr_02_name,
                                        self._get_info_from_element(child_02))
                                        
                                children_03 = child_02.getchildren()
                                if len(children_03) > 0:
                                    parent_03 = getattr(parent_02, attr_02_name)
                                    for child_03 in children_03:
                                        attr_03_name = self._get_attr_name(parent_03, 
                                                                     child_03.tag)
                                                                     
                                        setattr(parent_03,
                                                attr_03_name,
                                                self._get_info_from_element(child_03))
                        
        
        return child 
    
    @property
    def Z(self):
        """
        get z information
        """
        return self._Z
                 
    @Z.setter
    def Z(self, z_object):
        """
        set z object
        """
        
        if type(z_object) is not mtz.Z:
            raise MT_XML_Error('To set Z, input needs to be an mtpy.core.z.Z object')
        
        self._Z = z_object
        
    @property
    def Tipper(self):
        """
        get Tipper information
        """        
        return self._Tipper
                    
    @Tipper.setter
    def Tipper(self, t_object):
        """
        set z object
        """
        
        if type(t_object) is not mtz.Tipper:
            raise MT_XML_Error('To set Z, input needs to be an mtpy.core.z.Z object')
        
        self._Tipper = t_object

#==============================================================================
# exceptions
#==============================================================================
class MT_XML_Error(Exception):
    pass