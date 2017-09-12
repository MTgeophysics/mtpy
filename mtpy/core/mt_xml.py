# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:39:40 2017

@author: jpeacock
"""

import os
import mtpy.core.edi as mtedi
import xml.etree.cElementTree as ET
import datetime
from xml.dom import minidom

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
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    @property
    def value(self):
        return self._value
    @property
    def attr(self):
        return self._attr
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
                                 **{'Filename':XML_element('Filename', None, None),
                                    'GroupKey':XML_element('GroupKey', None, 0),
                                    'OrderKey':XML_element('OrderKey', None, 0)})
                                    
        
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
                                  'ConditionsOfUse':XML_element('ConditionsOfUse', None, conditions_of_use)})

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
                                    'Electrode':XML_element('Electrode', None, None,
                                                      **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None)}),
                                    'Magnetometer':XML_element('Magnetometer', None, None,
                                                        **{'Type':XML_element('Type', None, None),
                                                         'Manufacturer':XML_element('Manufacturer', None, None),
                                                         'Id':XML_element('Id', None, None)}),
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
                                                                     'Author':XML_element('Author', None, None)}),
                                        'SignConvention':XML_element('SignConvention', None, r'exp(+i\omega t)'),
                                        'RemoteRef':XML_element('RemoteRef', {'type':'Robust Remote Processing'}, None),
                                        'RemoteInfo':XML_element('RemoteInfo', None, None, 
                                                           **{'Project':XML_element('Project', None, None),
                                                              'Survey':XML_element('Survey', None, None),
                                                              'ID':XML_element('ID', None, None),
                                                              'Name':XML_element('Name', None, None),
                                                              'YearCollected':XML_element('YearCollected', None, None),
                                                              'Location':XML_element('Location', {'datum':'WGS84'}, None,
                                                                               **{'Latitude':XML_element('Latitude', None, None),
                                                                                  'Longitude':XML_element('Longitude', None, None),
                                                                                  'Elevation':XML_element('Elevation', {'units':'meters'}, None)})
                                                               })})
        self.InputChannels = XML_element('InputChannels', {'ref':'site', 'units':'m'}, None)
        self.OutputChannels = XML_element('OutputChannels', {'ref':'site', 'units':'m'}, None)
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
                                   
        
        for key in kwargs.keys():
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
                print 'Not including: {0}'.format(attr_00_name)                     

        # write the file
        with open(self.cfg_fn, 'w') as fid:
            fid.write('\n'.join(line_list))
            
        # show the user something happened
        print '-'*50
        print '    Wrote xml configuration file to {0}'.format(self.cfg_fn)
        print '-'*50
        
        

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
                                 for a_key in XML_element_obj._attr.keys()])
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
# Site
#==============================================================================
#class XML_Site(object):
#    """
#    Site
#    """
#
#    def __init__(self, **kwargs):
#        self.project = XML_element('Project', None, None)
#        self.end_date = XML_element('End', None, None)
#        self.id = XML_element('Id', None, None)
#        self.Location = XML_Location()
#        self.run_list = XML_element('RunList', None, None)
#        self.start_date = XML_element('Start', None, None)
#        self.survey = XML_element('Survey', None, None)
#        self.year_collected = XML_element('YearCollected', None, None)
#        self.acqby = XML_element('AcquiredBy', None, None)
#
#
#class XML_Location(object):
#    
#    def __init__(self, **kwargs):
#        self.latitude = XML_element('Latitude', None, None)
#        self.longitude = XML_element('Longitude', None, None)
#        self.elevation = XML_element('Elevation', {'units':'meters'}, None),
#        self.declination = XML_element('Declination', {'epoch':'1995'}, None)
#
#class XML
#==============================================================================
#  EDI to XML
#==============================================================================
class MT_XML(XML_Config):
    """
    convert an EDI file to XML format
    
    """
    
    def __init__(self, **kwargs):
        
        XML_Config.__init__(self, **kwargs )
        self.edi_fn = None
        self.xml_fn = None
        self.cfg_fn = None
        
        self.parent_element = None
        
        self.edi_obj = None
        
        self._order_list = ['Description',
                            'ProductId',
                            'SubType',
                            'Notes',
                            'Tags',
                            'ExternalUrl',
                            'TimeSeriesArchived',
                            'Image',
                            'Original',
                            'Attachment',
                            'Provenance',
                            'Copyright',
                            'Site',
                            'FieldNotes',
                            'ProcessingInfo',
                            'StatisticalEstimates',
                            'DataTypes',
                            'InputChannels',
                            'OutputChannels',
                            'Data',
                            'PeriodRange']
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            

        
    def read_edi(self, edi_fn=None):
        if edi_fn is not None:
            self.edi_fn = edi_fn
            
        self.edi_obj = mtedi.Edi(self.edi_fn)
        self.get_edi_info()
        
    def read_cfg(self, cfg_fn=None):
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn

        self.read_cfg_file(self.cfg_fn)
        
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
        

    
    def get_edi_info(self, cfg_fn=None, edi_fn=None):
        """
        get information from config file and edi file
        """
        
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn
        if edi_fn is not None:
            self.edi_fn = edi_fn
            
        if self.cfg_fn is not None:
            self.read_cfg_file()
        
        # --> extract information from EDI files
        # Site information
        self.Site.DateCollected = XML_element('DateCollected', None, 
                                        self.edi_obj.Header.acqdate)
        self.Site.Id._value = self.edi_obj.station
        self.Site.AcquiredBy._value = self.edi_obj.Header.acqby
        self.Site.Start._value = self.edi_obj.Header.acqdate
        self.Site.End._value = self.edi_obj.Header.acqdate
        self.Site.RunList._value = '1'
        self.Site.Location.Latitude._value = '{0:.6f}'.format(self.edi_obj.Header.lat)
        self.Site.Location.Longitude._value = '{0:.6f}'.format(self.edi_obj.Header.lon)
        self.Site.Location.Elevation._value = '{0:.3f}'.format(self.edi_obj.Header.elev)

       
        # processing information
        self.ProcessingInfo.RemoteInfo.Project._value = self.Project._value       
        self.ProcessingInfo.RemoteInfo.Survey._value = self.Survey._value
        self.ProcessingInfo.RemoteInfo.YearCollected._value = self.Site.DateCollected._value
        
       
        # Field Notes
        self.FieldNotes.Magnetometer.HX = XML_element('HX', None, 
                                                str(self.edi_obj.Define_measurement.meas_hx.id))
        self.FieldNotes.Magnetometer.HY = XML_element('HY', None, 
                                                str(self.edi_obj.Define_measurement.meas_hy.id))
        try:        
            self.cfg_obj.FieldNotes.Magnetometer.HZ = XML_element('HZ', None,
                                                            str(self.edi_obj.Define_measurement.meas_hz.id))
        except AttributeError:
            pass
        #TODO: need to fill in more information on dipoles and magnetometers        
        
        # Input Channels
        attr_dict = {'name':'hx', 
                     'z': '{0:.1f}'.format(0),
                     'y':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hx.y),
                     'x':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hx.x),
                     'orientation':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hx.azm)}
        self.InputChannels.Magnetic_HX = XML_element('Magnetic', attr_dict, None)
        
        attr_dict = {'name':'hy', 
                     'z': '{0:.1f}'.format(0),
                     'y':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hy.y),
                     'x':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hy.x),
                     'orientation':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hy.azm)}
        self.InputChannels.Magnetic_HY = XML_element('Magnetic', attr_dict, None)
        
        # Output Channels
        try:
            attr_dict = {'name':'hz', 
                         'z': '{0:.1f}'.format(0),
                         'y':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hz.y),
                         'x':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hz.x),
                         'orientation':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_hz.azm)}
            self.OutputChannels.Magnetic_HZ = XML_element('Magnetic', attr_dict, None)
        except AttributeError:
            print 'No HZ Information'
        
        attr_dict = {'name':'ex', 
                     'z': '{0:.1f}'.format(0),
                     'y':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ex.y),
                     'x':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ex.x),
                     'x2':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ex.y2),
                     'y2':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ex.x2)}
        self.OutputChannels.Electric_EX = XML_element('Electric', attr_dict, None)
                                                       
        attr_dict = {'name':'ey', 
                     'z': '{0:.1f}'.format(0),
                     'y':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ey.y),
                     'x':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ey.x),
                     'x2':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ey.y2),
                     'y2':'{0:.1f}'.format(self.edi_obj.Define_measurement.meas_ey.x2)}
        self.OutputChannels.Electric_EY = XML_element('Electric', attr_dict, None)
   
        self.PeriodRange._attr = {'min':'{0:.5g}'.format(1./self.edi_obj.Z.freq.max()),
                                  'max':'{0:.5g}'.format(1./self.edi_obj.Z.freq.min()),
                                  'units':'seconds'}

        self.format_data()
        
    def format_data(self):
        """
        format the Z and tipper data apporpriately
        """
        # --> useful variables
        comp_dict_z = {(0, 0):('Zxx', 'Hx', 'Ex'), 
                       (0, 1):('Zxy', 'Hy', 'Ex'), 
                       (1, 0):('Zyx', 'Hx', 'Ey'),
                       (1, 1):('Zyy', 'Hy', 'Ey')}
                            
        comp_dict_t = {(0, 0):('Tx', 'Hx', 'Hz'),
                       (0, 1):('Ty', 'Hy', 'Hz')}
                       
        header_dict = {}
        header_dict['Z'] = XML_element('Z',{'units':'[mV/km]/[nT]', 
                                      'type':'complex', 
                                      'size':'2 2'}, None)
        header_dict['Z.VAR'] = XML_element('Z.VAR', {'type':'real', 'size':'2 2'},
                                     None)
        header_dict['T'] = XML_element('T', {'units':'[]', 
                                       'type':'complex',
                                       'size':'1 2'}, None)
        header_dict['T.VAR'] = XML_element('T.VAR', {'type':'real', 
                                               'size':'1 2'}, None)
        
        attr_dict = {}
        attr_dict['z'] = 'z'
        attr_dict['z.var'] = 'z_err'
        attr_dict['t'] = 'tipper'
        attr_dict['t.var'] = 'tipper_err'

        nf = self.edi_obj.Z.freq.size 
        
        # determine whether to write tipper data or not
        if self.edi_obj.Tipper.tipper is not None:
            nz_tipper = mtedi.np.any(self.edi_obj.Tipper.tipper) == 0
            if nz_tipper == True:
                write_tipper = False
            else:
                write_tipper = True
        else: 
            write_tipper = False
            
        # set the estimates to write out
        if write_tipper == True:
            estimates = ['Z', 'Z.VAR',  'T', 'T.VAR']
        else:
            estimates = ['Z', 'Z.VAR']

        # make the data element
        self.Data = XML_element('Data', {'count':str(nf)}, None)  
        
        # loop through each period and add appropriate information
        for f_index, freq in enumerate(self.edi_obj.Z.freq):
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
                estimate_name = estimate.replace('.', '_')
                setattr(d_attr, estimate_name, header_dict[estimate])
                c_attr = getattr(d_attr, estimate_name)
                if 'z' in estimate.lower():
                    count = 0
                    for e_index in range(2):
                        for h_index in range(2):
                            c = comp_dict_z[(e_index, h_index)]
                            c_dict = {'name':c[0], 'input':c[1], 'output':c[2]}
                            
                            if estimate.lower() == 'z':
                                z_value = self.edi_obj.Z.z[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e} {1:<+.8e}'.format(z_value.real, 
                                                                         z_value.imag)
                            elif estimate.lower() == 'z.var':
                                z_value = self.edi_obj.Z.z_err[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e}'.format(z_value)
                            
                            setattr(c_attr, 
                                    'value_{0:02}'.format(count),
                                    XML_element('value', c_dict, c_value)) 
                            count += 1
                            
                if 't' in estimate.lower() and write_tipper == True:
                    count = 0
                    for e_index in range(1):
                        for h_index in range(2):
                            c = comp_dict_t[(e_index, h_index)]
                            c_dict = {'name':c[0], 'input':c[1], 'output':c[2]}
                            
                            if estimate.lower() == 't':
                                z_value = self.edi_obj.Tipper.tipper[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e} {1:<+.8e}'.format(z_value.real, 
                                                                         z_value.imag)
                            elif estimate.lower() == 't.var':
                                z_value = self.edi_obj.Tipper.tipper_err[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e}'.format(z_value)
                                                                      
                            setattr(c_attr, 
                                    'value_{0:02}'.format(count),
                                    XML_element('value', c_dict, c_value)) 
                            count += 1
        
        
    def write_element(self, parent_et, XML_element_obj):
        """
        make a new element 
        """
        if XML_element_obj._attr is None:
            XML_element_obj._attr = {}

            
        new_element = ET.SubElement(parent_et, XML_element_obj._name, XML_element_obj._attr)
        new_element.text = XML_element_obj._value
        return new_element
       
    def write_xml_from_edi(self, edi_fn=None, xml_fn=None, cfg_fn=None):
        """
        write xml from edi
        """
        if edi_fn is not None:
            self.edi_fn = edi_fn
        if xml_fn is not None:
            self.xml_fn = xml_fn
        if edi_fn is not None:
            self.cfg_fn = cfg_fn
        
        if self.xml_fn is None:
            self.xml_fn = '{0}.xml'.format(self.edi_fn[0:-4])
        
        if self.cfg_fn is not None:
            self.read_cfg_file()
            
        if self.edi_fn is not None:
            self.read_edi()
        else:
            raise MT_XML_Error('Need to input an EDI file to convert')

        # make the top of the tree element        
        emtf = ET.Element('EM_TF')
        
        # loop over the important information sections
        for element in self._order_list:
            # get the information for the given element
            d_00_obj = getattr(self, element)
            
            element_00 = self.write_element(emtf, d_00_obj)

            key_00_list = self._get_attr_keys(d_00_obj)
            if len(key_00_list) !=0:
                for key_00 in key_00_list:
                    d_01_obj = getattr(d_00_obj, key_00)
                    element_01 = self.write_element(element_00, d_01_obj)
                    
                    key_01_list = self._get_attr_keys(d_01_obj)
                    if len(key_01_list) !=0:
                        for key_01 in key_01_list:
                            d_02_obj = getattr(d_01_obj, key_01)
                            element_02 = self.write_element(element_01, d_02_obj)
                            
                            key_02_list = self._get_attr_keys(d_02_obj)
                            if len(key_02_list) !=0:
                                for key_02 in key_02_list:
                                    d_03_obj = getattr(d_02_obj, key_02)
                                    element_03 = self.write_element(element_02, d_03_obj)
                                    
                                    key_03_list = self._get_attr_keys(d_03_obj)
                                    if len(key_03_list) !=0:
                                        for key_03 in key_03_list:
                                            d_04_obj = getattr(d_03_obj, key_03)
                                            element_04 = self.write_element(element_03, d_04_obj)
                                    
                                            key_04_list = self._get_attr_keys(d_04_obj)
                                            if len(key_04_list) !=0:
                                                for key_04 in key_04_list:
                                                    d_05_obj = getattr(d_04_obj, key_04)
                                                    element_05 = self.write_element(element_04, d_05_obj)
                
  
        #--> write xml file
        with open(self.xml_fn, 'w') as fid:
            fid.write(ET.tostring(emtf))
        
        print '-'*50
        print '    Wrote xml file to: {0}'.format(self.xml_fn)
        print '-'*50
#        # make a nice print out
#        reparsed = minidom.parseString(ET.tostring(emtf, 'utf-8'))
#        print(reparsed.toprettyxml(indent='    ')) 
        
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
        if hasattr(parent, attr_name):
            for ii in range(1, 10):
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
        
#==============================================================================
# exceptions
#==============================================================================
class MT_XML_Error(Exception):
    pass