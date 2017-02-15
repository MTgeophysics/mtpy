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

dt_fmt = '%Y-%m-%dT%H:%M:%S'
#==============================================================================
# Inputs
#==============================================================================
edi_fn = r"c:\Users\jpeacock\Documents\Geothermal\Washington\MSHS\EDI_Files_birrp_mshs\Rotated_m18_deg\ms11.edi"
cfg_fn = r"C:\Users\jpeacock\Documents\PyScripts\xml_cfg_test.cfg"


class Dummy(object):
    def __init__(self, name, attr, text, **kwargs):
        self._name = name
        self._attr = attr
        self._value = text
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

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
                     
estimates = [Dummy('Estimate', {'type':'real', 'name':'VAR'}, None,
                   **{'Description':Dummy('Description', None, 'Variance'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Error Estimate'),
                      'Tag':Dummy('Tag', None, 'Variance')}),
            
             Dummy('Estimate', {'type':'complex', 'name':'COV'}, None,  
                   **{'Description':Dummy('Description', None, 
                                          'Full covariance between each two TF components'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Error Estimate'),
                      'Tag':Dummy('Tag', None, 'Covariance')}),
             
             Dummy('Estimate', {'type':'complex', 'name':'INVSIGCOV'}, None,  
                   **{'Description':Dummy('Description', None, 
                                          'Inverse Coherent Signal Power Matrix S'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Signal Power Estimate'),
                      'Tag':Dummy('Tag', None, 'inverse_signal_covariance')}),
             
             Dummy('Estimate',{'type':'complex', 'name':'RESIDCOV'}, None,  
                   **{'Description':Dummy('Description',None, 
                                          'Residual Covariance N'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Error Estimate'),
                      'Tag':Dummy('Tag', None, 'Coherence')}),
             
             Dummy('Estimate', {'type':'complex', 'name':'COH'}, None,  
                   **{'Description':Dummy('Description', None, 'Coherence'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Signal Coherence'),
                      'Tag':Dummy('Tag', None, 'Coherence')}),
             
             Dummy('Estimate', {'type':'complex', 'name':'PREDCOH'}, None,  
                   **{'Description':Dummy('Description', None, 'Multiple Coherence'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Signal Coherence'),
                      'Tag':Dummy('Tag', None, 'Multiple_Coherence')}),
             
             Dummy('Estimate', {'type':'complex', 'name':'SIGAMP'}, None,  
                   **{'Description':Dummy('Description', None, 'Signal Amplitude'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Signal Power Estimates'),
                      'Tag':Dummy('Tag', None, 'Signal_Amplitude')}),
             
             Dummy('Estimate', {'type':'complex', 'name':'SIGNOISE'}, None,  
                   **{'Description':Dummy('Description', None, 'Signal Noise'),
                      'ExternalUrl':Dummy('ExternalUrl', None, None),
                      'Intention':Dummy('Intention', None, 'Error Estimates'),
                      'Tag':Dummy('Tag', None, 'Signal_Noise')})]
                       
data_types = [Dummy('DataType',
                    {'units':'[mV/km]/[nT]', 
                     'name':'Z', 
                     'input':'H',
                     'output':'E'},
                     None, 
                     **{'Description':Dummy('Description', None, 'MT impedance'),
                       'ExternalUrl':Dummy('ExternalUrl', None, None),
                       'Intention':Dummy('Intention', None, 'primary data type'),
                       'Tag':Dummy('Tag', None, 'impedance')}),
              
              Dummy('DataType',
                    {'units':'[]', 
                     'name':'Z', 
                     'input':'H',
                     'output':'H'},
                     None, 
                     **{'Description':Dummy('Description', None, 
                                            'Tipper-Vertical Field Transfer Function'),
                       'ExternalUrl':Dummy('ExternalUrl', None, None),
                       'Intention':Dummy('Intention', None, 'primary data type'),
                       'Tag':Dummy('Tag', None, 'tipper')})]
#==============================================================================
# Useful Functions
#==============================================================================                
class XML_Config(object):
    """
    Class to deal with configuration files for xml.
    
    
    """    
    
    def __init__(self, **kwargs):
        self.cfg_fn = None
        
        # Initialize the default attributes and values
        self.Description = Dummy('Description', 
                                 None, 
                                 'Magnetotelluric Transfer Functions')
                                 
        self.ProductId = Dummy('ProductID', None, None)
        
        self.Project = Dummy('Project', None, None)
        
        self.Survey = Dummy('Survey', None, None)
        
        self.Country = Dummy('Country', None, None)
       
        self.SubType = Dummy('SubType', None, 'MT_FT')
       
        self.Notes = Dummy('Notes', None, None)
        
        self.Tags = Dummy('Tags', None, 'impedance, tipper')
        
        
        self.Image = Dummy('Image', None, None, 
                           **{'PrimaryData':Dummy('PrimaryData', None, None),
                              'Filename':Dummy('Filename', None, None)})
                              
        
        self.Original = Dummy('Original', None, None,
                              **{'Attachment':Dummy('Attachement', None, None),
                                 'Filename':Dummy('Filename', None, None)})
        
       
        self.TimeSeriesArchived = Dummy('TimeSeriesArchived', None, None, 
                                        **{'Value':Dummy('Value', None, 0), 
                                           'URL':Dummy('URL', None, None)})
        
        self.ExternalUrl = Dummy('ExternalUrl', None, None, 
                                 **{'Description':Dummy('Description', None, None),
                                    'Url':Dummy('Url', None, None)})
        
        self.PrimaryData = Dummy('PrimaryData', None, None, 
                                 **{'Filename':Dummy('Filename', None, None),
                                    'GroupKey':Dummy('GroupKey', None, 0),
                                    'OrderKey':Dummy('OrderKey', None, 0)})
                                    
        
        self.Attachment = Dummy('Attachement', None, None, 
                                **{'Filename':Dummy('Filename', None, None),
                                   'Description':Dummy('Description', None, 
                                                       'Original file use to produce XML')})
                                   
        
        self.Provenance = Dummy('Provenance', None, None,
                                **{'CreatTime':Dummy('CreationTime', None, 
                                                     datetime.datetime.strftime(
                                                     datetime.datetime.utcnow(), 
                                                     dt_fmt)),
                                    'CreatingApplication':Dummy('CreatingApplication', 
                                                                None, 
                                                                'MTpy.core.mtxml'),
                                    'Submitter':Dummy('Submitter', None, None,
                                                      **{'Name':Dummy('Name', None, None),
                                                         'Email':Dummy('Email', None, None),
                                                         'Org':Dummy('Org', None, None),
                                                         'OrgURL':Dummy('OrgURL', None, None)}),
                                    'Creator':Dummy('Creator', None, None,
                                                    **{'Name':Dummy('Name', None, None),
                                                       'Email':Dummy('Email', None, None),
                                                       'Org':Dummy('Org', None, None),
                                                       'OrgURL':Dummy('OrgURL', None, None)})})
                                                       
        self.Copyright = Dummy('Copyright', None, None,
                               **{'Citation':Dummy('Citation', None, None,
                                                   **{'Title':Dummy('Title', None, None),
                                                      'Authors':Dummy('Authors', None, None),
                                                      'Year':Dummy('Year', None, None),
                                                      'Journal':Dummy('Journal', None, None),
                                                      'Volume':Dummy('Volume', None, None),
                                                      'DOI':Dummy('DOI', None, None)}),
                                  'ReleaseStatus':Dummy('ReleaseStatus', None, 'Closed'),
                                  'ConditionsOfUse':Dummy('CondictionsOfUse', None, conditions_of_use)})

        self.Site = Dummy('Site', None, None,
                          **{'Project':Dummy('Project', None, None),
                             'Survey':Dummy('Survey', None, None),
                             'YearCollected':Dummy('YearCollected', None, None),
                             'Id':Dummy('Id', None, None),
                             'Location':Dummy('Location', None, None,
                                              **{'Latitude':Dummy('Latitude', None, None),
                                                 'Longitude':Dummy('Longitude', None, None),
                                                 'Elevation':Dummy('Elevation', {'units':'meters'}, None),
                                                 'Declination':Dummy('Declination', {'epoch':'1995'}, None)}),
                             'AcquiredBy':Dummy('AcquiredBy', None, None),
                             'Start':Dummy('Start', None, None),
                             'End':Dummy('End', None, None),
                             'RunList':Dummy('RunList', None, None)})
                             
        self.FieldNotes = Dummy('FieldNotes', None, None,
                                **{'Instrument':Dummy('Instrument', None, None,
                                                      **{'Type':Dummy('Type', None, None),
                                                         'Manufacturer':Dummy('Manufacturer', None, None),
                                                         'Id':Dummy('Id', None, None),
                                                         'Settings':Dummy('Settings', None, None)}),
                                    'Electrode':Dummy('Electrode', None, None,
                                                      **{'Type':Dummy('Type', None, None),
                                                         'Manufacturer':Dummy('Manufacturer', None, None),
                                                         'Id':Dummy('Id', None, None)}),
                                    'Magnetometer':Dummy('Magnetometer', None, None,
                                                        **{'Type':Dummy('Type', None, None),
                                                         'Manufacturer':Dummy('Manufacturer', None, None),
                                                         'Id':Dummy('Id', None, None)}),
                                    'DataQualityNotes':Dummy('DataQualityNotes', None, None,
                                                             **{'Rating':Dummy('Rating', None, None),
                                                                 'GoodFromPeriod':Dummy('GoodFromPeriod', None, None),
                                                                 'GoodToPeriod':Dummy('GoodToPeriod', None, None),
                                                                 'Comments':Dummy('Comments', None, None)}),
                                    'DataQualityWarnings':Dummy('DataQualityWarnings', None, None,
                                                                **{'Flag':Dummy('Flag', None, 0),
                                                                    'Comments':Dummy('Comments', None, None)})})
                                  
        self.ProcessingInfo = Dummy('ProcessingInfo', None, None,
                                    **{'ProcessedBy':Dummy('ProcessedBy', None, None),
                                       'ProcessingSoftware':Dummy('ProcessingSoftware', None, None,
                                                                  **{'Name':Dummy('Name', None, None),
                                                                     'LastMod':Dummy('LastMod', None, None),
                                                                     'Author':Dummy('Author', None, None)}),
                                        'SignConvention':Dummy('SignConvention', None, r'exp(+i\omega t)'),
                                        'RemoteRef':Dummy('RemoteRef', {'type':'Robust Remote Processing'}, None),
                                        'RemoteInfo':Dummy('RemoteInfo', None, None, 
                                                           **{'Project':Dummy('Project', None, None),
                                                              'Survey':Dummy('Survey', None, None),
                                                              'ID':Dummy('ID', None, None),
                                                              'Name':Dummy('Name', None, None),
                                                              'YearCollected':Dummy('YearCollected', None, None),
                                                              'Location':Dummy('Location', {'datum':'WGS84'}, None,
                                                                               **{'Latitude':Dummy('Latitude', None, None),
                                                                                  'Longitude':Dummy('Longitude', None, None),
                                                                                  'Elevation':Dummy('Elevation', {'units':'meters'}, None)})
                                                               })})
        self.InputChannels = Dummy('InputChannels', {'ref':'site', 'units':'m'}, None)
        self.OutputChannels = Dummy('OutputChannels', {'ref':'site', 'units':'m'}, None)
        self.Data = Dummy('Data', {'count':0}, None)
        self.PeriodRange = Dummy('PeriodRange', None, None)
                                        
        self.Datum = Dummy('Datum', None, 'WGS84')
        self.Declination = Dummy('Declination', None, None)
                 
        self.StatisticalEstimates = Dummy('StatisticalEstimates', None, None)                           
        for ii, estimate in enumerate(estimates):
            setattr(self.StatisticalEstimates, 
                    'Estimate_{0:02}'.format(ii),
                    estimate)
        
        self.DataTypes = Dummy('DataTypes', None, None)
        for ii, d_type in enumerate(data_types):
            setattr(self.DataTypes, 'DataType_{0:02}'.format(ii), d_type)
                                   
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def read_cfg_file(self, cfg_fn=None):
        """
        Read in a cfg file making all key = value pairs attribures of 
        XML_Config.  Being sure all new attributes are Dummy objects.
        
        The assumed structure of the xml.cfg file is similar to:
            ``# XML Configuration File MTpy
 
            Attachement.Description = Original file use to produce XML
            Attachement.Filename = None
             
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
            if attr_00_name == 'Data':
                continue
            
            # get the given attribute
            attr_00 = getattr(self, attr_00_name)
            
            # make sure it is of Dummy instance
            if isinstance(attr_00, Dummy):
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
        
    def _read_cfg_line(self, line):
        """
        read a configuration file line to make the appropriate attribute
        have the correct values and attributes.
        
        porbably should think of a better name for Dummy objects that are
        attributes of self.
        """
        
        # split the line by the =
        line_list = line.strip().split('=')
        # split the keys by . to get the different attributes
        key_list = line_list[0].strip().split('.')
        value = line_list[1].strip()
        
        # loop over the keys to set them appropriately 
        for ii, key in enumerate(key_list):
            # get the name of the key and any attributes it might have
            name, attr = self._read_cfg_key(key)
            
            # if its the first key, see if its been made an attribute yet
            if ii == 0: 
                if not hasattr(self, name):
                    setattr(self, name, Dummy(name, None, None))
                # for looping purposes we need to get the current Dummy object
                cfg_attr = getattr(self, name)
                # be sure to set any attributes, need to do this here because
                # the test for hasattr will only make a new one if there
                # isn't one already, but since most things in the cfg file
                # are already attributes of self, they already exist.
                cfg_attr._attr = attr
            else:
                if not hasattr(cfg_attr, name):
                    setattr(cfg_attr, name, Dummy(name, None, None))
                cfg_attr = getattr(cfg_attr, name) 
                cfg_attr._attr = attr
        
        # set the value of the current Dummy object
        cfg_attr._value = value

            
    def _read_cfg_key(self, key):
        """
        read a key from a cfg file and check to see if has any attributes
        in the form of:
            parent.name(attribute=0)(attribute=2) = value
        """
        attr = {}
        if '(' and ')' in key:
            key_list = key.split('(')
            for key_attr in key_list:
                k_list = key_attr.replace(')', '').split('=')
                attr[k_list[0].strip()] = k_list[1].strip()
                
        return key, attr
        

    def _write_cfg_line(self, dummy_obj, parent=None):
        """
        write a configuration file line in the format of:
        parent.attribute = value
        """        
        
        if parent is None:
            parent_str = ''
            
        elif type(parent) is list:
            parent_str = '.'.join([p._name for p in parent]+[''])
            
        elif isinstance(parent, Dummy):
            parent_str = '{0}.'.format(parent._name)
        
        if dummy_obj._attr is not None:
            attr_str = ''.join(['({0}={1})'.format(a_key, dummy_obj._attr[a_key]) 
                                 for a_key in dummy_obj._attr.keys()])
        else:
            attr_str = ''
        return '{0}{1}{2} = {3}'.format(parent_str,
                                        dummy_obj._name, 
                                        attr_str, 
                                        dummy_obj._value)
                                     
    def _get_attr_keys(self, attribute):
        return [a_key for a_key in sorted(attribute.__dict__.keys()) 
                if a_key not in ['_name', '_value', '_attr']]

#==============================================================================
#  EDI to XML
#==============================================================================
class EDI_to_XML(XML_Config):
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
        

    
    def get_info(self, cfg_fn=None, edi_fn=None):
        """
        get information from config file and edi file
        """
        
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn
        if edi_fn is not None:
            self.edi_fn = edi_fn
            
        self.read_edi()
        self.read_cfg()
        
        # --> extract information from EDI files
        # Site information
        self.Site.DateCollected = Dummy('DateCollected', None, 
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
        self.ProcessingInfo.RemoteInfo.Project._value = self.Project       
        self.ProcessingInfo.RemoteInfo.Survey._value = self.Survey
        self.ProcessingInfo.RemoteInfo.YearCollected._value = self.Site.DateCollected
        
       
        # Field Notes
        self.FieldNotes.Magnetometer.HX = Dummy('HX', None, 
                                                str(self.edi_obj.Define_measurement.meas_hx.acqchan))
        self.FieldNotes.Magnetometer.HY = Dummy('HY', None, 
                                                str(self.edi_obj.Define_measurement.meas_hy.acqchan))
        try:        
            self.cfg_obj.FieldNotes.Magnetometer.HZ = Dummy('HZ', None,
                                                            str(self.edi_obj.Define_measurement.meas_hz.acqchan))
        except AttributeError:
            pass
        #TODO: need to fill in more information on dipoles and magnetometers        
        
        # Input Channels
        attr_dict = {'name':'hx', 
                     'z': 0,
                     'y':self.edi_obj.Define_measurement.meas_hx.y,
                     'x':self.edi_obj.Define_measurement.meas_hx.x,
                     'orientation':self.edi_obj.Define_measurement.meas_hx.azm}
        self.InputChannels.Magnetic_HX = Dummy('Magnetic', attr_dict, None)
        
        attr_dict = {'name':'hy', 
                     'z': 0,
                     'y':self.edi_obj.Define_measurement.meas_hy.y,
                     'x':self.edi_obj.Define_measurement.meas_hy.x,
                     'orientation':self.edi_obj.Define_measurement.meas_hy.azm}
        self.InputChannels.Magnetic_HY = Dummy('Magnetic', attr_dict, None)
        
        # Output Channels
        try:
            attr_dict = {'name':'hz', 
                         'z': 0,
                         'y':self.edi_obj.Define_measurement.meas_hz.y,
                         'x':self.edi_obj.Define_measurement.meas_hz.x,
                         'orientation':self.edi_obj.Define_measurement.meas_hz.azm}
            self.InputChannels.Magnetic_HZ = Dummy('Magnetic', attr_dict, None)
        except AttributeError:
            print 'No HZ Information'
        
        attr_dict = {'name':'ex', 
                     'z': 0,
                     'y':self.edi_obj.Define_measurement.meas_ex.y,
                     'x':self.edi_obj.Define_measurement.meas_ex.x,
                     'x2':self.edi_obj.Define_measurement.meas_ex.y2,
                     'y2':self.edi_obj.Define_measurement.meas_ex.x2}
        self.InputChannels.Electric_EX = Dummy('Electric', attr_dict, None)
                                                       
        attr_dict = {'name':'ey', 
                     'z': 0,
                     'y':self.edi_obj.Define_measurement.meas_ey.y,
                     'x':self.edi_obj.Define_measurement.meas_ey.x,
                     'x2':self.edi_obj.Define_measurement.meas_ey.y2,
                     'y2':self.edi_obj.Define_measurement.meas_ey.x2}
        self.InputChannels.Electric_EY = Dummy('Electric', attr_dict, None)
   
        self.PeriodRange._attr = {'min':'{0:.5g}'.format(1./self.edi_obj.Z.freq.max()),
                                  'max':'{0:.5g}'.format(1./self.edi_obj.Z.freq.min())}

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
        header_dict['Z'] = Dummy('Z',{'units':'[mV/km]/[nT]', 
                                      'type':'complex', 
                                      'size':'2 2'}, None)
        header_dict['Z.VAR'] = Dummy('Z.VAR', {'type':'real', 'size':'2 2'},
                                     None)
        header_dict['T'] = Dummy('T', {'units':'[]', 
                                       'type':'complex',
                                       'size':'1 2'}, None)
        header_dict['T.VAR'] = Dummy('T.VAR', {'type':'real', 
                                               'size':'1 2'}, None)
        
        attr_dict = {}
        attr_dict['z'] = 'z'
        attr_dict['z.var'] = 'z_err'
        attr_dict['t'] = 'tipper'
        attr_dict['t.var'] = 'tipper_err'

        nf = self.edi_obj.Z.freq.size        

        # make the data element
        self.Data = Dummy('Data', {'count':str(nf)}, None)  
        
        # loop through each period and add appropriate information
        for f_index, freq in enumerate(self.edi_obj.Z.freq):
            f_name = 'Period_{0:02}'.format(f_index)
            # set attribute period name with the index value
            # we are setting _name to have the necessary information so
            # we can name the attribute whatever we want.
            setattr(self.Data, f_name,
                    Dummy('Period', {'value':'{0:.6g}'.format(1./freq)}, None))
                   
            # Get information from data
            for estimate in ['Z', 'Z.VAR', 'T', 'T.VAR']:
                d_attr = getattr(self.Data, f_name)
                setattr(d_attr, estimate, header_dict[estimate])
                c_attr = getattr(getattr(self.Data, f_name), estimate)
                if 'z' in estimate:
                    for e_index in range(2):
                        for h_index in range(2):
                            c = comp_dict_z[(e_index, h_index)]
                            c_dict = {'name':c[0], 'input':c[1], 'output':c[2]}
                            
                            if estimate == 'z':
                                z_value = self.edi_obj.Z.z[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e} {1:<+.8e}'.format(z_value.real, 
                                                                         z_value.imag)
                            elif estimate == 'z.var':
                                z_value = self.edi_obj.Z.z_err[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e}'.format(z_value)
                                                                      
                            setattr(c_attr, 'value_{0:02}'.format(e_index+h_index),
                                    Dummy('value', c_dict, c_value))
                if 't' in estimate and self.mt_obj.Tipper.tipper is not None:
                    for e_index in range(1):
                        for h_index in range(2):
                            c = comp_dict_t[(e_index, h_index)]
                            c_dict = {'name':c[0], 'input':c[1], 'output':c[2]}
                            
                            if estimate == 't':
                                z_value = self.edi_obj.Tipper.tipper[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e} {1:<+.8e}'.format(z_value.real, 
                                                                         z_value.imag)
                            elif estimate == 't.var':
                                z_value = self.edi_obj.Tipper.tipper_err[f_index, e_index, h_index]                            
                                c_value = '{0:<+.8e}'.format(z_value)
                                                                      
                            setattr(c_attr, 'value_{0:02}'.format(e_index+h_index),
                                    Dummy('value', c_dict, c_value))            
        
        
    def write_element(self, parent_et, dummy_obj):
#        try:
        if dummy_obj._attr is None:
            dummy_obj._attr = {}
#        except AttributeError:
#            raise ValueError(dummy_obj)
            
        new_element = ET.SubElement(parent_et, dummy_obj._name, dummy_obj._attr)
        new_element.text = dummy_obj._value
        return new_element
       
    def write_xml(self, edi_fn=None, xml_fn=None, cfg_fn=None):
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
        
        if self.edi_fn is not None:
            self.get_info()

        # make the top of the tree element        
        emtf = ET.Element('EM_TF')
        
        # loop over the important information sections
        for element in self._order_list[0:12]:
            print element
            if element == 'PeriodRange':
                break
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
        read in an xml file
        """
        pass
#    def get_info_from_element(element):
#        return_obj = mtxml.Dummy()
#        return_obj._name = element.tag
#        try:
#            return_obj._text = element.text.strip()
#        except AttributeError:
#            return_obj._text = None
#        return_obj._attr = element.attrib
#        
#        return return_obj
#        
#    def get_attr_name(parent, attr_name):
#        if hasattr(parent, attr_name):
#            for ii in range(1, 10):
#                new_attr_name = '{0}_{1:02}'.format(attr_name, ii)
#                if not hasattr(parent, new_attr_name):
#                    break
#    
#        else:
#            new_attr_name = attr_name
#        
#        return new_attr_name
#    
#    xml_fn = r"C:\Users\jpeacock\Documents\PyScripts\ORL09bc_J9.xml"
#    
#    e2xml = mtxml.EDI_to_XML()
#    
#    et_xml = ET.parse(xml_fn)
#    
#    
#    def read_element(element):
#        """
#        read a given element and return something useful
#        """
#        
#        
#        
#        child = get_info_from_element(element)
#        
#        children = element.getchildren()
#        if len(children) > 0:
#            for child_00 in children:
#                attr_name = get_attr_name(child, child_00.tag)
#                setattr(child, attr_name, get_info_from_element(child_00))
#                
#                children_01 = child_00.getchildren()
#                if len(children_01) > 0:
#                    parent_01 = getattr(child, attr_name)
#                    for child_01 in children_01:
#                        attr_01_name = get_attr_name(parent_01,
#                                                     child_01.tag)
#                                                     
#                        setattr(parent_01, 
#                                attr_01_name,
#                                get_info_from_element(child_01))
#                                
#                        children_02 = child_01.getchildren()
#                        if len(children_02) > 0:
#                            parent_02 = getattr(parent_01, attr_01_name)
#                            for child_02 in children_02:
#                                attr_02_name = get_attr_name(parent_02, 
#                                                             child_02.tag)
#                                                             
#                                setattr(parent_02,
#                                        attr_02_name,
#                                        get_info_from_element(child_02))
#                                        
#                                children_03 = child_02.getchildren()
#                                if len(children_03) > 0:
#                                    parent_03 = getattr(parent_02, attr_02_name)
#                                    for child_03 in children_03:
#                                        attr_03_name = get_attr_name(parent_03, 
#                                                                     child_03.tag)
#                                                                     
#                                        setattr(parent_03,
#                                                attr_03_name,
#                                                get_info_from_element(child_03))
#                        
#        
#        return child
#    
#    
#    k = mtxml.Dummy()
#    root = et_xml.getroot()
#    for element_00 in root.getchildren():
#        setattr(k, element_00.tag, read_element(element_00))
#        
##==============================================================================
## Do the dirty work
##==============================================================================
test = EDI_to_XML()
test.edi_fn = r"c:\Users\jpeacock\Documents\iMush\imush_edi_files_final\Interpolated\mshs86.edi"
test.cfg_fn = r"c:\Users\jpeacock\Documents\Test.cfg"
test.write_xml()
#
#test.cfg_obj.write_cfg_file(r"C:\Users\jpeacock\Documents\PyScripts\xml_cfg_test_out.cfg")


    

      