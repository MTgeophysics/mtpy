# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:39:40 2017

@author: jpeacock
"""

import os
import mtpy.core.mt as mt
import xml.etree.ElementTree as ET
import datetime
from xml.dom import minidom

dt_fmt = '%Y-%m-%dT%H:%M:%S'
#==============================================================================
# Inputs
#==============================================================================
edi_fn = r"c:\Users\jpeacock\Documents\Geothermal\Washington\MSHS\EDI_Files_birrp_mshs\Rotated_m18_deg\ms11.edi"
cfg_fn = r"C:\Users\jpeacock\Documents\PyScripts\xml_cfg_test.cfg"


class Dummy(object):
    def __init__(self, **kwargs):
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

conditions_of_use = "All data and metadata for this survey are available free \
                     of charge and may be copied freely, duplicated and further\
                     distributed provided this data set is cited as the\
                     reference. While the author(s) strive to provide data and \
                     metadata of best possible quality, neither the author(s) \
                     of this data set, not IRIS make any claims, promises, or \
                     guarantees about the accuracy, completeness, or adequacy \
                     of this information, and expressly disclaim liability for\
                     errors and omissions in the contents of this file. \
                     Guidelines about the quality or limitations of the data \
                     and metadata, as obtained from the author(s), are \
                     included for informational purposes only."
                     
estimates = [Dummy(**{'_name':'Estimate(type=real)(name=VAR)',
                      'Description':'Variance',
                       'ExternalURL':None,
                       'Intention':'error estimate',
                       'Tag':'variance'}),
             Dummy(**{'_name':'Estimate(type=complex)(name=COV)',
                      'Description':'Full covariance between each two TF components',
                      'ExternalURL':None,
                      'Intention':'error estimate',
                      'Tag':'covariance'}),
             Dummy(**{'_name':'Estimate(type=complex)(name=INVSIGCOV)',
                       'Description':'Inverse Coheren Signal Power Matrix S',
                       'ExternalURL':None,
                       'Intention':'signal power estimate',
                       'Tag':'inverse_signal_covariance'}),
              Dummy(**{'_name':'Estimate(type=complex)(name=RESIDCOV)',
                       'Description':'Residual Covariance N',
                       'ExternalURL':None,
                       'Intention':'error estimate',
                       'Tag':'residual_covariance'}),
              Dummy(**{'_name':'Estimate(type=complex)(name=COH)',
                       'Description':'Coherence',
                       'ExternalURL':None,
                       'Intention':'signal coherence',
                       'Tag':'coherence'}),
              Dummy(**{'_name':'Estimate(type=complex)(name=PREDCOH)',
                       'Description':'Multiple Coherence',
                       'ExternalURL':None,
                       'Intention':'signal coherence',
                       'Tag':'multiple_coherence'}),
              Dummy(**{'_name':'Estimate(type=complex)(name=SIGAMP)',
                       'Description':'Signal Amplitude',
                       'ExternalURL':None,
                       'Intention':'signal power estimate',
                       'Tag':'signal_amplitude'}),
              Dummy(**{'_name':'Estimate(type=complex)(name=SIGNOISE)',
                       'Description':'Signal Noise',
                       'ExternalURL':None,
                       'Intention':'error estimate',
                       'Tag':'signal_noise'})]
                       
data_types = [Dummy(**{'_name':'DataType(units=[mV/km]/[nT])(type=complex)(name=Z)(input=H)(output=E)',
                       'Description':'MT impedance',
                       'ExternalURL':None,
                       'Intention':'primary data type',
                       'Tag':'impedance'}),
              Dummy(**{'_name':'DataType(units=[])(type=complex)(name=T)(input=H)(output=H)',
                       'Description':'Tipper-Vertical Field Transfer Function',
                       'ExternalURL':None,
                       'Intention':'primary data type',
                       'Tag':'tipper'})]
#==============================================================================
# Useful Functions
#==============================================================================                
class XML_Config(object):
    """
    class to deal with configuration files for xml
    """    
    
    def __init__(self, **kwargs):
        self.cfg_fn = None
        self.cfg_dict = None
        
        # Initialize the default attributes and values
        self.Description = 'Magnetotelluric Transfer Functions'
        self.ProductID = None
        self.Project = None
        self.Survey = None
        self.Country = None
        self.SubType = 'MT_TF'
        self.Notes = None
        self.Tags = 'impedance, tipper'
        
        self.Image = Dummy(**{'PrimaryData':None,
                              'Filename':None,
                              '_name':'Image'})
                              
        self.Original = Dummy(**{'Attachment':None,
                                 'Filename':None,
                                 '_name':'Original'})
        
        self.TimeSeriesArchived = Dummy(**{'Value':0, 
                                           'URL':None,
                                           '_name':'TimeSeriesArchived'})
        self.ExternalURL = Dummy(**{'Description':None,
                                    'URL':None,
                                    '_name':'ExternalURL'})
        self.PrimaryData = Dummy(**{'Filename':None,
                                    'GroupKey':0,
                                    'OrderKey':0,
                                    '_name':'PrimaryData'})
                                    
        self.Attachment = Dummy(**{'Filename':None,
                                   'Description':'Original file use to produce XML',
                                   '_name':'Attachment'})
                                   
        
        self.Provenance = Dummy(**{'CreatTime':datetime.datetime.strftime(
                                               datetime.datetime.utcnow(), 
                                               dt_fmt),
                                    'CreatingApplication':'MTpy.core.xml',
                                    'Submitter':Dummy(**{'Name':None,
                                                         'Email':None,
                                                         'Org':None,
                                                         'OrgURL':None,
                                                         '_name':'Submitter'}),
                                    'Creator':Dummy(**{'Name':None,
                                                       'Email':None,
                                                       'Org':None,
                                                       'OrgURL':None,
                                                       '_name':'Creator'}),
                                    '_name':'Provenance'})
                                                       
        self.Copyright = Dummy(**{'Citation':Dummy(**{'Title':None,
                                                      'Authors':None,
                                                      'Year':None,
                                                      'Journal':None,
                                                      'Volume':None,
                                                      'DOI':None,
                                                      '_name':'Citation'}),
                                  'ReleaseStatus':'Closed',
                                  'ConditionsOfUse':conditions_of_use,
                                  '_name':'Copyright'})
                                  
        self.ProcessingInfo = Dummy(**{'ProcessedBy':None,
                                       'ProcessingSoftware':Dummy(**{'Name':None,
                                                                     'LastMod':None,
                                                                     'Author':None,
                                                                     '_name':'ProcessingSoftware'}),
                                        'SignConvention':r'exp(+i\omega t)',
                                        'RemoteRef(type=Robust Remote Processing)':None,
                                        'RemoteInfo':Dummy(**{'_name':'RemoteInfo',
                                                              'Project':None,
                                                              'Survey':None,
                                                              'ID':None,
                                                              'Name':None,
                                                              'YearCollected':None,
                                                              'Location':Dummy(**{'_name':'Location(datum=WGS84)',
                                                                                  'Latitude':None,
                                                                                  'Longitude':None,
                                                                                  'Elevation(units=meters)':None}),
                                                               }),
                                        '_name':'ProcessingInfo'})
                                        
        self.Datum = None
        self.Declination = None
        
        self.Instrument = Dummy(**{'Type':None,
                                   'Manufacturer':None,
                                   'Id':None,
                                   'Settings':None,
                                   '_name':'Instrument'})

        self.Electrode = Dummy(**{'Type':None,
                                  'Manufacturer':None,
                                  'Id':None,
                                  '_name':'Electrode'})

        self.Magnetometer = Dummy(**{'Type':None,
                                     'Manufacturer':None,
                                     'Id':None,
                                     '_name':'Magnetometer'})
                                     
        self.DataQualityNotes = Dummy(**{'Rating':None,
                                         'GoodFromPeriod':None,
                                         'GoodToPeriod':None,
                                         'Comments':None,
                                         '_name':'DataQualityNotes'})
                                         
        self.DataQualityWarnings = Dummy(**{'Flag':0,
                                            'Comments':None,
                                            '_name':'DataQualityWarnings'})
                                            
        self.StatisticalEstimates = estimates
        self.DataTypes = data_types
                                   
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def read_cfg_file(self, cfg_fn=None):
        if cfg_fn is not None:        
            self.cfg_fn = cfg_fn
            
        if not os.path.isfile(self.cfg_fn):
            raise NameError('Could not find {0}'.format(self.cfg_fn))
            
        with open(self.cfg_fn, 'r') as fid:
            lines = fid.readlines()
                        
        self.cfg_dict = {}    
        for line in lines:
            # skip comments
            if line[0] == '#':
                pass
            
            # skip blank lines
            elif line == '\n' or line == '\r':
                pass
            # else assume the line is metadata separated by = 
            else:
                line_list = line.strip().split('=')
                key = line_list[0].strip()
                value = line_list[1].strip()
                
                # if there is a dot then there is a tree of key words
                if key.find('.') > 0:
                    key_list = key.split('.')
                    if len(key_list) == 2:
                        setattr(getattr(self, key_list[0]), 
                                key_list[1], value)
                    elif len(key_list) == 3:
                        setattr(getattr(getattr(self, key_list[0]), key_list[1]),
                                key_list[2], value)

                else:
                    setattr(self, key, value)
                    
    def write_cfg_file(self, cfg_fn=None):
        """
        write out configuration file        
        """
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn
            
        line_list = []
        
        for attr_00 in sorted(self.__dict__.keys()):
            attr_00_value = getattr(self, attr_00)
            
            if attr_00 == 'Data':
                continue
            
            if type(attr_00_value) in [int, float, str]:
                line_list.append(self._write_cfg_line(attr_00, 
                                                      getattr(self, attr_00)))
            elif isinstance(attr_00_value, Dummy):
                for attr_01 in sorted(attr_00_value.__dict__.keys()):
                    attr_01_value = getattr(attr_00_value, attr_01)
                    if type(attr_01_value) in [float, int, str] or attr_01_value is None:
                        line_list.append(self._write_cfg_line('{0}.{1}'.format(attr_00,
                                                                               attr_01),
                                                              attr_01_value))
                    elif isinstance(attr_01_value, Dummy):
                        for attr_02 in sorted(attr_01_value.__dict__.keys()):
                            attr_02_value = getattr(attr_01_value, attr_02)
                            if type(attr_02_value) in [float, int, str] or attr_02_value is None:
                                line_list.append(self._write_cfg_line('{0}.{1}.{2}'.format(attr_00,
                                                                                           attr_01,
                                                                                           attr_02),
                                                                      attr_02_value))
                            elif isinstance(attr_02_value, Dummy):
                                for attr_03 in sorted(attr_02_value.__dict__.keys()):
                                    attr_03_value = getattr(attr_02_value, attr_03)
                                    if type(attr_03_value) in [float, int, str] or attr_03_value is None:
                                        line_list.append(self._write_cfg_line('{0}.{1}.{2}.{3}'.format(attr_00,
                                                                                                       attr_01,   
                                                                                                       attr_02,
                                                                                                       attr_03),
                                                                              attr_03_value))
                                    elif isinstance(attr_03_value, Dummy):
                                        for attr_04 in sorted(attr_03_value.__dict__.keys()):
                                            print attr_04
                                            attr_04_value = getattr(attr_03_value, attr_03)
                                            if type(attr_03_value) in [float, int, str] or attr_04_value is None:
                                                line_list.append(self._write_cfg_line('{0}.{1}.{2}.{3}.{4}'.format(attr_00,
                                                                                                                   attr_01,   
                                                                                                                   attr_02,
                                                                                                                   attr_03,
                                                                                                                   attr_04),
                                                                                      attr_04_value))
                line_list.append(' ')
        with open(self.cfg_fn, 'w') as fid:
            fid.write('\n'.join(line_list))
            
        print '-'*50
        print '    Wrote xml configuration file to {0}'.format(self.cfg_fn)
        print '-'*50

    def _write_cfg_line(self, key, value):
        return '{0} = {1}'.format(key, value)

#==============================================================================
#  EDI to XML
#==============================================================================
class EDI_to_XML(object):
    """
    convert an EDI file to XML format
    
    """
    
    def __init__(self, **kwargs):
        
        self.edi_fn = None
        self.xml_fn = None
        self.cfg_fn = None
        
        self.parent_element = None
        
        self.mt_obj = None
        self.cfg_obj = XML_Config()
        
        self._order_list = ['Description',
                            'ProductID',
                            'SubType',
                            'Notes',
                            'Tags',
                            'ExternalURL',
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
            
        self.mt_obj = mt.MT(self.edi_fn)
        
    def read_cfg(self, cfg_fn=None):
        if cfg_fn is not None:
            self.cfg_fn = cfg_fn
            
        self.cfg_obj = XML_Config()
        self.cfg_obj.read_cfg_file(self.cfg_fn)
        
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
        
    def make_element(self, parent, name):
        name, meta_dict = self._get_name(name)
        
        return ET.SubElement(parent, name, meta_dict)
        
    def write_element(self, parent_et, info_dict):
        
        for key in sorted(info_dict.keys()):
            if key == '_name':
                pass
            else:
                key_name, meta_dict = self._get_name(key)
                    
                k = ET.SubElement(parent_et, key_name, meta_dict)
                k.text = info_dict[key]
    
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
        self.cfg_obj.Site = Dummy()
        self.cfg_obj.Site._name = 'Site'
        self.cfg_obj.Site.Project = self.cfg_obj.Project
        self.cfg_obj.Site.Survey = self.cfg_obj.Survey
        self.cfg_obj.Site.DateCollected = self.mt_obj.edi_object.Header.acqdate
        self.cfg_obj.Site.Id = self.mt_obj.station
        self.cfg_obj.Site.AcquiredBy = self.mt_obj.edi_object.Header.acqby
        self.cfg_obj.Site.Start = self.mt_obj.edi_object.Header.acqdate
        self.cfg_obj.Site.End = self.mt_obj.edi_object.Header.acqdate
        self.cfg_obj.Site.RunList = '1'
        self.cfg_obj.Site.DataQualityNotes = self.cfg_obj.DataQualityNotes
        self.cfg_obj.Site.DataQualityWarnings = self.cfg_obj.DataQualityWarnings
        self.cfg_obj.Site.Location = Dummy(**{'Latitude':'{0:.6f}'.format(self.mt_obj.lat),
                                              'Longitude':'{0:.6f}'.format(self.mt_obj.lon),
                                              'Elevation(units=meters)':'{0:.3f}'.format(self.mt_obj.elev),
                                              'Declination(epoch=1995)':self.cfg_obj.Declination,
                                              '_name':'Location(datum={0})'.format(self.cfg_obj.Datum)})
        
       
        # processing information
        self.cfg_obj.ProcessingInfo.RemoteInfo.Project = self.cfg_obj.Project       
        self.cfg_obj.ProcessingInfo.RemoteInfo.Survey = self.cfg_obj.Survey
        self.cfg_obj.ProcessingInfo.RemoteInfo.YearCollected = self.cfg_obj.Site.DateCollected
        
       
        # Field Notes
        self.cfg_obj.FieldNotes = Dummy()
        self.cfg_obj.FieldNotes._name = 'FieldNotes(run=1)'
        self.cfg_obj.FieldNotes.Instrument = self.cfg_obj.Instrument
        self.cfg_obj.FieldNotes.Electrode = self.cfg_obj.Electrode
        self.cfg_obj.FieldNotes.Magnetometer = self.cfg_obj.Magnetometer
        self.cfg_obj.FieldNotes.Magnetometer.HX = str(self.mt_obj.edi_object.Define_measurement.meas_hx.acqchan)
        self.cfg_obj.FieldNotes.Magnetometer.HY = str(self.mt_obj.edi_object.Define_measurement.meas_hy.acqchan)
        try:        
            self.cfg_obj.FieldNotes.Magnetometer.HZ = str(self.mt_obj.edi_object.Define_measurement.meas_hz.acqchan)
        except AttributeError:
            pass
        #TODO: need to fill in more information on dipoles and magnetometers        
        
        # Input Channels
        self.cfg_obj.InputChannels = Dummy()
        self.cfg_obj.InputChannels._name = 'InputChannels(units=m)(ref=site)'

        self.cfg_obj.InputChannels.magnetic_hx = Dummy()
        hx = '(name=hx)(z=0)(y={0:.1f})(x={1:.1f})(orientation={2:.1f})'.format(
              self.mt_obj.edi_object.Define_measurement.meas_hx.y,
              self.mt_obj.edi_object.Define_measurement.meas_hx.x,
              self.mt_obj.edi_object.Define_measurement.meas_hx.azm)
        self.cfg_obj.InputChannels.magnetic_hx._name = 'Magnetic'+hx 
        
        self.cfg_obj.InputChannels.magnetic_hy = Dummy()
        hy = '(name=hy)(z=0)(y={0:.1f})(x={1:.1f})(orientation={2:.1f})'.format(
              self.mt_obj.edi_object.Define_measurement.meas_hy.y,
              self.mt_obj.edi_object.Define_measurement.meas_hy.x,
              self.mt_obj.edi_object.Define_measurement.meas_hy.azm)
        self.cfg_obj.InputChannels.magnetic_hy._name = 'Magnetic'+hy
        
        # Output Channels
        self.cfg_obj.OutputChannels = Dummy()
        self.cfg_obj.OutputChannels._name = 'OutputChannels(units=m)(ref=site)'
        try:
            hz = '(name=hz)(z=0)(y={0:.1f})(x={1:.1f})(orientation={2:.1f})'.format(
                  self.mt_obj.edi_object.Define_measurement.meas_hz.y,
                  self.mt_obj.edi_object.Define_measurement.meas_hz.x,
                  self.mt_obj.edi_object.Define_measurement.meas_hz.azm)
            self.cfg_obj.OutputChannels.magnetic_hz = Dummy()
            self.cfg_obj.OutputChannels.magnetic_hz._name = 'Magnetic'+hz 
        except AttributeError:
            print 'No HZ Information'
            
        ex = '(name=ex)(z=0)(y={0:.1f})(x={1:.1f})(y2={2:.1f})(x2={3:.1f})'.format(
                  self.mt_obj.edi_object.Define_measurement.meas_ex.y,
                  self.mt_obj.edi_object.Define_measurement.meas_ex.x,
                  self.mt_obj.edi_object.Define_measurement.meas_ex.y2,
                  self.mt_obj.edi_object.Define_measurement.meas_ex.x2)
        self.cfg_obj.OutputChannels.electric_ex = Dummy()
        self.cfg_obj.OutputChannels.electric_ex._name = 'Electric'+ex
            
        ey = '(name=ey)(z=0)(y={0:.1f})(x={1:.1f})(y2={2:.1f})(x2={3:.1f})'.format(
                  self.mt_obj.edi_object.Define_measurement.meas_ey.y,
                  self.mt_obj.edi_object.Define_measurement.meas_ey.x,
                  self.mt_obj.edi_object.Define_measurement.meas_ey.y2,
                  self.mt_obj.edi_object.Define_measurement.meas_ey.x2)
        self.cfg_obj.OutputChannels.electric_ey = Dummy()
        self.cfg_obj.OutputChannels.electric_ey._name = 'Electric'+ey
            
            
        self.cfg_obj.PeriodRange = Dummy()
        self.cfg_obj.PeriodRange._name = 'PeriodRange(min={0:.9f})(max={1:.9f})'.format(
                                         (1./self.mt_obj.Z.freq.max()),
                                         (1./self.mt_obj.Z.freq.min()))
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
        header_dict['z'] = 'Z(units=[mV/km]/[nT])(type=complex)(size=2 2)'
        header_dict['z.var'] = 'Z.VAR(type=real)(size=2 2)'
        header_dict['t'] = 'T(units=[])(type=complex)(size=1 2)'
        header_dict['t.var'] = 'T.VAR(type=real)(size=1 2)'
        
        attr_dict = {}
        attr_dict['z'] = 'z'
        attr_dict['z.var'] = 'z_err'
        attr_dict['t'] = 'tipper'
        attr_dict['t.var'] = 'tipper_err'

        nf = self.mt_obj.Z.freq.size        

        # make the data element
        self.cfg_obj.Data = Dummy(**{'_name':'Data(count={0})'.format(nf)})  
        
        # loop through each period and add appropriate information
        for f_index in range(nf):
            # set attribute period name with the index value
            # we are setting _name to have the necessary information so
            # we can name the attribute whatever we want.
            setattr(self.cfg_obj.Data, 
                    'Period_{0:02}'.format(f_index), 
                    Dummy(**{'_name':'Period(units=sec)(value={0:.6g})'.format(1./self.mt_obj.Z.freq[f_index])}))
            
            # Get information from data
            for estimate in ['z', 'z.var', 't', 't.var']:
                value_dict = {'_name':header_dict[estimate]}
                if 'z' in estimate:
                    for e_index in range(2):
                        for h_index in range(2):
                            c = comp_dict_z[(e_index, h_index)]
                            key_name = 'value(name={0})(input={1})(output={2})'.format(c[0], c[1], c[2])
                            if estimate == 'z':
                                z_value = getattr(self.mt_obj.Z, attr_dict[estimate])[f_index, e_index, h_index]                            
                                key_value = '{0:<+.8e} {1:<+.8e}'.format(z_value.real, 
                                                                         z_value.imag)
                            elif estimate == 'z.var':
                                z_value = getattr(self.mt_obj.Z, attr_dict[estimate])[f_index, e_index, h_index]                            
                                key_value = '{0:<+.8e}'.format(z_value)
                                                                      
                            value_dict[key_name] = key_value
                if 't' in estimate and self.mt_obj.Tipper.tipper is not None:
                    for e_index in range(1):
                        for h_index in range(2):
                            c = comp_dict_t[(e_index, h_index)]
                            key_name = 'value(name={0})(input={1})(output={2})'.format(c[0], c[1], c[2])
                            if estimate == 't':
                                z_value = getattr(self.mt_obj.Tipper, attr_dict[estimate])[f_index, e_index, h_index]                            
                                key_value = '{0:<+.8e} {1:<+.8e}'.format(z_value.real, 
                                                                         z_value.imag)
                            elif estimate == 't.var':
                                z_value = getattr(self.mt_obj.Tipper, attr_dict[estimate])[f_index, e_index, h_index]                            
                                key_value = '{0:<+.8e}'.format(z_value)
                                                                      
                            value_dict[key_name] = key_value
    
                # set the period attribute to have attributes for each
                # components of Z            
                setattr(getattr(self.cfg_obj.Data, 
                                'Period_{0:02}'.format(f_index)), 
                        estimate.capitalize(),
                        Dummy(**value_dict))
            
        
        
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
        
        self.get_info()

        # make the top of the tree element        
        emtf = ET.Element('EM_TF')
        
        # loop over the important information sections
        for element in self._order_list:
            # get the information for the given element
            value = getattr(self.cfg_obj, element)
            
            # check if it is just a value
            if type(value) in [float, int, str] or value is None:
                self.write_element(emtf, {element:value})
                
            # if its a class Dummy, then check for single values or another
            # class Dummy, probably a better way to code this with while loops
            elif isinstance(value, Dummy):
                print value._name
                # make a new tree limb
                new_element = self.make_element(emtf, value._name)
                
                # loop over attributes within the Dummy class, skipping _name
                for attr in sorted(value.__dict__.keys()):
                    if '_name' in attr:
                        pass
                    else:
                        attr_value = getattr(value, attr)
                        if isinstance(attr_value, Dummy):
                            # make a new tree limb
                            new_element_02 = self.make_element(new_element, 
                                                              attr_value._name)
                            
                            # loop over attributes within the Dummy class, skipping _name
                            for attr_02 in sorted(attr_value.__dict__.keys()):
                                if '_name' in attr_02:
                                    pass
                                else:
                                    attr_02_value = getattr(attr_value, attr_02)
                                    if isinstance(attr_02_value, Dummy):
                                        new_element_03 = self.make_element(new_element_02,
                                                                           attr_02_value._name)
                                        # loop over attributes within the Dummy class, skipping _name
                                        for attr_03 in sorted(attr_02_value.__dict__.keys()):
                                            if '_name' in attr_03:
                                                pass
                                            else:
                                                attr_03_value = getattr(attr_02_value, attr_03)
                                                if isinstance(attr_03_value, Dummy):
                                                    new_element_04 = self.make_element(new_element_03, 
                                                                                       attr_03_value._name)
                                                    # loop over attributes within the Dummy class, skipping _name
                                                    for attr_04 in sorted(attr_03_value.__dict__.keys()):
                                                        if '_name' in attr_04:
                                                            pass
                                                        else:
                                                            attr_04_value = getattr(attr_03_value, attr_04)
                                                            if isinstance(attr_04_value, Dummy):
                                                                new_element_05 = self.make_element(new_element_04, 
                                                                                                   attr_04_value._name)
                                                                for attr_05 in sorted(attr_04_value.__dict__.keys()):
                                                                    if '_name' in attr_05:
                                                                        pass
                                                                    else:
                                                                        attr_05_value = getattr(attr_04_value, attr_05)
                                                                        self.write_element(new_element_05,
                                                                                           {attr_05:attr_05_value})
                                                            else:
                                                                self.write_element(new_element_04, 
                                                                                   {attr_04:attr_04_value})
                                                else:
                                                    self.write_element(new_element_03,
                                                                       {attr_03:attr_03_value})
                                        
                                    elif type(attr_02_value) in [float, int, str] or attr_02_value is None:
                                        self.write_element(new_element_02, 
                                                           {attr_02:attr_02_value})

                        elif type(attr_value) in [float, int, str] or attr_value is None:
                            self.write_element(new_element, {attr:attr_value})
            
            # write out estimates
            elif type(value) in [list, tuple]:
                new_element = self.make_element(emtf, element)
                for element_ii in value:
                    newer_element = self.make_element(new_element, 
                                                      element_ii._name)
                    self.write_element(newer_element, element_ii.__dict__)
  
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
#        self.xml_fn = xml_fn
#        
#        self.cfg_obj = XML_Config()
#        
##==============================================================================
## Do the dirty work
##==============================================================================
#test = EDI_to_XML()
#test.edi_fn = r"c:\Users\jpeacock\Documents\Geothermal\Washington\MSHS\EDI_Files_birrp_mshs\Rotated_m18_deg\ms11.edi"
#test.cfg_fn = r"C:\Users\jpeacock\Documents\PyScripts\xml_cfg_test.cfg"
#test.write_xml()
#
#test.cfg_obj.write_cfg_file(r"C:\Users\jpeacock\Documents\PyScripts\xml_cfg_test_out.cfg")


    

      