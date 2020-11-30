# -*- coding: utf-8 -*-
"""
==================
metadata
==================

This module deals with metadata as defined by the MT metadata standards.
`metadata documentation 
<https://github.com/kujaku11/MTarchive/blob/tables/docs/mt_metadata_guide.pdf>`_.

There are multiple containers for each type of metadata, named appropriately.

Each container will be able to read and write:
    * dictionary
    * json
    * xml
    * csv?
    * pandas.Series
    * anything else?

Because a lot of the name words in the metadata are split by '.' there are some
issues we need to deal with.  I wrote in get and set attribute functions
to handle these types of names so the user shouldn't have to work about
splitting the names themselves.

These containers will be the building blocks for the metadata and how they are
interchanged between the HDF5 file and the user.  A lot of the metadata you
can get directly from the raw time series files, but the user will need to
input a decent amount on their own.  Dictionaries are the most fundamental
type we should be dealing with.

Each container has an attribute called _attr_dict which dictates if the
attribute is included in output objects, the data type, whether it is a
required parameter, and the style of output.  This should help down the road
with validation and keeping the data types consistent.  And if things change
you should only have to changes these dictionaries.

self._attr_dict = {'nameword':{'type': str, 'required': True, 'style': 'name'}}

Created on Sun Apr 24 20:50:41 2020

:copyright:
    Jared Peacock (jpeacock@usgs.gov)
    
:license: 
    MIT


"""
# =============================================================================
# Imports
# =============================================================================
import json
import pandas as pd
import numpy as np
import logging
import textwrap

from collections import OrderedDict
from collections.abc import Iterable
from operator import itemgetter

from mtpy.core.standards.schema import (
    Standards,
    validate_attribute,
    validate_type,
    MTSchemaError,
)
from mtpy.utils.mttime import MTime, MTTimeError
from mtpy.core.standards import helpers

ATTR_DICT = Standards().ATTR_DICT

# =============================================================================
# write doc strings
# =============================================================================
def wrap_description(description, column_width):
    """
    split a description into separate lines
    """
    d_lines = textwrap.wrap(description, column_width)
    if len(d_lines) < 9:
        d_lines += [""] * (9 - len(d_lines))

    return d_lines


def write_lines(attr_dict, c1=45, c2=45, c3=15):
    """
     write table lines
    :param lines_list: DESCRIPTION

    """
    line = "       | {0:<{1}}| {2:<{3}} | {4:<{5}}|"
    hline = "       +{0}+{1}+{2}+".format(
        "-" * (c1 + 1), "-" * (c2 + 2), "-" * (c3 + 1)
    )
    mline = "       +{0}+{1}+{2}+".format(
        "=" * (c1 + 1), "=" * (c2 + 2), "=" * (c3 + 1)
    )

    lines = [
        hline,
        line.format("**Metadata Key**", c1, "**Description**", c2, "**Example**", c3),
        mline,
    ]

    for key, entry in attr_dict.items():
        d_lines = wrap_description(entry["description"], c2)
        e_lines = wrap_description(entry["example"], c3)
        # line 1 is with the entry
        lines.append(line.format(f"**{key}**", c1, d_lines[0], c2, e_lines[0], c3))
        # line 2 skip an entry in the
        lines.append(line.format("", c1, d_lines[1], c2, e_lines[1], c3))
        # line 3 required
        lines.append(
            line.format(
                f"Required: {entry['required']}", c1, d_lines[2], c2, e_lines[2], c3
            )
        )
        # line 4 blank
        lines.append(line.format("", c1, d_lines[3], c2, e_lines[3], c3))

        # line 5 units
        lines.append(
            line.format(f"Units: {entry['units']}", c1, d_lines[4], c2, e_lines[4], c3)
        )

        # line 6 blank
        lines.append(line.format("", c1, d_lines[5], c2, e_lines[5], c3))

        # line 7 type
        lines.append(
            line.format(f"Type: {entry['type']}", c1, d_lines[6], c2, e_lines[6], c3)
        )

        # line 8 blank
        lines.append(line.format("", c1, d_lines[7], c2, e_lines[7], c3))

        # line 9 type
        lines.append(
            line.format(f"Style: {entry['style']}", c1, d_lines[8], c2, e_lines[8], c3)
        )

        # line 10 blank
        if len(d_lines) > 9:
            lines.append(line.format("", c1, d_lines[9], c2, "", c3))
            for d_line in d_lines[10:]:
                lines.append(line.format("", c1, d_line, c2, "", c3))

        lines.append(hline)

    return "\n".join(lines)


# =============================================================================
#  Base class that everything else will inherit
# =============================================================================
class Base:
    """
    A Base class that is common to most of the Metadata objects

    Includes:
        
        * to_json
        * from_json
        * to_dict
        * from_dict
        * to_series
        * from_series
        
    """

    def __init__(self, attr_dict={}, **kwargs):

        self._attr_dict = attr_dict

        self._class_name = validate_attribute(self.__class__.__name__)

        self.logger = logging.getLogger(f"{__name__}.{self._class_name}")

        for name, value in kwargs.items():
            self.set_attr_from_name(name, value)

    def __str__(self):
        meta_dict = self.to_dict()[self._class_name.lower()]
        lines = ["{0}:".format(self._class_name)]
        for name, value in meta_dict.items():
            lines.append("\t{0} = {1}".format(name, value))
        return "\n".join(lines)

    def __repr__(self):
        return self.to_json()

    def __eq__(self, other):
        # if other is None:
        #     self.logger.debug(f"Input is None, cannot compare with {self._class_name}")
        #     return False
        if isinstance(other, (type(self), Base, dict, str, pd.Series)):
            home_dict = self.to_dict()[self._class_name]
            if isinstance(other, Base):
                other_dict = other.to_dict()[self._class_name]
            elif isinstance(other, dict):
                other_dict = other
            elif isinstance(other, str):
                other_dict = OrderedDict(
                    sorted(json.loads(other).items(), key=itemgetter(0))
                )
            elif isinstance(other, pd.Series):
                other_dict = OrderedDict(
                    sorted(other.to_dict().items(), key=itemgetter(0))
                )
            if other_dict == home_dict:
                return True
            else:
                for key, value in home_dict.items():
                    try:
                        other_value = other_dict[key]
                        if value != other_value:
                            msg = f"{key}: {value} != {other_value}"
                            self.logger.info(msg)
                    except KeyError:
                        msg = "Cannot find {0} in other".format(key)
                        self.logger.info(msg)

                return False
        raise ValueError(f"Cannot compare {self._class_name} with {type(other)}")

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.get_attribute_list())

    def get_attribute_list(self):
        """
        return a list of the attributes 
        """

        return sorted(list(self._attr_dict.keys()))

    def attribute_information(self, name=None):
        """
        return a descriptive string of the attribute if none returns for all
    
        :param key: DESCRIPTION, defaults to None
        :type key: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if name:
            try:
                v_dict = OrderedDict(
                    sorted(self._attr_dict[name].items(), key=itemgetter(0))
                )
            except KeyError as error:
                msg = "{0} not attribute {1} found".format(error, name)
                self.logger.error(msg)
                raise MTSchemaError(msg)

            lines = ["{0}:".format(name)]
            for key, value in v_dict.items():
                lines.append("\t{0}: {1}".format(key, value))
        else:
            lines = []
            for name, v_dict in self._attr_dict.items():
                lines.append("{0}:".format(name))
                v_dict = OrderedDict(sorted(v_dict.items(), key=itemgetter(0)))
                for key, value in v_dict.items():
                    lines.append("\t{0}: {1}".format(key, value))
                lines.append("=" * 50)

        print("\n".join(lines))

    def _validate_name(self, name):
        """
        validate the name to conform to the standards
        name must be:
            * all lower case {a-z; 1-9}
            * must start with a letter
            * categories are separated by '.'
            * words separated by '_'

        {object}.{name_name}

        '/' will be replaced with '.'
        converted to all lower case

        :param name: name name
        :type name: string
        :return: valid name name
        :rtype: string

        """
        return validate_attribute(name)

    def _validate_type(self, value, v_type, style=None):
        """
        validate type from standards
        
        """
        # need to have this if the user sets a value with a class, there is not
        # a good way to validate the class object, but all elements within the
        # class object will be validated, so it seems fine to skip it.
        if isinstance(
            value,
            (
                Declination,
                Location,
                Fdsn,
                Rating,
                DataQuality,
                Citation,
                Copyright,
                Provenance,
                Person,
                Diagnostic,
                Battery,
                Electrode,
                TimingSystem,
                TimePeriod,
                Orientation,
                Software,
                Filtered,
                Filter,
                DataLogger,
                TransferFunction,
                Survey,
                Station,
                Run,
                Channel,
                Auxiliary,
                Electric,
                Magnetic,
            ),
        ):
            return value
        # return if the value is None, this may need to change in the future
        # if an empty list or something else should be returned
        if not isinstance(value, (list, tuple, np.ndarray)):
            if value in [None, "None", "none", "unknown"]:
                return None
        # hack to get around h5py reference types, in the future will need
        # a more robust test.
        if v_type == "h5py_reference":
            return value

        # return value if the value type is not defined.
        if v_type is None:
            msg = (
                "standards data type is unknown, if you want to "
                + "propogate this attribute using to_dict, to_json or "
                + "to_series, you need to add attribute description using "
                + "class function add_base_attribute."
                + "Example: \n\t>>> Run.add_base_attribute(new, 10, "
                + '{"type":float, "required": True, "units": None, '
                + '"style": number})'
            )
            self.logger.info(msg)
            return value

        # if not a python type but a string organize into a dictionary
        if not isinstance(v_type, type) and isinstance(v_type, str):
            type_dict = {"string": str, "integer": int, "float": float, "boolean": bool}
            v_type = type_dict[validate_type(v_type)]
        else:
            msg = "v_type must be a string or type not {0}".format(v_type)

        # check style for a list
        if isinstance(value, v_type):
            if style:
                if v_type is str and "list" in style:
                    value = value.replace("[", "").replace("]", "").split(",")
                    value = [ss.strip() for ss in value]
            return value

        # if value is not of v_type
        else:
            msg = "value={0} must be {1} not {2}"
            info = "converting {0} to {1}"
            # if the value is a string, convert to appropriate type
            if isinstance(value, str):
                if v_type is int:
                    try:
                        self.logger.debug(info.format(type(value), v_type))
                        return int(value)
                    except ValueError as error:
                        self.logger.exception(error)
                        raise MTSchemaError(msg.format(value, v_type, type(value)))
                elif v_type is float:
                    try:
                        self.logger.debug(info.format(type(value), v_type))
                        return float(value)
                    except ValueError as error:
                        self.logger.exception(error)
                        raise MTSchemaError(msg.format(value, v_type, type(value)))
                elif v_type is bool:
                    if value.lower() in ["false", "0"]:
                        self.logger.debug(info.format(value, False))
                        return False
                    elif value.lower() in ["true", "1"]:
                        self.logger.debug(info.format(value, True))
                        return True
                    else:
                        self.logger.exception(msg.format(value, v_type, type(value)))
                        raise MTSchemaError(msg.format(value, v_type, type(value)))
                elif v_type is str:
                    return value

            # if a number convert to appropriate type
            elif isinstance(value, (int, np.int_)):
                if v_type is float:
                    self.logger.debug(info.format(type(value), v_type))
                    return float(value)
                elif v_type is str:
                    self.logger.debug(info.format(type(value), v_type))
                    return "{0:.0f}".format(value)
                return int(value)

            # if a number convert to appropriate type
            elif isinstance(value, (float, np.float_)):
                if v_type is int:
                    self.logger.debug(info.format(type(value), v_type))
                    return int(value)
                elif v_type is str:
                    self.logger.debug(info.format(type(value), v_type))
                    return "{0}".format(value)
                return float(value)

            # if a list convert to appropriate entries to given type
            elif isinstance(value, Iterable):
                if v_type is str:
                    if isinstance(value, np.ndarray):
                        value = value.astype(np.unicode_)
                    value = [f"{v}" for v in value]
                elif v_type is int:
                    value = [int(float(v)) for v in value]
                elif v_type is float:
                    value = [float(v) for v in value]
                elif v_type is bool:
                    value_list = []
                    for v in value:
                        if v in [True, "true", "True", "TRUE"]:
                            value_list.append(True)
                        elif v in [False, "false", "False", "FALSE"]:
                            value_list.append(False)
                    value = value_list
                return value

            else:
                self.logger.exception(msg.format(value, v_type, type(value)))
                raise MTSchemaError(msg.format(value, v_type, type(value)))
        return None

    def _validate_option(self, name, option_list):
        """
        validate the given attribute name agains possible options and check
        for aliases
        
        :param name: DESCRIPTION
        :type name: TYPE
        :param option_list: DESCRIPTION
        :type option_list: TYPE
        :param alias_list: DESCRIPTION
        :type alias_list: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if name is None:
            return True, False, None

        options = [ss.lower() for ss in option_list]
        other_possible = False
        if "other" in options:
            other_possible = True
        if name.lower() in options:
            return True, other_possible, None
        elif name.lower() not in options and other_possible:
            msg = (
                "{0} not found in options list {1}, but other options"
                + " are allowed.  Allowing {2} to be set to {0}."
            )
            return True, other_possible, msg

        return False, other_possible, "{0} not found in options list {1}"

    def __setattr__(self, name, value):
        """
        set attribute based on metadata standards

        """
        # skip these attribute because they are validated in the property
        # setter.
        skip_list = [
            "latitude",
            "longitude",
            "elevation",
            "start_date",
            "end_date",
            "start",
            "end",
            "name",
            "applied",
            "logger",
            "Electric",
        ]

        if hasattr(self, "_attr_dict"):
            if name[0] != "_":
                if not name in skip_list:
                    self.logger.debug("Setting {0} to {1}".format(name, value))
                    v_dict = self._attr_dict[name]
                    v_type = self._get_standard_type(name)
                    try:
                        value = self._validate_type(value, v_type, v_dict["style"])
                    except (MTSchemaError, ValueError) as error:
                        msg = f"{name} failed because {error}"
                        self.logger.error(msg)
                        raise MTSchemaError(msg)
                    # check options
                    if v_dict["style"] == "controlled vocabulary":
                        options = v_dict["options"]
                        accept, other, msg = self._validate_option(value, options)
                        if not accept:
                            self.logger.error(msg.format(value, options))
                            raise MTSchemaError(msg.format(value, options))
                        if other and not accept:
                            self.logger.warning(msg.format(value, options, name))

        super().__setattr__(name, value)

    def _get_standard_type(self, name):
        """
        helper function to get the standard type for the given name
        """
        name = self._validate_name(name)
        try:
            standards = self._attr_dict[name]
            return standards["type"]
        except KeyError:
            if name[0] != "_":
                msg = (
                    "{0} is not defined in the standards. "
                    + " Should add attribute information with "
                    + "add_base_attribute if the attribute is going to "
                    + "propogate via to_dict, to_json, to_series"
                )
                self.logger.info(msg.format(name))
            return None

    def get_attr_from_name(self, name):
        """
        Access attribute from the given name.

        The name can contain the name of an object which must be separated
        by a '.' for  e.g. {object_name}.{name} --> location.latitude

        ..note:: this is a helper function for names with '.' in the name for
                 easier getting when reading from dictionary.

        :param name: name of attribute to get.
        :type name: string
        :return: attribute value
        :rtype: type is defined by the attribute name

        :Example:

        >>> b = Base(**{'category.test_attr':10})
        >>> b.get_attr_from_name('category.test_attr')
        10

        """
        name = self._validate_name(name)
        v_type = self._get_standard_type(name)

        if "." in name:
            value = helpers.recursive_split_getattr(self, name)
        else:
            value = getattr(self, name)

        try:
            return self._validate_type(value, v_type)
        except MTSchemaError as error:
            msg = f"cannot retrieve {name} because {error}"
            raise MTSchemaError(msg)

    def set_attr_from_name(self, name, value):
        """
        Helper function to set attribute from the given name.

        The name can contain the name of an object which must be separated
        by a '.' for  e.g. {object_name}.{name} --> location.latitude

        ..note:: this is a helper function for names with '.' in the name for
                 easier getting when reading from dictionary.

        :param name: name of attribute to get.
        :type name: string
        :param value: attribute value
        :type value: type is defined by the attribute name

        :Example: 

        >>> b = Base(**{'category.test_attr':10})
        >>> b.set_attr_from_name('category.test_attr', '10')
        >>> print(b.category.test_attr)
        '10'
        """
        if "." in name:
            try:
                helpers.recursive_split_setattr(self, name, value)
            except AttributeError as error:
                msg = (
                    "{0} is not in the current standards.  "
                    + "To properly add the attribute use "
                    + "add_base_attribute."
                )

                self.logger.error(msg.format(name))
                self.logger.exception(error)
                raise AttributeError(error)
        else:
            setattr(self, name, value)

    def add_base_attribute(self, name, value, value_dict):
        """
        Add an attribute to _attr_dict so it will be included in the
        output dictionary
        
        :param name: name of attribute
        :type name: string
        
        :param value: value of the new attribute
        :type value: described in value_dict
        
        :param value_dict: dictionary describing the attribute, must have keys
            ['type', 'required', 'style', 'units', 'alias', 'description',
             'options', 'example']
        :type name: string
    
        * type --> the data type [ str | int | float | bool ]
        * required --> required in the standards [ True | False ]
        * style --> style of the string
        * units --> units of the attribute, must be a string
        * alias --> other possible names for the attribute
        * options --> if only a few options are accepted, separated by | or 
          comma.b [ option_01 | option_02 | other ]. 'other' means other options 
          available but not yet defined.
        * example --> an example of the attribute
        
        :Example:
            
        >>> extra = {'type': str,
        >>> ...      'style': 'controlled vocabulary',
        >>> ...      'required': False,
        >>> ...      'units': celsius,
        >>> ...      'description': 'local temperature',
        >>> ...      'alias': ['temp'],
        >>> ...      'options': [ 'ambient', 'air', 'other'],
        >>> ...      'example': 'ambient'}
        >>> r = Run()
        >>> r.add_base_attribute('temperature', 'ambient', extra)

        """
        name = self._validate_name(name)
        self._attr_dict.update({name: value_dict})
        self.set_attr_from_name(name, value)
        self.logger.debug("Added {0} to _attr_dict with {1}".format(name, value_dict))
        self.logger.debug(
            "set {0} to {1} as type {2}".format(name, value, value_dict["type"])
        )

    def to_dict(self, nested=False, single=False, required=True):
        """
        make a dictionary from attributes, makes dictionary from _attr_list.
        
        :param nested: make the returned dictionary nested
        :type nested: [ True | False ] , default is False
        
        :param single: return just metadata dictionary -> meta_dict[class_name]
        :type single: [ True | False ], default is False
        
        :param required: return just the required elements and any elements with
                         non-None values
        
        """
        meta_dict = {}
        for name in list(self._attr_dict.keys()):
            try:
                value = self.get_attr_from_name(name)
            except AttributeError as error:
                msg = "{0}: setting {1} to None.  ".format(
                    error, name
                ) + "Try setting {0} to the desired value".format(name)
                self.logger.debug(msg)
                value = None

            if required:
                if value is not None or self._attr_dict[name]["required"]:
                    meta_dict[name] = value
            else:
                meta_dict[name] = value

        if nested:
            meta_dict = helpers.structure_dict(meta_dict)

        meta_dict = {
            self._class_name.lower(): OrderedDict(
                sorted(meta_dict.items(), key=itemgetter(0))
            )
        }

        if single:
            meta_dict = meta_dict[list(meta_dict.keys())[0]]

        return meta_dict

    def from_dict(self, meta_dict):
        """
        fill attributes from a dictionary
        
        :param meta_dict: dictionary with keys equal to metadata.
        :type meta_dict: dictionary
        
        """
        if not isinstance(meta_dict, (dict, OrderedDict)):
            msg = "Input must be a dictionary not {0}".format(type(meta_dict))
            self.logger.error(msg)
            raise MTSchemaError(msg)

        class_name = list(meta_dict.keys())[0]
        if class_name.lower() != self._class_name.lower():
            msg = (
                "name of input dictionary is not the same as class type "
                "input = {0}, class type = {1}".format(class_name, self._class_name)
            )
            self.logger.debug(msg)

        # be sure to flatten the dictionary first for easier transform
        meta_dict = helpers.flatten_dict(meta_dict[class_name])
        for name, value in meta_dict.items():
            self.set_attr_from_name(name, value)

    def to_json(self, nested=False, indent=" " * 4, required=True):
        """
        Write a json string from a given object, taking into account other
        class objects contained within the given object.
        
        :param nested: make the returned json nested
        :type nested: [ True | False ] , default is False
        
        """

        return json.dumps(
            self.to_dict(nested=nested, required=required),
            cls=helpers.NumpyEncoder,
            indent=indent,
        )

    def from_json(self, json_str):
        """
        read in a json string and update attributes of an object

        :param json_str: json string
        :type json_str: string

        """
        if not isinstance(json_str, str):
            msg = "Input must be valid JSON string not {0}".format(type(json_str))
            self.logger.error(msg)
            raise MTSchemaError(msg)

        self.from_dict(json.loads(json_str))

    def from_series(self, pd_series):
        """
        Fill attributes from a Pandas series
        
        .. note:: Currently, the series must be single layered with key names
                  separated by dots. (location.latitude)
        
        :param pd_series: Series containing metadata information
        :type pd_series: pandas.Series
        
        ..todo:: Force types in series
        
        """
        if not isinstance(pd_series, pd.Series):
            msg = "Input must be a Pandas.Series not type {0}".format(type(pd_series))
            self.logger.error(msg)
            MTSchemaError(msg)
        for key, value in pd_series.iteritems():
            self.set_attr_from_name(key, value)

    def to_series(self, required=True):
        """
        Convert attribute list to a pandas.Series
        
        .. note:: this is a flattened version of the metadata
        
        :return: pandas.Series
        :rtype: pandas.Series

        """

        return pd.Series(self.to_dict(single=True, required=required))

    def to_xml(self, string=False, required=True):
        """
        make an xml element for the attribute that will add types and 
        units.  
        
        :param string: output a string instead of an XML element
        :type string: [ True | False ], default is False
        
        :return: XML element or string

        """
        element = helpers.dict_to_xml(
            self.to_dict(nested=True, required=required), self._attr_dict
        )
        if not string:
            return element
        else:
            return helpers.element_to_string(element)

    def from_xml(self, xml_element):
        """
        
        :param xml_element: XML element
        :type xml_element: etree.Element
        
        :return: Fills attributes accordingly

        """

        self.from_dict(helpers.element_to_dict(xml_element))


# ============================================================================
# Location class, be sure to put locations in decimal degrees, and note datum
# ============================================================================
class Declination(Base):
    __doc__ = write_lines(ATTR_DICT["declination"])

    def __init__(self, **kwargs):

        self.value = None
        self.epoch = None
        self.model = None
        self.comments = None
        super(Declination, self).__init__(attr_dict=ATTR_DICT["declination"], **kwargs)


class Location(Base):
    __doc__ = write_lines(ATTR_DICT["location"])

    def __init__(self, **kwargs):

        self.comments = None
        self.datum = "WGS84"
        self.declination = Declination()

        self._elevation = 0.0
        self._latitude = 0.0
        self._longitude = 0.0
        self.x = 0.0
        self.y = 0.0
        self.z = 0
        self.x2 = 0.0
        self.y2 = 0.0
        super(Location, self).__init__(attr_dict=ATTR_DICT["location"], **kwargs)

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, lat):
        self._latitude = self._assert_lat_value(lat)

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, lon):
        self._longitude = self._assert_lon_value(lon)

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, elev):
        self._elevation = self._assert_elevation_value(elev)

    def _assert_lat_value(self, latitude):
        """
        Make sure the latitude value is in decimal degrees, if not change it.
        And that the latitude is within -90 < lat > 90.

        :param latitude: latitude in decimal degrees or other format
        :type latitude: float or string
        """
        if latitude in [None, "None", "none", "unknown"]:
            self.logger.debug("Latitude is None, setting to 0")
            return 0.0
        try:
            lat_value = float(latitude)

        except TypeError:
            self.logger.debug("Could not convert {0} setting to 0".format(latitude))
            return 0.0

        except ValueError:
            self.logger.debug("Latitude is a string {0}".format(latitude))
            lat_value = self._convert_position_str2float(latitude)

        if abs(lat_value) >= 90:
            msg = (
                "latitude value = {0} is unacceptable!".format(lat_value)
                + ".  Must be |Latitude| > 90"
            )
            self.logger.error(msg)
            raise ValueError(msg)

        return lat_value

    def _assert_lon_value(self, longitude):
        """
        Make sure the longitude value is in decimal degrees, if not change it.
        And that the latitude is within -180 < lat > 180.

        :param latitude: longitude in decimal degrees or other format
        :type latitude: float or string
        """
        if longitude in [None, "None", "none", "unknown"]:
            self.logger.debug("Longitude is None, setting to 0")
            return 0.0
        try:
            lon_value = float(longitude)

        except TypeError:
            self.logger.debug("Could not convert {0} setting to 0".format(longitude))
            return 0.0

        except ValueError:
            self.logger.debug("Longitude is a string {0}".format(longitude))
            lon_value = self._convert_position_str2float(longitude)

        if abs(lon_value) >= 180:
            msg = (
                "longitude value = {0} is unacceptable!".format(lon_value)
                + ".  Must be |longitude| > 180"
            )
            self.logger.error(msg)
            raise ValueError(msg)

        return lon_value

    def _assert_elevation_value(self, elevation):
        """
        make sure elevation is a floating point number

        :param elevation: elevation as a float or string that can convert
        :type elevation: float or str
        """

        try:
            elev_value = float(elevation)
        except (ValueError, TypeError):
            msg = "Could not convert {0} to a number setting to 0".format(elevation)
            self.logger.debug(msg)
            elev_value = 0.0

        return elev_value

    def _convert_position_float2str(self, position):
        """
        Convert position float to a string in the format of DD:MM:SS.

        :param position: decimal degrees of latitude or longitude
        :type position: float

        :returns: latitude or longitude in format of DD:MM:SS.ms
        """

        assert type(position) is float, "Given value is not a float"

        deg = int(position)
        sign = 1
        if deg < 0:
            sign = -1

        deg = abs(deg)
        minutes = (abs(position) - deg) * 60.0
        # need to round seconds to 4 decimal places otherwise machine precision
        # keeps the 60 second roll over and the string is incorrect.
        sec = np.round((minutes - int(minutes)) * 60.0, 4)
        if sec >= 60.0:
            minutes += 1
            sec = 0

        if int(minutes) == 60:
            deg += 1
            minutes = 0

        position_str = "{0}:{1:02.0f}:{2:05.2f}".format(
            sign * int(deg), int(minutes), sec
        )
        self.logger.debug("Converted {0} to {1}".format(position, position_str))

        return position_str

    def _convert_position_str2float(self, position_str):
        """
        Convert a position string in the format of DD:MM:SS to decimal degrees

        :param position: latitude or longitude om DD:MM:SS.ms
        :type position: float

        :returns: latitude or longitude as a float
        """

        if position_str in [None, "None"]:
            return None

        p_list = position_str.split(":")
        if len(p_list) != 3:
            msg = "{0} not correct format, should be DD:MM:SS".format(position_str)
            self.logger.error(msg)
            raise ValueError(msg)

        deg = float(p_list[0])
        minutes = self._assert_minutes(float(p_list[1]))
        sec = self._assert_seconds(float(p_list[2]))

        # get the sign of the position so that when all are added together the
        # position is in the correct place
        sign = 1
        if deg < 0:
            sign = -1

        position_value = sign * (abs(deg) + minutes / 60.0 + sec / 3600.0)

        self.logger.debug("Converted {0} to {1}".format(position_str, position_value))

        return position_value

    def _assert_minutes(self, minutes):
        if not 0 <= minutes < 60.0:
            msg = (
                "minutes should be 0 < > 60, currently {0:.0f}".format(minutes)
                + " conversion will account for non-uniform"
                + "timne. Be sure to check accuracy."
            )
            self.logger.warning(msg)

        return minutes

    def _assert_seconds(self, seconds):
        if not 0 <= seconds < 60.0:
            msg = (
                "seconds should be 0 < > 60, currently {0:.0f}".format(seconds)
                + " conversion will account for non-uniform"
                + "timne. Be sure to check accuracy."
            )
            self.logger.warning(msg)

        return seconds


# ==============================================================================
# Instrument
# ==============================================================================
class Instrument(Base):
    __doc__ = write_lines(ATTR_DICT["instrument"])

    def __init__(self, **kwargs):

        self.id = None
        self.manufacturer = None
        self.type = None
        self.model = None
        super().__init__(attr_dict=ATTR_DICT["instrument"], **kwargs)


# =============================================================================
# FDSN
# =============================================================================
class Fdsn(Base):
    __doc__ = write_lines(ATTR_DICT["fdsn"])

    def __init__(self, **kwargs):
        self.id = None
        self.network = None
        self.channel_code = None
        self.new_epoch = None

        super().__init__(attr_dict=ATTR_DICT["fdsn"], **kwargs)


# ==============================================================================
# Data Quality
# ==============================================================================
class Rating(Base):
    __doc__ = write_lines(ATTR_DICT["rating"])

    def __init__(self, **kwargs):
        self.author = None
        self.method = None
        self.value = 0.0

        super().__init__(attr_dict=ATTR_DICT["rating"], **kwargs)


class DataQuality(Base):
    __doc__ = write_lines(ATTR_DICT["data_quality"])

    def __init__(self, **kwargs):

        self.rating = Rating()
        self.warnings = None

        super().__init__(attr_dict=ATTR_DICT["data_quality"], **kwargs)


# ==============================================================================
# Citation
# ==============================================================================
class Citation(Base):
    __doc__ = write_lines(ATTR_DICT["citation"])

    def __init__(self, **kwargs):
        self.author = None
        self.title = None
        self.journal = None
        self.volume = None
        self.doi = None
        self.year = None
        super().__init__(attr_dict=ATTR_DICT["citation"], **kwargs)


# ==============================================================================
# Copyright
# ==============================================================================
class Copyright(Base):
    __doc__ = write_lines(ATTR_DICT["copyright"])

    def __init__(self, **kwargs):
        self.citation = Citation()
        self.conditions_of_use = "".join(
            [
                "All data and metadata for this survey are ",
                "available free of charge and may be copied ",
                "freely, duplicated and further distributed ",
                "provided this data set is cited as the ",
                "reference. While the author(s) strive to ",
                "provide data and metadata of best possible ",
                "quality, neither the author(s) of this data ",
                "set, not IRIS make any claims, promises, or ",
                "guarantees about the accuracy, completeness, ",
                "or adequacy of this information, and expressly ",
                "disclaim liability for errors and omissions in ",
                "the contents of this file. Guidelines about ",
                "the quality or limitations of the data and ",
                "metadata, as obtained from the author(s), are ",
                "included for informational purposes only.",
            ]
        )
        self.release_license = None
        self.comments = None
        super().__init__(attr_dict=ATTR_DICT["copyright"], **kwargs)


# ==============================================================================
# Provenance
# ==============================================================================
class Provenance(Base):
    __doc__ = write_lines(ATTR_DICT["provenance"])

    def __init__(self, **kwargs):

        self._creation_dt = MTime()
        self._creation_dt.now()
        self.creating_application = "MTH5"
        self.creator = Person()
        self.submitter = Person()
        self.software = Software()
        self.log = None
        self.comments = None
        super().__init__(attr_dict=ATTR_DICT["provenance"], **kwargs)

    @property
    def creation_time(self):
        return self._creation_dt.iso_str

    @creation_time.setter
    def creation_time(self, dt_str):
        self._creation_dt.from_str(dt_str)


# ==============================================================================
# Person
# ==============================================================================
class Person(Base):
    __doc__ = write_lines(ATTR_DICT["person"])

    def __init__(self, **kwargs):

        self.email = None
        self.author = None
        self.organization = None
        self.comments = None
        # self.url = None
        super().__init__(attr_dict=ATTR_DICT["person"], **kwargs)


# =============================================================================
# diagnostic
# =============================================================================
class Diagnostic(Base):
    __doc__ = write_lines(ATTR_DICT["diagnostic"])

    def __init__(self, **kwargs):
        self.units = None
        self.start = None
        self.end = None
        super().__init__(attr_dict=ATTR_DICT["diagnostic"], **kwargs)


# =============================================================================
# Battery
# =============================================================================
class Battery(Base):
    __doc__ = write_lines(ATTR_DICT["battery"])

    def __init__(self, **kwargs):

        self.type = None
        self.id = None
        self.voltage = Diagnostic()
        self.comments = None
        super().__init__(attr_dict=ATTR_DICT["battery"], **kwargs)


# =============================================================================
# Electrode
# =============================================================================
class Electrode(Base):
    __doc__ = write_lines(ATTR_DICT["electrode"])

    def __init__(self, **kwargs):

        self.id = None
        self.manufacturer = None
        self.type = None
        self.model = None
        self.latitude = 0.0
        self.longitude = 0.0
        self.elevation = 0.0
        self.x = 0.0
        self.x2 = 0.0
        self.y = 0.0
        self.y2 = 0.0
        super().__init__(attr_dict=ATTR_DICT["electrode"], **kwargs)


# =============================================================================
# Timing System
# =============================================================================
class TimingSystem(Base):
    __doc__ = write_lines(ATTR_DICT["timing_system"])

    def __init__(self, **kwargs):

        self.type = None
        self.drift = None
        self.drift_units = None
        self.uncertainty = None
        self.uncertainty_units = None
        self.comments = None
        super().__init__(attr_dict=ATTR_DICT["timing_system"], **kwargs)


class TimePeriod(Base):
    __doc__ = write_lines(ATTR_DICT["time_period"])

    def __init__(self, **kwargs):

        self._start_dt = MTime()
        self._end_dt = MTime()
        super().__init__(attr_dict=ATTR_DICT["time_period"], **kwargs)

    @property
    def start(self):
        return self._start_dt.iso_str

    @start.setter
    def start(self, start_date):
        self._start_dt.from_str(start_date)

    @property
    def end(self):
        return self._end_dt.iso_str

    @end.setter
    def end(self, stop_date):
        self._end_dt.from_str(stop_date)

    @property
    def start_date(self):
        return self._start_dt.date

    @start_date.setter
    def start_date(self, start_date):
        self._start_dt.from_str(start_date)

    @property
    def end_date(self):
        return self._end_dt.date

    @end_date.setter
    def end_date(self, stop_date):
        self._end_dt.from_str(stop_date)


class Orientation(Base):
    __doc__ = write_lines(ATTR_DICT["orientation"])

    def __init__(self, **kwargs):
        self.reference_frame = "geographic"
        self.method = None

        super().__init__(attr_dict=ATTR_DICT["orientation"], **kwargs)


# ==============================================================================
# Software
# ==============================================================================
class Software(Base):
    __doc__ = write_lines(ATTR_DICT["software"])

    def __init__(self, **kwargs):
        self.name = None
        self.version = None
        self._author = Person()

        super().__init__(attr_dict=ATTR_DICT["software"], **kwargs)

    @property
    def author(self):
        return self._author.author

    @author.setter
    def author(self, value):
        self._author.author = value


# =============================================================================
# filter
# =============================================================================
class Filtered(Base):
    __doc__ = write_lines(ATTR_DICT["filtered"])

    def __init__(self, **kwargs):
        self._name = []
        self._applied = []
        self.name = None
        self.applied = None
        self.comments = None
        super().__init__(attr_dict=ATTR_DICT["filtered"], **kwargs)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, names):
        if names is None:
            self._name = ["none"]
            return

        if isinstance(names, str):
            self._name = [ss.strip().lower() for ss in names.split(",")]
        elif isinstance(names, list):
            self._name = [ss.strip().lower() for ss in names]
        elif isinstance(names, np.ndarray):
            names = names.astype(np.unicode_)
            self._name = [ss.strip().lower() for ss in names]
        else:
            msg = "names must be a string or list of strings not {0}, type {1}"
            self.logger.error(msg.format(names, type(names)))
            raise MTSchemaError(msg.format(names, type(names)))

        check = self._check_consistency()
        if not check:
            self.logger.debug(
                "Filter names and applied lists are not the "
                + "same size. Be sure to check the inputs."
                + " names = {0}, applied = {1}".format(self._name, self._applied)
            )

    @property
    def applied(self):
        return self._applied

    @applied.setter
    def applied(self, applied):
        if not isinstance(applied, (list, tuple)):
            if applied in [None, "none", "None", "NONE", "null", 0, "0"]:
                self._applied = [False]
                return
        if applied == []:
            self.applied = [False]
            return

        if isinstance(applied, str):
            if applied.find("[") >= 0:
                applied = applied.replace("[", "").replace("]", "")

            if applied.count(",") > 0:
                applied_list = [ss.strip().lower() for ss in applied.split(",")]
            else:
                applied_list = [ss.lower() for ss in applied.split()]
        elif isinstance(applied, list):
            applied_list = applied
        elif isinstance(applied, bool):
            applied_list = [applied]
        # the returned type from a hdf5 dataset is a numpy array.
        elif isinstance(applied, np.ndarray):
            applied_list = list(applied)
            if applied_list == []:
                applied_list = [False]
        else:
            msg = "applied must be a string or list of strings not {0}"
            self.logger.error(msg.format(applied))
            raise MTSchemaError(msg.format(applied))

        bool_list = []
        for app_bool in applied_list:
            if app_bool is None:
                bool_list.append(False)
            if isinstance(app_bool, str):
                if app_bool.lower() in ["false", "0"]:
                    bool_list.append(False)
                elif app_bool.lower() in ["true", "1"]:
                    bool_list.append(True)
                else:
                    msg = "Filter.applied must be [ True | False ], not {0}"
                    self.logger.error(msg.format(app_bool))
                    raise MTSchemaError(msg.format(app_bool))
            elif isinstance(app_bool, bool):
                bool_list.append(app_bool)
            else:
                msg = "Filter.applied must be [True | False], not {0}"
                self.logger.error(msg.format(app_bool))
        self._applied = bool_list

        # check for consistency
        check = self._check_consistency()
        if not check:
            self.logger.debug(
                "Filter names and applied lists are not the "
                + "same size. Be sure to check the inputs."
                + ". name = {0}, applied = {1}".format(self._name, self._applied)
            )

    def _check_consistency(self):
        # check for consistency
        if self._name is not None:
            if self._applied is None:
                self.logger.warning("Need to input filter.applied")
                return False
            if len(self._name) == 1:
                if len(self._applied) == 1:
                    return True
            elif len(self._name) > 1:
                if len(self._applied) == 1:
                    self.logger.info(
                        "Assuming all filters have been "
                        + "applied as {0}".format(self._applied[0])
                    )
                    return True
                elif len(self._applied) > 1:
                    if len(self._applied) != len(self._name):
                        self.logger.waring(
                            "Applied and filter names "
                            + "should be the same length. "
                            + "Appied={0}, names={1}".format(
                                len(self._applied), len(self._name)
                            )
                        )
                        return False
        else:
            return False


class Filter(Base):
    __doc__ = write_lines(ATTR_DICT["filter"])

    def __init__(self, **kwargs):
        self.name = None
        self.type = None
        self.units_in = None
        self.units_out = None
        self._calibration_dt = MTime()
        self.operation = None
        self.normalization_frequency = None
        self.normalization_factor = None
        self.cutoff = None
        self.n_poles = None
        self.n_zeros = None
        self.comments = None
        self.conversion_factor = None

        super().__init__(attr_dict=ATTR_DICT["filter"], **kwargs)

    @property
    def calibration_date(self):
        return self._calibration_dt.date

    @calibration_date.setter
    def calibration_date(self, value):
        self._calibration_dt.from_str(value)


# =============================================================================
# Data logger
# =============================================================================
class DataLogger(Base):
    __doc__ = write_lines(ATTR_DICT["datalogger"])

    def __init__(self, **kwargs):
        self.id = None
        self.manufacturer = None
        self.type = None
        self.model = None
        self.timing_system = TimingSystem()
        self.firmware = Software()
        self.power_source = Battery()
        super().__init__(attr_dict=ATTR_DICT["datalogger"], **kwargs)


# =============================================================================
# transfer function
# =============================================================================
class TransferFunction(Base):
    __doc__ = write_lines(ATTR_DICT["transfer_function"])

    def __init__(self, **kwargs):

        self.processed_by = Person()
        self.software = Software()
        self.units = "millivolts_per_kilometer_per_nanotesla"
        self.sign_convention = "+"
        self.runs_processed = []
        self.remote_references = []
        self.processing_parameters = []
        self._processed_date = MTime()

        super().__init__(attr_dict=ATTR_DICT["transfer_function"], **kwargs)

    @property
    def processed_date(self):
        return self._processed_date.date

    @processed_date.setter
    def processed_date(self, value):
        self._processed_date = value


# ==============================================================================
# Site details
# ==============================================================================
class Survey(Base):
    __doc__ = write_lines(ATTR_DICT["survey"])

    def __init__(self, **kwargs):

        self.acquired_by = Person()
        self.fdsn = Fdsn()
        self.citation_dataset = Citation()
        self.citation_journal = Citation()
        self.comments = None
        self.country = None
        self.datum = None
        self.geographic_name = None
        self.name = None
        self.northwest_corner = Location()
        self.project = None
        self.project_lead = Person()
        self.release_license = "CC-0"
        self.southeast_corner = Location()
        self.summary = None
        self.survey_id = None
        self.time_period = TimePeriod()

        super().__init__(attr_dict=ATTR_DICT["survey"], **kwargs)


# =============================================================================
# Station Class
# =============================================================================
class Station(Base):
    __doc__ = write_lines(ATTR_DICT["station"])

    def __init__(self, **kwargs):
        self.id = None
        self.fdsn = Fdsn()
        self.geographic_name = None
        self.datum = None
        self.num_channels = None
        self.channels_recorded = []
        self.run_list = []
        self.channel_layout = None
        self.comments = None
        self.data_type = None
        self.orientation = Orientation()
        self.acquired_by = Person()
        self.provenance = Provenance()
        self.location = Location()
        self.time_period = TimePeriod()
        self.transfer_function = TransferFunction()

        super().__init__(attr_dict=ATTR_DICT["station"], **kwargs)

    @property
    def run_names(self):
        runs = []
        for rr in self.run_list:
            if isinstance(rr, Run):
                runs.append(rr.id)
            else:
                runs.append(rr)
        return runs


# =============================================================================
# Run
# =============================================================================
class Run(Base):
    __doc__ = write_lines(ATTR_DICT["run"])

    def __init__(self, **kwargs):
        self.id = None
        self.sample_rate = None
        # self.channels_recorded_auxiliary = []
        # self.channels_recorded_electric = []
        # self.channels_recorded_magnetic = []
        self.comments = None
        self._n_chan = None
        self.data_type = None
        self.acquired_by = Person()
        self.provenance = Provenance()
        self.time_period = TimePeriod()
        self.data_logger = DataLogger()
        self.metadata_by = Person()
        self.fdsn = Fdsn()
        self._ex = Electric()
        self._ey = Electric()
        self._hx = Magnetic()
        self._hy = Magnetic()
        self._hz = Magnetic()
        self._rrhx = Magnetic()
        self._rrhy = Magnetic()
        self._temperature = Auxiliary()

        super().__init__(attr_dict=ATTR_DICT["run"], **kwargs)

    @property
    def n_channels(self):
        number = 0
        for channel in ["auxiliary", "electric", "magnetic"]:
            channel_list = getattr(self, "channels_recorded_{0}".format(channel))
            if channel_list is not None:
                number += len(channel_list)
        return number

    @property
    def channels_recorded_all(self):
        """
        
        :return: a list of all channels recorded
        :rtype: TYPE

        """

        all_channels = []
        for recorded in ["electric", "magnetic", "auxiliary"]:
            rec_list = getattr(self, f"channels_recorded_{recorded}")
            if rec_list is None:
                continue
            else:
                all_channels += rec_list

        return all_channels

    @property
    def channels_recorded_electric(self):
        rchannels = []
        for comp in ["ex", "ey"]:
            obj = getattr(self, comp)
            if obj.component is None:
                continue
            if obj.component.lower() in [comp]:
                rchannels.append(comp)
        return rchannels

    @property
    def channels_recorded_magnetic(self):
        rchannels = []
        for comp in ["hx", "hy", "hz"]:
            obj = getattr(self, comp)
            if obj.component is None:
                continue
            if obj.component.lower() in [comp]:
                rchannels.append(comp)
        return rchannels

    @property
    def channels_recorded_auxiliary(self):
        rchannels = []
        for comp in ["temperature"]:
            obj = getattr(self, comp)
            if obj.component is None:
                continue
            if obj.component.lower() in [comp]:
                rchannels.append(comp)
        return rchannels

    @property
    def ex(self):
        return self._ex

    @ex.setter
    def ex(self, value):
        if not isinstance(value, Electric):
            msg = f"Input must be metadata.Electric not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Electric object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["ex"]:
            msg = f"Input Electric.component must be ex not {value.component}"
            self.logger.error(msg)
            raise ValueError(msg)
        self._ex.from_dict(value.to_dict())

    @property
    def ey(self):
        return self._ey

    @ey.setter
    def ey(self, value):
        if not isinstance(value, Electric):
            msg = f"Input must be metadata.Electric not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Electric object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["ey"]:
            msg = f"Input Electric.component must be ey not {value.component}"
            self.logger.error(msg)
            raise ValueError(msg)
        self._ey.from_dict(value.to_dict())

    @property
    def hx(self):
        return self._hx

    @hx.setter
    def hx(self, value):
        if not isinstance(value, Magnetic):
            msg = f"Input must be metadata.Magnetic not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Magnetic object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["hx"]:
            msg = f"Input Magnetic.component must be hx not {value.component}"
            self.logger.error(ValueError)
            raise ValueError(msg)
        self._hx.from_dict(value.to_dict())

    @property
    def hy(self):
        return self._hy

    @hy.setter
    def hy(self, value):
        if not isinstance(value, Magnetic):
            msg = f"Input must be metadata.Magnetic not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Magnetic object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["hy"]:
            msg = f"Input Magnetic.component must be hy not {value.component}"
            self.logger.error(ValueError)
            raise ValueError(msg)
        self._hy.from_dict(value.to_dict())

    @property
    def hz(self):
        return self._hz

    @hz.setter
    def hz(self, value):
        if not isinstance(value, Magnetic):
            msg = f"Input must be metadata.Magnetic not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Magnetic object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["hz"]:
            msg = f"Input Magnetic.component must be hz not {value.component}"
            self.logger.error(ValueError)
            raise ValueError(msg)
        self._hz.from_dict(value.to_dict())

    @property
    def rrhx(self):
        return self._rrhx

    @rrhx.setter
    def rrhx(self, value):
        if not isinstance(value, Magnetic):
            msg = f"Input must be metadata.Magnetic not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Magnetic object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["rrhx"]:
            msg = f"Input Magnetic.component must be rrhx not {value.component}"
            self.logger.error(ValueError)
            raise ValueError(msg)
        self._rrhx.from_dict(value.to_dict())

    @property
    def rrhy(self):
        return self._rrhy

    @rrhy.setter
    def rrhy(self, value):
        if not isinstance(value, Magnetic):
            msg = f"Input must be metadata.Magnetic not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component is None:
            msg = "assuming initial empty Magnetic object"
            self.logger.debug(msg)
        elif value.component.lower() not in ["rrhy"]:
            msg = f"Input Magnetic.component must be rrhy not {value.component}"
            self.logger.error(ValueError)
            raise ValueError(msg)
        self._rrhy.from_dict(value.to_dict())

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        if not isinstance(value, Auxiliary):
            msg = f"Input must be metadata.Magnetic not {type(value)}"
            self.logger.error(msg)
            raise ValueError(msg)
        if value.component.lower() not in ["temperature"]:
            msg = f"Input Auxiliary.component must be temperature not {value.component}"
            self.logger.error(ValueError)
            raise ValueError(msg)
        self._temperature.from_dict(value.to_dict())


# =============================================================================
# Base Channel
# =============================================================================
class Channel(Base):
    __doc__ = write_lines(ATTR_DICT["channel"])

    def __init__(self, **kwargs):
        self.type = "auxiliary"
        self.units = None
        self.channel_number = None
        self.channel_id = None
        self.comments = None
        self._component = None
        self.sample_rate = 0.0
        self.measurement_azimuth = 0.0
        self.measurement_tilt = 0.0
        self.data_quality = DataQuality()
        self.filter = Filtered()
        self.location = Location()
        self.time_period = TimePeriod()
        self.translated_azimuth = None
        self.translated_tilt = None
        self.sensor = Instrument()
        self.fdsn = Fdsn()

        super().__init__(attr_dict=ATTR_DICT["channel"], **kwargs)

    @property
    def component(self):
        return self._component

    @component.setter
    def component(self, value):
        if value is not None:
            self._component = value.lower()


# =============================================================================
# auxiliary channel
# =============================================================================
class Auxiliary(Channel):
    __doc__ = write_lines(ATTR_DICT["channel"])

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


# =============================================================================
# Electric Channel
# =============================================================================
class Electric(Channel):
    __doc__ = write_lines(ATTR_DICT["electric"])

    def __init__(self, **kwargs):
        self.dipole_length = 0.0
        self.positive = Electrode()
        self.negative = Electrode()
        self.contact_resistance = Diagnostic()
        self.ac = Diagnostic()
        self.dc = Diagnostic()
        self.units_s = None

        Channel.__init__(self, **kwargs)
        self.type = "electric"

        self._attr_dict = ATTR_DICT["electric"]


# =============================================================================
# Magnetic Channel
# =============================================================================
class Magnetic(Channel):
    __doc__ = write_lines(ATTR_DICT["magnetic"])

    def __init__(self, **kwargs):
        self.sensor = Instrument()
        self.h_field_min = Diagnostic()
        self.h_field_max = Diagnostic()

        Channel.__init__(self, **kwargs)
        self.type = "magnetic"

        self._attr_dict = ATTR_DICT["magnetic"]
