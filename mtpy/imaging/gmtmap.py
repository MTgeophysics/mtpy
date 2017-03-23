#!/usr/bin/env python

"""
mtpy/utils/gmtmap.py

Functions for the generation of geographical maps, based on a local GMT installation.

Functionality provided by the gmtpy.py package - Copyright 2009, S. Heimann (GFZ Potsdam, Germany)


@Geoscience Australia, 2017
(Alison Kirkby)

"""

#=================================================================


import numpy as np

import os.path as op

import mtpy.modeling.modem_new as mtmn
import mtpy.analysis.pt as mtpt

from mtpy.utils.exceptions import *

#=================================================================


class GMTEllipse():
    """
    class to write gmt data to plot in GMT
    """

    def __init__(self, **kwargs):
        self.workdir = kwargs.pop('workdir', None)
        self.filename = kwargs.pop('filename', None)
        self.data_type = None  # 'data', 'response', or 'residual'
        self.data_array = kwargs.pop('data_array', None)
        self.colorby = kwargs.pop('colorby', 'phimin')
        self.mt_dict = None

        self.period = None
        self.plot_period = kwargs.pop('plot_frequency', None)
        self.plot_period_index = kwargs.pop('plot_frequency_index', None)

    def _check_data_type(self, datafile, respfile):
        if ((datafile is None) and (self.data_type in ['data', 'residual'])):
            self.data_type = None
            print "Warning - provided input files and specified data type are incompatible: updating data_type"
        if ((respfile is None) and (
                self.data_type in ['response', 'residual'])):
            self.data_type = None
            print "Warning - provided input files and specified data type are incompatible: updating data_type"

    def _construct_filename(self):

    def get_plot_period_index(self):

        if self.plot_period_index is None:
            if self.plot_period is not None:
                period_diff = np.abs(
                    np.log10(
                        self.period) -
                    np.log10(
                        self.plot_period))
                self.plot_period_index = list(np.argsort(period_diff)).index(0)
                return
            else:
                self.plot_period_index = 0

        self.plot_period = self.period[self.plot_period_index]

    def read_ModEM(self, datafile=None, respfile=None):

        if ((datafile is None) and (respfile is None)):
            print "Please provide data and/or response file"

        # check data type compatible with provided input files
        self._check_data_type(datafile, respfile)

        # determine whether to automatically choose data type
        if self.data_type is None:
            find_datatype = True
        else:
            find_datatype = False

        # read files
        if datafile is not None:
            mdObj = mtmn.Data()
            mdObj.read_data_file(datafile)
            if self.workdir is None:
                self.workdir = op.dirname(datafile)
            if find_datatype:
                self.data_type = 'data'
            self.period = mdObj.period_list
        if respfile is not None:
            mrObj = mtmn.Data()
            mrObj.read_data_file(respfile)
            if self.workdir is None:
                self.workdir = op.dirname(respfile)
            if find_datatype:
                if self.data_type == 'data':
                    self.data_type = 'residual'
                else:
                    self.data_type = 'response'
            self.period = mrObj.period_list

        # get period index and period for plotting
        self.get_plot_period_index()

        # get mt_dict containing data, responses, or residual depending on
        # data_type
        if self.data_type == 'data':
            self.mt_dict = mdObj.mt_dict
        elif self.data_type == 'response':
            self.mt_dict = mrObj.mt_dict
        elif self.data_type == 'residual':
            self.mt_dict = {}
            for key in mdObj.mt_dict.keys():
                self.mt_dict[key] = mtpt.ResidualPhaseTensor(pt_object1=mdObj.mt_dict[key],
                                                             pt_object2=mrObj.mt_dict[key])

    def build_data_array(self):

        if self.mt_dict is None:
            print "Cannot save GMT, please read a ModEM data and/or response file first"

        self.data_array = np.zeros((len(self.mt_dict), 6))

        for i in range(len(self.mt_dict)):
            for ii, att in enumerate(
                    ['lon', 'lat', self.colorby, 'azimuth', 'phimin', 'phimax', 'skew']):
                self.data_array[i, ii] = getattr(self.mt_dict[i], att)

    def write_data(self):

        if self.data_array is None:
            self.build_data_array()

        # write to text file in correct format
        fmt = ['%+11.6f', '%+10.6f'] + ['%+9.4f'] * 2 + ['%8.4f'] * 2
        np.savetxt(op.join(savepath, filename), self.data_array, fmt=fmt)


class GMTScript():
    """
    class to write a template gmt script for plotting data in gmt

    """

    def __init__(self, workdir, **kwargs):

        self.workdir = workdir
        self.xlim = None
        self.ylim = None
        self.plotdata_dict = []  # list containing GMTData objects
        self.cmap = 'polar'
        self.clim = None
        self.mapsize = '18c'
