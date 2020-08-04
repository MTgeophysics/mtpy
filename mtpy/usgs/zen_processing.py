# -*- coding: utf-8 -*-
"""
ZEN PROCESSING TOOLS
=======================
    * Interface with BIRRP
    * TODO: interface with EMTF
    * Prepare files for Zonge Processing codes
Created on Fri Sep 16 14:29:43 2016
@author: jpeacock
"""
#==============================================================================
import numpy as np
import datetime
import os
import sys
from io import StringIO
from pathlib import Path
import pandas as pd

import mtpy.processing.birrp as birrp
import mtpy.utils.configfile as mtcfg
import mtpy.utils.exceptions as mtex
import mtpy.imaging.plotnresponses as plotnresponses
import mtpy.imaging.plotresponse as plotresponse
import mtpy.usgs.zen as zen
import mtpy.core.edi as mtedi
import mtpy.core.ts as mtts
from mtpy.usgs import z3d_collection as zc

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator

#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'
datetime_sec = '%Y-%m-%d %H:%M:%S'
#==============================================================================

#==============================================================================
# Survey configuration file
#==============================================================================
class SurveyConfig(object):
    """
    survey config class
    will setup a survey configuration file that has the form of:
    [station]
        b_instrument_amplification = 1
        b_instrument_type = coil
        b_logger_gain = 1
        b_logger_type = zen
        b_xaxis_azimuth = 0
        b_yaxis_azimuth = 90
        box = 26
        date = 2015/06/09
        e_instrument_amplification = 1
        e_instrument_type = Ag-Agcl electrodes
        e_logger_gain = 1
        e_logger_type = zen
        e_xaxis_azimuth = 0
        e_xaxis_length = 100
        e_yaxis_azimuth = 90
        e_yaxis_length = 100
        elevation = 2113.2
        hx = 2274
        hy = 2284
        hz = 2254
        lat = 37.7074236995
        location = Earth
        lon = -118.999542099
        network = USGS
        notes = Generic config file
        rr_box = 25
        rr_date = 2015/06/09
        rr_hx = 2334
        rr_hy = 2324
        rr_lat = 37.6909139779
        rr_lon = -119.028707542
        rr_station = 302
        sampling_interval = all
        save_path = \home\mtdata\survey_01\mt_01
        station = 300
        station_type = mt
    """
    def __init__(self, **kwargs):
        self.b_instrument_amplification = 1
        self.b_instrument_type = 'induction coil'
        self.b_logger_gain = 1
        self.b_logger_type = 'zen'
        self.b_xaxis_azimuth = 0
        self.b_yaxis_azimuth = 90
        self.box = 24
        self.date = '01/01/00'
        self.e_instrument_amplification = 1
        self.e_instrument_type = 'Ag-Agcl electrodes'
        self.e_logger_gain = 1
        self.e_logger_type = 'zen'
        self.e_xaxis_azimuth = 0
        self.e_xaxis_length = 100
        self.e_yaxis_azimuth = 90
        self.e_yaxis_length = 100
        self.elevation = 2113.2
        self.hx = 2324
        self.hy = 2314
        self.hz = 2334
        self.lat = 0.0
        self.location = 'Earth'
        self.lon = 0.0
        self.network = 'USGS'
        self.notes = 'Generic config file'
        self.sampling_interval = 'all'
        self.station = 'mb000'
        self.station_type = 'mt'
        self.save_path = None

        self.rr_lat = None
        self.rr_lon = None
        self.rr_station = None
        self.rr_date = None
        self.rr_box = None

        for key in kwargs:
            setattr(self, key, kwargs[key])

    def from_df(self, z3d_df):
        """

        :param z3d_df: DESCRIPTION
        :type z3d_df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        z3d_df.remote = z3d_df.remote.astype(str)
        s_df = z3d_df[z3d_df.remote == 'False']
        s_df.start = pd.to_datetime(s_df.start)

        self.b_xaxis_azimuth = s_df[s_df.component == 'hx'].azimuth.mode()[0]
        self.b_yaxis_azimuth = s_df[s_df.component == 'hy'].azimuth.mode()[0]
        self.box = s_df.zen_num.mode()[0]
        self.date = s_df.start.min().isoformat()
        self.e_instrument_type = 'Ag-Agcl electrodes'
        self.e_logger_type = 'zen'
        self.e_xaxis_azimuth = s_df[s_df.component == 'ex'].azimuth.mode()[0]
        self.e_xaxis_length = s_df[s_df.component == 'ex'].dipole_length.mode()[0]
        self.e_yaxis_azimuth = s_df[s_df.component == 'ey'].azimuth.mode()[0]
        self.e_yaxis_length = s_df[s_df.component == 'hx'].dipole_length.mode()[0]
        self.elevation = s_df.elevation.median()
        self.hx = s_df[s_df.component == 'hx'].coil_number.mode()[0]
        self.hy = s_df[s_df.component == 'hy'].coil_number.mode()[0]
        try:
            self.hz = s_df[s_df.component == 'hz'].coil_number.mode()[0]
        except IndexError:
            self.hz = ''
        self.lat = s_df.latitude.median()
        self.location = 'Earth'
        self.lon = s_df.longitude.median()
        self.network = 'USGS'
        self.notes = 'Generic config file'
        self.sampling_interval = 'all'
        self.station = s_df.station.mode()[0]
        self.station_type = 'mt'

        rr_df = z3d_df[z3d_df.remote == 'True']
        rr_df.start = pd.to_datetime(rr_df.start)
        if len(rr_df) > 0:
            self.rr_lat = []
            self.rr_lon = []
            self.rr_station = []
            self.rr_date = []
            self.rr_box = []
            for station in rr_df.station.unique():
                rdf = rr_df[rr_df.station == station]
                self.rr_lat.append(rdf.latitude.median())
                self.rr_lon.append(rdf.longitude.median())
                self.rr_station.append(rdf.station.mode()[0])
                self.rr_date.append(rdf.start.min().isoformat())
                self.rr_box.append(rdf.zen_num.mode()[0])

    def write_survey_config_file(self, save_path=None):
        """
        write a survey config file to save path
        """

        if save_path is not None:
            self.save_path = Path(save_path)
        fn = self.save_path.joinpath('{0}.cfg'.format(self.station))
        mtcfg.write_dict_to_configfile({self.station:self.__dict__}, fn)

        #('Wrote survey config file to {0}'.format(fn))

        return fn

    def read_survey_config_file(self, survey_cfg_fn, station):
        """
        Parameters
        ----------
        survey_dfg_fn : TYPE
            DESCRIPTION.
        station : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        survey_cfg_dict = mtcfg.read_survey_configfile(survey_cfg_fn)[station]

        for key, value in survey_cfg_dict.items():
            setattr(self, key, value)

#==============================================================================
# Z3D files to EDI using BIRRP
#==============================================================================
class Z3D2EDI(object):
    """
    go from z3d files to .edi using BIRRP as the processing code
    Arguments
    ---------------
        **station_z3d_dir** : string
                          full path to station directory where the Z3D files
                          to be processed are.
        **rr_station_z3d_dir** : string
                            full path to remote reference directory where
                            Z3D files exist.
        **birrp_exe** : string
                        full path to birrp executable
        **coil_cal_path** : string
    Methods
    ----------
        **get_z3d_fn_blocks** : get z3d files in blocks of sampling rate and
                                date
        **make_mtpy_ascii_files** : make mtpy ascii files from given blocks
                                    of file names
        **get_schedules_fn_from_dir** : get mtpy ascii file names in blocks
                                        by sampling rate and date from a TS
                                        directory.
        **get_schedules_fn** : get mtpy ascii file names in blocks
                               by sampling rate and date from an array.
        **write_script_files** : write birrp script files for each sampling
                                 rate schedule block
        **run_birrp** : run birrp from a script file and write an .edi file
        **write_edi_file** : write edi file from birrp outputs and given
                            survey parameters.
        **plot_responses** : plots all edi files output from each sampling rate
        **process_data** : a convinience function to go from Z3D files to
                           edi files.
        **combine_edi_files** : combine all edi files from each sampling rate
                                given a frequency range for each.
    Example
    ------------
        >>> import mtpy.usgs.zen_processing as zp
        >>> zp_obj = zp.Z3D2EDI()
        >>> zp_obj.station_z3d_dir = r"/home/data/mt01"
        >>> zp_obj.rr_station_z3d_dir = r"/home/data/mt02"
        >>> zp_obj.birrp_exe = r"/home/bin/birrp52"
        >>> zp_obj.coil_cal_path = r"/home/data/ant_calibration"
        >>> plot_obj, comb_edi = zp_obj.process_data(df_list=[4096, 256, 16])
    """

    def __init__(self, station_z3d_dir=None, **kwargs):

        self._station_z3d_dir = None
        self._station_ts_dir = None
        self._rr_station_z3d_dir = None
        self._rr_station_ts_dir = None

        self.station_z3d_dir = station_z3d_dir

        self.survey_config = SurveyConfig(save_path=self.station_z3d_dir)
        self.survey_config_fn = None
        self.birrp_config_fn = None
        self.birrp_exe = r"/home/peacock/Documents/birrp52/SourceCode/birrp52_big"
        self.calibration_path = None
        self.calibration_dict = {}
        self.num_comp = 5
        self.df_list = [4096, 256, 4]
        self.max_blocks = 4
        self._max_nread = 20000000
        # number of lines for birrp to skip in file ts header
        self._header_len = 23
        self._tol_dict = {4096: {'s_diff': 5 * 60 * 4096,
                                      'min_points': 2**18},
                               256: {'s_diff': 3 * 3600 * 256,
                                     'min_points': 2**19},
                               4: {'s_diff': 4 * 3600 * 4,
                                   'min_points': 2**14}}

        # data types for different aspects of getting information
        if sys.version_info[0] == 2:
            self._ts_fn_dtype = np.dtype([(u'station','S6'),
                                          (u'npts', np.int),
                                          (u'df', np.int),
                                          (u'start_dt', 'S22'),
                                          (u'end_dt', 'S22'),
                                          (u'comp', 'S2'),
                                          (u'fn', 'S100'),
                                          (u'calibration_fn', 'S100')])

            self._birrp_fn_dtype = np.dtype([('fn', 'S100'),
                                             ('nread', np.int),
                                             ('nskip', np.int),
                                             ('comp', 'S2'),
                                             ('calibration_fn', 'S100'),
                                             ('rr', np.bool),
                                             ('rr_num', np.int),
                                             ('start_dt', 'S22'),
                                             ('end_dt', 'S22')])
        elif sys.version_info[0] == 3:
            self._ts_fn_dtype = np.dtype([('station','U6'),
                                          ('npts', np.int),
                                          ('df', np.int),
                                          ('start_dt', 'U22'),
                                          ('end_dt', 'U22'),
                                          ('comp', 'U2'),
                                          ('fn', 'U100'),
                                          ('calibration_fn', 'U100')])

            self._birrp_fn_dtype = np.dtype([('fn', 'U100'),
                                             ('nread', np.int),
                                             ('nskip', np.int),
                                             ('comp', 'U2'),
                                             ('calibration_fn', 'U100'),
                                             ('rr', np.bool),
                                             ('rr_num', np.int),
                                             ('start_dt', 'U22'),
                                             ('end_dt', 'U22')])

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    @property
    def station_z3d_dir(self):
        return self._station_z3d_dir

    @station_z3d_dir.setter
    def station_z3d_dir(self, station_z3d_dir):
        if station_z3d_dir in [None, 'None']:
            self._station_z3d_dir = None
        else:
            self._station_z3d_dir = Path(station_z3d_dir)

    @property
    def rr_station_z3d_dir(self):
        return self._rr_station_z3d_dir

    @rr_station_z3d_dir.setter
    def rr_station_z3d_dir(self, rr_station_z3d_dir):
        if rr_station_z3d_dir in [None, 'None']:
            self._rr_station_z3d_dir = None
        if isinstance(rr_station_z3d_dir, (str, Path)):
            self._rr_station_z3d_dir = [Path(rr_station_z3d_dir)]
        elif isinstance(rr_station_z3d_dir, list):
            self._rr_station_z3d_dir = [Path(r) for r in rr_station_z3d_dir]

    @property
    def station_ts_dir(self):
        if self._station_ts_dir is None:
            return self.station_z3d_dir.joinpath('TS')
        else:
            return self._station_ts_dir

    @station_ts_dir.setter
    def station_ts_dir(self, station_ts_dir):
        if station_ts_dir in [None, 'None']:
            self._station_ts_dir = None
        else:
            self._station_ts_dir = Path(station_ts_dir)

    @property
    def rr_station_ts_dir(self):
        if self._rr_station_ts_dir is None and \
            self._rr_station_z3d_dir is not None:
            return [p.joinpath('TS') for p in self.rr_station_z3d_dir]
        else:
            return self._rr_station_ts_dir

    @rr_station_ts_dir.setter
    def rr_station_ts_dir(self, rr_station_ts_dir):
        if rr_station_ts_dir in [None, 'None']:
            self._rr_station_ts_dir = None
        if isinstance(rr_station_ts_dir, (str, Path)):
            self._rr_station_ts_dir = [Path(rr_station_ts_dir)]
        elif isinstance(rr_station_ts_dir, list):
            self._rr_station_ts_dir = [Path(r) for r in rr_station_ts_dir]

    def get_calibrations(self, calibration_path):
        """
        Get coil calibration file names from a given directory

        :param calibration_path: path to calibration files
        :type calibration_path: string or Path

        :return: sets calibration_dict internally for use later
        """
        self.calibration_path = calibration_path
        if self.calibration_path is None:
            print('ERROR: Calibration path is None')
            self.calibration_dict = {}

        if not isinstance(self.calibration_path, Path):
            self.calibration_path = Path(self.calibration_path)

        if not self.calibration_path.exists():
            print('WARNING: could not find calibration path: '
                  '{0}'.format(self.calibration_path))
            self.calibration_dict = {}

        self.calibration_dict = {}
        for cal_fn in self.calibration_path.glob('*.csv'):
            cal_num = cal_fn.stem
            self.calibration_dict[cal_num] = cal_fn

    def convert_z3d_to_mtts(self, station_z3d_dir, rr_station_z3d_dir=None,
                            use_blocks_dict=None, overwrite=False,
                            combine=True, notch_dict=None,
                            combine_sampling_rate=4, calibration_path=None):
        """
        Convert Z3D files into MTTS objects and write ascii files for input
        into BIRRP.  Will write a survey configuration file that can be read
        in when making EDI files. Writes the DataFrame to a .csv file in the
        Z3D directory.

        :param station_z3d_dir: path to station z3d folder.  Will be
                                converted to a Path object on setting.
        :type station_z3d_dir: string or Path object

        :param rr_station_z3d_dir: path to remote station z3d folder.  Will be
                                converted to a Path object on setting.
        :type rr_station_z3d_dir: string or Path object, optional

        :param use_blocks_dict: Dictionary of blocks to use.  If None is input
                                then will use all blocks. Keys are sampling
                                rates, values are a list of block numbers
                                --> {256:[0, 1, 3], 4:[0]}
        :type use_blocks_dict: dictionary , optional

        :param overwrite: Overwrite existing ascii files, defaults to False
        :type overwrite: [ True | False ], optional

        :param combine: combine all Z3D files for each component into a single
                        file samples at combine_sampling_rate.
                        defaults to True
        :type combine: [ True | False ], optional

        :param notch_dict: dictionary of notches to filter out. Can be useful
                           for data with 60 Hz. Keys are sampling rates,
                           values are notch frequency and width. If an empty
                           dictionary is entered 60 Hz and harmonics will be
                           filtered out. defaults to None
        :type notch_dict: TYPE, optional

        :param combine_sampling_rate: sampling rate for combined file,
                                      defaults to 4
        :type combine_sampling_rate: int, optional

        :param calibration_path: path to calibration files, defaults to None
        :type calibration_path: string or Path, optional

        :return: dataframe containing information on Z3D files to be used
                 later
        :rtype: pandas.DataFrame

        :return: path to csv file
        :rtype: string

        :Example: ::

            >>> zp_obj = zp.Z3D2EDI()
            >>> rr_z3d_station_dir = ['/home/mt02', '/home/mt03']
            >>> use_blocks_dict = {4096: [0, 1, 2], 256: [1, 2], 4:[0]}
            >>> notch_dict = {4096: {}, 256: None, 4: None}
            >>> kw_dict = {'rr_station_z3d_dir': rr_station_z3d_dir,
                           'use_blocks_dict': use_blocks_dict,
                           'notch_dict': notch_dict}
            >>> z3d_df, z3d_csv = zp_obj.convert_z3d_to_mtts(r"/home/mt01",
                                                             **kw_dict)

        """
        if calibration_path is not None:
            self.calibration_path = calibration_path

        self.station_z3d_dir = Path(station_z3d_dir)
        if rr_station_z3d_dir is not None:
            self.rr_station_z3d_dir = rr_station_z3d_dir

        kw_dict = {'block_dict': use_blocks_dict,
                   'notch_dict': notch_dict,
                   'overwrite': overwrite,
                   'combine': combine,
                   'combine_sampling_rate': combine_sampling_rate,
                   'calibration_path': self.calibration_path}

        zc_obj = zc.Z3DCollection()
        station_df, station_csv = zc_obj.from_dir_to_mtts(self.station_z3d_dir,
                                                          **kw_dict)
        if self.rr_station_z3d_dir is not None:
            kw_dict['remote'] = True
            if not isinstance(self.rr_station_z3d_dir, list):
                rr_z3d_dir = [self.rr_station_z3d_dir]
            self.rr_station_ts_dir = []
            for rr_path in self.rr_station_z3d_dir:
                rr_df, rr_csv = zc_obj.from_dir_to_mtts(rr_path,**kw_dict)
                station_df = station_df.append(rr_df)
                self.rr_station_ts_dir.append(Path(rr_path).joinpath('TS'))
        processing_csv = Path(station_z3d_dir).joinpath(
                    '{0}_processing_df.csv'.format(self.station_z3d_dir.name))
        station_df.to_csv(processing_csv)

        self.station_ts_dir = Path(station_z3d_dir).joinpath('TS')

        # write configuration file for edi, this should be deprecated later
        self.survey_config.from_df(station_df)
        self.survey_config.save_path = self.station_ts_dir
        self.survey_config_fn = self.survey_config.write_survey_config_file()

        return station_df, processing_csv

    def make_block_entry(self, entry, nskip, nread, rr_num=0):
        """
        Make a block entry from a given entry and parameters

        :param entry: dataframe row entry
        :type entry: pandas.DataFrame

        :param nkip: number of points to skip
        :type nkip: int

        :param nread: number of points to read
        :type nread: int

        :param rr_num: remote reference number, defaults to 0
        :type rr_num: int, optional

        :return: dictionary of the required values
        :rtype: dictionary

        """

        r_dict = dict([('fn', entry.fn_ascii),
                       ('nread', nread),
                       ('nskip', nskip),
                       ('comp', entry.component),
                       ('calibration_fn', entry.cal_fn),
                       ('rr', entry.remote),
                       ('rr_num', rr_num),
                       ('start', entry.start),
                       ('stop', entry.stop),
                       ('sampling_rate', entry.sampling_rate),
                       ('station', entry.station)])
        return r_dict

    def compare_times(self, entry, start, stop):
        """
        Compare start and stop time from a given row entry.  Calculates the
        number of points to skip and read.

        :param entry: row entry from z3d DataFrame
        :type entry: pandas.DataFrame

        :param start: maximum start time for a given block
        :type start: pandas.Timestamp

        :param stop: minimum stop time for a given block
        :type stop: pandas.Timestamp

        :return: dictionary that includes the number of points to skip, read
                 and time difference at the beginning and end
        :rtype: dictionary.

        """

        sr = entry.sampling_rate
        info_dict = {'nskip': 0, 'nread': 0, 'start_diff': 0, 'end_diff': 0}

        # estimate time difference at beginning
        start_time_diff = sr * (start - entry.start).total_seconds()
        info_dict['start_diff'] = start_time_diff
        # if difference is positive entry starts before station
        if start_time_diff > 0:
            info_dict['nskip'] = self._header_len + start_time_diff
            info_dict['nread'] = entry.nread - start_time_diff
        else:
            info_dict['nskip'] = self._header_len
            info_dict['nread'] = entry.nread

        # check the end times
        end_time_diff = sr * (entry.stop - stop).total_seconds()
        info_dict['end_diff'] = end_time_diff
        # if end diff is positive entry ends after station
        if end_time_diff > 0:
            info_dict['nread'] -= end_time_diff

        return info_dict

    def make_block_df(self, block_list):
        """
        convert a given block list into a DataFrame with the appropriate
        data types specifically for the start and stop.

        :param block_list: list of entries from make_block_entry for a given
                           schedule block
        :type block_list: list

        :return: DataFrame with the appropriate data types
        :rtype: pandas.DataFrame

        """

        block_df = pd.DataFrame(block_list)
        block_df = block_df.infer_objects()
        block_df.start = pd.to_datetime(block_df.start)
        block_df.stop = pd.to_datetime(block_df.stop)

        return block_df

    def align_block(self, block_df):
        """
        Align a given block DataFrame such that the start and stop times are
        the same for each component.  Uses compare_times to align.

        :param block_df: block data frame from make_block_df
        :type block_df: pandas.DataFrame

        :return: aligned data frame
        :rtype: pandas.DataFrame

        """

        b_start = block_df.start.max()
        b_stop = block_df.stop.min()

        for entry in block_df.itertuples():
            diff_dict = self.compare_times(entry, b_start, b_stop)
            for key in ['nskip', 'nread']:
                block_df.at[entry.Index, key] = diff_dict[key]

        block_df.nread = block_df.nread.min()
        return block_df

    def get_birrp_dict(self, df, df_list=[4096, 256, 4],
                       use_blocks_dict=None):
        """
        Make a dictionary that can be used to write script files for birrp.
        This will align each schedule block.

        :param df: z3d data frame
        :type df: pandas.DataFrame

        :param df_list: list of sampling rates to process,
                        defaults to [4096, 256, 4]
        :type df_list: list

        :param use_blocks_dict: dictionary of blocks to use for each sampling
                                rate.  defaults to None which uses all blocks
                                possible.
        :type use_blocks_dict: dictionary

        :return: an array of rec arrays for each block
        :rtype: dictionary{sampling_rate:numpy.recarray}

        :Example: ::

            >>> b_dict = zp_obj.get_birrp_dict(df, df_list=[256, 4],
                                               use_blocks_dict={256:[0, 2],
                                                                '4':[0]})

        """
        zc_obj = zc.Z3DCollection()
        use_blocks_dict = zc_obj._validate_block_dict(df, use_blocks_dict)

        birrp_dict = {}
        for sr in df.sampling_rate.unique():
            if sr not in df_list:
                continue
            sr_list = []
            sr_df = df[(df.sampling_rate == sr) & (df.remote == 'False')]
            rr_df = df[(df.sampling_rate == sr) & (df.remote == 'True')]

            rr_stations = dict([(rr, ii) for ii, rr in
                                enumerate(rr_df.station.unique())])

            # sort through station blocks first
            block_count = 0
            blocks_read_total = 0
            for block in sr_df.block.unique():
                if block not in use_blocks_dict[sr]:
                    continue
                if block_count > self.max_blocks:
                    print('WARNING: Max blocks of {0} reached for {1}'.format(
                          self.max_blocks, sr))
                    break
                block_list = []
                block_df = sr_df[sr_df.block == block]
                # find the latest start time
                start = block_df.start.max()
                stop = block_df.stop.min()
                if str(stop) == 'NaT':
                    print('WARNING: Skipping block {0} for {1}.'.format(block,
                                                                        sr) +
                          '\n\tReason: no end time')
                    continue

                # get information for station block and align
                for entry in block_df.itertuples():
                    block_list.append(self.make_block_entry(entry, 0,
                                                            entry.n_samples))

                # make the block into a dataframe
                block_birrp_df = self.make_block_df(block_list)

                # --> get remote reference blocks
                # check for start time
                rr_block_list = []
                for rr_entry in rr_df.itertuples():
                    if rr_entry.start > stop:
                        print('INFO: Skipping {0} starts after station'.format(
                              rr_entry.station))
                        continue
                    t_diff = abs((rr_entry.start - start).total_seconds()) * sr
                    # check to see if the difference is within given tolerance
                    if t_diff <= self._tol_dict[sr]['s_diff']:
                        # check number of samples
                        rr_samples = rr_entry.n_samples - t_diff
                        if rr_samples < self._tol_dict[sr]['min_points']:
                            print('WARNING: skipping {0} block {1} df {2} at {3}'.format(
                                  rr_entry.station, rr_entry.block,
                                  rr_entry.sampling_rate,
                                  rr_entry.start) +
                                  '\n\tNot enough points {0}'.format(rr_samples))
                        # make a block entry and append
                        else:
                            rr_block_list.append(self.make_block_entry(rr_entry,
                                                 0,
                                                 rr_entry.n_samples,
                                                 rr_stations[rr_entry.station]))
                            
                # check to make sure there are remote references
                if len(rr_block_list) > 1:
                    rr_block_birrp_df = self.make_block_df(rr_block_list)
                    block_birrp_df = block_birrp_df.append(rr_block_birrp_df)

                # align block and append
                sr_list.append(self.align_block(block_birrp_df.reset_index()))
                block_count += 1
                blocks_read_total += block_birrp_df.nread.mean()
                if blocks_read_total > self._max_nread:
                    dn = blocks_read_total - self.max_nread
                    block_birrp_df.nread = block_birrp_df.nread.mean() - dn
                    print('WARNING: reached maximum points to read' +
                          ' {0}. \nCutting last block to {1}.'.format(
                          self._max_nread, block_birrp_df.nread.mean()))

            birrp_dict[sr] = np.array([a_df.to_records(index=False) for a_df
                                       in sr_list])

        return birrp_dict

    def write_script_files(self, birrp_arr_dict, save_path=None,
                           birrp_params_dict={}, **kwargs):
        """
        Write BIRRP script files for each sampling rate in the given
        dictionary.  This will sort the appropriate parameters from the given
        data.

        :param birrp_arr_dict: dictionary with keys as sampling rates and
                               values as arrays of recarrays for each schedule
                               block.
        :type birrp_arr_dict: dictionary{sampling_rate:numpy.ndarray(
                                                       numpy.recarray)}

        :param save_path: path to save script file, defaults to
                          station_ts_dir/BF/sampling_rate/station.script
        :type save_path: string or Path

        :param birrp_params_dict: dictionary of birrp parameters, default is
                                  an empty dictionary that will use default
                                  parameters that work well in most cases.

        .. seealso:: mtpy.processing.birrp

        :return: list of script file paths
        :rtype: list

        :Example: ::

            >>> sfn_list = zp_obj.write_script_files(birrp_dict,
                                                     birrp_params={'tbw':2})

        """

        # make save path
        if save_path is None:
            save_path = self.station_ts_dir.joinpath('BF')
        elif not isinstance(save_path, Path):
            save_path = Path(save_path)
        if not save_path.exists():
            save_path.mkdir()

        script_fn_list = []
        # loop through by keys, which should be sampling rates
        for df_key, fn_arr in birrp_arr_dict.items():
            # make a path unique to the sampling rate
            bf_path = save_path.joinpath('{0:.0f}'.format(df_key))
            if not bf_path.exists():
                bf_path.mkdir()

            # get the fn_array, and make sure that it is a ndarray type
            birrp_fn_arr = birrp_arr_dict[df_key]

            # make a script object passing on the desired birrp parameters
            try:
                birrp_script_obj = birrp.ScriptFile(fn_arr=fn_arr)
            except birrp.ScriptFileError as error:
                print('ERROR: {0}'.format(error))
                print('WARNING: Skipping script file for {0}'.format(df_key))
                continue

            # get station name
            remotes = np.where(fn_arr['rr'] == False)
            station = str(np.unique(fn_arr[remotes]['station'])[0])

            # add parameters to birrp_params_dict
            birrp_params_dict['ofil'] = str(bf_path.joinpath(station))
            if df_key == 16:
                birrp_params_dict['nfft'] = 2**16
                birrp_params_dict['nsctmax'] = 11
            if df_key == 4:
                if birrp_fn_arr['nread'].sum(axis=0)[0] / 2**16 < 6:
                    birrp_params_dict['nfft'] = 2**15
                    birrp_params_dict['nsctmax'] = 11
                else:
                    birrp_params_dict['nfft'] = 2**16
                    birrp_params_dict['nsctmax'] = 12

            birrp_script_obj.from_dict(birrp_params_dict)

            # write script file
            b_script = bf_path.joinpath('{0}.script'.format(station))
            try:
                birrp_script_obj.write_script_file(script_fn=b_script)
            except ValueError as error:
                print(fn_arr['fn'])
                raise ValueError(error)
                
            # write a birrp parameter configuration file
            birrp_script_obj.write_config_file(bf_path.joinpath(station))
            script_fn_list.append(birrp_script_obj.script_fn)

        return script_fn_list

    def run_birrp(self, script_fn_list=None, birrp_exe=None):
        """
        run birrp given the specified files

        :param script_fn_list: list of script file paths
        :type script_fn_list: list

        :param birrp_exe: path to BIRRP executable
        :type birrp_exe: string
        """

        if script_fn_list is None:
            raise IOError('Need to input a script file or list of script files')

        if birrp_exe is not None:
            self.birrp_exe = birrp_exe


        if type(script_fn_list) is list:
            self.edi_fn = []
            for script_fn in script_fn_list:
                out_str = birrp.run(self.birrp_exe, script_fn)
                print('INFO: BIRRP Processing \n {0}'.format(out_str))

                output_path = os.path.dirname(script_fn)
                try:
                    self.edi_fn.append(self.write_edi_file(output_path,
                                       survey_config_fn=self.survey_config_fn,
                                       birrp_config_fn=self.birrp_config_fn))
                except Exception as error:
                    print('ERROR: {0} did not run properly'.format(script_fn))
                    print('ERROR: {0}'.format(error))

        elif type(script_fn_list) is str:
            out_str = birrp.run(self.birrp_exe, script_fn_list)

            output_path = os.path.dirname(script_fn_list)
            self.edi_fn = self.write_edi_file(output_path,
                                      survey_config_fn=self.survey_config_fn,
                                      birrp_config_fn=self.birrp_config_fn)

    def write_edi_file(self, birrp_output_path, survey_config_fn=None,
                       birrp_config_fn=None):
        """
        Write an edi file from outputs of BIRRP

        :param birrp_output_path: full path to directory with BIRRP outputs
        :type birrp_output_path: string

        :param survey_config_fn: path to survey configuration file
        :type survey_config_fn: string

        :param birrp_config_fn: path to BIRRP configuration file
        :type birrp_config_fn: string

        :return: path to edi file
        :rtype: string

        """
        if self.survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn

        if self.survey_config_fn is None:
            ts_find = birrp_output_path.find('TS')
            if ts_find > 0:
                ts_dir = birrp_output_path[0:ts_find+2]
                for fn in os.listdir(ts_dir):
                    if fn[-4:] == '.cfg':
                        self.survey_config_fn = os.path.join(ts_dir, fn)
                        self.survey_config.read_survey_config_file(
                            self.survey_config_fn, 
                            self.station_z3d_dir.name)

        j2edi_obj = birrp.J2Edi(station=self.survey_config.station,
                                survey_config_fn=self.survey_config_fn,
                                birrp_dir=birrp_output_path,
                                birrp_config_fn=self.birrp_config_fn)

        edi_fn = j2edi_obj.write_edi_file()

        return edi_fn

    def plot_responses(self, edi_fn_list=None):
        """
        plot all the edi files that were created.

        :param edi_fn_list: list of EDI file paths to plot
        :type edi_fn_list: list

        :return: plot object
        :rtype: mtpy.imaging.plotnresponse.PlotMultipleResponses

        :Example: ::

            >>> rp = zp_obj.plot_responses(["/home/mt01/TS/BF/4096/mt01edi",
                                            "/home/mt01/TS/BF/256/mt01.edi",
                                            "/home/mt01/TS/BF/4/mt01.edi"])
            >>> rp.tipper_limits = (-.5, .5)
            >>> rp.res_limits = (.1, 1000)
            >>> rp.redraw_plot()
        """
        if edi_fn_list is not None:
            self.edi_fn = edi_fn_list

        if type(self.edi_fn) is list:
            # check file lengths to make sure there are no zero length files
            for edi_fn in self.edi_fn:
                fn_size = os.path.getsize(edi_fn)
                if fn_size < 3000:
                    self.edi_fn.remove(edi_fn)
                if len(self.edi_fn) == 0:
                    raise ValueError('No good .edi files where produced')
            resp_plot = plotnresponses.PlotMultipleResponses(fn_list=self.edi_fn,
                                                         plot_style='compare',
                                                         plot_tipper='yri')
        elif type(self.edi_fn) is str:
            if os.path.getsize(self.edi_fn) < 3000:
                raise ValueError('No good .edi files where produced')
            resp_plot = plotresponse.PlotResponse(fn=self.edi_fn,
                                                  plot_tipper='yri')

        return resp_plot

    def process_data(self, df_fn=None, df_list=[4096, 256, 4], plot=True,
                     notch_dict={}, use_blocks_dict=None,  overwrite=False,
                     sr_dict={4096:(1000., 4),
                              1024:(3.99, 1.),
                              256:(3.99, .126),
                              4:(.125, .0001)},
                     birrp_param_dict={}, **kwargs):
        """
        process_data is a convinience function that will process Z3D files
        and output an .edi file.  The workflow is to convert Z3D files to
        MTpy ascii format, apply a notch filter if desired, decimate if
        desired.  Do the same for a remote reference if specified.  Make a
        BIRRP script file for each sampling rate, process data, take outputs
        and write an edi file and plot result.  The whole process will take
        a few minutes so be patient.

        :param df_fn: path to existing dataframe that was used for processing
                      useful if changing parameters to a previously processed
                      station.
        :type df_fn: string

        :param df_list: list of sampling rates
        :type df_list: list, defaults to [4096, 256, 4]

        :param notch_dict: dict(sampling_rate: notches to filter)
                           dictionary of notches to filter for each sampling
                           rate.  Keys are sampling rates (int) and values
                           are dictionary of notches to filter
                           ex. {4096:{60, 120, 180}} to filter out 60 Hz noise
                           *default* is {} which filters out 60 Hz noise and
                           harmonics.
        :type notch_dict: dictionary

        :param sr_dict: dict(sampling_rate: (max_freq, max_freq))
                        dictionary of min and max frequencies to use for
                        each sampling rate when making an .edi file.  The
                        defaults usually work well, but check the plot
                        to see if there is a better frequency range for
                        each sampling rate.
                        *default* is
                          {4096:(1000., 4),
                          1024:(3.99, 1.),
                          256:(3.99, .126),
                          16:(.125, .0001)}
        :type sr_dict: dictionary

        :param use_blocks_dict: dict(sampling_rate: list of blocks to use)
                               dictionary with sampling rate as keys and
                               list of which blocks to use as values
                               ex. {4096: [0, 2]} to skip the second block
                               use value 'all' to use all blocks
                               *default* is None which will use all blocks
        :type use_blocks_dict: dictionary

        :param birrp_param_dict: dict(birrp_param: value)
                                 dictionary of birrp parameters where the
                                 key is the birrp parameter and value is
                                 what you want that parameter to be.  See
                                 mtpy.processing.birrp.write_script_file
                                 or the BIRRP manual for details on birrp
                                 parameters
                                 ex. {'nar':5}
                                 *default* is {}
        :type birrp_param_dict: dictionary

        :return: plot_response object
        :rtype: mtpy.imaging.plotnresponse.PlotMultipleResponses

        :return: full path to combined EDI file
        :rtype: string

        :return: processing dataframe
        :rtype: pandas.DataFrame

        :Example: ::

            >>> b_param_dict = {'ilev': ilev,
                                'c2threshb': .35,
                                'c2threshe': .35,
                                'c2thresh1': .35,
                                'ainuin': .95,
                                'ainlin': .05,
                                'nar': 11,
                                'tbw': 2}
            >>> zp_obj = zp.Z3D2EDI(station_z3d_dir='/home/mt01')
            >>> zp_obj.rr_station_z3d_dir = ['/home/mt02',  '/home/mt03']
            >>> zp_obj.birrp_exe = '/home/bin/birrp52.exe'
            >>> zp_obj.calibration_path = '/home/mt/calibrations'
            >>> plot_obj, comb_edi_fn, zdf = zp_obj.process_data(
                                                 birrp_param_dict=b_param_dict)

        """
        if df_list is not None:
            self.df_list = df_list

        # get start time
        st = datetime.datetime.now()

        # make files into mtpy files
        if df_fn is not None:
            zc_obj = zc.Z3DCollection()
            z3d_df = zc_obj.from_csv(df_fn)
        else:
            # skip the block dict, want to look through all the files to get the
            # data frame.
            kw_dict = {'use_blocks_dict': None,
                       'overwrite': overwrite}
            z3d_df, cfn = self.convert_z3d_to_mtts(self.station_z3d_dir,
                                                   self.rr_station_z3d_dir,
                                                   **kw_dict)
        # make birrp dictionary
        birrp_dict = self.get_birrp_dict(z3d_df,
                                         df_list=df_list,
                                         use_blocks_dict=use_blocks_dict)

        # write script files for birrp
        sfn_list = self.write_script_files(birrp_dict,
                                           birrp_params_dict=birrp_param_dict,
                                           **kwargs)

        # run birrp
        self.run_birrp(sfn_list)

        # combine edi files
        comb_edi_fn = self.combine_edi_files(self.edi_fn, sr_dict)

        if comb_edi_fn is not None:
            self.edi_fn.append(comb_edi_fn)

        # plot the output
        if plot:
            r_plot = self.plot_responses()
        else:
            r_plot = None

        et = datetime.datetime.now()
        t_diff = (et - st).total_seconds()
        print('INFO: All processing took {0:02.0f}:{1:02.0f} minutes'.format(
              t_diff // 60, t_diff % 60))

        return r_plot, comb_edi_fn, z3d_df

    def combine_edi_files(self, edi_fn_list,
                          sr_dict={4096:(1000., 4),
                                   256:(3.99, .126),
                                   4:(.125, .00001)}):
        """
        combine the different edi files that are computed for each sampling
        rate. For now just a simple cutoff

        :param edi_fn_list: list of EDI files to combine
        :dtype edi_fn_list: list

        :param sr_dict: dictionary of frequency limits for each sampling rate
                        {4096:(1000., 4), 256:(3.99, .126), 4:(.125, .00001)}
        :dtype sr_dict: dictionary

        :return: full path to combined edi file
        :rtype: string

        :example: ::

            >>> zp_obj.combine_edi_files(["/home/mt01/TS/BF/4096/mt01edi",
                                          "/home/mt01/TS/BF/256/mt01.edi",
                                          "/home/mt01/TS/BF/4/mt01.edi"],
                                         sr_dict={4096:(1000., 50),
                                                  256:(49.99, .0625),
                                                  4:(.0624, .00001)})

        """

        if isinstance(edi_fn_list, str):
            print('WARNING: Only one edi file, skipping combining')
            return edi_fn_list

        if len(edi_fn_list) == 1:
            print('WARNING: Only one edi file, skipping combining')
            return edi_fn_list[0]

        data_arr = np.zeros(100,
                            dtype=[('freq', np.float),
                                   ('z', (np.complex, (2, 2))),
                                   ('z_err', (np.float, (2, 2))),
                                   ('tipper', (np.complex, (2, 2))),
                                   ('tipper_err', (np.float, (2, 2)))])

        count = 0
        for edi_fn in edi_fn_list:
            if not isinstance(edi_fn, Path):
                edi_fn = Path(edi_fn)
            # get sampling rate from file name
            sr_key = int(edi_fn.parts[-2])
            if sr_key in list(sr_dict.keys()):
                try:
                    edi_obj = mtedi.Edi(edi_fn)
                    # locate frequency range
                    f_index = np.where((edi_obj.Z.freq >= sr_dict[sr_key][1]) &
                                       (edi_obj.Z.freq <= sr_dict[sr_key][0]))


                    data_arr['freq'][count:count+len(f_index[0])] = edi_obj.Z.freq[f_index]
                    data_arr['z'][count:count+len(f_index[0])] = edi_obj.Z.z[f_index]
                    data_arr['z_err'][count:count+len(f_index[0])] = edi_obj.Z.z_err[f_index]
                    if edi_obj.Tipper.tipper is not None:
                        data_arr['tipper'][count:count+len(f_index[0])] = edi_obj.Tipper.tipper[f_index]
                        data_arr['tipper_err'][count:count+len(f_index[0])] = edi_obj.Tipper.tipper_err[f_index]

                    count += len(f_index[0])
                except IndexError:
                    pass
                    print('ERROR: Something went wrong with processing {0}'.format(edi_fn))

            else:
                pass
                print('WARNING: {0} was not in combining dictionary'.format(sr_key))

        # now replace
        data_arr = data_arr[np.nonzero(data_arr['freq'])]
        sort_index = np.argsort(data_arr['freq'])

        # check to see if the sorted indexes are descending or ascending,
        # make sure that frequency is descending
        if data_arr['freq'][0] > data_arr['freq'][1]:
            sort_index = sort_index[::-1]

        data_arr = data_arr[sort_index]
        new_z = mtedi.MTz.Z(data_arr['z'],
                            data_arr['z_err'],
                            data_arr['freq'])

        # check for all zeros in tipper, meaning there is only
        # one unique value
        if np.unique(data_arr['tipper']).size > 1:
            new_t = mtedi.MTz.Tipper(data_arr['tipper'],
                                     data_arr['tipper_err'],
                                     data_arr['freq'])

        else:
            new_t = mtedi.MTz.Tipper()

        edi_obj = mtedi.Edi(edi_fn_list[0])
        edi_obj.Z = new_z
        edi_obj.Tipper = new_t
        edi_obj.Data_sect.nfreq = new_z.z.shape[0]

        n_edi_fn = Path.joinpath(Path(self.station_ts_dir),
                                 '{0}_comb.edi'.format(Path(self.station_ts_dir).name))
        # n_edi_fn = os.path.join(self.station_z3d_dir,
        #                         '{0}_comb.edi'.format(os.path.basename(self.station_z3d_dir)))
        n_edi_fn = edi_obj.write_edi_file(new_edi_fn=n_edi_fn)

        return n_edi_fn

#==============================================================================

# this should capture all the print statements from zen
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = self._stringio = StringIO()
        # sys.stdout = self._stringio = BytesIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


#==============================================================================
def compute_mt_response(survey_dir, station='mt000', copy_date=None,
                        birrp_exe=r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts.exe",
                        ant_calibrations=r"c:\MT\Ant_calibrations",
                        process_df_list=[256],
                        use_blocks_dict={256:[0, 1]},
                        notch_dict={256:None}):
    """
    This code will down load Z3D files from a Zen that is in SD Mode,
    convert the Z3D files to ascii format, then process them for each
    sampling rate using Alan Chave's BIRRP code.  The outputs are then
    converted to .edi files and plotted.
    You need 2 things to run this code:
        * mtpy --> a Python package for MT and can be found at
                   https://github.com/geophysics/mtpy
        * BIRRP executable --> you can get this from Alan Chave at WHOI
                               if you are using it for non-commercial projects.
    ..note:: This code is quite specific to my setup, so let me know what
             doesn't work for you so I can generalize it.
    Arguments
    ----------------
        **survey_dir** : string
                         full path to the directory where you are storing
                         the station data.  ie. /home/MT/Survey_00
        **station** : string
                      name of the station you are down loading.
                      *default* is 'mt000'
        **copy_date** : string
                        copy all files on and after this date
                        format is YYYY-MM-DD
                        *default* is None, which copies all files on the SD
                        cards.
        **birrp_exe** : string
                        full path to the BIRRP executable on your machine
                        *default* is the location on my machine
        **ant_calibrations** : string
                               full path to a folder that contains the coil
                               calibration data.  These must be in seperate
                               .csv files for each coil named by corresponding
                               coil name. If you're coil is 2884, then you
                               need a calibration file named Ant2884_cal.csv
                               in which the data is freq,real,imaginary
        **process_df_list** : list
                              list of sampling rates to process
    Returns
    -----------------
        **rp_plot** : mtpy.imaging.plotnresponses object
                      ploting object of data, if you want to change how the
                      output is plot change the attributes of rp_plot
    Outputs
    -----------------
        **copy_from_sd.log** : file
                               contains information on how files were copied
                               from the SD cards.
        **processing.log** : file
                             a log file of how the program ran
        **survey_dir/station/TS** : directory
                                   contains the time series data in .ascii
                                   format
        **survey_dir/station/TS/BF** : directory
                                       contains the processing results from
                                       BIRRP for each sampling rate in the
                                       data in subsequent directories
        **survey_dir/station/TS/station.cfg** : file
                                                configuration file of the
                                                station parameters
    Example
    ------------------------
        >>> import zen_processing_data as zpd
        >>> zpd.compute_mt_response(r"/home/mt/survey_00",
                                    station='mt010',
                                    copy_date='2015-05-22',
                                    birrp_exe=r"/home/bin/birrp52.exe",
                                    ant_calibrations=r"/home/ant_calibrations",
                                    process_df_list=[1024, 256])
    """

    station_z3d_dir = os.path.join(survey_dir, station)

    st = datetime.datetime.now()
    #--> Copy data from files
    try:
        if copy_date is None:
            zen.copy_from_sd(station, save_path=survey_dir)
        else:
            zen.copy_from_sd(station, save_path=survey_dir,
                             copy_date=copy_date, copy_type='after')
    except IOError:
        pass
        print('ERROR: No files copied from SD cards')
        print('INFO: Looking in  {0} for Z3D files'.format(station_z3d_dir))

    #--> process data

    with Capturing() as output:
        z2edi = Z3D2EDI(station_z3d_dir)
        z2edi.birrp_exe = birrp_exe
        z2edi.calibration_path = ant_calibrations
        try:
            rp = z2edi.process_data(df_list=process_df_list,
                                    use_blocks_dict=use_blocks_dict)
        except mtex.MTpyError_inputarguments:
            print('WARNING: Data not good!! Did not produce a proper .edi file')
            et = datetime.datetime.now()
            print('--> took {0} seconds'.format((et-st).total_seconds()))
            rp = None

    #--> write log file
    log_fid = open(os.path.join(station_z3d_dir, 'Processing.log'), 'w')
    log_fid.write('\n'.join(output))
    log_fid.close()

    return rp