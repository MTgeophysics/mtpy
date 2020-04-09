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
        self.hz = s_df[s_df.component == 'hz'].coil_number.mode()[0]
        self.lat = s_df.latitude.median()
        self.location = 'Earth'
        self.lon = s_df.longitude.median()
        self.network = 'USGS'
        self.notes = 'Generic config file'
        self.sampling_interval = 'all'
        self.station = s_df.station.mode()[0]
        self.station_type = 'mt'

        rr_df = z3d_df[z3d_df.remote == 'True']
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
                               4: {'s_diff': 3 * 3600 * 4,
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
        return self._station_ts_dir
    
    @station_ts_dir.setter
    def station_ts_dir(self, station_ts_dir):
        if station_ts_dir in [None, 'None']:
            self._station_ts_dir = None
        else:
            self._station_ts_dir = Path(station_ts_dir)
            
    @property
    def rr_station_ts_dir(self):
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
        get coil calibrations
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

        :param station_z3d_dir: DESCRIPTION
        :type station_z3d_dir: TYPE

        :param rr_station_z3d_dir: DESCRIPTION, defaults to None
        :type rr_station_z3d_dir: TYPE, optional

        :param use_blocks_dict: DESCRIPTION, defaults to None
        :type use_blocks_dict: Dictionary of blocks to use
                               ex. {4096:[0, 2], 256:[2, 3]}
        :param overwrite: DESCRIPTION, defaults to False
        :type overwrite: TYPE, optional

        :param notch_dict: DESCRIPTION, defaults to None
        :type notch_dict: TYPE, optional

        :param combine: DESCRIPTION, defaults to True
        :type combine: TYPE, optional

        :param combine_sampling_rate: DESCRIPTION, defaults to 4
        :type combine_sampling_rate: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE

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

        :param entry: DESCRIPTION
        :type entry: TYPE
        :param nkip: DESCRIPTION
        :type nkip: TYPE
        :param nread: DESCRIPTION
        :type nread: TYPE
        :param rr_num: DESCRIPTION, defaults to 1
        :type rr_num: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

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

        :param entry_01: DESCRIPTION
        :type entry_01: TYPE
        :param entry_02: DESCRIPTION
        :type entry_02: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

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

        :param block_list: DESCRIPTION
        :type block_list: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        block_df = pd.DataFrame(block_list)
        block_df = block_df.infer_objects()
        block_df.start = pd.to_datetime(block_df.start)
        block_df.stop = pd.to_datetime(block_df.stop)

        return block_df

    def align_block(self, block_df):
        """

        :param block_df: DESCRIPTION
        :type block_df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

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

        :param station_df: DESCRIPTION
        :type station_df: TYPE
        :param remote_df: DESCRIPTION
        :type remote_df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

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
                        continue
                    t_diff = abs((rr_entry.start - start).total_seconds()) * sr
                    # check to see if the difference is within given tolerance
                    if t_diff < self._tol_dict[sr]['s_diff']:
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
        write script files in the new method
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
            station = np.unique(fn_arr[remotes]['station'])[0]

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
            birrp_script_obj.write_script_file(script_fn=b_script)

            # write a birrp parameter configuration file
            birrp_script_obj.write_config_file(bf_path.joinpath(station))
            script_fn_list.append(birrp_script_obj.script_fn)

        return script_fn_list

    def run_birrp(self, script_fn_list=None, birrp_exe=None):
        """
        run birrp given the specified files
        """

        if script_fn_list is None:
            raise IOError('Need to input a script file or list of script files')

        if birrp_exe is not None:
            self.birrp_exe = birrp_exe


        if type(script_fn_list) is list:
            self.edi_fn = []
            for script_fn in script_fn_list:
                birrp.run(self.birrp_exe, script_fn)

                output_path = os.path.dirname(script_fn)
                self.edi_fn.append(self.write_edi_file(output_path,
                                   survey_config_fn=self.survey_config_fn,
                                   birrp_config_fn=self.birrp_config_fn))
        elif type(script_fn_list) is str:
            birrp.run(self.birrp_exe, script_fn_list)

            output_path = os.path.dirname(script_fn_list)
            self.edi_fn = self.write_edi_file(output_path,
                                      survey_config_fn=self.survey_config_fn,
                                      birrp_config_fn=self.birrp_config_fn)

    def write_edi_file(self, birrp_output_path, survey_config_fn=None,
                       birrp_config_fn=None):
        """
        write an edi file from outputs of birrp
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
                        self.survey_config.read_survey_config_file()

        j2edi_obj = birrp.J2Edi(station=self.survey_config.station,
                                survey_config_fn=self.survey_config_fn,
                                birrp_dir=birrp_output_path,
                                birrp_config_fn=self.birrp_config_fn)

        edi_fn = j2edi_obj.write_edi_file()

        return edi_fn

    def plot_responses(self, edi_fn_list=None):
        """
        plot all the edi files that were created.
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
                              16:(.125, .0001)},
                     birrp_param_dict={}, **kwargs):
        """
        process_data is a convinience function that will process Z3D files
        and output an .edi file.  The workflow is to convert Z3D files to
        MTpy ascii format, apply a notch filter if desired, decimate if
        desired.  Do the same for a remote reference if specified.  Make a
        BIRRP script file for each sampling rate, process data, take outputs
        and write an edi file and plot resutl.  The whole process will take
        a few minutes so be patient.
        Arguments
        ---------------
            **df_list** : list of sampling rates
                          list of sampling rates (int) to process, options are
                          - 4096
                          - 1024
                          - 256
                          - 16
            **max_blocks** : int
                            maximum number of blocks to process, this cannot
                            exceed the number of blocks that BIRRP was
                            compiled with.  *default* is 3
            **notch_dict** : dict(sampling_rate: notches to filter)
                            dictionary of notches to filter for each sampling
                            rate.  Keys are sampling rates (int) and values
                            are dictionary of notches to filter
                            ex. {4096:{60, 120, 180}} to filter out 60 Hz noise
                            *default* is {} which filters out 60 Hz noise and
                            harmonics.
            **sr_dict** : dict(sampling_rate: (max_freq, max_freq))
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
            **use_blocks_dict** : dict(sampling_rate: list of blocks to use)
                                  dictionary with sampling rate as keys and
                                  list of which blocks to use as values
                                  ex. {4096: [0, 2]} to skip the second block
                                  use value 'all' to use all blocks
                                  *default* is
                                      {4096:'all', 256:'all', 16:'all'}
            **birrp_param_dict** : dict(birrp_param: value)
                                   dictionary of birrp parameters where the
                                   key is the birrp parameter and value is
                                   what you want that parameter to be.  See
                                   mtpy.processing.birrp.write_script_file
                                   or the BIRRP manual for details on birrp
                                   parameters
                                   ex. {'nar':5}
                                   *default* is {}
        Returns
        -------------
            **plot_obj** : plot_response object
            **comb_edi_fn** : string
                              full path to combined edi file
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
        rate.
        for now just a simple cutoff
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

def get_remote_reference_schedule(survey_path, plot=True):
    """
    Get a detailed list of which stations recorded at the same time, just from
    the raw Z3D files.
    Arguments
    ----------------
        **survey_path** : string
                          full path to the survey directory
        **plot** : [ True | False]
                   True to plot the schedules as line bars.
                   *default* is True
    Outputs
    -----------------
        **station_path/Remote_Reference_List.txt** : file with a list of
                                                     dates with station
                                                     and sampling rate
        **plot** if plot is True.
    :Example: ..
        >>> import mtpy.usgs.zen_processing as zp
        >>> zp.get_remote_reference_schedule(r"\home\MT_Data\Survey_01")
    """
    date_dict = {}
    rr_dict = {}
    for station in os.listdir(survey_path):
        station_path = os.path.join(survey_path, station)
        if os.path.isdir(station_path) is True:
            fn_list = [os.path.join(station_path, fn)
                       for fn in os.listdir(station_path)
                       if fn.lower().endswith('ex.z3d')]

            if len(fn_list) == 0:
                print('ERROR: No Z3D files found in folder: {0} '.format(station))
                continue

            station_date_arr = np.zeros(len(fn_list), dtype=[('df', np.float),
                                                             ('start_dt', '|S20')])
            for f_index, fn in enumerate(fn_list):
                zd = zen.Zen3D(fn)
                zd.read_all_info()
                station_date_arr[f_index]['df'] = zd.df
                station_date_arr[f_index]['start_dt'] = zd.zen_schedule
                try:
                    rr_dict[zd.zen_schedule]
                except KeyError:
                    rr_dict[zd.zen_schedule] = []

                rr_dict[zd.zen_schedule].append((station, zd.df))

            if len(np.nonzero(station_date_arr)[0]) != 0:
                date_dict[station] = station_date_arr[np.nonzero(station_date_arr['df'])]

    #---------------------------------------------
    # print out in a useful way

    lines = []
    for key in sorted(rr_dict.keys()):
        k_list = key.split(',')
        k_date = k_list[0]
        k_time = k_list[1]
        lines.append('Date: {0}, Time: {1}'.format(k_date, k_time))
        lines.append('-'*60)
        for k_tuple in rr_dict[key]:
            lines.append('\tStation: {0}, Sampling Rate: {1:.0f}'.format(k_tuple[0],
                         k_tuple[1]))

        lines.append('='*60)
    with open(os.path.join(survey_path, 'Remote_Reference_List.txt'), 'w') as fid:
        fid.write('\n'.join(lines))



    #-------------------------------------

    if plot is True:
        plt.rcParams['font.size'] = 14

        datetime_display = '%m-%d, %H:%M:%S'

        df_dict = {4096:(.7, .1, 0), 1024:(.5, .5, 0), 256:(0, .2, .8)}
        # plot the results is a compeling graph
        fig = plt.figure(2, [12, 10])
        ax = fig.add_subplot(1, 1, 1)
        y_labels = ['', '']

        for k_index, key in enumerate(sorted(date_dict.keys())):
            y_labels.append(key)
            for ii, k_arr in enumerate(date_dict[key]):
                x_date_0 = datetime.datetime.strptime(k_arr['start_dt'], 
                                                      datetime_fmt)
                try:
                    x_date_1 = datetime.datetime.strptime(date_dict[key]['start_dt'][ii+1],
                                                          datetime_fmt)
                except IndexError:
                    x_date_1 = x_date_0

                y_values = [k_index, k_index]

                l1, = ax.plot([x_date_0, x_date_1], y_values,
                              lw=8, color=df_dict[k_arr['df']])

        fig.autofmt_xdate(rotation=60)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(datetime_display))
        ax.xaxis.set_tick_params(width=2, size=5)

        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_ticklabels(y_labels)
        ax.yaxis.set_tick_params(width=2, size=5)
        ax.set_ylim(-1, len(y_labels)-2)

        ax.grid(which='major', linestyle='--', color=(.7, .7, .7))
        ax.set_axisbelow(True)

        l_4096 = plt.Line2D([0, 1], [0, 0], lw=8, color=df_dict[4096])
        l_1024 = plt.Line2D([0, 1], [0, 0], lw=8, color=df_dict[1024])
        l_256 = plt.Line2D([0, 1], [0, 0], lw=8, color=df_dict[256])

        fig.legend([l_4096, l_1024, l_256], ['4096', '1024', '256'], ncol=3,
                   loc='upper center')

        plt.show()

#==============================================================================

# this should capture all the print statements from zen
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
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
                        max_blocks=2,
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
        z2edi.coil_cal_path = ant_calibrations
        try:
            rp, comb_edi = z2edi.process_data(df_list=process_df_list,
                                              max_blocks=max_blocks,
                                              notch_dict=notch_dict)
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

def get_z3d_info(z3d_path):
    """
    get information on z3d files
    """
    if not isinstance(z3d_path, Path):
        z3d_path = Path(z3d_path)
    ### need to get all the files for one channel
    fn_dict = dict([(key, []) for key in ['ex', 'ey', 'hx', 'hy', 'hz']])
    ### get all z3d files within a given folder, will look through recursively
    fn_list = [fn_path for fn_path in z3d_path.rglob('*')
               if fn_path.suffix in ['.z3d', '.Z3D']]
    ### loop over files, read just the metadata and get important information
    for fn in fn_list:
        z_obj = zen.Zen3D(fn)
        z_obj.read_all_info()
        fn_dict[z_obj.component].append({'start':z_obj.zen_schedule.isoformat(),
                                         'df':z_obj.df,
                                         'fn':z_obj.fn})

    return fn_dict

def combine_z3d_files(z3d_path, new_sampling_rate=4, t_buffer=8*3600):
    """
    Combine all z3d files for a given station and given component for
    processing and getting the long period estimations.

    :param str z3d_path: full path to z3d files
    :param str component: component to combine
    :param int new_sampling_rate: new sampling rate of the data
    :param int t_buffer: buffer for the last time series, should be length
                         of longest schedule chunk
    """
    st = datetime.datetime.now()
    attr_list = ['station', 'channel_number', 'component', 'coordinate_system',
                 'dipole_length', 'azimuth', 'units', 'lat', 'lon', 'elev',
                 'datum', 'data_logger', 'instrument_id', 'calibration_fn',
                 'declination',  'fn', 'conversion', 'gain']

    fn_df = get_z3d_info(z3d_path)
    sv_path = Path.joinpath(Path(z3d_path), 'TS')
    if not sv_path.exists():
        sv_path.mkdir()

    return_fn_list = []
    for comp in ['ex', 'ey', 'hx', 'hy', 'hz']:
        if len(list(sv_path.glob('*combined_4.{0}'.format(comp)))) == 1:
            comp_fn = list(sv_path.glob('*combined_4.{0}'.format(comp)))[0]
            print('INFO: skipping {0} already exists'.format(comp_fn))
            return_fn_list.append(comp_fn)
            continue
        if len(fn_df[comp]) == 0:
            print('WARNING: Skipping {0} because no Z3D files found.'.format(comp))
            continue
        comp_df = pd.DataFrame(fn_df[comp])
        ### sort the data frame by date
        comp_df = comp_df.sort_values('start')

        ### get start date and end at last start date, get time difference
        start_dt = datetime.datetime.fromisoformat(comp_df.start.min())
        end_dt = datetime.datetime.fromisoformat(comp_df.start.max())
        t_diff = (end_dt - start_dt).total_seconds()

        ### make a new MTTS object that will have a length that is buffered
        ### at the end to make sure there is room for the data, will trimmed
        new_ts = mtts.MTTS()
        new_ts.ts = np.zeros(int((t_diff + t_buffer) * new_sampling_rate))
        new_ts.sampling_rate = new_sampling_rate
        new_ts.start_time_utc = start_dt

        ### make an attribute dictionary that can be used to fill in the new
        ### MTTS object
        attr_dict = dict([(key, []) for key in attr_list])
        ### loop over each z3d file for the given component
        for index, row in comp_df.iterrows():
            z_obj = zen.Zen3D(row['fn'])
            z_obj.read_z3d()
            t_obj = z_obj.ts_obj
            if t_obj.component in ['ex', 'ey']:
                t_obj.ts.data /= (t_obj.dipole_length/1000)
                t_obj.units = 'mV/km'
                print('INFO: Using scales {0} = {1} m'.format(t_obj.component,
                                                        t_obj.dipole_length))
            ### decimate to the required sampling rate
            t_obj.decimate(int(z_obj.df/new_sampling_rate))
            ### fill the new time series with the data at the appropriate times
            new_ts.ts.data[(new_ts.ts.index >= t_obj.ts.index[0]) &
                            (new_ts.ts.index <= t_obj.ts.index[-1])] = t_obj.ts.data
            ### get the end date as the last z3d file
            end_date = z_obj.ts_obj.ts.index[-1]
            ### fill attribute data frame
            for attr in attr_list:
                attr_dict[attr].append(getattr(t_obj, attr))

        ### need to trim the data
        new_ts.ts = new_ts.ts.data[(new_ts.ts.index >= start_dt) &
                                   (new_ts.ts.index <= end_date)].to_frame()

        ### fill gaps with forwards or backwards values, this seems to work
        ### better than interpolation and is faster than regression.
        ### The gaps should be max 13 seconds if everything went well
        new_ts.ts.data[new_ts.ts.data == 0] = np.nan
        new_ts.ts.data.fillna(method='ffill', inplace=True)

        ### fill the new MTTS with the appropriate metadata
        attr_df = pd.DataFrame(attr_dict)
        for attr in attr_list:
            try:
                attr_series = attr_df[attr][attr_df[attr] != 0]
                try:
                    setattr(new_ts, attr, attr_series.median())
                except TypeError:
                    setattr(new_ts, attr, attr_series.mode()[0])
            except ValueError:
                print('WARNING:  could not set {0}'.format(attr))


        ascii_fn = '{0}_combined_{1}.{2}'.format(new_ts.station,
                                                 int(new_ts.sampling_rate),
                                                 new_ts.component)

        sv_fn_ascii = sv_path.joinpath(ascii_fn)
        new_ts.write_ascii_file(sv_fn_ascii.as_posix())

        return_fn_list.append(sv_fn_ascii)

    et = datetime.datetime.now()
    compute_time = (et - st).total_seconds()
    print('INFO: Combining took {0:.2f} seconds'.format(compute_time))
    return return_fn_list
