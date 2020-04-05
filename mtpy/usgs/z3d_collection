#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z3DCollection
=================

An object to hold Z3D file information to make processing easier.


Created on Sat Apr  4 12:40:40 2020

@author: peacock
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd
from pathlib import Path

from mtpy.usgs import zen
from mtpy.core import ts as mtts

# =============================================================================
# Collection of Z3D Files
# =============================================================================
class Z3DCollection(object):
    """
    An object to deal with a collection of Z3D files. Metadata and information
    are contained with in Pandas DataFrames for easy searching.

    The main data frame should have keys and data types as follows:

        ================ ========= ===========================================
        Key              Type      Description
        ================ ========= ===========================================
        'station'        string    station name
        'start'          string    start time in isoformat (converted to
                                   pandas.Timestamp internally)
        'stop'           string    stop time in isoformat (converted to
                                   pandas.Timestamp internally)
        'sampling_rate'  float     samping rate in samples/second
        'component'      string    [ ex | ey | hx | hy | hz ]
        'fn_z3d'         string    full path to z3d file
        'azimuth'        float     azimuth of sensor (degrees)
        'dipole_length'  float     dipole length (meters)
        'coil_number'    string    mag sensor serial number
        'latitude'       float     latitude (decimal degrees)
        'longitude'      float     longitude (decimal degrees)
        'elevation'      float     elevation (meters)
        'n_samples'      integer   number of samples in time series
        'fn_ascii'       string    full path to ascii file
        ================ ========= ===========================================

    :Example:

        >>> from mtpy.usgs import z3d_collections as zc
        >>> zc_obj = zc.Z3DCollection(r"/home/z3d_files")
        >>> z3d_fn_list = zc.get_z3d_fn_list()
        >>> z3d_df = zc.get_z3d_info(z3d_fn_list)
        >>> z3d_df.to_csv(r"/home/z3d_files/z3d_info.csv")
        >>> z3d_df_final = zc.convert_to_mtts(z3d_df)
        >>> z3d_df_final.to_csv(r"/home/z3d_files/station_info.csv")
    """

    def __init__(self, z3d_path=None):
        self._z3d_path = None
        self.z3d_path = z3d_path
        self.ts_path = None

        self._keys_dict = {'station': 'station',
                           'start': 'start',
                           'stop': 'stop',
                           'sampling_rate': 'df',
                           'component': 'component',
                           'fn_z3d':'fn',
                           'azimuth': 'azimuth',
                           'dipole_length': 'dipole_len',
                           'coil_number': 'coil_num',
                           'latitude': 'lat',
                           'longitude': 'lon',
                           'elevation': 'elev',
                           'n_samples': 'n_samples',
                           'fn_ascii': 'fn_ascii'}

        self._dtypes = {'station': str,
                       'start': str,
                       'stop': str,
                       'sampling_rate': float,
                       'component': str,
                       'fn_z3d': str,
                       'azimuth': float,
                       'dipole_length': float,
                       'coil_number': str,
                       'latitude': float,
                       'longitude': float,
                       'elevation': float,
                       'n_samples': int,
                       'fn_ascii': str}
    @property
    def z3d_path(self):
        """
        Path object to z3d directory
        """
        return self._z3d_path

    @z3d_path.setter
    def z3d_path(self, z3d_path):
        """
        :param z3d_path: path to z3d files
        :type z3d_path: string or Path object

        sets z3d_path as a Path object
        """
        if z3d_path is None:
            return
        if not isinstance(z3d_path, Path):
            z3d_path = Path(z3d_path)
        self._z3d_path = z3d_path

    def get_z3d_fn_list(self, z3d_path=None):
        """
        Get a list of z3d files in a given directory

        :param z3d_path: Path to z3d files
        :type z3d_path: [ str | pathlib.Path object]
        :return: list of z3d files
        :rtype: list

        :Example: ::

            >>> zc = Z3DCollection()
            >>> z3d_fn_list = zc.get_z3d_fn_list(z3d_path=r"/home/z3d_files")
        """
        if z3d_path is not None:
            self.z3d_path = z3d_path
        if not z3d_path.exists():
            raise ValueError('Error: Directory {0} does not exist'.format(
                             self.z3d_path))

        z3d_list = [fn_path for fn_path in z3d_path.rglob('*')
                    if fn_path.suffix in ['.z3d', '.Z3D']]
        return z3d_list

    def get_z3d_info(self, z3d_fn_list):
        """
        Get general z3d information and put information in a dataframe

        :param z3d_fn_list: List of files Paths to z3d files
        :type z3d_fn_list: list

        :return: Dataframe of z3d information
        :rtype: Pandas.DataFrame

        :Example: ::

            >>> zc_obj = zc.Z3DCollection(r"/home/z3d_files")
            >>> z3d_fn_list = zc.get_z3d_fn_list()
            >>> z3d_df = zc.get_z3d_info(z3d_fn_list)
            >>> # write dataframe to a file to use later
            >>> z3d_df.to_csv(r"/home/z3d_files/z3d_info.csv")

        """

        z3d_info_list = []
        for z3d_fn in z3d_fn_list:
            z3d_obj = zen.Zen3D(z3d_fn)
            z3d_obj.read_all_info()
            z3d_obj.start = z3d_obj.zen_schedule.isoformat()
            # set some attributes to null to fill later
            z3d_obj.stop = None
            z3d_obj.n_samples = 0
            z3d_obj.fn_ascii = None
            # make a dictionary of values to put into data frame
            entry = dict([(key, getattr(z3d_obj, value)) for key, value in
                          self._keys_dict.items()])
            # add a block entry that will be filled later
            entry['block'] = 0
            z3d_info_list.append(entry)

        # make pandas dataframe and set data types
        z3d_df = pd.DataFrame(z3d_info_list)
        z3d_df = z3d_df.astype(self._dtypes)
        z3d_df.start = pd.to_datetime(z3d_df.start, errors='coerce')
        z3d_df.stop = pd.to_datetime(z3d_df.stop, errors='coerce')

        # assign block numbers
        for sr in z3d_df.sampling_rate.unique():
            starts = sorted(z3d_df[z3d_df.sampling_rate == sr].start.unique())
            for block_num, start in enumerate(starts):
                z3d_df.loc[(z3d_df.start == start), 'block'] = block_num

        return z3d_df

    def to_csv(self, z3d_df, fn_basename=None):
        """
        write data frame to a csv file, sort of redundant, a helper function

        :param z3d_df: DataFrame holding info on Z3D Files
        :type z3d_df: pandas.DataFrame
        :param save_basename: basename of file,
                              defaults to 'staion_z3d_info.csv'
        :type save_basename: string, optional
        :return: full path to saved file
        :rtype: string

        """
        sv_path = Path(z3d_df.fn_z3d[0]).parent
        if fn_basename is None:
            fn_basename = '{0}_z3d_info.csv'.format(z3d_df.station.unique()[0])
        sv_fn = sv_path.joinpath(fn_basename)
        z3d_df.to_csv(sv_fn)

        return sv_fn

    def from_csv(self, df_fn):
        """
        read in a csv file holding information on z3d files, sets data types

        :param df_fn: full path to csv file
        :type df_fn: string
        :return: DataFrame of the csv file
        :rtype: pandas.DataFrame

        """

        z3d_df = pd.read_csv(df_fn)
        z3d_df = z3d_df.astype(self._dtypes)
        z3d_df.start = pd.to_datetime(z3d_df.start, errors='coerce')
        z3d_df.stop = pd.to_datetime(z3d_df.stop, errors='coerce')

        return z3d_df

    def _make_block_dict(self, z3d_df):
        """

        """
        block_dict = {}
        for sr in z3d_df.sampling_rate.unique():
            block_dict[sr] = list(z3d_df[z3d_df.sampling_rate == sr].block.unique())

        return block_dict

    def convert_to_mtts(self, z3d_df, block_dict=None, notch_dict=None,
                        dec=1, overwrite=False, combine=True,
                        combine_sampling_rate=4):
        """
        Convert z3d files to MTTS objects and write ascii files if they do
        not already exist.

        :param z3d_df: dataframe holding information about z3d files see help
                       for more information on data frame structure
        :type z3d_df: pandas.DataFrame
        :param block_dict: dictionary of blocks to use. Has keys of sample
                           rate and values of a list of blocks to use,
                           defaults to None
        :type block_dict: dictionary, optional
        :param notch_dict: dictionary of notches to apply for each sample rate
                           keys are sampling rates, values are notch dice
                           defaults to None, if an empy dictionary is used
                           then notches at 60 Hz and harmonics is applied
        :type notch_dict: dictionary, optional
        :param dec: decimation frequency, defaults to 1
        :type dec: int, optional
        :return: dataframe filled with timeseries information
        :rtype: pandas.DataFrame

        .. todo:: Add examples of notch dict, block dict

        """
        # if the block dictionary is empty make one that covers all files
        if block_dict is None:
            block_dict = self._make_block_dict(z3d_df)

        # loop over each entry in the data frame
        for entry in z3d_df.itertuples():
            # check to see if the file already exists
            if entry.fn_ascii is not None:
                # if the file exists and no overwrite get information and skip
                if Path(entry.fn_ascii).exists() and overwrite is False:
                    print('INFO: Skipping {0}'.format(entry.fn_ascii))
                    ts_obj = mtts.MTTS()
                    ts_obj.read_ascii_header(entry.fn_ascii)
                    z3d_df.at[entry.Index, 'stop'] = pd.Timestamp(ts_obj.stop_time_utc)
                    z3d_df.at[entry.Index, 'n_samples'] = ts_obj.n_samples
                    z3d_df.at[entry.Index, 'start'] = pd.Timestamp(ts_obj.start_time_utc)
                    z3d_df.at[entry.Index, 'fn_ascii'] = ts_obj.fn
                    continue

            # check if the block should be used
            if entry.block in block_dict[entry.sampling_rate]:
                z3d_obj = zen.Zen3D(entry.fn_z3d)
                z3d_obj.read_z3d()
                ts_obj = z3d_obj.ts_obj

                # write mtpy mt file
                z3d_obj.write_ascii_mt_file(notch_dict=notch_dict, dec=dec)

                # get information from time series and fill data frame
                z3d_df.at[entry.Index, 'stop'] = pd.Timestamp(ts_obj.stop_time_utc)
                z3d_df.at[entry.Index, 'n_samples'] = ts_obj.n_samples
                z3d_df.at[entry.Index, 'start'] = pd.Timestamp(ts_obj.start_time_utc)
                z3d_df.at[entry.Index, 'fn_ascii'] = z3d_obj.fn_mt_ascii

        if combine:
            csr = combine_sampling_rate
            z3d_df = self.combine_z3d_files(z3d_df, new_sampling_rate=csr)

        return z3d_df

    def combine_z3d_files(self, z3d_df, new_sampling_rate=4, t_buffer=3600):
        """
        Combine all z3d files for a given station and given component for
        processing to get long period estimations.

        :param str z3d_path: full path to z3d files
        :param str component: component to combine
        :param int new_sampling_rate: new sampling rate of the data
        :param int t_buffer: buffer for the last time series, should be length
                             of longest schedule chunk
        """
        attr_list = ['station', 'channel_number', 'component',
                     'coordinate_system', 'dipole_length', 'azimuth', 'units',
                     'lat', 'lon', 'elev', 'datum', 'data_logger',
                     'instrument_id', 'calibration_fn', 'declination',
                     'fn', 'conversion', 'gain']
        sv_path = Path(z3d_df.fn_ascii[0]).parent

        combined_entries = []
        for comp in ['ex', 'ey', 'hx', 'hy', 'hz']:
            if len(list(sv_path.glob('*combined_4.{0}'.format(comp)))) == 1:
                comp_fn = list(sv_path.glob('*combined_4.{0}'.format(comp)))[0]
                print('INFO: skipping {0} already exists'.format(comp_fn))
                continue
            # sort out files for the given component
            comp_df = z3d_df[z3d_df.component == comp].copy()
            if len(comp_df) == 0:
                print('Warning: Skipping {0} because no Z3D files found.'.format(comp))
                continue

            # sort the data frame by date
            comp_df = comp_df.sort_values('start')

            # get start date and end at last start date, get time difference
            start_dt = comp_df.start.min()
            end_dt = comp_df.stop.max()
            t_diff = int((end_dt - start_dt).total_seconds())

            # make a new MTTS object that will have a length that is buffered
            # at the end to make sure there is room for the data, will trimmed
            new_ts = mtts.MTTS()
            new_ts.ts = np.zeros(int((t_diff + t_buffer) * new_sampling_rate))
            new_ts.sampling_rate = new_sampling_rate
            new_ts.start_time_utc = start_dt.isoformat()

            # make an attribute dictionary that can be used to fill in the new
            # MTTS object
            attr_dict = dict([(key, []) for key in attr_list])
            # loop over each z3d file for the given component
            for row in comp_df.itertuples():
                z_obj = zen.Zen3D(row.fn_z3d)
                z_obj.read_z3d()
                t_obj = z_obj.ts_obj
                if row.component in ['ex', 'ey']:
                    t_obj.ts.data /= (row.dipole_length/1000)
                    t_obj.units = 'mV/km'
                    print('Using scales {0} = {1} m'.format(row.component,
                                                            row.dipole_length))
                # decimate to the required sampling rate
                t_obj.decimate(int(z_obj.df/new_sampling_rate))
                # fill the new time series with the data at the appropriate times
                new_ts.ts.data[(new_ts.ts.index >= t_obj.ts.index[0]) &
                                (new_ts.ts.index <= t_obj.ts.index[-1])] = t_obj.ts.data
                # get the end date as the last z3d file
                end_date = z_obj.ts_obj.ts.index[-1]
                # fill attribute data frame
                for attr in attr_list:
                    attr_dict[attr].append(getattr(t_obj, attr))

            # need to trim the data
            new_ts.ts = new_ts.ts.data[(new_ts.ts.index >= start_dt) &
                                       (new_ts.ts.index <= end_date)].to_frame()

            # fill gaps with forwards or backwards values, this seems to work
            # better than interpolation and is faster than regression.
            # The gaps should be max 13 seconds if everything went well
            new_ts.ts.data[new_ts.ts.data == 0] = np.nan
            new_ts.ts.data.fillna(method='ffill', inplace=True)

            # fill the new MTTS with the appropriate metadata
            attr_df = pd.DataFrame(attr_dict)
            for attr in attr_list:
                try:
                    attr_series = attr_df[attr][attr_df[attr] != 0]
                    try:
                        setattr(new_ts, attr, attr_series.median())
                    except TypeError:
                        setattr(new_ts, attr, attr_series.mode()[0])
                except ValueError:
                    print('Warning: could not set {0}'.format(attr))

            ascii_fn = '{0}_combined_{1}.{2}'.format(new_ts.station,
                                                     int(new_ts.sampling_rate),
                                                     new_ts.component)

            sv_fn_ascii = sv_path.joinpath(ascii_fn)
            new_ts.write_ascii_file(sv_fn_ascii.as_posix())

            entry = {'station': new_ts.station,
                     'start': new_ts.start_time_utc,
                     'stop': new_ts.stop_time_utc,
                     'sampling_rate': new_ts.sampling_rate,
                     'component': new_ts.component,
                     'fn_z3d': None,
                     'aximuth': new_ts.azimuth,
                     'dipole_length': new_ts.dipole_length,
                     'coil_number': new_ts.instrument_id,
                     'latitude': new_ts.lat,
                     'longitude': new_ts.lon,
                     'elevation': new_ts.elev,
                     'n_samples': new_ts.n_samples,
                     'fn_ascii': sv_fn_ascii}

            combined_entries.append(entry)
        
        # make data frame of combined information and append to existing
        # data frame
        combined_df = pd.DataFrame(combined_entries)
        full_df = z3d_df.append(combined_df)

        return full_df
