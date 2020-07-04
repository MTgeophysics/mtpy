# -*- coding: utf-8 -*-
"""
.. module:: TS
   :synopsis: Deal with MT time series

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

#==============================================================================
# Imports
#==============================================================================
import os
import datetime
import dateutil

import numpy as np
import pandas as pd
import scipy.signal as signal

import mtpy.utils.gis_tools as gis_tools
import mtpy.processing.filter as mtfilter

import matplotlib.pyplot as plt

#==============================================================================

#==============================================================================
class MTTS(object):
    """
    MT time series object that will read/write data in different formats
    including hdf5, txt, miniseed.

    The foundations are based on Pandas Python package.

    The data are store in the variable ts, which is a pandas dataframe with
    the data in the column 'data'.  This way the data can be indexed as a
    numpy array:

        >>> MTTS.ts['data'][0:256]

        or

        >>> MTTS.ts.data[0:256]

    Also, the data can be indexed by time (note needs to be exact time):

        >>> MTTS.ts['2017-05-04 12:32:00.0078125':'2017-05-05 12:35:00]

    Input ts as a numpy.ndarray or Pandas DataFrame

    ==================== ==================================================
    Metadata              Description
    ==================== ==================================================
    azimuth              clockwise angle from coordinate system N (deg)
    calibration_fn       file name for calibration data
    component            component name [ 'ex' | 'ey' | 'hx' | 'hy' | 'hz']
    coordinate_system    [ geographic | geomagnetic ]
    datum                datum of geographic location ex. WGS84
    declination          geomagnetic declination (deg)
    dipole_length        length of dipole (m)
    data_logger          data logger type
    instrument_id        ID number of instrument for calibration
    lat                  latitude of station in decimal degrees
    lon                  longitude of station in decimal degrees
    n_samples            number of samples in time series
    sampling_rate        sampling rate in samples/second
    start_time_epoch_sec start time in epoch seconds
    start_time_utc       start time in UTC
    station              station name
    units                units of time series
    ==================== ==================================================

    .. note:: Currently only supports hdf5 and text files

    ======================= ===============================================
    Method                  Description
    ======================= ===============================================
    read_hdf5               read an hdf5 file
    write_hdf5              write an hdf5 file
    write_ascii_file        write an ascii file
    read_ascii_file         read an ascii file
    ======================= ===============================================


    :Example: ::

        >>> import mtpy.core.ts as ts
        >>> import numpy as np
        >>> MTTS = ts.MTTS()
        >>> MTTS.ts = np.random.randn(1024)
        >>> MTTS.station = 'test'
        >>> MTTS.lon = 30.00
        >>> MTTS.lat = -122.00
        >>> MTTS.component = 'HX'
        >>> MTTS.units = 'counts'
        >>> MTTS.write_hdf5(r"/home/test.h5")


    """

    def __init__(self, **kwargs):

        self.station = 'mt00'
        self.channel_number = 1
        self.component = None
        self.coordinate_system = 'geomagnetic'
        self.dipole_length = 0
        self.azimuth = 0
        self.units = 'mV'
        self._lat = 0.0
        self._lon = 0.0
        self._elev = 0.0
        self._n_samples = 0
        self._sampling_rate = 1
        self.datum = 'WGS84'
        self.data_logger = 'Zonge Zen'
        self.instrument_id = None
        self.calibration_fn = None
        self.declination = 0.0
        self._ts = pd.DataFrame({'data':[0]})
        self.fn = None
        self.conversion = None
        self.gain = None
        self._end_header_line = 0

        self._date_time_fmt = '%Y-%m-%d %H:%M:%S.%f'
        self._attr_list = ['station',
                           'sampling_rate',
                           'start_time_utc',
                           'stop_time_utc',
                           'n_samples',
                           'component',
                           'channel_number',
                           'coordinate_system',
                           'dipole_length',
                           'elev',
                           'azimuth',
                           'units',
                           'lat',
                           'lon',
                           'datum',
                           'data_logger',
                           'instrument_id',
                           'calibration_fn',
                           'declination',
                           'gain',
                           'conversion']

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    ###-------------------------------------------------------------
    ## make sure some attributes have the correct data type
    # make sure that the time series is a pandas data frame
    @property
    def ts(self):
        return self._ts

    @ts.setter
    def ts(self, ts_arr):
        """
        if setting ts with a pandas data frame, make sure the data is in a
        column name 'data'
        """
        if isinstance(ts_arr, np.ndarray):
            self._ts = pd.DataFrame({'data':ts_arr})
            self._set_dt_index(self.start_time_utc, self.sampling_rate)

        elif isinstance(ts_arr, pd.core.frame.DataFrame):
            try:
                ts_arr['data']
                self._ts = ts_arr
                self._set_dt_index(self.start_time_utc, self.sampling_rate)

            except AttributeError:
                raise MTTSError('Data frame needs to have a column named "data" '+\
                                   'where the time series data is stored')
        else:
            raise MTTSError('Data type {0} not supported'.format(type(ts_arr))+\
                              ', ts needs to be a numpy.ndarray or pandas DataFrame')

        self._n_samples = self.ts.data.size

    ##--> Latitude
    @property
    def lat(self):
        """Latitude in decimal degrees"""
        return self._lat

    @lat.setter
    def lat(self, latitude):
        """
        latitude in either decimal degrees or hh:mm:ss
        """
        self._lat = gis_tools.assert_lat_value(latitude)

    ##--> Longitude
    @property
    def lon(self):
        """Longitude in decimal degrees"""
        return self._lon

    @lon.setter
    def lon(self, longitude):
        """
        longitude in either decimal degrees or hh:mm:ss
        """
        self._lon = gis_tools.assert_lon_value(longitude)

    ##--> elevation
    @property
    def elev(self):
        """elevation in elevation units"""
        return self._elev

    @elev.setter
    def elev(self, elevation):
        """elevation in elevation units"""
        self._elev = gis_tools.assert_elevation_value(elevation)

    #--> number of samples just to make sure there is consistency
    @property
    def n_samples(self):
        """number of samples"""
        return int(self._n_samples)

    @n_samples.setter
    def n_samples(self, num_samples):
        """number of samples (int)"""
        self._n_samples = int(num_samples)

    def _check_for_index(self):
        """
        check to see if there is an index in the time series
        """
        if len(self._ts) > 0:
            return True
        else:
            return False

    #--> sampling rate
    @property
    def sampling_rate(self):
        """sampling rate in samples/second"""
        if self._check_for_index():
            if isinstance(self._ts.index[0], int):
                sr = self._sampling_rate
            else:
                sr = 1E9/self._ts.index[0].freq.nanos
        else:
            sr = self._sampling_rate
        return np.round(sr, 0)

    @sampling_rate.setter
    def sampling_rate(self, sampling_rate):
        """
        sampling rate in samples/second

        type float
        """
        try:
            sr = float(sampling_rate)
        except (ValueError):
            raise MTTSError("Input sampling rate should be a float not {0}".format(type(sampling_rate)))
        self._sampling_rate = sr
        if self._check_for_index():
            if isinstance(self._ts.index[0], int):
                return
            else:
                if 1E9/self._ts.index[0].freq.nanos == self._sampling_rate:
                    return
                else:
                    if self.start_time_utc is not None:
                        self._set_dt_index(self.start_time_utc,
                                           self._sampling_rate)


    ## set time and set index
    @property
    def start_time_utc(self):
        """start time in UTC given in time format"""
        if self._check_for_index():
            if isinstance(self._ts.index[0], int):
                return None
            else:
                return self._ts.index[0].isoformat()

    @start_time_utc.setter
    def start_time_utc(self, start_time):
        """
        start time of time series in UTC given in some format or a datetime
        object.

        Resets epoch seconds if the new value is not equivalent to previous
        value.

        Resets how the ts data frame is indexed, setting the starting time to
        the new start time.
        """

        if not isinstance(start_time, datetime.datetime):
            start_time = dateutil.parser.parse(start_time)

        if self._check_for_index():
            if isinstance(self._ts.index[0], int):
                self._set_dt_index(start_time.isoformat(),
                                   self._sampling_rate)
            else:
                if start_time.isoformat() == self.ts.index[0].isofromat():
                    return
                else:
                    self._set_dt_index(start_time.isoformat(),
                                       self._sampling_rate)

        # make a time series that the data can be indexed by
        else:
            raise MTTSError('No Data to set start time for, set data first')

    ## epoch seconds
    @property
    def start_time_epoch_sec(self):
        """start time in epoch seconds"""
        if self._check_for_index():
            if isinstance(self._ts.index[0], int):
                return None
            else:
                return self.ts.index[0].timestamp()
        else:
            return None

    @start_time_epoch_sec.setter
    def start_time_epoch_sec(self, epoch_sec):
        """
        start time in epoch seconds

        Resets start_time_utc if different

        Resets how ts data frame is indexed.
        """
        try:
            epoch_sec = float(epoch_sec)
        except ValueError:
            raise MTTSError("Need to input epoch_sec as a float not {0} {1".format(type(epoch_sec), self.fn_ascii))

        dt_struct = datetime.datetime.utcfromtimestamp(epoch_sec)
        # these should be self cosistent
        try:
            if self.ts.index[0] != dt_struct:
                self.start_time_utc = dt_struct
        except IndexError:
            print('setting time')
            self.start_time_utc = dt_struct

    @property
    def stop_time_epoch_sec(self):
        """
        End time in epoch seconds
        """
        if self._check_for_index():
            if isinstance(self._ts.index[-1], int):
                return None
            else:
                return self.ts.index[-1].timestamp()
        else:
            return None

    @property
    def stop_time_utc(self):
        """
        End time in UTC
        """
        if self._check_for_index():
            if isinstance(self._ts.index[-1], int):
                return None
            else:
                return self._ts.index[-1].isoformat()

    def _set_dt_index(self, start_time, sampling_rate):
        """
        get the date time index from the data

        :param start_time: start time in time format
        :type start_time: string
        """
        if len(self.ts) == 0:
            return

        if start_time is None:
            print('Start time is None, skipping calculating index')
            return
        dt_freq = '{0:.0f}N'.format(1./(sampling_rate)*1E9)

        dt_index = pd.date_range(start=start_time,
                                 periods=self.ts.data.size,
                                 freq=dt_freq)

        self.ts.index = dt_index
        print("   * Reset time seies index to start at {0}".format(start_time))

    def apply_addaptive_notch_filter(self, notches=None, notch_radius=0.5,
                                     freq_rad=0.5, rp=0.1):
        """
        apply notch filter to the data that finds the peak around each
        frequency.

        see mtpy.processing.filter.adaptive_notch_filter

        :param notch_dict: dictionary of filter parameters.
                           if an empty dictionary is input the filter looks
                           for 60 Hz and harmonics to filter out.
        :type notch_dict: dictionary

        """
        if notches is None:
            notches = list(np.arange(60, 1860, 120))

        kwargs = {'df':self.sampling_rate,
                  'notches':notches,
                  'notchradius':notch_radius,
                  'freqrad':freq_rad,
                  'rp':rp}

        ts, filt_list = mtfilter.adaptive_notch_filter(self.ts.data, **kwargs)

        self.ts.data = ts

        print('\t Filtered frequency with bandstop:')
        for ff in filt_list:
            try:
                print('\t\t{0:>6.5g} Hz  {1:>6.2f} db'.format(np.nan_to_num(ff[0]),
                                                             np.nan_to_num(ff[1])))
            except ValueError:
                pass

    # decimate data
    def decimate(self, dec_factor=1):
        """
        decimate the data by using scipy.signal.decimate

        :param dec_factor: decimation factor
        :type dec_factor: int

        * refills ts.data with decimated data and replaces sampling_rate

        """
        # be sure the decimation factor is an integer
        dec_factor = int(dec_factor)

        if dec_factor > 1:
            if dec_factor > 8:
                n_dec = np.log2(dec_factor)/np.log2(8)
                dec_list = [8] * int(n_dec) + [int(2**(3 * n_dec % 1))]
                decimated_data = signal.decimate(self.ts.data, 8, n=8)
                for dec in dec_list[1:]:
                    if dec == 0:
                        break
                    decimated_data = signal.decimate(decimated_data,
                                                     dec,
                                                     n=8)
            else:
                decimated_data = signal.decimate(self.ts.data, dec_factor, n=8)
            start_time = str(self.start_time_utc)
            self.ts = decimated_data
            self.sampling_rate /= float(dec_factor)
            self._set_dt_index(start_time, self.sampling_rate)

    def low_pass_filter(self, low_pass_freq=15, cutoff_freq=55):
        """
        low pass the data

        :param low_pass_freq: low pass corner in Hz
        :type low_pass_freq: float

        :param cutoff_freq: cut off frequency in Hz
        :type cutoff_freq: float

        * filters ts.data
        """

        self.ts = mtfilter.low_pass(self.ts.data,
                                    low_pass_freq,
                                    cutoff_freq,
                                    self.sampling_rate)

    ###------------------------------------------------------------------
    ### read and write file types
    def write_hdf5(self, fn_hdf5, compression_level=0, compression_lib='blosc'):
        """
        Write an hdf5 file with metadata using pandas to write the file.

        :param fn_hdf5: full path to hdf5 file, has .h5 extension
        :type fn_hdf5: string

        :param compression_level: compression level of file [ 0-9 ]
        :type compression_level: int

        :param compression_lib: compression library *default* is blosc
        :type compression_lib: string

        :returns: fn_hdf5

        .. seealso:: Pandas.HDf5Store
        """

        hdf5_store = pd.HDFStore(fn_hdf5, 'w',
                                 complevel=compression_level,
                                 complib=compression_lib)

        # might want to re index because the time string takes up a lot of
        # storage
        #df = pd.DataFrame({'data':self.ts.data})
        #df.index = np.arange(df.data.size)
        hdf5_store['time_series'] = self.ts

        # add in attributes
        for attr in self._attr_list:
            setattr(hdf5_store.get_storer('time_series').attrs,
                    attr,
                    getattr(self, attr))

        hdf5_store.flush()
        hdf5_store.close()

        return fn_hdf5

    def read_hdf5(self, fn_hdf5, compression_level=0, compression_lib='blosc'):
        """
        Read an hdf5 file with metadata using Pandas.

        :param fn_hdf5: full path to hdf5 file, has .h5 extension
        :type fn_hdf5: string

        :param compression_level: compression level of file [ 0-9 ]
        :type compression_level: int

        :param compression_lib: compression library *default* is blosc
        :type compression_lib: string

        :returns: fn_hdf5

        .. seealso:: Pandas.HDf5Store
        """
        self.fn = fn_hdf5

        hdf5_store = pd.HDFStore(fn_hdf5, 'r', complib=compression_lib)

        self.ts = hdf5_store['time_series']

        for attr in self._attr_list:
            value = getattr(hdf5_store.get_storer('time_series').attrs, attr)
            setattr(self, attr, value)

        hdf5_store.close()

    def write_ascii_file(self, fn_ascii, chunk_size=4096):
        """
        Write an ascii format file with metadata

        :param fn_ascii: full path to ascii file
        :type fn_ascii: string

        :param chunk_size: read in file by chunks for efficiency
        :type chunk_size: int

        :Example: ::

            >>> ts_obj.write_ascii_file(r"/home/ts/mt01.EX")

        """

        st = datetime.datetime.utcnow()

        # get the number of chunks to write
        chunks = int(self.ts.shape[0]/chunk_size)

        # make header lines
        header_lines = ['# *** MT time series text file for {0} ***'.format(self.station)]
        header_lines += ['# {0} = {1}'.format(attr, getattr(self, attr))
                        for attr in sorted(self._attr_list)]

        # write to file in chunks
        with open(fn_ascii, 'w') as fid:
            # write header lines first
            fid.write('\n'.join(header_lines))

            # write time series indicator
            fid.write('\n# *** time_series ***\n')

            # write in chunks
            for cc in range(chunks):
                # changing the dtype of the array is faster than making
                # a list of strings with 22 places to incorporate exponential
                # form
                ts_lines = np.array(self.ts.data[cc*chunk_size:(cc+1)*chunk_size],
                                    dtype='U22')

                fid.write('\n'.join(list(ts_lines)))
                # be sure to write a new line after each chunk otherwise
                # they run together
                fid.write('\n')

            # be sure to write the last little bit
            fid.write('\n'.join(list(np.array(self.ts.data[(cc+1)*chunk_size:],
                                              dtype='U22'))))


        # get an estimation of how long it took to write the file
        et = datetime.datetime.utcnow()
        time_diff = et-st

        print('--> Wrote {0}'.format(fn_ascii))
        print('    Took {0:.2f} seconds'.format(time_diff.seconds+time_diff.microseconds*1E-6))

    def read_ascii_header(self, fn_ascii):
        """
        Read an ascii metadata

        :param fn_ascii: full path to ascii file
        :type fn_ascii: string

        :Example: ::

            >>> ts_obj.read_ascii_header(r"/home/ts/mt01.EX")
        """
        if not os.path.isfile(fn_ascii):
            raise MTTSError('Could not find {0}, check path'.format(fn_ascii))
        self.fn = fn_ascii

        with open(self.fn, 'r') as fid:
            line = fid.readline()
            count = 0
            attr_dict = {}
            find_old = False
            while line.find('#') == 0:
                line_list = line[1:].strip().split('=')
                if len(line_list) == 2:
                    key = line_list[0].strip()
                    try:
                        value = float(line_list[1].strip())
                    except ValueError:
                        value = line_list[1].strip()
                    attr_dict[key] = value

                # skip the header lines
                elif line.find('***') > 0:
                    pass
                else:
                    line_list = line[1:].strip().split()
                    if len(line_list) == 9:
                        print('Reading old MT TS format')
                        find_old = True
                        self.station = line_list[0]
                        self.component = line_list[1].lower()
                        self.sampling_rate = float(line_list[2])
                        self.start_time_epoch_sec = float(line_list[3])
                        # skip setting number of samples
                        self.units = line_list[5]
                        self.lat = float(line_list[6])
                        self.lon = float(line_list[7])
                        self.elev = float(line_list[8])
                count += 1
                line = fid.readline()
        self._end_header_line = count

        if find_old:
            return

        # make a dummy time series to get end time etc
        self.ts = np.zeros(int(attr_dict['n_samples']))
        for key, value in attr_dict.items():
            try:
                setattr(self, key, value)
            except AttributeError:
                if key not in ['n_samples',
                               'start_time_epoch_sec',
                               'start_time_utc',
                               'stop_time_utc']:
                    print('Could not set {0} to {1}'.format(key, value))

    def read_ascii(self, fn_ascii):
        """
        Read an ascii format file with metadata

        :param fn_ascii: full path to ascii file
        :type fn_ascii: string

        :Example: ::

            >>> ts_obj.read_ascii(r"/home/ts/mt01.EX")
        """

        self.read_ascii_header(fn_ascii)

        start_time = self.start_time_utc

        self.ts = pd.read_csv(self.fn,
                              sep='\n',
                              skiprows=self._end_header_line,
                              memory_map=True,
                              names=['data'])
        self._set_dt_index(start_time, self.sampling_rate)
        print(self.start_time_utc)
        print('Read in {0}'.format(self.fn))

    def plot_spectra(self, spectra_type='welch', **kwargs):
        """
        Plot spectra using the spectral type

        .. note:: Only spectral type supported is welch

        :param spectra_type: [ 'welch' ]
        :type spectral_type: string

        :Example: ::

            >>> ts_obj = mtts.MTTS()
            >>> ts_obj.read_hdf5(r"/home/MT/mt01.h5")
            >>> ts_obj.plot_spectra()

        """

        s = Spectra()
        param_dict = {}
        if spectra_type == 'welch':
            param_dict['fs'] = kwargs.pop('sampling_rate',
                                                       self.sampling_rate)
            param_dict['nperseg'] = kwargs.pop('nperseg', 2**12)
            s.compute_spectra(self.ts.data, spectra_type, **param_dict)

#==============================================================================
# Error classes
#==============================================================================
class MTTSError(Exception):
    pass

#==============================================================================
#  spectra
#==============================================================================
class Spectra(object):
    """
    compute spectra of time series
    """

    def __init__(self, **kwargs):
        self.spectra_type = 'welch'

    def compute_spectra(self, data, spectra_type, **kwargs):
        """
        compute spectra according to input type
        """

        if spectra_type.lower() == 'welch':
            self.welch_method(data, **kwargs)

    def welch_method(self, data, plot=True, **kwargs):
        """
        Compute the spectra using the Welch method, which is an average
        spectra of the data.  Computes short time window of length nperseg and
        averages them to reduce noise.

        Arguments
        ------------

        """

        f, p = signal.welch(data, **kwargs)

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.loglog(f, p, lw=1.5)
            ax.set_xlabel('Frequency (Hz)',
                          fontdict={'size':10, 'weight':'bold'})
            ax.set_ylabel('Power (dB)',
                          fontdict={'size':10, 'weight':'bold'})
            ax.axis('tight')
            ax.grid(which='both')

            plt.show()

        return f, p




