# -*- coding: utf-8 -*-
"""
Created on Wed May 13 19:10:46 2020

@author: jpeacock
"""

import datetime
import logging
import numpy as np
from copy import deepcopy

from dateutil import parser as dtparser
from dateutil.tz.tz import tzutc

from mtpy.utils.exceptions import MTTimeError

# =============================================================================
#  Get leap seconds
# =============================================================================
leap_second_dict = {
    0: {"min": datetime.date(1980, 1, 1), "max": datetime.date(1981, 7, 1)},
    1: {"min": datetime.date(1981, 7, 1), "max": datetime.date(1982, 7, 1)},
    2: {"min": datetime.date(1982, 7, 1), "max": datetime.date(1983, 7, 1)},
    3: {"min": datetime.date(1983, 7, 1), "max": datetime.date(1985, 7, 1)},
    4: {"min": datetime.date(1985, 7, 1), "max": datetime.date(1988, 1, 1)},
    5: {"min": datetime.date(1988, 1, 1), "max": datetime.date(1990, 1, 1)},
    6: {"min": datetime.date(1990, 1, 1), "max": datetime.date(1991, 1, 1)},
    7: {"min": datetime.date(1991, 1, 1), "max": datetime.date(1992, 7, 1)},
    8: {"min": datetime.date(1992, 7, 1), "max": datetime.date(1993, 7, 1)},
    9: {"min": datetime.date(1993, 7, 1), "max": datetime.date(1994, 7, 1)},
    10: {"min": datetime.date(1994, 7, 1), "max": datetime.date(1996, 1, 1)},
    11: {"min": datetime.date(1996, 1, 1), "max": datetime.date(1997, 7, 1)},
    12: {"min": datetime.date(1997, 7, 1), "max": datetime.date(1999, 1, 1)},
    13: {"min": datetime.date(1999, 1, 1), "max": datetime.date(2006, 1, 1)},
    14: {"min": datetime.date(2006, 1, 1), "max": datetime.date(2009, 1, 1)},
    15: {"min": datetime.date(2009, 1, 1), "max": datetime.date(2012, 6, 30)},
    16: {"min": datetime.date(2012, 7, 1), "max": datetime.date(2015, 7, 1)},
    17: {"min": datetime.date(2015, 7, 1), "max": datetime.date(2017, 1, 1)},
    18: {"min": datetime.date(2017, 1, 1), "max": datetime.date(2021, 7, 1)},
}


def calculate_leap_seconds(year, month, day):
    """
    get the leap seconds for the given year to convert GPS time to UTC time
    
    .. note:: GPS time started in 1980
    
    .. note:: GPS time is leap seconds ahead of UTC time, therefore you
              should subtract leap seconds from GPS time to get UTC time.
              
    =========================== ===============================================
    Date Range                  Leap Seconds
    =========================== ===============================================
    1981-07-01 - 1982-07-01     1
    1982-07-01 - 1983-07-01     2
    1983-07-01 - 1985-07-01     3
    1985-07-01 - 1988-01-01     4
    1988-01-01 - 1990-01-01     5
    1990-01-01 - 1991-01-01     6
    1991-01-01 - 1992-07-01     7
    1992-07-01 - 1993-07-01     8
    1993-07-01 - 1994-07-01     9
    1994-07-01 - 1996-01-01     10
    1996-01-01 - 1997-07-01     11
    1997-07-01 - 1999-01-01     12
    1999-01-01 - 2006-01-01     13
    2006-01-01 - 2009-01-01     14
    2009-01-01 - 2012-07-01     15
    2012-07-01 - 2015-07-01     16
    2015-07-01 - 2017-01-01     17
    2017-01-01 - ????-??-??     18
    =========================== ===============================================
    
    """

    # make the date a datetime object, easier to test
    given_date = datetime.date(int(year), int(month), int(day))

    # made an executive decision that the date can be equal to the min, but
    # not the max, otherwise get an error.
    for leap_key in sorted(leap_second_dict.keys()):
        if (
            given_date < leap_second_dict[leap_key]["max"]
            and given_date >= leap_second_dict[leap_key]["min"]
        ):
            return int(leap_key)

    return None


# ==============================================================================
# convenience date-time container
# ==============================================================================
class MTime:
    """
    Date and Time container based on datetime and dateutil.parsers

    Will read in a string or a epoch seconds into a datetime.datetime object
    assuming the time zone is UTC.  If UTC is not the timezone you need to
    correct the time before inputing.  Use datetime.timezone to shift time.

    Outputs can be an ISO formatted string YYYY-MM-DDThh:mm:ss.ssssss+00:00:

        >>> t = MTtime()
        >>> t.iso_str
        '1980-01-01T00:00:00+00:00'

    .. note:: if microseconds are 0 they are omitted.

    or Epoch seconds (float):

        >>> t.epoch_seconds
        315532800.0


    Convenience getters/setters are provided as properties for the different
    parts of time.

        >>> t = MTtime()
        >>> t.year = 2020
        >>> t.year
        2020

    """

    def __init__(self, time=None, gps_time=False):

        self.logger = logging.getLogger(
            "{0}.{1}".format(__name__, self.__class__.__name__)
        )
        self.dt_object = self.now()

        if time is not None:
            if isinstance(time, str):
                self.logger.debug("Input time is a string, will be parsed")
                self.from_str(time)

            elif isinstance(time, (int, float)):
                self.logger.debug(
                    "Input time is a number, assuming epoch " + "seconds in UTC"
                )
                self.epoch_seconds = time
            elif isinstance(time, (np.datetime64)):
                self.logger.debug(
                    "Input time is a np.datetime64 "
                    + "dt_object set to datetime64.tolist()."
                )
                self.dt_object = self.validate_tzinfo(time.tolist())

            elif isinstance(time, (datetime.datetime)):
                self.logger.debug(
                    "Input time is a np.datetime64 "
                    + "dt_object set to datetime64.tolist()."
                )
                self.dt_object = self.validate_tzinfo(time)

            else:
                msg = "input time must be a string, float, or int, not {0}"
                self.logger.error(msg.format(type(time)))

        else:
            self.logger.debug(
                "Initiated with None, dt_object is set to "
                + "default time 1980-01-01 00:00:00"
            )
            self.from_str("1980-01-01 00:00:00")

        if gps_time:
            leap_seconds = calculate_leap_seconds(self.year, self.month, self.day)
            self.logger.debug(
                f"Converting GPS time to UTC with {leap_seconds} leap seconds"
            )
            self.dt_object -= datetime.timedelta(seconds=leap_seconds)

    def __str__(self):
        return self.iso_str

    def __repr__(self):
        return self.iso_str

    def __eq__(self, other):
        if other is None:
            self.logger.warning("Other is None, cannot compare")
            return False
        if isinstance(other, datetime.datetime):
            return bool(self.dt_object == other)

        elif isinstance(other, MTime):
            return bool(self.dt_object == other.dt_object)

        elif isinstance(other, str):
            return bool(self.iso_str == other)

        elif isinstance(other, (int, float)):
            return bool(self.epoch_seconds == float(other))
        else:
            msg = "Cannot compare {0} of type {1} with MTime Object".format(
                other, type(other)
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, datetime.datetime):
            return bool(self.dt_object < other)

        elif isinstance(other, MTime):
            return bool(self.dt_object < other.dt_object)

        elif isinstance(other, str):
            return bool(self.iso_str < other)

        elif isinstance(other, (int, float)):
            return bool(self.epoch_seconds < float(other))

        else:
            msg = "Cannot compare {0} of type {1} with MTime Object".format(
                other, type(other)
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

    def __le__(self, other):
        if isinstance(other, datetime.datetime):
            return bool(self.dt_object <= other)

        elif isinstance(other, MTime):
            return bool(self.dt_object <= other.dt_object)

        elif isinstance(other, str):
            return bool(self.iso_str <= other)

        elif isinstance(other, (int, float)):
            return bool(self.epoch_seconds <= float(other))
        else:
            msg = "Cannot compare {0} of type {1} with MTime Object".format(
                other, type(other)
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

    def __gt__(self, other):
        return not self.__lt__(other)

    def __ge__(self, other):
        if isinstance(other, datetime.datetime):
            return bool(self.dt_object >= other)

        elif isinstance(other, MTime):
            return bool(self.dt_object >= other.dt_object)

        elif isinstance(other, str):
            return bool(self.iso_str >= other)

        elif isinstance(other, (int, float)):
            return bool(self.epoch_seconds >= float(other))
        else:
            msg = "Cannot compare {0} of type {1} with MTime Object".format(
                other, type(other)
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

    def __add__(self, other):
        """
        add time only using datetime.timedelta, otherwise it does not make 
        sense to at 2 times together.  
        
        """
        if isinstance(other, (int, float)):
            other = datetime.timedelta(seconds=other)
            self.logger.debug("Assuming other time is in seconds")

        if not isinstance(other, (datetime.timedelta)):
            msg = (
                "Adding times does not make sense, must use "
                + "datetime.timedelta to add time. \n"
                + "\t>>> add_time = datetime.timedelta(seconds=10) \n"
                + "\t>>> mtime_obj + add_time"
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

        return MTime(self.dt_object + other)

    def __sub__(self, other):
        """
        Get the time difference between to times in seconds.
        
        :param other: other time value
        :type other: [ str | float | int | datetime.datetime | np.datetime64 ]
        :return: time difference in seconds
        :rtype: float

        """

        if isinstance(other, type(self)):
            other_seconds = other.epoch_seconds
        elif isinstance(other, (str)):
            other_dt = self.validate_tzinfo(dtparser.parse(other))
            other_seconds = other_dt.timestamp()
        elif isinstance(other, (float, int)):
            other_seconds = float(other)
        elif isinstance(other, (np.datetime64)):
            other_seconds = other.astype(np.float)
        elif isinstance(other, (datetime.datetime)):
            other_seconds = other.timestamp()

        else:
            msg = "Cannot compare {0} of type {1} with MTime Object".format(
                other, type(other)
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

        return self.epoch_seconds - other_seconds

    @property
    def iso_str(self):
        return self.dt_object.isoformat()

    @property
    def iso_no_tz(self):
        return self.dt_object.isoformat().split("+", 1)[0]

    @property
    def epoch_seconds(self):
        return self.dt_object.timestamp()

    @epoch_seconds.setter
    def epoch_seconds(self, seconds):
        self.logger.debug(
            "reading time from epoch seconds, assuming UTC " + "time zone"
        )
        dt = datetime.datetime.utcfromtimestamp(seconds)
        dt = dt.replace(tzinfo=datetime.timezone.utc)
        self.dt_object = dt

    def from_str(self, dt_str):
        try:
            self.dt_object = self.validate_tzinfo(dtparser.parse(dt_str))
        except dtparser.ParserError as error:
            msg = (
                f"{error}. "
                + "Input must be a valid datetime string, see "
                + "https://docs.python.org/3.8/library/datetime.html"
            )
            self.logger.error(msg)
            raise MTTimeError(msg)

    def validate_tzinfo(self, dt_object):
        """
        make sure the timezone is UTC
        """

        if dt_object.tzinfo == datetime.timezone.utc:
            return dt_object

        elif isinstance(dt_object.tzinfo, tzutc):
            return dt_object.replace(tzinfo=datetime.timezone.utc)

        elif dt_object.tzinfo is None:
            return dt_object.replace(tzinfo=datetime.timezone.utc)

        elif dt_object.tzinfo != datetime.timezone.utc:
            raise ValueError("Time zone must be UTC")

    @property
    def date(self):
        return self.dt_object.date().isoformat()

    @property
    def year(self):
        return self.dt_object.year

    @year.setter
    def year(self, value):
        self.dt_object = self.dt_object.replace(year=value)

    @property
    def month(self):
        return self.dt_object.month

    @month.setter
    def month(self, value):
        self.dt_object = self.dt_object.replace(month=value)

    @property
    def day(self):
        return self.dt_object.day

    @day.setter
    def day(self, value):
        self.dt_object = self.dt_object.replace(day=value)

    @property
    def hour(self):
        return self.dt_object.hour

    @hour.setter
    def hour(self, value):
        self.dt_object = self.dt_object.replace(hour=value)

    @property
    def minutes(self):
        return self.dt_object.minute

    @minutes.setter
    def minutes(self, value):
        self.dt_object = self.dt_object.replace(minute=value)

    @property
    def seconds(self):
        return self.dt_object.second

    @seconds.setter
    def seconds(self, value):
        self.dt_object = self.dt_object.replace(second=value)

    @property
    def microseconds(self):
        return self.dt_object.microsecond

    @microseconds.setter
    def microseconds(self, value):
        self.dt_object = self.dt_object.replace(microsecond=value)

    def now(self):
        """
        set date time to now

        :return: current UTC time
        :rtype: datetime with UTC timezone

        """
        self.dt_object = self.validate_tzinfo(datetime.datetime.utcnow())

    def copy(self):
        """ make a copy of the time """
        return deepcopy(self)


def get_now_utc(fmt="iso"):
    """
    Get the current time in UTC format
    
    :param str fmt: return format can be iso, date, obj
    
    :return: format specified by fmt
    :rtype: string, :class:`mtpy.utils.mttime.MTime`
    
    """

    m_obj = MTime()
    m_obj.now()
    if fmt == "iso":
        return m_obj.iso_str
    elif fmt == "date":
        return m_obj.date
    elif fmt == "obj":
        return m_obj
