import pyasdf
import obspy.core.trace as trace
from obspy.core.utcdatetime import UTCDateTime

import numpy as np
import time

from datetime import datetime

import re


class TSData():
    def __init__(self, filename: str = None, numofsamples: int = 1000):
        self.wavelist = {}
        self.wavemeta = {}

        if filename is not None:
            self.loadFile(filename)

        self.wavecache = {}
        self.numofsamples = numofsamples


    def loadFile(self, filename: str):
        rawdata = pyasdf.ASDFDataSet(filename, mode="r")

        for stationname in rawdata.waveforms.list():
            for network in rawdata.waveforms[stationname].StationXML:
                if network.code not in self.wavelist:
                    self.wavelist[network.code] = {}
                for station in network:
                    if station.code not in self.wavelist[network.code]:
                        self.wavelist[network.code][station.code] = {}
                    for channel in station:
                        wavename = network.code + '.' + station.code + '.' + channel.location_code + '.' + channel.code
                        if wavename not in self.wavelist[network.code][station.code]:
                            self.wavelist[network.code][station.code][wavename] = [str(channel)]
                        else:
                            self.wavelist[network.code][station.code][wavename].append(str(channel))
                        self.wavemeta[str(channel)] = (rawdata, channel, wavename)


    def getwaveform(self, waveform: str, starttime: datetime=None, endtime: datetime=None):

        if starttime is None or endtime is None:
            timewindow = None
        else:
            timewindow = endtime-starttime

        #print('!!!!!!timewindow=',timewindow)

        if waveform in self.wavecache and (timewindow is None or abs(self.wavecache[waveform][0] - timewindow)/timewindow<0.01):
            return self.readcache(waveform, starttime, endtime)
        elif self.writecache(waveform, timewindow):
            return self.readcache(waveform, starttime, endtime)
        else:
            return self.readdisc(rawdata, channel, wavename, starttime, endtime)


    def writecache(self, waveform, timewindow):
        # print('writecache', timewindow)
        _, channel, _ = self.wavemeta[waveform]
        if (channel.end_date - channel.start_date) / channel.sample_rate < 1e+9:
            outwave, wavename, start_date, end_date = \
                self.readdisc(waveform, channel.start_date, channel.end_date, False)

            if timewindow is None:
                timewindow = end_date - start_date

            # print('rate=',timewindow, channel.sample_rate, self.numofsamples)
            # print('length=', len(outwave.data))
            rate = timewindow*channel.sample_rate/self.numofsamples
            rate = int(rate)

            if rate == 0:
                pass
            else:
                if rate <= 1:
                    pass
                elif rate > 16:
                    tmp = trace.Trace()
                    tmp.data = outwave.data[::rate].copy()
                    tmp.meta['delta'] = outwave.meta['delta'] * rate
                    tmp.meta['starttime'] = outwave.meta['starttime']
                    outwave = tmp  # .decimate(1, True)
                elif rate >= 1:
                    outwave.decimate(rate)

            # print('writecache length=',len(outwave.data))

            #self.wavecache[waveform] = timewindow, outwave, wavename, start_date, end_date
            self.wavecache[waveform] = timewindow, np.vstack((outwave.times()+outwave.meta['starttime'].timestamp, outwave.data)) , wavename, start_date, end_date

            # print(type(outwave.data),'!'*10)

            print(np.vstack((outwave.data, outwave.times())).shape,'!'*10,'writecache')
            print(np.vstack((outwave.times()+outwave.meta['starttime'].timestamp, outwave.data)).shape,'writecache')
            print(waveform)

            return True
        else:
            return False

    def readcache(self, waveform: str, starttime: datetime, endtime: datetime):
        # print('readcache', starttime, endtime)
        timewindow, wave, wavename, start_date, end_date = self.wavecache[waveform]
        # print(timewindow, wave, wavename, start_date, end_date)

        if starttime is None:
            starttime = start_date

        if endtime is None:
            endtime = end_date


        #outwave = wave.slice(starttime, endtime)
        print(wave[0,0],wave[0,-1],starttime, endtime)

        head = np.argmax(wave[0,:] > starttime)
        tail = np.argmax(wave[0,:] > endtime)
        if tail == 0 and wave[0,0]<endtime:
            tail = len(wave[0,:])
        outwave = wave[:,head: tail]

        print(outwave.shape,'readcache',head, tail)

        return outwave, wavename, starttime, endtime

    def readdisc(self, waveform: str, starttime: datetime, endtime: datetime, resample: bool=True):
        print('readdisc', starttime, endtime)
        rawdata, channel, wavename = self.wavemeta[waveform]
        ntwk = re.sub('([^.]+)(.*)','\\1', wavename)
        sttn = re.sub('([^.]+\.)([^.]+)(.*)','\\2', wavename)

        if starttime is None:
            starttime = channel.start_date

        if endtime is None:
            endtime = channel.end_date

        outwave = rawdata.get_waveforms(network=ntwk, station=sttn, location=channel.location_code, \
                            channel=channel.code, starttime=starttime, endtime=endtime, tag="raw_recording")
        #print(ntwk, sttn, channel.location_code, channel.code, starttime, endtime, "raw_recording")
        #print(outwave)

        if len(outwave)>0:
            outwave = outwave[0]

            rate = round(float(len(outwave.data)) / self.numofsamples)
            if resample == False or rate<=1:
                pass
            elif rate>16:
                tmp = trace.Trace()
                tmp.data = outwave.data[::rate].copy()
                tmp.meta['delta'] = outwave.meta['delta'] * rate
                tmp.meta['starttime'] = outwave.meta['starttime']
                outwave = tmp  # .decimate(1, True)
            elif rate>=1:
                outwave.decimate(rate)
        else:
            outwave = trace.Trace()
            outwave.data = np.array([np.nan]*self.numofsamples)
            outwave.meta['starttime'] = starttime
            outwave.meta['delta'] = (endtime-starttime)/self.numofsamples


        return outwave, wavename, channel.start_date, channel.end_date



    def getlist(self):
        return self.wavelist

