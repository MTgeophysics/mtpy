import pyasdf
import obspy.core.trace as trace
from obspy.core.utcdatetime import UTCDateTime

import numpy as np
import time

from datetime import datetime

import re


class TSData():
    def __init__(self, filename: str = None, numofsamples: int = 400, cachesize = 1e+9):
        self.wavelist = {}
        self.wavemeta = {}
        print("ini","!"*10)

        if filename is not None:
            self.loadFile(filename)

        self.wavecache = {}
        self.numofsamples = numofsamples
        self.cachesize = cachesize


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
        print(len(self.wavemeta))

    def getwaveform(self, waveform: str, starttime: datetime=None, endtime: datetime=None):

        if starttime is None or endtime is None:
            timewindow = None
        else:
            timewindow = endtime-starttime

        if waveform in self.wavecache and (timewindow is None or abs(self.wavecache[waveform][0] - timewindow)/timewindow<0.1):
            return self.readcache(waveform, starttime, endtime)
        elif self.writecache(waveform, timewindow):
            return self.readcache(waveform, starttime, endtime)
        else:
            outwave, wavename, start_date, end_date, gaps= self.readdisc(waveform, starttime, endtime)
            wave = np.vstack((outwave.times() + outwave.meta['starttime'].timestamp, outwave.data))

            wave = wave[:, wave[0, :] != np.nan]
            wave = wave[:, wave[1, :] != np.nan]
            return wave, wavename, start_date, end_date, gaps


    def writecache(self, waveform, timewindow):
        #print('writecache', timewindow)
        _, channel, _ = self.wavemeta[waveform]
        if (channel.end_date - channel.start_date) / channel.sample_rate < self.cachesize:
            outwave, wavename, start_date, end_date, gaps = \
                self.readdisc(waveform, channel.start_date, channel.end_date, False)

            if timewindow is None:
                timewindow = end_date - start_date


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

            wave = np.vstack((outwave.times() + outwave.meta['starttime'].timestamp, outwave.data))
            wave = np.array(wave).copy()
            self.wavecache[waveform] = timewindow, wave , wavename, start_date, end_date, gaps


            return True
        else:
            return False

    def readcache(self, waveform: str, starttime: datetime, endtime: datetime):
        #print('readcache', starttime, endtime)
        timewindow, wave, wavename, start_date, end_date, gaps = self.wavecache[waveform]

        if starttime is None:
            starttime = start_date

        if endtime is None:
            endtime = end_date



        head = int((starttime - start_date)/(wave[0,1]-wave[0,0]))
        tail = int((endtime - start_date) / (wave[0,1] - wave[0,0]))

        if tail >= wave.shape[1]:
            tail = wave.shape[1]-1
        outwave = wave[:,head: tail]

        return outwave, wavename, start_date, end_date, gaps

    def readdisc(self, waveform: str, starttime: datetime, endtime: datetime, resample: bool=True, fill_value:str='latest'):
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



        if len(outwave)>0:
            gaps = outwave.get_gaps()

            mergewave = outwave[0]
            for w in outwave[1:]:
                mergewave = mergewave.__add__(w,fill_value=fill_value)

            outwave = mergewave





            rate = round(float(len(outwave.data)) / self.numofsamples)
            if resample == False or rate<=1:
                pass
            elif rate>16:
                tmp = trace.Trace()
                tmp.data = outwave.data[::rate].copy()
                tmp.meta['delta'] = outwave.meta['delta'] * rate
                tmp.meta['starttime'] = outwave.meta['starttime']
                outwave = tmp  # .decimate(1, True)
                #print(outwave.meta['endtime'], 'new endtime  ')
            elif rate>=1:
                outwave.decimate(rate)
        else:
            outwave = trace.Trace()
            outwave.data = np.array([np.nan]*self.numofsamples)
            outwave.meta['starttime'] = starttime
            outwave.meta['delta'] = (endtime-starttime)/self.numofsamples

        print(channel.start_date, channel.end_date,'==================')

        return outwave, wavename, channel.start_date, channel.end_date, gaps



    def getlist(self):
        return self.wavelist

