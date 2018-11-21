import pyasdf
import obspy.core.trace as trace
import numpy as np
import time

from datetime import datetime

import re


class TSData():
    def __init__(self, filename: str = None):
        self.wavelist = {}
        self.wavemeta = {}

        if filename is not None:
            self.loadFile(filename)


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


    def getwaveform(self, waveform: str, starttime: datetime=None, endtime: datetime=None, numofsamples: int=0):
        rawdata, channel, wavename = self.wavemeta[waveform]

        if starttime is None:
            starttime = channel.start_date

        if endtime is None:
            endtime = starttime+1000

        ntwk = re.sub('([^.]+)(.*)','\\1', wavename)
        sttn = re.sub('([^.]+\.)([^.]+)(.*)','\\2', wavename)

        outwave = rawdata.get_waveforms(network=ntwk, station=sttn, location=channel.location_code, \
                            channel=channel.code, starttime=starttime, endtime=endtime, tag="raw_recording")

        if len(outwave)>0:
            outwave = outwave[0]

            if numofsamples==0:
                pass
            else:
                rate = round(float(len(outwave.data)) / numofsamples)
                if rate<=1:
                    pass
                elif rate>16:
                    tmp = trace.Trace()
                    tmp.data = outwave.data[::rate].copy()
                    tmp.meta['delta'] = outwave.meta['delta'] * rate
                    tmp.meta['starttime'] = outwave.meta['starttime']
                    outwave = tmp.decimate(1, True)
                elif rate>=1:
                    outwave.decimate(rate)
        else:
            outwave = trace.Trace()
            outwave.data = np.array([np.nan]*numofsamples)
            outwave.meta['starttime'] = starttime
            outwave.meta['delta'] = (endtime-starttime)/numofsamples


        return outwave, wavename, channel.start_date, channel.end_date



    def getlist(self):
        return self.wavelist

