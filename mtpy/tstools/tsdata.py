import pyasdf
import obspy.core.trace as trace
import numpy as np
import time

from tswave import TSWave
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


    def getwaveform(self, wave, starttime=None, endtime=None):
        rawdata, channel, wavename = self.wavemeta[wave]

        if starttime is None:
            starttime = channel.start_date
            endtime = starttime+1000

        print(starttime,endtime)

        ntwk = re.sub('([^.]+)(.*)','\\1', wavename)
        sttn = re.sub('([^.]+\.)([^.]+)(.*)','\\2', wavename)
        outwave = rawdata.get_waveforms(network=ntwk, station=sttn, location=channel.location_code, \
                            channel=channel.code, starttime=starttime, endtime=endtime, tag="raw_recording")

        return outwave, wavename, channel.start_date, channel.end_date



    def getlist(self):
        return self.wavelist

