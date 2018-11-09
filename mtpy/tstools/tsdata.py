import pyasdf
import obspy.core.trace as trace
import numpy as np
import time

from tswave import TSWave
import re


class TSData():
    def __init__(self, filename: str = None):
        self.wavelist = {}

        if filename is not None:
            self.loadFile(filename)


    def loadFile(self, filename: str):
        rawdata = pyasdf.ASDFDataSet(filename, mode="r")

        count = 0
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
                            self.wavelist[network.code][station.code][wavename] = [(rawdata, channel)]
                        else:
                            self.wavelist[network.code][station.code][wavename].append((rawdata, channel))


    def getwaveform(self, wave, starttime=None, endtime=None):
        if starttime is None:
            starttime = wave.channelitem.start_date
            endtime = starttime+1000

        print(starttime,"starttime")
        print(endtime,"endtime")

        ntwk = re.sub('([^.]+)(.*)','\\1',wave.wavename)
        sttn = re.sub('([^.]+\.)([^.]+)(.*)','\\2',wave.wavename)
        outwave = self.rawdata.get_waveforms(network=ntwk, station=sttn, location=wave.channelitem.location_code, \
                            channel=wave.channelitem.code, starttime=starttime, endtime=endtime, tag="raw_recording")
        print("""self.rawdata.get_waveforms(network="""+ntwk+""", station="""+sttn+""", 
        location="""+wave.channelitem.location_code+""", channel="""+wave.channelitem.code+""", 
        starttime="""+str(wave.channelitem.start_date)+""", endtime="""+str(wave.channelitem.end_date)+""", tag="raw_recording")""")


        return outwave



    def getlist(self):
        return self.wavelist

