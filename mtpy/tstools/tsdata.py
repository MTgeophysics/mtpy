
import pyasdf
import obspy.core.trace as trace
import numpy as np
import time

from tswave import TSWave
import re
class TSData():
    def __init__(self, filename=None):

        self.wavelist = {}

        if filename is not None:
            self.loadFile(filename)



    def loadFile(self, filename):
        self.rawdata = pyasdf.ASDFDataSet(filename, mode="r")

        count = 0
        for stationname in self.rawdata.waveforms.list():
            for network in self.rawdata.waveforms[stationname].StationXML:
                if network.code not in self.wavelist:
                    self.wavelist[network.code] = {}
                for station in network:
                    if station.code not in self.wavelist[network.code]:
                        self.wavelist[network.code][station.code] = {}
                    for channel in station:
                        wavename = network.code+'.'+station.code+'.'+channel.location_code + '.' + channel.code
                        if wavename not in self.wavelist[network.code][station.code]:
                            self.wavelist[network.code][station.code][wavename] = [channel]
                            #print(wavename, self.wavelist[network.code][station.code][wavename])
                        else:
                            self.wavelist[network.code][station.code][wavename].append(channel)







    def getwaveform(self, wave, starttime=None, endtime=None):
        print("getwave")
        print(type(wave), 'type1')
        print(type(wave.channelitem),'type2')
        print(dir(wave.channelitem))
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



    def getwaveformresample(self, shift):

        if self.decirate==1 and shift<0:
            pass
        elif self.decirate==16 and shift>0:
            pass
        elif shift==0:
            pass
        else:
            decirate = int(self.decirate + shift)

            if decirate < 1:
                decirate = 1
            if decirate > 16:
                decirate = 16




            self.currenttimes = trace.Trace(
                self.rawtimes).decimate(decirate,True).copy()


            self.currentdata = trace.Trace(
                self.rawdata).decimate(decirate,True).copy()



            gap = self.end-self.start
            self.start = int(float(self.start) * self.decirate / decirate)
            self.end = self.start+gap
            if self.end<1:
                self.end=1
            if self.end>len(self.currentdata):
                self.end=len(self.currentdata)



            self.decirate = decirate

            self.times = np.array(self.currenttimes[self.start: self.end])
            self.data = np.array(self.currentdata[self.start: self.end])

        return [self.times,
                self.data]


    def recommendaspect(self):
        self.aspect = float(self.times.max()-self.times.min())/(self.data.max()-self.data.min())/10

        return self.aspect

    def getList(self):
        return self.wavelist

