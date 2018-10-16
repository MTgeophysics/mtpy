
import pyasdf
import obspy.core.trace as trace
import numpy as np
import time

from tswave import TSWave

class TSData():
    def __init__(self, filename=None):

        self.dict = {}

        if filename is not None:
            self.loadFile(filename)



    def loadFile(self, filename):
        rawdata = pyasdf.ASDFDataSet(filename, mode="r")


        for name in rawdata.waveforms.list():
            for w  in rawdata.waveforms[name].raw_recording:
                fullname = name+'.'+w.meta['location']+'.'+w.meta['channel']

                if fullname in self.dict:
                    self.dict[fullname].expandWave(w)
                else:
                    self.dict[fullname] = TSWave(w)






    def getwaveform(self, waveformname, starttime=None, endtime=None):
        wave = self.dict[waveformname].getWave(starttime, endtime)
        return wave


    # def getwaveformshift(self, shift):
    #
    #     shift = (self.endtime-self.starttime+1)*shift
    #     self.starttime = self.starttime + shift
    #     self.endtime = self.endtime + shift
    #
    #     return [wave for wave in self.vis]
    #
    # def getwaveform
    #
    #     if self.start==0 and shift<0:
    #         pass
    #     elif self.end == len(self.currentdata) and shift >0:
    #         pass
    #     else:
    #         if self.start+shift<0:
    #             shift = -self.start
    #         if self.end+shift>len(self.currentdata):
    #             shift = len(self.currentdata) -self.end
    #
    #         self.start = int(self.start+shift)
    #         self.end = int(self.end+shift)
    #
    #         self.times = np.array(self.currenttimes[self.start: self.end])
    #         self.data = np.array(self.currentdata[self.start: self.end])
    #
    #     return [self.times,
    #             self.data]

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
        return list(self.dict)