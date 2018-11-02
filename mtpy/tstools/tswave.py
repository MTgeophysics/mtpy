import numpy as np
import obspy.core.trace as trace
from obspy.core.utcdatetime import UTCDateTime

class TSWave():
    def __init__(self, trace):
        self.times = []
        self.data = []
        self.decirate = None
        self.sampled = False
        self.sampledTrace = None

        self.setWave(trace)

    def setWave(self, trace):
        self.trace = trace

    def expandWave(self, trace):
        self.trace = self.trace.__add__(trace)

    def starttime(self):
        return self.trace.meta['starttime']

    def endtime(self):
        return self.trace.meta['endtime']

    def getWave(self, starttime=None, endtime=None, decirate=None):

        if starttime is None:
            starttime = self.trace.meta['starttime']
        else:
            starttime = UTCDateTime(starttime)

        if endtime is None:
            endtime = self.trace.meta['endtime']
        else:
            endtime = UTCDateTime(endtime)

        if decirate is None:
            decirate = round((endtime - starttime) / 1000)

        if decirate<1:
            decirate = 1

        if decirate != self.decirate:
            self.decirate = decirate
            self.sampled = False

        if self.sampled == True:
            pass
        elif decirate>16:
            print('resampled', decirate)
            self.sampledTrace = trace.Trace()
            self.sampledTrace.data = self.trace.data[::decirate].copy()
            self.sampledTrace.meta['delta'] = self.trace.meta['delta']*decirate
            self.sampledTrace.meta['starttime'] = self.trace.meta['starttime']
            self.sampledTrace = self.sampledTrace.decimate(1, True)
        else:
            print('resampled', decirate)
            self.sampledTrace = self.trace.copy().decimate(self.decirate, True)


        self.sampled = True

        return self.sampledTrace.slice(starttime, endtime)