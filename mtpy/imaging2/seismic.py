#!/bin/env python
"""
Description:
    Classes for reading 2D segy data and for time-to-depth conversion of seismic data
    based on a given velocity model

References:
 
CreationDate:   11/28/17
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     11/28/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import logging, traceback

from collections import defaultdict
from scipy.interpolate import interp1d
from scipy.integrate import quad

import numpy as np
from obspy.io.segy import segy

class Segy:
    def __init__(self, segy_fn, pick_every=1):
        '''
        Class for reading 2D segy files. This class assumes that the traces are
        CDP-ordered.

        :param segy_fn: segy file name
        :param pick_every: by default, every trace is read from the segy file, but
                           pick_every>1 allows skipping traces.
        '''
        self._sfn = segy_fn
        self._sf = segy._read_segy(self._sfn)
        self._pick_every = pick_every

        # Read traces
        self._mtraces = []  # trace samples
        self._xs = []  # x coordinates
        self._ys = []  # y coordinates
        self._ns = []  # num samples
        self._cdps = []  # ensemble number (cdps)
        self._ls = []  # sample length in seconds
        count = 0
        for itr, tr in enumerate(self._sf.traces):
            if (itr % self._pick_every == 0):
                self._mtraces.append(tr.data)
                self._xs.append(tr.header.source_coordinate_x)
                self._ys.append(tr.header.source_coordinate_y)
                self._ns.append(tr.header.number_of_samples_in_this_trace)
                self._ls.append(tr.header.sample_interval_in_ms_for_this_trace / 1e3) # convert ms to s
                self._cdps.append(tr.header.ensemble_number)
                count += 1
        # end for

        self._ntraces = len(self._sf.traces)

        self._ns = np.array(self._ns)
        self._ls = np.array(self._ls)
        self._xs = np.array(self._xs)
        self._ys = np.array(self._ys)
        self._cdps = np.array(self._cdps)
        self._dist = np.sqrt(self._xs ** 2 + self._ys ** 2) # compute euclidean distance along profile
        self._ts = np.linspace(0, np.max(self._ls), np.max(self._ns))

        self._mtraces = np.array(self._mtraces)
        self._mt, self._md = np.meshgrid(self._ts, self._dist)  # mesh grid of distance and time
        self._mt, self._mc = np.meshgrid(self._ts, self._cdps)  # mesh grid of cdp and time

        self._dist -= np.min(self._dist)  # distance along profile starts from 0
    # end func

    def getAttribute(self, key, dist):
        '''
        Returns attribute value -- for the given key -- of the closest trace at a given
        distance along the seismic profile.

        :param key: attribute key; see key-value pairs below:
                    'trace' : trace samples
                    'x': x-coordinate of trace
                    'y': y-coordinate of trace
                    'cdp': CDP of trace
                    'ts': time samples
        :param dist: distance along profile
        :return: attribute value of trace samples for a trace at a given distance along the
                 seismic profile.
        '''

        if(key not in ['trace', 'x', 'y', 'cdp', 'ts']):
            assert 0, "Invalid key; should be one of ['trace', 'x', 'y', 'cdp', 'ts']"
        if (dist < np.max(self._dist)):
            idx = np.argmin(np.fabs(self._dist - dist))

            if (key == 'trace'): return self._mtraces[idx]
            elif (key == 'x'): return self._xs[idx]
            elif (key == 'y'): return self._ys[idx]
            elif (key == 'cdp'): return self._cdps[idx]
            elif (key == 'ts'): return self._ts
        else:
            return None
        # end func
# end class

class VelocityModel:
    def __init__(self, stacking_velocity_fn, ni=50):
        '''

        :param velocity_fn: stacking velocity file name
        :param ni: number of time intervals in each velocity profile; a larger value
            produces a smoother model at the expense of increased computational
            cost. This value should be ideally ~2x the average number of intervals
            in the stacking velocity file.
        '''

        f = None
        try:
            f = open(stacking_velocity_fn)
        except Exception as err:
            print 'Failed to read %s' % (stacking_velocity_fn)
            logging.error(traceback.format_exc())
            exit(-1)

        self._ni = ni

        # read stacking velocity file and extract time intervals and
        # corresponding velocities for each CDP location
        self._cdps = []
        self._times = []
        self._vels = []
        tempList = []
        for line in f:
            if ('HANDVEL' in line):
                cdp = line.split()[1]
                self._cdps.append(cdp)
                if (len(tempList)):
                    tempList = np.int_(np.array(tempList))
                    self._times.append(tempList[::2] / 1e3)
                    self._vels.append(tempList[1::2])
                    tempList = []
            elif ('END' in line):
                tempList = np.int_(np.array(tempList))
                self._times.append(tempList[::2] / 1e3)
                self._vels.append(tempList[1::2])
            else:
                if ('*' not in line):
                    items = line.split()
                    for i in items:
                        tempList.append(i)
        # end for
        f.close()
        self._cdps = np.int_(self._cdps)

        # generate depth model
        self._generateDepthModel()
    # end func

    def _getIntervalVelocity(self, cdp, t):
        '''
        function to retrieve velocity within correct inverval

        :param cdp: cdp number
        :param t: time sample
        :return: returns corresponding interval velocity
        '''

        idx = np.argmax(np.logical_and(t >= self._intStarts[cdp], t < self._intEnds[cdp]))
        return self._intVels[cdp][idx]
    # end func

    def _generateDepthModel(self):
        # Generate interval velocities using Dix' formula
        self._intVels = defaultdict(list)
        self._intStarts = defaultdict(list)
        self._intEnds = defaultdict(list)

        for icdp, cdp in enumerate(self._cdps):
            for i in range(len(self._times[icdp]) - 1):
                self._intVels[cdp].append(((self._times[icdp][i + 1] * self._vels[icdp][i + 1] ** 2 -
                                            self._times[icdp][i] * self._vels[icdp][i] ** 2) /
                                           (self._times[icdp][i + 1] - self._times[icdp][i])) ** 0.5)
                self._intStarts[cdp].append(self._times[icdp][i])
                self._intEnds[cdp].append(self._times[icdp][i + 1])
            # end for
            self._intVels[cdp] = np.array(self._intVels[cdp])
            self._intStarts[cdp] = np.array(self._intStarts[cdp])
            self._intEnds[cdp] = np.array(self._intEnds[cdp])
        # end for

        # Integrate each velocity profile to compute depth profile
        self._ts = np.linspace(self._times[0][0], self._times[0][-1], self._ni)
        self._depth_ios = []

        depths = np.zeros(self._ts.shape)
        for icdp, cdp in enumerate(self._cdps):
            vs = np.array([self._getIntervalVelocity(cdp, t) for t in self._ts])
            io = interp1d(self._ts, vs)
            for it, t in enumerate(self._ts):
                # Doing much more work here than is needed. Should be
                # done using a recursive function; TODO.
                depths[it] = quad(io, 0, t)[0]
            # end for
            self._depth_ios.append(interp1d(self._ts, depths))
        # end for

        # Compute mean depth profile
        self._mean_depth_profile = np.zeros(self._ts.shape)
        for icdp, cdp in enumerate(self._cdps):
            io = self._depth_ios[icdp]

            for it, t in enumerate(self._ts):
                self._mean_depth_profile[it] += io(t)
        # end for
        self._mean_depth_profile /= float(len(self._cdps))
        self._mean_depth_profile_io = interp1d(self._ts,
                                               self._mean_depth_profile)
    # end func

    def getDepth(self, cdp, ts, nn=1):
        '''
        cdp: cdp number. If cdp is None, use mean depth profile; nn is
             ignored
        ts: time samples in seconds
        nn: number of neighbours to be used for calculating depth
            value using inverse distance weighted interpolation;

        '''
        if (cdp == None):
            return self._mean_depth_profile_io(ts)
        else:
            idx = np.argsort(np.fabs(self._cdps - cdp))[:nn]
            ds = np.zeros(ts.shape)

            for inn in range(nn):
                io = self._depth_ios[idx[inn]]
                for it, t in enumerate(ts):
                    ds[it] += io(t)
            return ds / float(nn)
        # end if
    # end func
# end class

def main():
    """
    define main function
    :return:
    """

    utils = os.path.dirname(__file__)
    mtpy = os.path.dirname(utils)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'seismic')
    s_fn = os.path.join(ModEM_files, 'seismic.sgy')
    v_fn = os.path.join(ModEM_files, 'stacking_velocities.txt')

    s = Segy(segy_fn=s_fn)
    v = VelocityModel(v_fn, ni=20)
    return
# end


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()