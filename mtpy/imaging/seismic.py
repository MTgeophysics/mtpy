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
from scipy.spatial import cKDTree
import numpy

try:
    from obspy.io.segy import segy
except ImportError:
    print("Could not find Obspy, check installation")

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
        self._si = []  # sample interval
        count = 0
        for itr, tr in enumerate(self._sf.traces):
            if (itr % self._pick_every == 0):
                self._mtraces.append(tr.data)
                self._xs.append(tr.header.source_coordinate_x)
                self._ys.append(tr.header.source_coordinate_y)
                self._ns.append(tr.header.number_of_samples_in_this_trace)
                self._si.append(tr.header.sample_interval_in_ms_for_this_trace/1e6) # convert from micro seconds to s
                self._cdps.append(tr.header.ensemble_number)
                count += 1
        # end for

        self._ntraces = len(self._sf.traces)

        self._ns = numpy.array(self._ns)
        self._si = numpy.array(self._si)
        self._xs = numpy.array(self._xs)
        self._ys = numpy.array(self._ys)
        self._cdps = numpy.array(self._cdps)
        self._station_space = numpy.sqrt((self._xs[:-1] - self._xs[1:])**2 + (self._ys[:-1] - self._ys[1:])**2)
        self._dist = numpy.array([numpy.sum(self._station_space[:i]) for i in range(len(self._station_space))]) # compute euclidean distance along profile
        self._ts = numpy.linspace(0, (numpy.max(self._ns)-1)*numpy.max(self._si),
                                  numpy.max(self._ns))
        self._mtraces = numpy.array(self._mtraces)
        self._mt, self._md = numpy.meshgrid(self._ts, self._dist)  # mesh grid of distance and time
        self._mt, self._mc = numpy.meshgrid(self._ts, self._cdps)  # mesh grid of cdp and time

        self._dist -= numpy.min(self._dist)  # distance along profile starts from 0
    # end func

    def getDistances(self):
        '''
        Distance of each trace from the start of the 2D line
        :return: numpy array of distances of shape (nt), where nt is the number of
                 traces.
        '''

        return self._dist
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
        if (dist <= numpy.max(self._dist)):
            idx = numpy.argmin(numpy.fabs(self._dist - dist))

            if (key == 'trace'): return self._mtraces[idx, :]
            elif (key == 'x'): return self._xs[idx]
            elif (key == 'y'): return self._ys[idx]
            elif (key == 'cdp'): return self._cdps[idx]
            elif (key == 'ts'): return self._ts
            else:
                raise ValueError('Attribute "%s" not found' % key)
        else:
            return None
    # end func

    def getMigratedProfile(self, velocity_model, ntraces=-1, ndepths=-1, max_depth=60e3, time_shift=0, nn=1):
        '''
        Computes and returns a depth-migrated image, defined by distances
        and depths along the profile, based on the given velocity model.
        Note that the velocity_model may not necessarily contain depth-profiles
        for all CDP values in the segy file, in which case, depth-profiles for
        closest CDP values are used.

        :param velocity_model: either an instance of class VelocityModel, containing
                                stacking velocities for the given 2D seismic
                                line, or a floating point value representing
                                a constant velocity (m/s)
        :param ntraces: number of equally spaced interpolated traces to compute along
                        a seismic line; when set to -1, ntraces is internally reset to
                        the number of traces in the SegY file
        :param ndepths: number of depth values to compute along a seismic line,
                        starting from 0 to 'max_depth'. When set to -1, ndepth
                        is internally reset to the number of samples in a trace
        :param max_depth: the maximum depth to which depth values are computed.
        :param time_shift: Time shift in ms to remove topography
        :param nn: While computing the depth-migration, this parameter dictates
                   whether to use the mean of the closest 'nn' or the mean of all
                   available cdp depth-profiles available in the velocity_model,
                   if set to -1. In the case of the latter, the maximum depth
                   computed for each trace will be the same.
        :return: mdepths: 2D array of depths of shape (depths, distances),
                 mdistances: 2D array of distances of shape (depths, distances)
                 mvals: 2D array of amplitudes of shape (depths, distances),
                        by default, amplitudes for depths outside of the
                        depth-ranges available are set to -9999
                 xy_list: 2D array of spatial coordinates (x, y) of traces along
                          the seismic line of shape (ntraces, 2)
        '''

        constant_veolocity_model = False
        if(not isinstance(velocity_model, VelocityModel)):
            try:
                velocity_model = numpy.float_(velocity_model)
                constant_veolocity_model = True
            except:
                raise ValueError('Invalid velocity_model')
        # end if

        if (ntraces == -1): ntraces = self._ntraces
        if (ndepths == -1): ndepths = self._mtraces.shape[1]

        gdepth    = numpy.linspace(0, max_depth, ndepths)
        gdistance = numpy.linspace(0, numpy.max(self._dist), ntraces)

        mdepths, mdistances = numpy.meshgrid(gdepth, gdistance)
        mvals = numpy.zeros(mdepths.shape)

        shiftIdx = None
        if(time_shift>0):
            time_shift /= 1e3  # convert to seconds
            if (time_shift >= self._ts[-1]):
                raise RuntimeError
            shiftIdx = numpy.ceil(time_shift / numpy.max(self._si))
            shiftIdx = numpy.int_(shiftIdx)
        # end if

        gx = []
        gy = []
        for idist, dist in enumerate(gdistance):

            gx.append(self.getAttribute('x', dist))
            gy.append(self.getAttribute('y', dist))

            cdp = None
            if(nn != -1): cdp = self.getAttribute('cdp', dist)

            ts = None
            amps = None
            if (time_shift == 0):
                ts = self._ts
                amps = self.getAttribute('trace', dist)
            else:
                ts = self._ts[shiftIdx:] - time_shift
                amps = self.getAttribute('trace', dist)[shiftIdx:]
            # end if

            if(constant_veolocity_model):
                ds = ts*velocity_model
            else:
                ds = velocity_model.getDepth(cdp, ts, nn)
            # end if

            io = interp1d(ds, amps,
                          bounds_error=False, fill_value=-9999)
            mvals[idist, :] = io(gdepth)
        # end for

        xy_list = numpy.array([[x, y] for x, y in zip(gx, gy)])

        return mdepths, mdistances, mvals, xy_list
    # end func
# end class

class VelocityModel:
    def __init__(self, stacking_velocity_fn, ni=50):
        '''
        Class for computing interval velocities using Dix formula, based on provided
        stacking velocities.

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
            print(('Failed to read %s' % (stacking_velocity_fn)))
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
                    tempList = numpy.int_(numpy.array(tempList))
                    ts = tempList[::2] / 1e3
                    vs = tempList[1::2]
                    # extend beyond given time range
                    ts = numpy.append(ts, 1e6)
                    vs = numpy.append(vs, vs[-1])

                    self._times.append(ts)
                    self._vels.append(vs)
                    tempList = []
            elif ('END' in line):
                tempList = numpy.int_(numpy.array(tempList))
                ts = tempList[::2] / 1e3
                vs = tempList[1::2]
                # extend beyond given time range
                ts = numpy.append(ts, 1e6)
                vs = numpy.append(vs, vs[-1])
                
                self._times.append(ts)
                self._vels.append(vs)
            else:
                if ('*' not in line):
                    items = line.split()
                    for i in items:
                        tempList.append(i)
        # end for
        f.close()
        self._cdps = numpy.int_(self._cdps)

        # Create Kd-Tree for CDP queries
        cdpArray = numpy.expand_dims(self._cdps, 1)
        self._cdp_tree = cKDTree(cdpArray)

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

        idx = numpy.argmax(numpy.logical_and(t >= self._intStarts[cdp], t <= self._intEnds[cdp]))
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
            self._intVels[cdp] = numpy.array(self._intVels[cdp])
            self._intStarts[cdp] = numpy.array(self._intStarts[cdp])
            self._intEnds[cdp] = numpy.array(self._intEnds[cdp])
        # end for

        # Generate uniformly spaced time intervals up to the latest time value in the stacking-velocities
        # file; finally add the entry that extends the time-range to 1e6, as done in __init__.
        self._ts = numpy.linspace(self._times[0][0], self._times[0][-2], self._ni-1)
        self._ts = numpy.append(self._ts, self._times[0][-1])

        # Integrate each velocity profile to compute depth profile
        self._depth_ios = []
        self._cdp_mean_interval_velocity = []
        depths = numpy.zeros(self._ts.shape)
        for icdp, cdp in enumerate(self._cdps):
            vs = numpy.array([self._getIntervalVelocity(cdp, t) for t in self._ts])
            self._cdp_mean_interval_velocity.append(numpy.mean(vs))

            io = interp1d(self._ts, vs)
            for it, t in enumerate(self._ts):
                # progressively build up the depth-profile
                if (it == 0):
                    depths[it] = quad(io, 0, t)[0]
                else:
                    depths[it] = depths[it - 1] + quad(io, self._ts[it - 1], t)[0]
            # end for
            depths /= 2.0 # two-way-traveltime
            self._depth_ios.append(interp1d(self._ts, depths))
        # end for
        self._cdp_mean_interval_velocity = numpy.array(self._cdp_mean_interval_velocity)

        # Compute mean depth profile
        self._mean_depth_profile = numpy.zeros(self._ts.shape)
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
        nn: number of closest cdp depth-profiles to be used for calculating depth
            values; must be >= 1

        '''
        if (cdp == None):
            return self._mean_depth_profile_io(ts)
        else:
            _, idx = self._cdp_tree.query([cdp], nn)
            ds = numpy.zeros(ts.shape)

            if (type(idx) == int): idx = [idx]
            for inn in range(nn):
                io = self._depth_ios[idx[inn]]
                ds += io(ts)
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
