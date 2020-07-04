# -*- coding: utf-8 -*-
"""
Created on Thu Jun 05 13:55:13 2014

@author: Alison Kirkby

"""

import mtpy.core.edi as mtedi
import os
import os.path as op
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import mtpy.utils.exceptions as MTex
import mtpy.utils.calculator as MTcc
import mtpy.analysis.geometry as MTg
import cmath
import math


class Control():

    def __init__(self, **input_parameters):

        self.run_input = [1, 0, 0.1, 40, 1.05, 1, 0]
        # define control file parameters
        self.iteration_max = 100  # max number of iterations
        self.penalty_type_structure = 6
        self.penalty_type_anisotropy = 2  # type of structure and anisotropy penalties
        # values for the structure penalty weights
        self.penalty_weight_structure = [0.1, 1.0, 10.0]
        # values for the anisotropy penalty weights
        self.penalty_weight_anisotropy = [0.1, 1.0, 10.0]
        self.working_directory = '.'

        for key in list(input_parameters.keys()):
            setattr(self, key, input_parameters[key])

        if not os.path.exists(self.working_directory):
            os.mkdir(self.working_directory)

    def write_ctlfile(self):
        """
        write control file
        """

        # create control file
        # control file name is hardcoded into software!
        ctlfile = open(os.path.join(
            self.working_directory, 'inregulm.dat'), 'wb')

        # define number of weights
        nw_struct = len(self.penalty_weight_structure)
        nw_aniso = len(self.penalty_weight_anisotropy)

        for thing in [(2, self.iteration_max), (nw_struct, nw_aniso), (self.penalty_type_structure, self.penalty_type_anisotropy)]:
            ctlfile.write('%1i%6i\n' % thing)
        for thing in [self.penalty_weight_structure, self.penalty_weight_anisotropy]:
            ctlfile.write('  '.join([str(i) for i in thing]) + '\n')
        ctlfile.close()
        print("written control file to {}".format(self.working_directory))

    inmodel_kwds = ['inmodel_dictionary']


class Inmodel():
    """
    **inmodel_
    """

    def __init__(self, **input_parameters):
        self.working_directory = '.'
        self.inmodel_modeldir = None
        self.inmodelfile = 'inmodel.dat'
        # dictionary containing values for
        self.inmodel_dictionary = {0: [100, 100, 0]}
        # inmodel file, in format topdepth: [minres,maxres,strike]

        for key in list(input_parameters.keys()):
            setattr(self, key, input_parameters[key])

    def build_inmodel(self):
        """
        build an inmodel file to be used as a constraint
        need to give it a dictionary containing values (list of rmin,rmax and strike) and bottom depths
        depths are the keys, resistivities are the values
        and a modeldir - needs to have the same steps as
        the model planned to run.
        """

        modelf = open(os.path.join(self.inmodel_modeldir, 'ai1mod.dat'))
        modelf.readline()

        flag = True
        model = []

        ii = 1
        while flag:
            try:
                m = [float(i) for i in modelf.readline().strip().split()]
                if len(m) > 0:
                    if ii % 2 == 0:
                        model.append(m)
                    ii += 1
            except:
                flag = False

        model = np.array(model)
        model[:, 2:] = 0.
        mvals = model[:, 2:]
        mi = model[:, 0]
        mdepths = [0.] + list(model[:, 1])
        mthick = np.array([mdepths[i + 1] - mdepths[i]
                           for i in range(len(mi))])

        keys = list(self.inmodel_dictionary.keys())

        keys.sort()
        for key in keys:
            cond = model[:, 1] >= key
            mvals[cond] = np.array(self.inmodel_dictionary[key])

        self.inmodel = np.vstack([mi, mthick, mvals.T]).T

    def write_inmodel(self, wd=None):
        """
        """

        if wd is not None:
            self.working_directory = wd

        if not hasattr(self, 'inmodel'):
            self.build_inmodel()

        np.savetxt(os.path.join(self.working_directory, 'inmodel.dat'),
                   self.inmodel,
                   fmt=['%5i', '%11.4e', '%11.4e', '%11.4e', '%11.4e'])
        print("written inmodel file to {}".format(self.working_directory))

    def read_inmodel(self):
        """
        read the inmodel file
        """

        # read in file
        inmodel = np.loadtxt(os.path.join(
            self.working_directory, self.inmodelfile))

        # convert layer thicknesses to depths
        depths = np.array([[sum(inmodel[:i, 1]), sum(inmodel[:i + 1, 1])]
                           for i in range(len(inmodel))]).flatten()
        values = np.zeros((len(inmodel) * 2, 5))
        ii = 0
        for val in inmodel:
            for i in range(2):
                values[ii] = val
                ii += 1

        self.inmodel = np.vstack(
            [values[:, 0], depths, values[:, 2], values[:, 3], values[:, 4]]).T

    def get_boundaries(self):
        """
        get points at which the resistivity changes in the inmodel file

        """

        if not hasattr(self, 'inmodel'):
            try:
                self.read_inmodel()
            except IOError:
                print("please define working directory")
                return

        data = self.inmodel

        bd = []

        for i in range(len(data) - 1):
            if data[i, 2] != data[i + 1, 2]:
                bd.append(data[i, 1])
            elif data[i, 2] != data[i + 1, 2]:
                bd.append(data[i, 1])

        self.boundary_depths = bd


class Data():
    """
    deals with input data from 1d inversions, including creating a data file
    and reading a data file afterwards to compare with inversion responses

    """

    def __init__(self, working_directory, **input_parameters):
        self.working_directory = working_directory
        self.respfile = 'ai1dat.dat'
        self.datafile = None
        self.errorfloor = np.ones([2, 2]) * 0.1
        self.errorfloor_type = 'relative'  # relative, absolute or offdiagonals
        self.edipath = None
        self.mode = 'I'

        for key in list(input_parameters.keys()):
            if hasattr(self, key):
                setattr(self, key, input_parameters[key])

        # default working directory is epath if it is specified, otherwise
        # current directory
        if self.working_directory is None:
            if self.edipath is not None:
                self.working_directory = os.path.dirname(self.edipath)
            else:
                self.working_directory = '.'

    def build_data(self):
        """
        create data to write to an input file

        """

        # read edi file to edi object
        self.edi_object = mtedi.Edi(self.edipath)

        # define z
        zr = np.real(self.edi_object.Z.z)
        # sign of imaginary component needs to be reversed for the pek1d
        # inversion code
        zi = -np.imag(self.edi_object.Z.z)
        ze = self.edi_object.Z.z_err
        z = zr + 1j * zi

        # set errorfloors
        if type(self.errorfloor) in [int, float]:
            self.errorfloor = np.ones([2, 2]) * self.errorfloor

        if self.errorfloor_type in ['relative', 'offdiagonals']:
            zer = ze / np.abs(z)
            for i in range(2):
                for j in range(2):
                    zer[:, i, j][(zer[:, i, j] < self.errorfloor[
                                  i, j])] = self.errorfloor[i, j]
            ze = np.abs(z) * zer

            if self.errorfloor_type == 'offdiagonals':
                for i in range(2):
                    for iz in range(len(z)):
                        if ze[iz, i, i] < ze[iz, i, 1 - i]:
                            ze[iz, i, i] = ze[iz, i, 1 - i]

        elif self.errorfloor_type == 'absolute':
            for i in range(2):
                for j in range(2):
                    ze[:, i, j][(ze[:, i, j] < self.errorfloor[
                                 i, j])] = self.errorfloor[i, j]

        # define header info for data file
        header = '{:>5}\n{:>5}'.format(self.mode, len(self.edi_object.Z.resistivity))

        # create data array
        data_list = [1. / self.edi_object.Z.freq]
        for i in range(2):
            for j in range(2):
                if self.mode == 'I':
                    dd = [zr, ze, zi, ze]

                for d in dd:
                    data_list.append(d[:, i, j])

        self.header = header
        self.data = np.vstack(data_list).T
        self.z = zr + 1j * zi
        self.z_err = ze

    def write_datafile(self, wd=None):
        """
        write data to file

        """

        if wd is not None:
            self.working_directory = wd

        self.build_data()

        # define format list for writing data file
        fmt = ['%14.5f'] + ['%12.5e'] * 16

        # define file name and save data file
        fname_bas = self.edi_object.station.split('_')[0]
        self.datafile = fname_bas + '.dat'
        fname = os.path.join(self.working_directory, self.datafile)

        np.savetxt(fname, self.data, fmt=fmt, header=self.header, comments='')

    def read_datafile(self):
        """
        read data file into the data object.
        calculate resistivity and phase

        """
        if self.datafile is None:

            default_files = ['ai1dat.dat', 'ai1mod.dat', 'ai1fit.dat',
                             'inmodel.dat', 'inregulm.dat']
            dlst = [i for i in os.listdir(self.working_directory) if
                    (i[-4:] == '.dat') and (i not in default_files)]
            if len(dlst) == 1:
                self.datafile = dlst[0]
            else:
                print("please define datafile")
                return

        # define path to file
        datafpath = os.path.join(self.working_directory, self.datafile)
        self.mode = open(datafpath).readline().strip().split()[0]
        data = np.loadtxt(datafpath, skiprows=2)
        self.freq = 1. / data[:, 0]

        if self.mode == 'I':
            zr = np.vstack([data[:, i]
                            for i in range(len(data[0])) if (i - 1) % 4 == 0])
            ze = np.vstack([data[:, i]
                            for i in range(len(data[0])) if (i - 2) % 4 == 0])
            zi = -np.vstack([data[:, i]
                             for i in range(len(data[0])) if (i - 3) % 4 == 0])
            z = zr + 1j * zi
            self.z = z.T.reshape(len(z[0]), 2, 2)
            self.z_err = ze.T.reshape(len(z[0]), 2, 2)

            # make a frequency array that has the same shape as z
            freq2 = np.zeros(np.shape(self.z))
            for i in range(len(freq2)):
                freq2[i, :, :] = 1. / data[:, 0][i]

#           calculate resistivity
            self.resistivity = 0.2 * np.abs(self.z)**2 / freq2

            q = np.zeros(np.shape(self.resistivity))
#            q[(zr<0)&(zi<0)] = np.pi
#            q[(zr<0)&(zi>0)] = -np.pi
            phase = np.zeros([len(self.z), 2, 2])
            res = np.zeros([len(self.z), 2, 2])
            self.resistivity_err = np.zeros([len(self.z), 2, 2])
            self.phase_err = np.zeros([len(self.z), 2, 2])

            self.q = q
            for iz in range(len(self.z)):
                for i in range(2):
                    for j in range(2):
                        phase[iz, i, j] = np.rad2deg(
                            cmath.phase(self.z[iz, i, j]))
                        res[iz, i, j] = 0.2 * \
                            np.abs(self.z[iz, i, j])**2 / self.freq[iz]
                        r_err, phi_err = MTcc.z_error2r_phi_error(
                            np.real(self.z[iz, i, j]),
                            self.z_err[iz, i, j],
                            np.imag(self.z[iz, i, j]),
                            self.z_err[iz, i, j])

                        self.resistivity_err[iz, i, j] = \
                            0.4 * np.abs(self.z[iz, i, j]) /\
                            self.freq[iz] * r_err
                        self.phase_err[iz, i, j] = phi_err

            phase[phase < -180] += 360
            self.phase = phase
            self.resistivity = res

        elif self.mode == 'R':
            res = np.vstack([data[:, i]
                             for i in range(len(data[0])) if (i - 1) % 4 == 0])
            self.resistivity = res.T.reshape(len(res[0]), 2, 2)
            res_err = np.vstack([data[:, i]
                                 for i in range(len(data[0])) if (i - 2) % 4 == 0])
            self.resistivity_err = res_err.T.reshape(len(res_err[0]), 2, 2)

            phs = np.vstack([data[:, i]
                             for i in range(len(data[0])) if (i - 3) % 4 == 0])
            self.phase = phs.T.reshape(len(phs[0]), 2, 2)
            phs_err = np.vstack([data[:, i]
                                 for i in range(len(data[0])) if (i - 4) % 4 == 0])
            self.phase_err = phs_err.T.reshape(len(phs_err[0]), 2, 2)

    def rotate(self, rotation_angle):
        """
        use mtpy.analysis.geometry to rotate a z array and recalculate res and phase

        """
        from . import pek1dclasses as pek1dc

        if not hasattr(self, 'z'):
            self.read_datafile()

        new_z = np.zeros_like(self.z)
        new_ze = np.zeros_like(self.z_err, dtype=float)

#        for iz,zarray in enumerate(self.z):
        new_z, new_ze = MTg.MTz.rotate_z(
            self.z, rotation_angle, z_err_array=self.z_err)

        self.z = new_z
        self.z_err = new_ze

        self.resistivity, self.resistivity_err, self.phase, self.phase_err = \
            pek1dc._compute_res_phase(self.z, self.z_err, self.freq)

        self.rotation_angle = rotation_angle


class Response():
    """
    deals with responses from 1d inversions
    """

    def __init__(self, wkdir, **input_parameters):

        self.working_directory = wkdir
        self.respfile = 'ai1dat.dat'
        self.misfit_threshold = 1.1
        self.station = None

        for key in list(input_parameters.keys()):
            if hasattr(self, key):
                setattr(self, key, input_parameters[key])

        self.read_respfile()

    def read_respfile(self):
        """
        read respfile into a data object
        """

        # define path to file
        respfpath = os.path.join(self.working_directory, self.respfile)
        respf = open(respfpath)

        # find out number of models
        n = 0
        for line in respf.readlines():
            if 'REG' in line:
                n += 1

        # load model responses into an array
        resp = np.genfromtxt(respfpath, skiprows=1, invalid_raise=False)
        resmod = np.vstack([resp[:, i]
                            for i in range(len(resp[0])) if (i - 1) % 2 == 0])
        phsmod = np.vstack([resp[:, i] for i in range(
            len(resp[0])) if i != 0 and (i - 2) % 2 == 0])
        period = resp[:len(resp) / n, 0]

        self.resistivity = resmod.T.reshape(n, len(resp) / n, 2, 2)
        self._phase = phsmod.T.reshape(n, len(resp) / n, 2, 2)
        self.freq = 1. / period
        zabs = np.zeros((n, len(resp) / n, 2, 2))
        for m in range(n):
            for f in range(len(self.freq)):
                zabs[m, f] = (self.resistivity[m, f] * self.freq[f] / 0.2)**0.5
        zr = zabs * np.cos(np.deg2rad(self._phase))
        zi = -zabs * np.sin(np.deg2rad(self._phase))
        self.z = zr + 1j * zi
        self.phase = -self._phase

    def rotate(self, rotation_angle):
        """
        use mtpy.analysis.geometry to rotate a z array and recalculate res and phase

        """
        from . import pek1dclasses as pek1dc

        if not hasattr(self, 'z'):
            self.read_respfile()

        new_z = np.zeros_like(self.z)
        z_err = np.zeros_like(self.z, dtype=float)

        for iz, zarray in enumerate(self.z):
            new_z[iz], ze = MTg.MTz.rotate_z(zarray, rotation_angle)

        self.z = new_z

        self.resistivity = np.zeros_like(self.z, dtype=float)
        self.phase = np.zeros_like(self.z, dtype=float)

        for iz in range(len(self.z)):
            r, re, p, pe = pek1dc._compute_res_phase(
                self.z[iz], z_err[iz], self.freq)
            self.resistivity[iz] = r
#            self.resistivity_err[iz] = re
            self.phase[iz] = p
#            self.phase_err[iz] = pe

        self.rotation_angle = rotation_angle


class Fit():
    """
    deals with outputs from 1d inversions
    """

    def __init__(self, wkdir, **input_parameters):

        self.working_directory = wkdir
        self.fitfile = 'ai1fit.dat'
        self.respfile = 'ai1dat.dat'
        self.misfit_threshold = 1.1
        self.station = None

        for key in list(input_parameters.keys()):
            if hasattr(self, key):
                setattr(self, key, input_parameters[key])

        self.read_fit()

    def find_nperiods(self):
        """
        find number of periods used in inversion
        """

        # find out number of periods
        respfpath = os.path.join(self.working_directory, self.respfile)
        respf = open(respfpath)
        respf.readline()

        n = 0
        line = respf.readline()
        while 'REG' not in line:
            line = respf.readline()
            n += 1
        self.n_periods = n - 1

    def read_fit(self):
        """
        read fit file to give structure and anisotropy penalties and penalty weights
        """
        # load the file with fit values in it
        fit = np.loadtxt(os.path.join(self.working_directory, self.fitfile))
#        print os.path.join(self.working_directory,self.fitfile)
#        print np.shape(fit)
        # find number of periods
        self.find_nperiods()

        # total misfit
        self.misfit_mean = (fit[:, 5] / (self.n_periods * 8.))**0.5

        # structure and anisotropy penalty
        self.penalty_structure = fit[:, 6]
        self.penalty_anisotropy = fit[:, 7]
        self.weight_structure = fit[:, 2]
        self.weight_anisotropy = fit[:, 4]
        self.modelno = fit[:, 0]
        self.fit = fit

    def find_bestmodel(self):
        """
        find the smoothest model that fits the data within self.misfit_threshold
        """

        self.read_fit()
        fit = self.fit

        # define parameters
        mis = self.misfit_mean
        s = self.penalty_structure / np.median(self.penalty_structure)
        a = self.penalty_anisotropy / np.median(self.penalty_anisotropy)

        # define function to minimise
        f = a * s * np.abs(a - s) / (a + s)

        # define the parameters relating to the best model
        self.params_bestmodel = fit[
            f == min(f[mis < min(mis) * self.misfit_threshold])][0]
        self.params_fittingmodels = fit[mis < min(mis) * self.misfit_threshold]


class Model():
    """
    deals with outputs from 1d inversions
    """

    def __init__(self, wkdir, **input_parameters):

        self.working_directory = wkdir
        self.modelfile = 'ai1mod.dat'
        self.respfile = 'ai1dat.dat'
        self.fitfile = 'ai1fit.dat'
        self.inmodelfile = 'inmodel.dat'
        self.datafile = None
        self.modelno = 1
        self.models = None
        self.misfit_threshold = 1.1
        self.station = None
        self.Fit = None
        self.Resp = None
        self.Data = None
        self.x = 0.
        self.y = 0.
        self.input_parameters = input_parameters

        for key in list(input_parameters.keys()):
            if hasattr(self, key):
                setattr(self, key, input_parameters[key])

        if self.station is None:
            self.station = os.path.basename(
                self.working_directory).split('_')[0]

        self.read_model()
        self.read_fit()
        self.read_response()
        self.read_datafile()
        self._calculate_fit_vs_freq()

    def read_model(self):
        """
        read all models into an array
        """

        fpath = os.path.join(self.working_directory, self.modelfile)
#        print fpath
        nlayers = 0
        flag = True
        modelf = open(fpath)
        modelf.readline()

        while flag:
            try:
                nlayers = int(modelf.readline().strip().split()[0])
            except:
                flag = False

        models = np.genfromtxt(fpath, skiprows=1, invalid_raise=False)
        self.models = models.reshape(
            0.5 * len(models) / nlayers, 2 * nlayers, 5)

    def read_fit(self):
        if self.Fit is None:
            self.Fit = Fit(self.working_directory, **self.input_parameters)

    def read_response(self):
        if self.Resp is None:
            self.Resp = Response(self.working_directory,
                                 **self.input_parameters)

    def read_datafile(self):
        if self.Data is None:

            self.Data = Data(working_directory=self.working_directory,
                             **self.input_parameters)
            self.Data.read_datafile()

    def _calculate_fit_vs_freq(self):

        misfit_real = ((np.real(
            self.Resp.z[self.modelno - 1]) - np.real(self.Data.z)) / self.Data.z_err)**2
        misfit_imag = ((np.imag(
            self.Resp.z[self.modelno - 1]) - np.imag(self.Data.z)) / self.Data.z_err)**2

        self.Fit.misfit = misfit_real + 1j * misfit_imag

    def check_consistent_strike(self, depth,
                                window=5,
                                threshold=15.):
        """
        check if a particular depth point corresponds to a consistent 
        strike direction

        """

        if self.models is None:
            self.read_model()

        # get model of interest
        model = self.models[self.modelno - 1]

        #
        depths = model[:, 1]
        closest_depth = depths[
            np.abs(depths - depth) == np.amin(np.abs(depths - depth))][0]
        cdi = list(depths).index(closest_depth)
        i1 = max(0, cdi - int(window / 2) * 2 - 1)
        i2 = min(len(model) - 2, cdi + int(window / 2) * 2 + 1)

        strikes = model[:, -1][i1:i2]

        return np.std(strikes) < threshold

    def find_max_anisotropy(self, min_depth=0.,
                            max_depth=None,
                            strike_window=5,
                            strike_threshold=10.):
        """
        find the point of maximum anisotropy in a model result within a given
        depth range. Check that the strike is stable below defined threshold

        """
        if self.models is None:
            self.read_model()
        print(self.station)
        # get model of interest
        model = self.models[self.modelno - 1]

        if max_depth is None:
            max_depth = np.amax(model[:, 1])

        # get values only between min and max depth
        model_filt = model[(model[:, 1] > min_depth) &
                           (model[:, 1] < max_depth)]

        aniso = 1. * model_filt[:, 3] / model_filt[:, 2]
        aniso_max = np.amax(aniso)
        # define an initial aniso max depth
        depth_aniso_max = model_filt[:, 1][aniso == aniso_max][0]

        i = 0
        while not self.check_consistent_strike(depth_aniso_max,
                                               window=strike_window,
                                               threshold=strike_threshold):

            aniso[aniso == aniso_max] = 1.
            aniso_max = np.amax(aniso)
            depth_aniso_max = model_filt[:, 1][aniso == aniso_max][0]
            i += 1
            if i > len(model_filt):
                print("can't get stable strike")
                break

        params = model_filt[aniso == aniso_max][0]
#        params[-1] = params[-1]%180

        self.anisotropy_max_parameters = params

    def update_location_from_file(self, xyfile, indices=[0, 999]):
        """
        updates x and y location from an xy file with format
        station x y
        can give indices to search on if the station name in the file
        is not exactly the same as defined in the model.

        """
        return


class Model_suite():
    """
    """

    def __init__(self,
                 working_directory,
                 **input_parameters):
        self.working_directory = working_directory
        self.model_list = []
        self.inmodel_list = []
        self.modelfile = 'ai1mod.dat'
        self.respfile = 'ai1dat.dat'
        self.fitfile = 'ai1fit.dat'
        self.inmodelfile = 'inmodel.dat'
        self.rotation_angle = 0
        self.modelno = 1
        self.station_list = []
        self.station_listfile = None
        self.station_search_indices = [0, 999]
        self.station_xyfile = None
        self.anisotropy_surface_file = 'model%03i_aniso_depth.dat'

        for key in list(input_parameters.keys()):
            setattr(self, key, input_parameters[key])

        if self.station_listfile is not None:
            try:
                self.station_list = [i.strip() for i in open(
                    self.station_listfile).readlines()]
            except:
                print("can't open station list file")

        if self.model_list == []:
            self.inmodel_list = []
            wd = self.working_directory
            folder_list = [os.path.join(wd, f) for f in os.listdir(
                wd) if os.path.isdir(os.path.join(wd, f))]
            if len(self.station_list) > 0:
                i1, i2 = self.station_search_indices
                folder_list2 = []
                for s in self.station_list:
                    for ff in folder_list:
                        if str.lower(os.path.basename(ff).split('_')[0][i1:i2]) == str.lower(s):
                            folder_list2.append(ff)
#                            print s
                folder_list = folder_list2
        for folder in folder_list:
            try:
                model = Model(folder)
                model.read_model()
                self.model_list.append(model)
            except IOError:
                print(folder, "model file not found")
            try:
                inmodel = Inmodel(working_directory=folder)
                inmodel.read_inmodel()
                self.inmodel_list.append(inmodel)
            except IOError:
                print("inmodel file not found")

        if self.station_xyfile is not None:
            self.update_multiple_locations_from_file()

    def get_aniso_peak_depth(self,
                             min_depth=0,
                             max_depth=None,
                             strike_threshold=10.,
                             strike_window=5):
        """
        get the min and max resistivities, depth and strike at point of maximum
        anisotropy between min and max depth.

        min and max depth can be float, integer or numpy array

        the depth is only selected if the strike is stable within parameters
        given by strike threshold and strike window.

        """

        model_params = np.zeros([len(self.model_list), 6])

        if type(min_depth) in [float, int]:
            min_depth = np.zeros(len(self.model_list)) + min_depth
        if type(max_depth) in [float, int]:
            max_depth = np.zeros(len(self.model_list)) + max_depth

        for i, model in enumerate(self.model_list):
            model.modelno = self.modelno
            model.find_max_anisotropy(min_depth=min_depth[i],
                                      max_depth=max_depth[i],
                                      strike_window=strike_window,
                                      strike_threshold=strike_threshold)
            x, y = model.x, model.y
            depth, te, tm, strike = model.anisotropy_max_parameters[1:]
            strike = strike + self.rotation_angle

            model_params[i] = x, y, depth, te, tm, strike

        self.anisotropy_max_parameters = model_params

        if '%' in self.anisotropy_surface_file:
            self.anisotropy_surface_file = self.anisotropy_surface_file % self.modelno

        np.savetxt(os.path.join(self.working_directory,
                                self.anisotropy_surface_file),
                   model_params,
                   header=' '.join(
                       ['x', 'y', 'z', 'resmin', 'resmax', 'strike']),
                   fmt=['%14.6f', '%14.6f', '%8.2f', '%8.2f', '%8.2f', '%8.2f'])

    def update_multiple_locations_from_file(self):
        """
        updates multiple x and y locations from an xy file with format
        station x y
        can give indices to search on if the station name in the file
        is not exactly the same as defined in the model.

        """
        xy = {}
        i1, i2 = self.station_search_indices

        for line in open(self.station_xyfile):
            line = line.strip().split()
            xy[str.lower(line[0])] = [float(line[1]), float(line[2])]

        for model in self.model_list:
            model.x, model.y = xy[str.lower(model.station[i1:i2])]

        self.x = np.array([m.x for m in self.model_list])
        self.y = np.array([m.y for m in self.model_list])

    def get_median_misfit(self):
        """
        """

        n = len(self.model_list)
        model_misfits = np.zeros(n)
        for m, model in enumerate(self.model_list):
            fit = Fit(model.working_directory,
                      fitfile=self.fitfile,
                      respfile=self.respfile)
            fit.read_fit()
            model_misfits[m] = fit.misfit[self.modelno - 1]

        self.model_misfits = model_misfits
        self.median_misfit = np.median(model_misfits)


def _compute_res_phase(z, z_err, freq):
    """
    calculates *resistivity*, *phase*, *resistivity_err*, *phase_err*

    values for resistivity are in in Ohm m and phase in degrees.

    """

    resistivity_err = np.zeros_like(z_err)
    phase_err = np.zeros_like(z_err)

    resistivity = np.zeros_like(z, dtype='float')
    phase = np.zeros_like(z, dtype='float')

    # calculate resistivity and phase
    for idx_f in range(len(z)):
        for i in range(2):
            for j in range(2):
                resistivity[idx_f, i, j] = np.abs(z[idx_f, i, j])**2 /\
                    freq[idx_f] * 0.2
                phase[idx_f, i, j] = math.degrees(cmath.phase(
                    z[idx_f, i, j]))

                if z_err is not None:

                    r_err, phi_err = MTcc.z_error2r_phi_error(
                        np.real(z[idx_f, i, j]),
                        z_err[idx_f, i, j],
                        np.imag(z[idx_f, i, j]),
                        z_err[idx_f, i, j])

                    resistivity_err[idx_f, i, j] = \
                        0.4 * np.abs(z[idx_f, i, j]) /\
                        freq[idx_f] * r_err
                    phase_err[idx_f, i, j] = phi_err
    return resistivity, resistivity_err, phase, phase_err
