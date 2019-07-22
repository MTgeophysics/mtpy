# -*- coding: utf-8 -*-
"""
Created on Fri Aug 01 15:52:47 2014

@author: Alison Kirkby
"""
import mtpy.modeling.occam2d as o2d
import mtpy.modeling.pek1dclasses as p1dc
import numpy as np
import os
import os.path as op
from . import pek2dforward as p2d
import string
import scipy.interpolate as si
import mtpy.utils.filehandling as fh
import mtpy.core.edi as mtedi


class Model():
    """
    class for creating and reading model files

    """

    def __init__(self, working_directory, **input_parameters):
        self.working_directory = working_directory
        self.edi_directory = None
        self.occam_configfile = None
        self.parameters_model = {}
        self.parameters_model['no_sideblockelements'] = 5
        self.parameters_model['no_bottomlayerelements'] = 4
        self.parameters_model['firstlayer_thickness'] = 100

        # model depth is in km!
        self.parameters_model['model_depth'] = 100
        self.parameters_model['no_layers'] = 25
        self.parameters_model['max_blockwidth'] = 1000

        self.parameters_data = {}
        self.parameters_data['strike'] = 0.
        self.n_airlayers = 5

        self.mesh = None
        self.meshlocations_x = None
        self.meshlocations_z = None
        self.meshblockwidths_x = None
        self.meshblockthicknesses_z = None
        self.profile_easts = None
        self.profile_norths = None
        self.inversion1d_dirdict = {}
        self.inversion1d_masterdir = '.'
        self.inversion1d_modelno = 0
        self.inversion1d_imethod = 'nearest'
        self.binsize_resistivitylog10 = 1.
        self.binsize_strike = 20.
        self.build_from_1d = False
        self.rotation = 0.
        self.modelfile = 'model'

        self.anisotropy_min_depth = 0.

        self.edifiles = []

        self.Data = None

        self.modelfile = 'model'

        update_dict = {}

        # correcting dictionary for upper case keys
        input_parameters_nocase = {}
        for key in list(input_parameters.keys()):
            input_parameters_nocase[key.lower()] = input_parameters[key]

        update_dict.update(input_parameters_nocase)

        for dictionary in [self.parameters_model]:
            for key in list(dictionary.keys()):
                if key in update_dict:
                    # check if entry exists:
                    try:
                        value = float(update_dict[key])
                        dictionary[key] = value
                    except:
                        value = update_dict[key]
                        dictionary[key] = value
                        if type(value) in [str]:
                            if value.strip().lower() == 'none':
                                dictionary[key] = None

        for key in update_dict:
            try:
                value = getattr(self, key)
                if update_dict[key] is not None:
                    try:
                        value = float(update_dict[key])
                        setattr(self, key, value)
                    except:
                        value = update_dict[key]
                        setattr(self, key, value)
                        if type(value) in [str]:
                            if value.strip().lower() == 'none':
                                setattr(self, key, None)
            except:
                continue

        self.input_parameters = update_dict

        if self.edifiles == []:
            if self.edi_directory is not None:
                try:
                    self.edifiles = [op.join(self.edi_directory,
                                             f) for f in os.listdir(self.edi_directory)]
                except IOError:
                    print("failed to find edi directory")
                    pass

    def build_model(self):
        """
        build model file string
        """
        # build the mesh
        self.build_mesh()
        self.build_aircells()
        if self.build_from_1d:
            self.get_1d_results()
            self.interpolate_1d_results()
            self.bin_resistivity_values()
        else:
            self.resistivity = np.ones([len(self.meshblockthicknesses_z),
                                        len(self.meshblockwidths_x), 3]) * 100

        self.get_station_meshblock_numbers()

    def write_modelfile(self):
        self.build_modelfilestring()
        outfile = open(op.join(self.working_directory, self.modelfile), 'w')
        outfile.write(self.modelfilestring)
        outfile.close()

    def read_model(self):

        # read model file
        modelf = open(op.join(self.working_directory, self.modelfile))

        # get nx, nz, and number of air layers n_airlayers
        for i in range(3):
            modelf.readline()
        nx, nz, self.n_airlayers = [
            int(n) for n in modelf.readline().strip().split()]
        self.n_airlayers -= 1

        # get mesh cell sizes
        meshx, meshz = [], []
        while len(meshx) < nx - 1:
            meshx += [float(n) for n in modelf.readline().strip().split()]
        while len(meshz) < nz - 1:
            meshz += [float(n) for n in modelf.readline().strip().split()]
        self.meshblockwidths_x = np.array(meshx)
        self.meshblockthicknesses_z = np.array(meshz)
        self.meshlocations_x = np.array([sum(self.meshblockwidths_x[:i])
                                         for i in range(len(self.meshblockwidths_x) + 1)])
        self.meshlocations_z = np.array([sum(self.meshblockthicknesses_z[:i])
                                         for i in range(len(self.meshblockthicknesses_z) + 1)])
        self.meshlocations_x -= self.meshlocations_x[
            self.parameters_model['no_sideblockelements']]
        self.meshlocations_z -= self.meshlocations_z[self.n_airlayers + 1]
        # get model block numbers
        modelblocks = []
        while len(modelblocks) < nz - 1:
            modelblocks.append(modelf.readline().strip().split())
        modelf.readline()
        self.modelblocknums = np.array(modelblocks)

        # get resistivitiy values for each model block number
        modelfile_reslines = []
        while len(modelfile_reslines) < len(np.unique(self.modelblocknums)):
            resline = modelf.readline().strip().split()
            for r in range(1, 14):
                if (r == 1) or (r >= 8):
                    resline[r] = int(resline[r])
                else:
                    resline[r] = float(resline[r])
            modelfile_reslines.append(resline)

        self.modelfile_reslines = np.array(modelfile_reslines)

        # assign resistivities to model blocks
        data = np.zeros(list(np.shape(self.modelblocknums)) + [6])
#        print np.shape(data)
        for rv in self.modelfile_reslines:
            data[self.modelblocknums == rv[0]] = rv[2:8]
        self.resistivity = data[:, :, :3]
        self.sds = data[:, :, 3:]

        # get column numbers of stations
        modelf.readline()
        nstations = int(modelf.readline().strip())
        modelf.readline()
        station_indices = []

        while len(station_indices) <= nstations:
            station_indices += [int(n) -
                                1 for n in modelf.readline().strip().split()]
        self.station_indices = station_indices

    def build_modelfilestring(self):
        # initialise a list containing info for model file
        modelfilestring = []
        # number of periods
        period = np.unique(np.around(1. / self.Data.frequencies, 2))
        period = list(period[period > 0.])
        modelfilestring.append('%5i' % len(period))
        # add frequency string
        modelfilestring.append(p2d.create_multiple_line_string(period,
                                                               5, '%10.2f'))
        # add string giving number of cells:
        modelfilestring.append(''.join(['%5i' % i for i in [len(self.meshlocations_x),
                                                            len(self.meshlocations_z) +
                                                            self.n_airlayers,
                                                            self.n_airlayers + 1]]))

        # add strings giving horizontal and vertical mesh steps

        meshz = list(self.meshblockthicknesses_zair) + \
            list(self.meshblockthicknesses_z)
        for meshstep in [self.meshblockwidths_x, meshz]:
            modelfilestring.append\
                (p2d.create_multiple_line_string(meshstep,
                                                 10, '%5.2f'))

        # add resistivity map
        rmap = ('0' * len(self.resistivity_map[0]) + '\n') * self.n_airlayers
        rmap += '\n'.join([''.join(i) for i in self.resistivity_map])
        modelfilestring.append(rmap)

        # add number of resistivity domains (+1 to include air)
        modelfilestring.append('%5i' % (len(list(self.resistivity_dict.keys())) + 1))

        # add dictionary contents, assuming rvertical = rmax, slant and dip zero
        # first, air layer, properties always the same
        modelfilestring.append(
            '0   0     -1.00      0.00      0.00      0.00      0.00      0.00')
        # second, dictionary contents
        for key in list(self.resistivity_dict.keys()):
            rlist = self.resistivity_dict[key]
            rlist.insert(2, rlist[1])
            rlist += [0., 0.]
            if rlist[1] / rlist[0] == 1:
                aniso = '   0'
            else:
                aniso = '   1'
            modelfilestring.append(
                ''.join([key, aniso] + ['%10.2f' % i for i in rlist]))

        # append bathymetry index, at this stage only 0 allowed:
        modelfilestring.append('%5i' % 0)

        # append number of calculation points (stations):
        modelfilestring.append('%5i' % len(self.stationblocknums))

        # append rotation
        modelfilestring.append('%10.2f' % self.rotation)

        # append station blocknums
        modelfilestring.append(p2d.create_multiple_line_string(self.stationblocknums,
                                                               5, '  %03i'))
        self.modelfilestring = '\n'.join(modelfilestring)

    def build_mesh(self):
        """
        create a mesh using occam2d
        """

        # create an occam2d setup object
        so = o2d.Setup(wd=self.working_directory,
                       edi_directory=self.edi_directory,
                       edifiles=self.edifiles,
                       configfile=self.occam_configfile,
                       strike=self.parameters_data['strike'],
                       **self.parameters_model)

        so.read_edifiles(edi_dir=self.edi_directory)
        # create an occam2d data object
        so.Data = o2d.Data(edilist=self.edifiles,
                           wd=so.wd, **so.parameters_data)
        # set up meshlocations
        so.setup_mesh_and_model()
        self.stationlocations = np.array(so.Data.stationlocations) / 1000.

        # set occam mesh attributes to pek2d object
        for attribute in ['meshlocations_x', 'meshlocations_z',
                          'meshblockwidths_x', 'meshblockthicknesses_z',
                          'profile_easts', 'profile_norths', 'Data']:
            if 'mesh' in attribute:
                attvalue = np.array(getattr(so, attribute)) / 1000.
            else:
                attvalue = getattr(so, attribute)
            setattr(self, attribute, attvalue)

        for attribute in ['firstlayer_thickness', 'model_depth',
                          'no_sideblockelements', 'no_bottomlayerelements',
                          'no_layers', 'max_blockwidth']:
            self.parameters_model[attribute] = so.parameters_inmodel[attribute]

        self.meshlocations_z = np.array([0.] + list(self.meshlocations_z))

        # get block centres
        self.blockcentres_x = [np.mean(self.meshlocations_x[i:i + 2]) for i in
                               range(len(self.meshlocations_x) - 1)]
        self.blockcentres_z = [np.mean(self.meshlocations_z[k:k + 2]) for k in
                               range(len(self.meshlocations_z) - 1)]

        # remove small cells and construct meshlocations
#        self.meshblockthicknesses_z = self.meshblockthicknesses_z\
#        [self.meshblockthicknesses_z>0.005]
#        self.meshblockwidths_x = self.meshblockwidths_x\
#        [self.meshblockwidths_x>0.005]
#        mbt = self.meshblockthicknesses_z
#        self.meshlocations_z = np.array([sum(mbt[:i]) for i in range(len(mbt)+1)])
#        mbw = self.meshblockwidths_x
#        self.meshlocations_x = np.array([sum(mbw[:i]) for i in range(len(mbw)+1)])

    def build_aircells(self):
        flt = np.log10(self.parameters_model['firstlayer_thickness'] / 1000.)
        md = np.log10(self.parameters_model['model_depth'])
        mz = np.logspace(flt, md, self.n_airlayers + 1)[::-1]
        self.meshblockthicknesses_zair = [
            mz[i] - mz[i + 1] for i in range(len(mz) - 1)]
        self.meshlocations_zair = mz

    def get_station_meshblock_numbers(self):
        """

        """

#        try:
        ivals = []
        ii = 2
        for j in range(len(self.blockcentres_x[:-1])):
            for sl in self.stationlocations:
                if (sl > self.blockcentres_x[j]) & (sl <= self.blockcentres_x[j + 1]):
                    ivals.append(ii)
            ii += 1
        self.stationblocknums = ivals
#        except AttributeError:
#            print "no stationlocations, please build mesh first"

    def get_1d_results(self, split='_'):
        """
        get 1d inversion results to apply to inputs of 2d model

        """
        if not hasattr(self, 'Data'):
            self.build_mesh()
        elif self.Data is None:
            self.build_mesh()
            fh.get_pathlist()

        self.inversion1d_dirdict = \
            fh.get_pathlist(self.inversion1d_masterdir,
                            search_stringlist=self.Data.stations,
                            start_dict=self.inversion1d_dirdict,
                            split=split,
                            folder=True)
        models1d = {}

        for key in list(self.inversion1d_dirdict.keys()):
            idir = self.inversion1d_dirdict[key]
            mod = p1dc.Model(idir)
            mod.read_model()
            models1d[key] = mod.models[self.inversion1d_modelno - 1]

        self.models1d = models1d

    def interpolate_1d_results(self):
        """
        interpolate 1d inversion results onto grid

        """

        if not hasattr(self, 'models1d'):
            self.get_1d_results()

        self.resistivity = np.zeros([len(self.meshblockthicknesses_z),
                                     len(self.meshblockwidths_x), 3])
        xvals = self.stationlocations
        model_list = []

        for key in list(self.models1d.keys()):
            model_list.append(self.models1d[key][:, 2:])

        yvals = self.models1d[key][:, 1]
        points = np.array(np.meshgrid(xvals, yvals)).T.reshape(
            len(xvals) * len(yvals), 2)

        # xi needs to be block centres
        xi = np.array(np.meshgrid(self.blockcentres_x, self.blockcentres_z)).T
#        xishape = np.shape(xi)
#        xi.reshape(np.product(xishape[:-1]),xishape[-1])

        # values are a 1d array constructed by stacking individual model
        # results
        for n in range(3):
            values = np.hstack([model[:, n] for model in model_list])
            if n < 2:
                values = np.log10(values)
# print si.griddata(points,values,xi,method=self.inversion1d_imethod).T
                self.resistivity[
                    :, :, n] = 10**(si.griddata(points, values, xi, method=self.inversion1d_imethod).T)

            else:
                #                f = si.interp2d(points[:,0],points[:,1],values)
                #                self.resistivity[:,:,n] = f(xi[:,:,0],xi[:,:,1])
                self.resistivity[:, :, n] = si.griddata(
                    points, values, xi, method=self.inversion1d_imethod).T
#                rvals = f(xi[:,:,0],xi[:,:,1])
#                self.resistivity[:,:,n] = rvals.reshape(xishape[:-1])

        # index in x direction before which resistivities are null
        n = int(len(self.resistivity) / 2)
        i1 = list(np.isfinite(self.resistivity[n, :, 0])).index(True)
        # index in x direction after which resistivities are null
        i2 = -1 - \
            list(np.isfinite(self.resistivity[n, :, 0][::-1])).index(True)

        n = int(len(self.resistivity[0]) / 2)
        k1 = list(np.isfinite(self.resistivity[:, n, 0])).index(True)
        # index in z direction after which resistivities are null
        k2 = -1 - \
            list(np.isfinite(self.resistivity[:, n, 0][::-1])).index(True)

        for i in range(0, i1):
            self.resistivity[:, i] = self.resistivity[:, i1]
        for i in range(i2, 0):
            self.resistivity[:, i] = self.resistivity[:, i2]
        for k in range(0, k1):
            self.resistivity[k] = self.resistivity[k1]
        for k in range(k2, 0):
            self.resistivity[k] = self.resistivity[k2]

        self.blockcentres_x = np.array(self.blockcentres_x)
        self.blockcentres_z = np.array(self.blockcentres_z)
        self.force_isotropy()
        self.resistivity[:, :, -1] = self.resistivity[:, :, -1] % 180

    def force_isotropy(self):
        """
        force isotropy at depths shallower than anisotropy_min_depth. Clears
        up some bins for resistivity - these are limited

        """

        rnew = 1. * self.resistivity
        rmedian = 1. * np.median(self.resistivity[:, :, :2], axis=2)

        for j in range(len(rnew)):
            for k in range(len(rnew[j])):
                for i in range(2):
                    if self.blockcentres_z[j] < self.anisotropy_min_depth:
                        rnew[j, k, i] = rmedian[j, k]

        self.resistivity = rnew

    def bin_resistivity_values(self):

        rdict, rmap, rbinned = p2d.bin_results(np.log10(self.resistivity[:, :, :-1]),
                                               self.binsize_resistivitylog10)
        self.strike_std = []
        for key in list(rdict.keys()):
            rdict[key] = [10**r for r in rdict[key]]
            rdict[key].append(
                np.median(self.resistivity[:, :, -1][rmap == key]))
            self.strike_std.append(
                np.std(self.resistivity[:, :, -1][rmap == key]))

        self.resistivity_map = rmap
        self.resistivity_dict = rdict
        self.resistivity_binned = 10**rbinned


def bin_results(in_array, binsize):
    """
    binsize can be a float or numpy array of length shape(in_array)[-1]
    """

    paramdict = {}
    keys = []

    inshape = [int(i) for i in np.shape(in_array)]

    rmmvals = in_array.reshape(np.product(inshape[:-1]), inshape[-1])
    rmm_rounded = np.around(rmmvals / binsize) * binsize

    letters = '123456789:;<=>?@' + string.uppercase + string.lowercase
    i = 0

    for r, rmm in enumerate(rmm_rounded):
        if list(rmm) not in list(paramdict.values()):
            try:
                paramdict[letters[i]] = list(rmm)
                keys.append(letters[i])
                i += 1
            except IndexError:
                print("Cannot assign any more indices, try a larger binsize")
                return
        else:
            for key in list(paramdict.keys()):
                if list(rmm) == paramdict[key]:
                    keys.append(key)

    keys = np.array(keys).reshape(inshape[:-1])

    return paramdict, keys, rmm_rounded.reshape(inshape)


def create_multiple_line_string(inlist, linelength, sformat):
    linestring = ''
    for i in range(len(inlist)):
        if (i % linelength == 0) and (i > 0):
            linestring += '\n'
        linestring += sformat % inlist[i]
    return linestring


class Response():
    """
    class to contain outputs of forward modelling

    """

    def __init__(self, working_directory, **input_parameters):
        self.working_directory = working_directory
        self.edi_directory = None
        self.datafile = 'MT_TAB_ROT.DAT'
        # dict with station names as keys and number in modelfile as values
        self.station_dict = {}

        for key in list(input_parameters.keys()):
            setattr(self, key, input_parameters[key])

    def read_response(self):
        """
        """

        output = np.genfromtxt(op.join(self.working_directory,
                                       self.datafile))
        self.station_numbers = np.unique(output[:, 0])

        self.period = np.unique(output[:, 3])
        self.nperiods = len(self.period)

        z = np.zeros([len(self.station_numbers), self.nperiods, 2, 2])
        res = np.ones_like(z)
        z = z + 1j * z

        for s in range(len(self.station_numbers)):
            ii = 4
            for i in range(2):
                for j in range(2):
                    o = output[output[:, 0] == self.station_numbers[s]]
                    z[s, :, i, j] = o[:, ii] + 1j * o[:, ii + 1]
                    res[s, :, i, j] = self.period * \
                        (4e5 / np.pi) * np.abs(z[s, :, i, j])**2
                    ii += 2

        self.z = z
        self.resistivity = res
        self.phase = np.rad2deg(np.arctan(np.imag(z) / np.real(z)))

    def find_edifiles(self):
        """
        """
        search_stringlist = list(self.station_dict.keys())
        self.edifiles = fh.get_pathlist(self.edi_directory,
                                        search_stringlist=search_stringlist,
                                        split='_',
                                        extension='.edi')
        print(list(self.edifiles.keys()))

    def read_edifiles(self):
        """
        """

        self.find_edifiles()

        if len(self.edifiles) > 0:
            self.edi_objects = [mtedi.Edi(filename=efile) for efile in
                                list(self.edifiles.values())]
