# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 10:21:05 2014

@author: Alison Kirkby

"""
from __future__ import print_function
import os
import os.path as op
import mtpy.utils.filehandling as fh
#import mtpy.utils.elevation_data as mted
import pek1dclasses as pek1dc
from sys import argv
from subprocess import call
import time
import numpy as np


def parse_arguments(arguments):
    """
    takes list of command line arguments obtained by passing in sys.argv
    reads these and returns a parser object
    """

    import argparse

    parser = argparse.ArgumentParser(
        description='Set up and run a set of 1d anisotropic model runs')
    parser.add_argument('-l', '--program_location',
                        help='path to the inversion program',
                        type=str, default=r'/home/547/alk547/aniso1d/ai1oz_ak')
    parser.add_argument('-r', '--run_input', nargs=7,
                        help='command line input for the inversion program',
                        type=float, default=[1, 0, 0.1, 40, 1.05, 1, 0])
    efhelp = 'error floor for impedence tensor or resisitivity values, provide 1,2 or 4 values.\n'
    efhelp += '1 value: same errorfloor applied to all 4 components of impedance tensor\n'
    efhelp += '2 values: first value applied to diagonals (xx and yy), second value applied to off diagonals\n'
    efhelp += '4 values: values applied in order to xx, xy, yx, yy'
    efhelp += 'if 3 values are provided then first 2 are taken. If > 4 are provided then first 4 are taken\n'

    parser.add_argument('-ef', '--errorfloor',
                        help=efhelp, nargs='*',
                        type=float, default=0.1)
    parser.add_argument('-eft', '--errorfloor_type',
                        help='type of error floor, absolute, relative or offdiagonals',
                        type=str, default='relative')
    parser.add_argument('-wd', '--working_directory',
                        help='working directory',
                        type=str, default='.')
    parser.add_argument('-el', '--edifolder_list', nargs='*',
                        help='list of folders containing edi files to use, full path or relative to working directory',
                        type=str, default=None)
    parser.add_argument('-ei', '--edifolder_identifier',
                        help='identifying string contained in folders of interest',
                        type=str, default='')
    parser.add_argument('-m', '--mode',
                        help='mode to put in data file, impedence (I) or resistivity (R) and phase',
                        type=str, default='I')
    parser.add_argument('-ps', '--penalty_type_structure',
                        help='number describing type of structure penalty',
                        type=int, default=6)
    parser.add_argument('-pa', '--penalty_type_anisotropy',
                        help='number describing type of anisotropy penalty',
                        type=int, default=2)
    parser.add_argument('-pws', '--penalty_weight_structure', nargs=3,
                        help='structure penalty weights to apply in the inversion, provide log10(minimum penalty weight), log10(maximum), number of values',
                        type=float, default=[0, 2, 3])
    parser.add_argument('-pwa', '--penalty_weight_anisotropy', nargs=3,
                        help='anisotropy penalty weights to apply in the inversion, provide log10(minimum penalty weight), log10(maximum), number of values',
                        type=float, default=[0, 2, 3])
    parser.add_argument('-imax', '--iteration_max',
                        help='maximum number of iterations',
                        type=int, default=100)
    parser.add_argument('-i', '--build_inmodel',
                        help='build inmodel, True or False',
                        type=bool, default=False)
    parser.add_argument('-ip', '--inmodel_parameters_file',
                        help='full path (or path relative to working directory) to file containing inmodel parameters',
                        type=str)
    parser.add_argument('-id', '--inmodel_modeldir',
                        help='full path to an output model file from previous run containing layer depths',
                        type=str)
    parser.add_argument('-s', '--master_savepath',
                        help='master directory to save suite of runs into',
                        default='inversion_suite')

    args = parser.parse_args(arguments)
    args.working_directory = os.path.abspath(args.working_directory)
    #args.run_input = args.run_input[0]
    for i in [0, 1, 3, 5, 6]:
        args.run_input[i] = int(args.run_input[i])

    if np.iterable(args.errorfloor):
        if len(args.errorfloor) == 1:
            args.errorfloor = args.errorfloor[0]
        elif (len(args.errorfloor) == 2) or (len(args.errorfloor) == 3):
            ef = args.errorfloor[:2]
            args.errorfloor = np.array([ef, ef[::-1]])
        elif len(args.errorfloor) == 4:
            args.errorfloor = np.reshape(args.errorfloor, [2, 2])

    return args


def create_inmodel_dictionary_from_file(input_file,
                                        x, y,
                                        working_directory=None):
    """
    update inmodel dictionary to get elevation details from file

    ------------------------------Parameters-----------------------------------    
    **input_file** full path to a csv file containing list of following parameters:
    elevation filename,offset,resmin,resmax,strike
    where:
    elevation filename = Full path to x y z file containing elevations of the 
                         constraining layer to put into the inversions, put
                         none if providing a constant elevation.
                         Numbers are converted to absolute values internally.
    offset = Constant depth of constraining layer, if provided in addition to
             elevation filename the offset is added/subtracted from depth 
             from file, positive down.
    resmin, resmax, strike = minimum and maximum resistivity values and strike 
                             of minimum resistivity for constraining layer
    **x** x position of station in same coordinate system as elevation file
    **y** y position of station in same coordinate system as elevation file

    """

    inmodel_dict = {}
    inmodel_list = []

    if working_directory is None:
        working_directory = os.path.abspath('.')

    for line in open(input_file).readlines()[1:]:
        line = line.strip().split(',')
        if str.lower(line[0]) != 'none':
            elevfn = os.path.join(working_directory, line[0])
    #        print "elevfn",elevfn
            try:
                elev = np.abs(mted.get_elevation(x, y, elevfn) / 1000.)
            except IOError:
                print("File not found, set elevation to zero instead")
                elev = 0.0
        else:
            elev = 0.0
        params = [float(pp) for pp in line[1:]]
   #     print elev,params
        inmodel_list.append([round(elev + params[0], 2), params[1:]])
    #print(x, y, "inmodel_list", inmodel_list, end=' ')  # end syntx error
    print(x, y, "inmodel_list", inmodel_list)
    i = 0
    while i < len(inmodel_list) - 1:
        print(i, end=' ')
        if inmodel_list[i][0] > inmodel_list[i + 1][0]:
            print("remove")
            inmodel_list.remove(inmodel_list[i])
        i += 1
    print("inmodel_list", inmodel_list, end=' ')
    for item in inmodel_list:
        try:
            print("item[0],item[1]", item[0], item[1], end=' ')
            inmodel_dict[item[0]] = item[1]
        except:
            print("couldn't assign value to dictionary")

    # print "inmodel_dict",inmodel_dict
    return inmodel_dict


def create_filelist(wd, subfolder_list=None, subfolder_identifier=None):
    """
    create a list of full paths to edi files    

    """

    edi_list = []

    if subfolder_list is None:
        subfolder_list = [folder for folder, sf,
                          f in os.walk(wd) if folder != wd]
    if subfolder_identifier is not None:
        subfolder_list = [
            f for f in subfolder_list if subfolder_identifier == op.basename(f)]

    for subfolder in subfolder_list:
        # print subfolder
        epath = os.path.join(wd, subfolder)
        edi_list += [os.path.join(epath, ff)
                     for ff in os.listdir(epath) if ff[-4:] == '.edi']

    return edi_list


def update_inputs():
    """
    update input parameters from command line

    """

    args = parse_arguments(argv[1:])
    cline_inputs = {}
    cline_keys = [i for i in dir(args) if i[0] != '_']

    for key in cline_keys:
        cline_inputs[key] = getattr(args, key)

    return cline_inputs


def generate_inputfiles(epath, **input_parameters):
    """
    generate input files for a model. 

    -----------------------Compulsory parameter--------------------------------
    **epath** the full path to the edi file.

    -----------------------Recommended parameter-------------------------------
    **wd** working directory, default is the edi directory. A new directory
           is created under this directory to put all the input files into

    ------------------------Optional Parameters--------------------------------
    **datafile** name for the input file, if not specified, name is taken from
                 the edi file
    **errorfloor_z** error floor for the input z values, can be an absolute
                     value or relative (e.g. 0.1 means 10%)
                     default is 0.1
    **errorfloor_type** type of error floor, either 'relative' or 'absolute'
                        default is relative.
    **type_struct** type of structure penalty, default is 6
    **type_aniso** type of anisotropy penalty, default is 2
    **value_struct** structural penalty weights to apply, default is [1,10,100]
    **value_aniso** anisotropy penalty weights to apply, default is [1,10,100]
    **imax** maximum number of iterations to run, default is 100


    to generate an a priori (inmodel) file, need to put keyword
    **build_inmodel** = True, default is False

    also need to specify the following parameters:
    **inmodel_vals**



    inmodel_modeldir = string, folder containing previous model run with same 
    resolution, necessary for constructing the layer depths in the inmodel file.
    inmodel_vals = dictionary structured as follows:
    {layer top depth:[minimum_resistivity, maximum_resistivity, strike]}

    """
    from . import pek1d

    data_kwds = ['working_directory', 'datafile', 'errorfloor',
                 'errorfloor_type', 'edipath', 'mode']
    control_kwds = ['penalty_type_structure', 'penalty_type_anisotropy',
                    'penalty_weight_structure', 'penalty_weight_anisotropy',
                    'iteration_max']
    inmodel_kwds = ['inmodel_dictionary']

    data_inputs = {'edipath': epath}
    control_inputs = {}
    inmodel_inputs = {}

    build_inmodel = False
    for key in list(input_parameters.keys()):
        if key in data_kwds:
            data_inputs[key] = input_parameters[key]
        if key in control_kwds:
            control_inputs[key] = input_parameters[key]
        if key in inmodel_kwds:
            inmodel_inputs[key] = input_parameters[key]
        if key == 'build_inmodel':
            build_inmodel = input_parameters[key]

    for pw in ['penalty_weight_structure', 'penalty_weight_anisotropy']:
        min, max, n = control_inputs[pw]
        control_inputs[pw] = np.logspace(min, max, n)
    Data = pek1dc.Data(**data_inputs)
    Data.build_data()

    # make a save path to match the edi file
    wd = input_parameters['working_directory']
    sp = input_parameters['master_savepath']
    savepath = fh.make_unique_folder(os.path.join(wd, sp),
                                     os.path.basename(Data.edipath).split('_')[0] + Data.mode)
    os.mkdir(savepath)
    Data.write_datafile(wd=savepath)

    # update the working directory to the new savepath
    control_inputs['working_directory'] = savepath
    inmodel_inputs['working_directory'] = savepath

    Ctl = pek1dc.Control(**control_inputs)
    Ctl.write_ctlfile()
    print(os.path.basename(Data.working_directory), build_inmodel)
    if build_inmodel:
        if 'inmodel_modeldir' in list(input_parameters.keys()):
            print(os.path.basename(Data.working_directory), end=' ')
            inmodel_dict = pek1d.create_inmodel_dictionary_from_file(input_parameters['inmodel_parameters_file'],
                                                                     Data.edi_object.lon, Data.edi_object.lat,
                                                                     working_directory=data_inputs['working_directory'])
            print(inmodel_dict)
            Inmodel = pek1dc.Inmodel(inmodel_modeldir=input_parameters['inmodel_modeldir'],
                                     inmodel_dictionary=inmodel_dict,
                                     **inmodel_inputs)
            Inmodel.write_inmodel()

    return Data


def build_run():
    """
    build input files and run a suite of models
    runs one inversion per processor, make sure you have enough processors!

    """
    try:
        from mpi4py import MPI
        mpi_import = True
    except:
        mpi_import = False

    # get command line arguments as a dictionary
    input_parameters = update_inputs()

    # categorise inputs
    build_parameters = ['working_directory', 'datafile', 'errorfloor',
                        'errorfloor_type', 'mode',
                        'penalty_type_structure', 'penalty_type_anisotropy',
                        'penalty_weight_structure', 'penalty_weight_anisotropy',
                        'iteration_max', 'inmodel_parameters_file', 'build_inmodel',
                        'inmodel_modeldir']

    # establish the rank of the computer
    if mpi_import:
        rank = MPI.COMM_WORLD.Get_rank()
    else:
        rank = 0

    # create a list of edi files to model
    edi_list = create_filelist(input_parameters['working_directory'],
                               subfolder_list=input_parameters[
                                   'edifolder_list'],
                               subfolder_identifier=input_parameters['edifolder_identifier'])
    # print "edi_list",edi_list
    # print 'working_directory',input_parameters['working_directory']
    # print 'edifolder_list',input_parameters['edifolder_list']
    # print 'edifolder_identifier',input_parameters['edifolder_identifier']
    # update input parameters for building of model
    build_inputs = {}
    for key in build_parameters:
        try:
            build_inputs[key] = input_parameters[key]
        except:
            pass

    # make a master directory under the working directory to save all runs into
    master_directory = os.path.join(
        input_parameters['working_directory'], input_parameters['master_savepath'])
    if rank == 0:
        if not os.path.exists(master_directory):
            os.mkdir(master_directory)

    build_inputs['master_savepath'] = master_directory
    # wait til master directory is made until progressing
    print("waiting for directory")
    while not os.path.isdir(master_directory):
        time.sleep(1)
        print('.', end=' ')

    # build a model
    time.sleep(rank)
    Data = generate_inputfiles(edi_list[rank], **build_inputs)
    time.sleep(5)
    os.chdir(Data.working_directory)

    # run the model
    print("running model on cpu number {} from directory {}".format(rank, Data.working_directory))
    print("current directory, {}".format(os.getcwd()))
    call([input_parameters['program_location']] + [Data.datafile] + [str(n)
                                                                     for n in input_parameters['run_input']])


if __name__ == '__main__':
    build_run()
