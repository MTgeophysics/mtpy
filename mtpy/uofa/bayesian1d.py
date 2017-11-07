#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 31.07.2013
@author: LK@UofA

mtpy/uofa/bayesian1d.py

Module for handling the UofA Bayesian 1D inversion/modelling code.

"""


import os
import sys
import os.path as op
import mtpy.utils.filehandling as MTfh
import mtpy.core.edi as EDI
import mtpy.utils.exceptions as MTex
import numpy as np


def generate_input_file(edifilename, outputdir=None):

    eo = EDI.Edi()
    eo.readfile(edifilename)
    filebase = op.splitext(op.split(edifilename)[-1])[0]

    outfilename1 = '{0}_bayesian1d_z.in'.format(filebase)
    outfilename2 = '{0}_bayesian1d_zvar.in'.format(filebase)
    outdir = op.split(edifilename)[0]

    if outputdir is not None:
        try:
            if not op.isdir(outputdir):
                os.makedirs(outputdir)
                outdir = outputdir
        except:
            pass

    outfn1 = op.join(outdir, outfilename1)
    outfn2 = op.join(outdir, outfilename2)

    outfn1 = MTfh.make_unique_filename(outfn1)
    outfn2 = MTfh.make_unique_filename(outfn2)

    freqs = eo.freq

    z_array = eo.Z.z
    z_err_array = eo.Z.z_err

    if len(freqs) != len(z_array):
        raise MTex.MTpyError_edi_file('ERROR in Edi file {0} - number of '
                                      'freqs different from length of Z array'.format(eo.filename))

    sorting = np.argsort(freqs)

    outstring1 = ''
    outstring2 = ''

    for idx in sorting:
        z = z_array[idx]
        z_err = z_err_array[idx]
        f = freqs[idx]
        outstring1 += '{0}\t'.format(f)
        outstring2 += '{0}\t'.format(f)
        for i in np.arange(2):
            for j in np.arange(2):
                if np.imag(z[i % 2, (j + 1) / 2]) < 0:
                    z_string = '{0}-{1}i'.format(np.real(z[i % 2, (j + 1) / 2]),
                                                 np.abs(np.imag(z[i % 2, (j + 1) / 2])))
                else:
                    z_string = '{0}+{1}i'.format(np.real(z[i % 2, (j + 1) / 2]),
                                                 np.imag(z[i % 2, (j + 1) / 2]))

                z_err_string = '{0}'.format(z_err[i % 2, (j + 1) / 2])

                outstring1 += '{0}\t'.format(z_string)
                outstring2 += '{0}\t'.format(z_err_string)

        outstring1 = outstring1.rstrip() + '\n'
        outstring2 = outstring2.rstrip() + '\n'

    Fout1 = open(outfn1, 'w')
    Fout2 = open(outfn2, 'w')
    Fout1.write(outstring1.expandtabs(4))
    Fout2.write(outstring2.expandtabs(4))
    Fout1.close()
    Fout2.close()

    return outfn1, outfn2
