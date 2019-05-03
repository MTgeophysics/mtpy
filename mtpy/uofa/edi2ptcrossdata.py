#!/usr/bin/env python
"""
    MTpy Script module

    edi2crossdata.py


    Read EDI file(s) and extract information needed to generate a PhaseTensorCross plot
    representation of data.

    Errors/uncertainties are either calculated by theoretical propagation of errors, or by
    statistical evaluation: a number of PTs is generated from the values of Z and its uncertainties. For each value of Z
    random numbers for its replacement are drawn from a normal distribution. The value of Zerr determines the sigma of the
    distribution.

    From these realisations the PT parameters are taken, and their respective means and standard deviations
    are interpreted as final outputs.

    This file contains a function with the aforementioned functionality as well as the required
    wrapper to run as a shell script.

    Actual PT Cross Plot is then generated using GMT
    (to be included later)

    Script usage
    ------------
    Input:
        - <EDI file> : file name
        - <output directory> : relative or absolute path to put the output file
        - <#iterations> : positive integer number of iterations;
                            value 0 results in straight error propagation (no statistics!)
        Optional arguments:
            - <sigma scaling> : the factor by which to multiply the Zerr value to obtain sigma
                                of the normal distribution; default = 1
            - batch process flag "-b"


    Output:
        - 1 file containing evaluated PT data


    @MTpy 2014, UofA (LK)

"""


import os
import sys
import os.path as op

import numpy as np
import mtpy.core.edi as MTedi
import mtpy.analysis.pt as MTpt

# for debugging:
import ipdb


def main():

    if len(sys.argv) < 4:
        print('\nNeed at least 2 arguments: <EDI file> '\
            '<output directory> <#iterations> \n\n'\
            'Optional arguments: \n [sigma scaling]\n'\
            ' [batch process flag "-b"] \n\n')
        return

    try:
        fn_in = sys.argv[1]
        fn_in = op.join(op.abspath(os.curdir), fn_in)
        edi_object = MTedi.Edi(filename=fn_in)
    except:
        print('\n\tERROR - File is not a valid EDI file: {0}\n'.format(fn_in))
        sys.exit()

    try:
        outdir = sys.argv[2]
        outdir = op.join(op.abspath(os.curdir), outdir)
        if not op.isdir(outdir):
            os.makedirs(outdir)
    except:
        print('\n\tERROR - Output directory does not exist and cannot be'\
            ' generated: {0}\n'.format(outdir))
        sys.exit()

    try:
        n_iterations = int(float(sys.argv[3]))
        if n_iterations < 0:
            raise
    except:
        print('\n\t ERROR - number of iterations must be a positive integer (incl. 0)\n')
        sys.exit()

    fn_out = None
    sigma_scaling = 1
    if len(sys.argv) > 4:
        try:
            sigma_scaling = float(sys.argv[4])
            if sigma_scaling <= 0:
                raise
        except:
            sigma_scaling = 1
            print('\nWARNING - Invalid sigma scale ..using 1 instead\n')

    try:
        print()
        print('Generating PTcross input data file')
        if n_iterations != 0:
            print('(evaluating {0} realisations of Z)'.format(n_iterations))
        print('\t...')

        generate_ptcrossdata_file(
            edi_object,
            n_iterations,
            sigma_scaling,
            outdir,
            fn_out)
    except:
        raise
        print('\n\tERROR - could not generate PT Cross Data file - check EDI file!\n')


def generate_ptcrossdata_file(
        edi_object, n_iterations, sigma_scaling, outdir, outfn=None):

    freqs = edi_object.freq

    station = edi_object.station
    # no spaces in file names:
    if len(station.split()) > 1:
        station = '_'.join(station.split())

    # Define and check validity of output file
    if outfn is None:
        fn = '{0}_PTcrossdata'.format(station)
        outfn = op.join(outdir, fn)

    outfn = op.realpath(outfn)

    try:
        Fout = open(outfn, 'w')
    except:
        print('\n\tERROR - Cannot generate output file!\n')
        raise
    if n_iterations == 0:
        Fout.write('# {0}   {1:+010.6f}   {2:+011.6f}\n'.format(station,
                                                                edi_object.lat, edi_object.lon))
    else:
        Fout.write('# {0}   {1:+010.6f}   {2:+011.6f} \t\t statistical evaluation of {3} realisations\n'.format(
            station, edi_object.lat, edi_object.lon, abs(int(n_iterations))))
    headerstring = '# lat \t\t lon \t\t freq \t\t Pmin  sigma \t Pmax  sigma \t alpha  '\
        'sigma \t beta  sigma \t ellipticity  \n'
    Fout.write(headerstring)

    if n_iterations == 0:

        pt = MTpt.PhaseTensor(z_object=edi_object.Z, freq=freqs)

        a = pt.alpha
        b = pt.beta

        phimin = pt.phimin[0]
        phiminerr = pt.phimin[1]
        phimax = pt.phimax[0]
        phimaxerr = pt.phimax[1]

        #e = pt.ellipticity
        #e = (pmax-pmin)/(pmax+pmin)

        for i, freq in enumerate(edi_object.freq):
            try:
                e = (phimax[i] - phimin[i]) / (phimax[i] + phimin[i])
                vals = '{10:.4f}\t{11:.4f}\t{0:.4e}\t{1: 3.2f}\t{2:3.2f}\t{3: 3.2f}\t{4:3.2f}\t{5: 3.2f}\t{6:3.2f}'\
                    '\t{7: 3.2f}\t{8:3.2f}\t{9:.3f}\n'.format(
                        freq, phimin[i], phiminerr[i], phimax[
                            i], phimaxerr[i], a[0][i] % 90, a[1][i] % 90,
                        b[0][i], b[1][i], e, edi_object.lat, edi_object.lon)
                Fout.write(vals)
            except:
                raise
                continue

        Fout.close()
        print('\n\t Done - Written data to file: {0}\n'.format(outfn))
        return

    # for all values n_iterations !=0 loop over abs(n_iterations) #
    # loops individual per frequency

    for idx, f in enumerate(freqs):

        z = edi_object.Z.z
        z_err = edi_object.Z.z_err

        lo_pts = []
        lo_ptserr = []

        cur_z = z[idx]
        cur_z_err = z_err[idx]

        # crude check for 'bad' values in z_err:
        for i in np.arange(4):
            a = cur_z_err[i / 2, i % 2]
            val = cur_z[i / 2, i % 2]

            try:
                dummy = float(a)
                if np.isnan(a) or np.isinf(a):
                    raise
            except:
                # correct by nearest neighbour value
                print(f, '\ncorrecting error value', a, '...', end=' ')
                rel_errs = []
                if idx + 1 != len(freqs):
                    next_z = z[idx + 1][i / 2, i % 2]
                    next_z_err = z_err[idx + 1][i / 2, i % 2]
                    rel_err = next_z_err / abs(next_z)
                    rel_errs.append(rel_err)
                if idx != 0:
                    last_z = z[idx - 1][i / 2, i % 2]
                    last_z_err = z_err[idx - 1][i / 2, i % 2]
                    rel_err = last_z_err / abs(last_z)
                    rel_errs.append(rel_err)
                rel_err_final = np.mean(rel_errs)
                err = rel_err_final * abs(val)
                cur_z_err[i / 2, i % 2] = err
                # print err

        # calculate random numbers:
        lo_rands = []
        for k in np.arange(4):
            randnums = sigma_scaling * \
                cur_z_err[k / 2, k %
                          2] * np.random.randn(2 * abs(int(n_iterations)))
            lo_rands.append(randnums)

        # Loop over |n_iterations| random realisations:

        lo_pmin = []
        lo_pmax = []

        lo_alphas = []
        lo_betas = []
        #lo_ellipticities = []

        # print 'running {0} iterations for {1} Hz'.format(n_iterations,f)
        for run in np.arange(abs(int(n_iterations))):
            tmp_z = np.array(
                [[complex(np.real(cur_z[0, 0]) + lo_rands[0][run],
                          np.imag(cur_z[0, 0]) + lo_rands[0][-run - 1]),
                  complex(np.real(cur_z[0, 1]) + lo_rands[1][run],
                          np.imag(cur_z[0, 1]) + lo_rands[1][-run - 1])],
                 [complex(np.real(cur_z[1, 0]) + lo_rands[2][run],
                          np.imag(cur_z[1, 0]) + lo_rands[2][-run - 1]),
                  complex(np.real(cur_z[1, 1]) + lo_rands[3][run],
                          np.imag(cur_z[1, 1]) + lo_rands[3][-run - 1])
                  ]]
            )
            tmp_pt = MTpt.PhaseTensor(
                z_array=tmp_z.reshape(
                    1, 2, 2), freq=f.reshape(1))
            pi1 = tmp_pt._pi1()[0]
            pi2 = tmp_pt._pi2()[0]
            lo_pmin.append(pi2 - pi1)
            lo_pmax.append(pi2 + pi1)

            alpha = tmp_pt.alpha[0][0]
            if alpha < 0 and alpha % 90 < 10:
                lo_alphas.append(alpha % 90 + 90)
            else:
                lo_alphas.append(alpha % 90)

            # if idx%10==0:
            #     print alpha,lo_alphas[-1]

            lo_betas.append(tmp_pt.beta[0][0])
            # in percent:
            # lo_ellipticities.append(100*tmp_pt.ellipticity[0][0])

        lo_alphas = np.array(lo_alphas)
        a = np.median(lo_alphas)
        aerr = np.median(np.abs(lo_alphas - a))
        #aerr = np.std(lo_alphas)
        # if idx%10==0:
        #     print '\t',a,aerr
        #     print

        b = np.mean(lo_betas)
        berr = np.std(lo_betas)

        # ipdb.set_trace()

        # convert to angles:
        phimin = np.mean([np.degrees(np.arctan(i)) for i in lo_pmin])
        phiminerr = np.std([np.degrees(np.arctan(i)) for i in lo_pmin])

        phimax = np.mean([np.degrees(np.arctan(i)) for i in lo_pmax])
        phimaxerr = np.std([np.degrees(np.arctan(i)) for i in lo_pmax])

        #e = np.mean(lo_ellipticities)
        #eerr = np.std(lo_ellipticities)
        e = (phimax - phimin) / (phimax + phimin)

        try:
            vals = '{10:.4f}\t{11:.4f}\t{0:.4e}\t{1: 3.2f}\t{2:3.2f}\t{3: 3.2f}\t{4:3.2f}\t{5: 3.2f}\t{6:3.2f}'\
                '\t{7: 3.2f}\t{8:3.2f}\t{9:.3f}\n'.format(
                    f, phimin, phiminerr, phimax, phimaxerr, a, aerr, b, berr, e, edi_object.lat, edi_object.lon)
            Fout.write(vals)
        except:
            raise
            continue

    Fout.close()
    print('\n\t Done - Written data to file: {0}\n'.format(outfn))
    return


if __name__ == '__main__':
    main()
