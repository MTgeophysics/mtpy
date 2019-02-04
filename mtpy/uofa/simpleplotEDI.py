#!/usr/bin/env python


import os
import sys
import os.path as op

import mtpy.core.edi as MTedi


def main():
    fn = sys.argv[1]

    if not op.isfile(fn):
        print('\n\tFile does not exist: {0}\n'.format(fn))
        sys.exit()
    saveplot = False
    if len(sys.argv) > 2:
        arg2 = sys.argv[2]
        if 's' in arg2.lower():
            saveplot = True

    fn = plotedi(fn, saveplot)


def plotedi(fn, saveplot=False, component=None):
    edi = MTedi.Edi()
    try:
        edi.readfile(fn)
    except:
        print('\n\tERROR - not a valid EDI file: {0}\n'.format(fn))
        sys.exit()

    # if saveplot is True:
    #     import matplotlib
        # matplotlib.use('Agg')

    import pylab

    lo_comps = []
    if component is not None:
        'n' in component.lower()
        try:
            if 'n' in component.lower():
                lo_comps.append('n')
        except:
            pass
        try:
            if 'e' in component.lower():
                lo_comps.append('e')
        except:
            pass
    if len(lo_comps) == 0:
        lo_comps = ['n', 'e']

    res_te = []
    res_tm = []
    phi_te = []
    phi_tm = []
    reserr_te = []
    reserr_tm = []
    phierr_te = []
    phierr_tm = []

    for r in edi.Z.resistivity:
        res_te.append(r[0, 1])
        res_tm.append(r[1, 0])

    for p in edi.Z.phase:
        phi_te.append(p[0, 1] % 90)
        phi_tm.append(p[1, 0] % 90)

    if pylab.np.mean(phi_te) > 90 and pylab.np.mean(phi_tm) > 90:
        phi_te = [i % 90 for i in phi_te]
        phi_tm = [i % 90 for i in phi_tm]

    for r in edi.Z.resistivity_err:
        reserr_te.append(r[0, 1])
        reserr_tm.append(r[1, 0])
    for p in edi.Z.phase_err:
        phierr_te.append(p[0, 1])
        phierr_tm.append(p[1, 0])

    periods = 1. / edi.freq

    resplotelement_xy = None
    resplotelement_yx = None

    axes = pylab.figure('EDI ' + fn)
    ax1 = pylab.subplot(211)

    if 'n' in lo_comps:
        resplotelement_xy = pylab.errorbar(
            periods, res_te, reserr_te, marker='x', c='b', fmt='x')
    if 'e' in lo_comps:
        resplotelement_yx = pylab.errorbar(
            periods, res_tm, reserr_tm, marker='x', c='r', fmt='x')
    pylab.xscale('log', nonposx='clip')
    pylab.yscale('log', nonposy='clip')
    minval = pylab.min(pylab.min(res_te, res_tm))
    maxval = pylab.max(pylab.max(res_te, res_tm))
    pylab.xlim(0.5 * pylab.min(periods), 2 * pylab.max(periods))

    # ylim([0.1,100])
    pylab.ylim([minval / 10, maxval * 10])

    pylab.autoscale(False)

    pylab.ylabel(r' $\rho$ (in $\Omega m$)')
    pylab.setp(ax1.get_xticklabels(), visible=False)
    # share x only
    ax2 = pylab.subplot(212, sharex=ax1)
    pylab.autoscale(False)

    # ylim(-45,135)
    if 'n' in lo_comps:
        pylab.errorbar(periods, phi_te, phierr_te, marker='x', c='b', fmt='x')
    if 'e' in lo_comps:
        pylab.errorbar(periods, phi_tm, phierr_tm, marker='x', c='r', fmt='x')
    pylab.ylabel('Phase angle ($\degree$)')
    pylab.xlabel('Period (in s)')
    pylab.plot([pylab.xlim()[0], pylab.xlim()[1]], [45, 45], '-.', c='0.7')
    pylab.ylim([-0, 90])

    ax1.legend([resplotelement_xy, resplotelement_yx], ['$E_{X}/B_Y$', '$E_Y/B_X$'], loc=2, ncol=1,
               numpoints=1, markerscale=0.8, frameon=True, labelspacing=0.3,
               prop={'size': 8}, fancybox=True, shadow=False)

    pylab.tight_layout()
    if saveplot is True:

        pylab.ioff()
        outfn = op.splitext(fn)[0] + '.png'
        pylab.savefig(outfn, bbox_inches='tight')
        pylab.close('all')
        pylab.ion()
        return outfn

    else:
        pylab.ion()
        pylab.show(block=True)

        return None


if __name__ == '__main__':
    main()
