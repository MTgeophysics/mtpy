# -*- coding: utf-8 -*-

"""Program to process a single station using information from files.

If you use the optimized birrp.exe and you get an error saying: didn't find all
the .bf files, something a miss.  Try running using the old version.  Something
in the optimization compiler changes something in the variable useage in the
Fortran program.

There is a beep when it finishes and can be loud depending on the computer
volume.

At the end I've added a plot section, if you want to save the plots change
the n to y, and hopefully that will work, not sure about windows 7.  If that
doesn't work you can leave it at n and save the plot when it comes up by
clicking on the save icon on the plot and save manually.

It the rotation angles are not correct or you want to change birrp parameters
change it in the processinginfofile and rerun the program.

Good luck

Jared Peacock 2011"""


#=========================================================================
# Import necessary packages
#=========================================================================

import os
import os.path as op
import sys
import mtpy.processing.birrptools as brp
import pp
import pickle
import mtpy.core.mttools as mt
import mtpy.imaging.mtplottools as mtplot
import shutil


def main():

    arglist = sys.srgv[1:]

    if len(arglist) < 5:
        sys.exit('ERROR -- provide 5 arguments: <station folders directory> <processing parameter file> <station info file> <BIRRP executable> <EDI files directory>')

    #=========================================================================
    # Input files
    #=========================================================================
    # directory where station folders are
    dirpath = op.abspath(arglist[0])

    # dirpath=r'g:\ParalanaSept2011'
    # dirpath=r'G:\University dos\Monash\Processing'
    # dirpath=r"c:\Sept2011"

    # file where all the processing parameters are, ie day, start time, end time
    # and birrp parameters like tbw, thetae, etc
    # processinginfofile=r'F:\InjectionJuly2011\Day199.txt'
    # processinginfofile=r'c:\Sept2011\Sept2011pf.txt'
    # processinginfofile=r'g:\University dos\Monash\Processing\sashapro.txt'
    # processinginfofile=r"/wolle/InjectionJuly2011/AdvPro24Hrs100Hz.txt"
    # processinginfofile=r"/wolle/InjectionJuly2011/InjectionHours.csv"
    processinginfofile = op.abspath(arglist[1])

    # file where the station info is, ie lat, long, ex, ey, notes
    # stationinfofile=r'c:\Sept2011\Sept2011Info.txt'
    # stationinfofile=r'c:\InjectionJuly2011\InjectionJuly2011Info.txt'
    stationinfofile = op.abspath(arglist[2])
    # stationinfofile=r'g:\University dos\Monash\Processing\SashaInfo.txt'

    # the location of birrp5.exe on your computer, can be the full path to the
    # executable like r"c:\BIRRP\birrp5Optimized.exe"
    # birrploc=r"c:\Peacock\PHD\BIRRP\birrp5_3pcs20E9ptsOptimized.exe"
    birrploc = op.abspath(arglist[3])
    # birrploc=r"c:\Peacock\PHD\BIRRP\birrp51lp.exe"

    # edipath=r"c:\Sept2011\EDIfiles"
    edipath = op.abspath(arglist[4])

    pstart = 0
    #=========================================================================
    # #get information from the processing file and put into a list of dictionaries
    #=========================================================================

    plst = brp.readProDict(processinginfofile, dirpath)

    #=========================================================================
    # Run in parallel
    #=========================================================================
    #
    # jobserver=pp.Server()
    #
    # if len(plst[pstart:])<jobserver.get_ncpus():
    #    jobserver.set_ncpus(len(plst[pstart:]))
    #
    # jobserver.set_ncpus(1)
    # print "Running ",jobserver.get_ncpus()," processors for runBIRRPpp"
    #
    # jobs=[]
    # for prodict in plst[pstart:]:
    #    jobs.append(jobserver.submit(brp.runBIRRPpp,
    #                                 args=(dirpath,prodict,
    #                                       stationinfofile,
    #                                       birrploc,1),
    #                                 depfuncs=(brp.writeedi,
    #                                           brp.read2c2,
    #                                           brp.bbcalfunc,
    #                                           brp.readrf,
    #                                           brp.bbconvz,
    #                                           brp.sigfigs,
    #                                           brp.scriptfilePrep,
    #                                           brp.writeScriptfile,
    #                                           brp.callBIRRP,
    #                                           brp.convertBIRRPoutputs,
    #                                           brp.writeLogfile,
    #                                           brp.writecoh,
    #                                           brp.writedat,
    #                                           brp.writeimp,
    #                                           brp.lpconvz,),
    #                                 modules=("MTtools as mt",
    #                                          "numpy as np",
    #                                          "os",
    #                                          "subprocess",
    #                                          "time",
    #                                          "datetime",
    #                                          "BIRRPTools as brp",
    #                                          "from scipy import interpolate",
    #                                          "fnmatch",
    #                                          "shutil",
    #                                          "MTPlotTools as mtplot",)))
    #
    # take outputs for each station and put them into a dictionary for next step
    # filelst=[]
    # for job in jobs:
    #    filelst.append(job())
    # jobserver.print_stats()
    #
    # jobserver.destroy()
    #
    # pfid=open(os.path.join(dirpath,'ppBFFiles.pkl'),'w')
    # pickle.dump(filelst,pfid)
    # pfid.close()

    #=========================================================================
    # Combine files, make script file, run birrp
    #=========================================================================
    # if you find that your responses are not scaled correctly, change the parameter
    # ffactor which multiplies the responses by that number. This might happen if the
    # gains are not quite right or the dipole lengths are not quite right.

    #
    flstall = []
    for ii, pdict in enumerate(plst[235:]):
        try:
            flst = brp.runBIRRPpp(dirpath, pdict, stationinfofile, birrploc,
                                  ffactor=1, edipath=edipath)
            flstall.append(flst)
        except ValueError:
            print('Did not run ', pdict['station'], pdict['day'], pdict['start'])
    #    brp.plotBFfiles(flst['edifile'],
    #                    cohfile=flst['cohfile'],
    #                    cfilelst=flst['cfilelst'],
    #                    save='y',show='n')
    #    shutil.copy(flst['edifile'],os.path.join(edipath,
    #                os.path.basename(flst['edifile'][:-4])+
    #                pdict['day']+pdict['start'][0:2])+'.edi')
    #=========================================================================
    # Plot files
    #=========================================================================

    # change save='n' to save='y' if want to save the plots, will save in a folder
    # called dirpath\plots

    # if you don't want to use the save icon in the plots you can type in the
    # interpreter plt.savefig(FullPathSaveName,fmt='pdf)
    # note that fmt can be jpg, eps, or svg

    # brp.plotBFfiles(flst['edifile'],cohfile=flst['cohfile'],save='y',show='n')

    # if this doesn't work try:
    # mtplot.plotResPhase(flst['edifile'],plotnum=2,fignum=1)
    #

if __name__ == '__main__':
    main()
