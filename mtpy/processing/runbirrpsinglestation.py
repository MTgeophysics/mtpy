# -*- coding: utf-8 -*-

"""Program to process a single station using information from files. 

If you use the optimized birrp.exe and you get an error saying: didn't find all
the .bf files, something a miss.  Try running using the old version.  Something
in the optimization compiler changes something in the variable useage in the
Fortran program.  

At the end I've added a plot section, if you want to save the plots change
the n to y, and hopefully that will work, not sure about windows 7.  If that 
doesn't work you can leave it at n and save the plot when it comes up by
clicking on the save icon on the plot and save manually.

It the rotation angles are not correct or you want to change birrp parameters
change it in the processinginfofile and rerun the program.

--if you want to change parameters in plst[ii] from the command line use:

    plst[ii][parameter]=newvalue
    
    where parameter can be anything in plst[ii].keys() (type that in to the 
    command line to see what parameters are available)
    
    then run flst=brp.runBIRRPpp(dirpath,plst[stationindex],stationinfofile,
                                 birrploc,ffactor=1)
                        
    from the command line.  If this makes it run or makes the parameters better
    than change it in the processinginfofile for later use.

--if you want to look at the time frequency plot:
    your should decimate the data down to something smaller than 50Hz
    import TFtools as tf
    
    bx=np.loadtxt(filename)
    bx=tf.decimatef(bx,decimation factor)
    psm,tsm,fsm,pst=tf.smethod(bx,L=11,nh=2**8,nfbins=2**9,tstep=2**6,
                               df=decimation frequency)
    
    where decimation frequency the frequency to which you decimated 
    so if you sampled at 500 Hz and you decimated to 50 Hz then df=50
    
    to plot the spectrogram
    tf.plottf(psm,tsm,fsm)
    
    there are a whole lot of parameters that you can change for plotting, 
    type in tf.plottf? on the command line and it should give you all the 
    parameters to change.
        

Good luck

Jared Peacock 2011"""


#===============================================================================
# Import necessary packages
#===============================================================================

import os
import os.path as op
import sys
import mtpy.processing.birrptools as brp 
import mtpy.core.mttools as mt
import mtpy.imaging.mtplottools as mtplot
import mtpy.core.z as Z


def main():
    #===============================================================================
    # Input files
    #===============================================================================


    arglist = sys.argv[1:]
    if len(arglist) < 4:
    sys.exit('ERROR - provide 4 arguments:<station folder directory> <processingparameter file> <stationparameter file> <BIRRP executable>')    
    
    #directory where station folders are
    dirpath = op.abspath(arglist[0])

    #file where all the processing parameters are, ie day, start time, end time
    #and birrp parameters like tbw, thetae, etc    
            
    processinginfofile = op.abspath(arglist[1])

    #file where the station info is, ie lat, long, ex, ey, notes

    stationinfofile = op.abspath(arglist[2])

    #the location of birrp5.exe on your computer, can be the full path to the 
    #executable like r"c:\BIRRP\birrp5Optimized.exe"
    #birrploc=r"c:\Peacock\PHD\BIRRP\birrp5_3pcs20E9ptsOptimized.exe"

    birrploc = op.abspath(arglist[3])

    #this is the index of which station to process which corresponds to the
    #line number in Notepad++ minus 2 of the processinginfofile.  So if you want 
    #to process the first station in processinginfofile which is line 2 in the 
    #notepad file, the statinindex will be 0.  

    stationindex=0
    #===============================================================================
    # #get information from the processing file and put into a list of dictionaries
    #===============================================================================

    plst=brp.readProDict(processinginfofile,dirpath)


    #===============================================================================
    # Combine files, make script file, run birrp
    #===============================================================================
    #if you find that your responses are not scaled correctly, change the parameter
    #ffactor which multiplies the responses by that number. This might happen if the
    #gains are not quite right or the dipole lengths are not quite right.

    #flst=brp.runBIRRPpp(dirpath,plst[stationindex],stationinfofile,birrploc,
    #                    ffactor=1)
    #                    
    #if you want to run multiple stations, one after the other uncomment the
    #following loop.  This will processes the station then plot the apparent
    #resistivity and phase of all 4 components, then plot the phase tensor 
    #components.  If you want to start plst from a different index, because you 
    #keep adding to the processinginfofile for each day, which I suggest doing so 
    #when you come back from the field all the info is one place, just change
    #the plst in enumrate(plst,1) to plst[start:stop] or plst[start:] for all 
    #stations after start.

    flstall=[]
    for ii,pdict in enumerate(plst,1):
        try:    
            flst=brp.runBIRRPpp(dirpath,pdict,stationinfofile,birrploc,
                            ffactor=1)
            flstall.append(flst)
            brp.plotBFfiles(flst['edifile'],cohfile=flst['cohfile'],save='y',
                            show='n')
    #        z1=Z.Z(flst['edifile'])
    #        z1.plotResPhase(fignum=ii,plottype=2)
    #        z1.plotPTAll(fignum=ii+len(plst))
        except TypeError:
            print 'Did not process ',pdict['station']
        except IOError:
            print 'Did not process ',pdict['station']
        except IndexError:
            print 'Did not process ',pdict['station']
        except ValueError:
            print 'Did not process ',pdict['station']

    #===============================================================================
    # Plot files 
    #===============================================================================

    #change save='n' to save='y' if want to save the plots, will save in a folder
    #called dirpath\plots

    #if you don't want to use the save icon in the plots you can type in the 
    #interpreter plt.savefig(FullPathSaveName,fmt='pdf)
    #note that fmt can be jpg, eps, or svg

    #brp.plotBFfiles(flst['edifile'],cohfile=flst['cohfile'],save='n',show='y')

    #if this doesn't work try:
    #mtplot.plotResPhase(flst['edifile'],plotnum=2,fignum=1)
    #or
    #z1=Z.Z(flst['edifile'])
    #z1.plotResPhase(fignum=1,plottype=2)
    #z1.plotPTAll(fignum=2

if __name__ == '__main__':
    main()
