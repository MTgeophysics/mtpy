#!/usr/bin/env python

#-----------------------------------------------------------------------
#	decimates 650Hz data to 1Hz.
#	need to change outdirname, indir_base
#	in main() need to change profile_prefix, station_idx,stationname
#-----------------------------------------------------------------------


import re
import os
import sys
import os.path as op
import shutil
import calendar


Pak2Asc_exe = '/stash/Working_Scripts/PakScriptsMay14/Pak2AscDeci650'
#fileroot = 'Pak2Asc-500-Data.as4'
temp_out_fn = 'Pak2Asc-500-Data.fw4.deci650'
outdirname = '/stash/roma_decimated/Magnetics_Processed_Phase_3/'

indir_base = '/media/Elements2/KEXRR_Apr/'

profile_prefix = 'L1'

starting_station = 0

number_of_stations = 10


def main():

    # if len(sys.argv) < 3:
    #     sys.exit('\n\tERROR - need 2 arguments as input: '\
    #         '<stationname> <station data location> \n')

    # stationname = sys.argv[1].upper()
    # stationdatafolder = sys.argv[2]

    #indir =   stationdatafolder
    print()
    profile = 1

    #profile_prefix = '{0:02d}'.format(profile)
    #profile_prefix = 'RRB'

    for i in range(number_of_stations):

        station_idx = i + starting_station
        stationname = '{0:02d}'.format(station_idx)
        #stationname = 'RRB'
        indir = op.abspath(op.join(indir_base, stationname))

        if not op.isdir(indir):
            print('WARNING - no folder found for station {0} ({1})\n'.format(
                stationname, indir))
            continue

        #outdir = op.abspath(op.join(outdirname,'{0}'.format(profile_prefix+'_'+stationname)))
        outdir = op.abspath(op.join(outdirname, '{0}'.format(stationname)))

        if not op.isdir(outdir):
            os.makedirs(outdir)

        try:
            unpack1station(stationname, indir, outdir)
        except:
            continue

        print()


def unpack1station(stationname, indir, outdir):

    cwd = op.abspath(os.curdir)

    # print
    all_subdirs = os.listdir(indir)
    all_subdirs = sorted([op.abspath(op.join(indir, i))
                          for i in all_subdirs if op.isdir(op.join(indir, i))])

    # print all_subdirs
    print('\t=====================\n\tUnpacking station {0}:\n\t====================='.format(stationname))
    counter6hrblocks = 0
    for subdir in all_subdirs:

        try:
            os.chdir(subdir)
            print('\n...switched to folder {0}...'.format(op.abspath(os.curdir)))

            lo_unpackables = os.listdir('.')
            lo_unpackables = [
                i for i in lo_unpackables if i.lower().endswith('.bz2')]
            if len(lo_unpackables) > 0:
                print('...unpacking...')
                for bz2 in lo_unpackables:
                    os.system('bunzip2 -f -k {0}'.format(bz2))

            lo_paks = os.listdir('.')
            lo_paks = [i for i in lo_paks if i.lower().endswith('.pak')]

            if len(lo_paks) == 0:
                print('no data - skipping subfolder {0}'.format(subdir))
                continue

            print('...conversion Pak2Asc...')
            try:
                os.system('{0} fw4 00000000.pak'.format(Pak2Asc_exe))
            except:
                continue

            counter6hrblocks += 1
            new_filename = 'sta{0}_decimated_{1:03d}'.format(
                stationname, counter6hrblocks)
            shutil.copy2(temp_out_fn, new_filename)
            print('...remove unpacked data...')
            os.system('rm -f *.pak')

            try:
                print('...moving data to {0}...'.format(outdir))
                shutil.copy2(new_filename, outdir)
            except:
                print("couldn't move data")

            print('...remove all the rest of temporary files...\n')
            os.system('rm -f Pak2Asc-*')
            os.system('rm -f {0}'.format(new_filename))
        except:
            continue

    os.chdir(cwd)

    print('\n\t {0} blocks done for station {1}!\n'.format(counter6hrblocks, stationname))
    print('Finished station {0}!\n====================='.format(stationname))

    return 0


if __name__ == '__main__':
    main()
