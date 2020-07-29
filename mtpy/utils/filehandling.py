#!/usr/bin/env python

"""
Helper functions for file handling. 

The various functions deal with renaming, sorting, 
concatenation of time series, extraction of names and times from filenames,
reading configuration files, ....


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import sys
import os
import os.path as op
import calendar
import time
import fnmatch
import shutil

import mtpy.utils.calculator as MTcc
import mtpy.utils.exceptions as MTex
import mtpy.utils.configfile as MTcf

#=================================================================

#define uncertainty for differences between time steps
epsilon = 1e-9

#=================================================================

lo_headerelements = ['station', 'channel','samplingrate','t_min',
                    'nsamples','unit','lat','lon','elev']

#=================================================================



def read_surface_ascii(ascii_fn):
    """
    read in surface which is ascii format ()
    unlike original function, returns numpy array of lon, lat, elev (no projections)

    The ascii format is assumed to be:
    ncols        2743
    nrows        2019
    xllcorner    111.791666666667 (lon of lower left)
    yllcorner    -45.341666666667 (lat of lower left)
    cellsize     0.016666666667
    NODATA_value  -9999
    elevation data origin (0,0) is NW upper left.
    NW --------------> E
    |
    |
    |
    |
    S
    """
    with open(ascii_fn, 'r') as dfid:
        d_dict = {}
        skiprows = 0
        for ii in range(6):
            dline = dfid.readline()
            dline = dline.strip().split()
            key = dline[0].strip().lower()
            value = float(dline[1].strip())
            d_dict[key] = value
            # check if key is an integer
            try:
                int(key)
            except:
                skiprows += 1
    # not required dfid.close()

    x0 = d_dict['xllcorner']
    y0 = d_dict['yllcorner']
    nx = int(d_dict['ncols'])
    ny = int(d_dict['nrows'])
    cs = d_dict['cellsize']

    elevation = np.loadtxt(ascii_fn, skiprows=skiprows)[::-1]  # ::-1 reverse an axis to put the southern line first

    # create lat and lon arrays from the dem file
    lon = np.arange(x0, x0 + cs * (nx), cs)
    lat = np.arange(y0, y0 + cs * (ny), cs)
    lon = np.linspace(x0, x0 + cs * (nx - 1), nx)
    lat = np.linspace(y0, y0 + cs * (ny - 1), ny)

    return lon, lat, elevation   # this appears correct


def read1columntext(textfile):
    """
    read a list from a one column text file
    """
    return [ff.strip() for ff in open(textfile).readlines()]

def read_stationdatafile(textfile,read_duplicates = True):
    """
    read a space delimited file containing station info of any sort - 
    3 columns: station x, y, ... - to a dictionary - station:[x,y,...]
    textfile = full path to text file
    read_duplicates = True/False - if stations are listed more than once do you 
                                   want to read all information or just the 
                                   first occurrence, default True
    
    example:
    import mtpy.utils.filehandling as fh
    stationdict = fh.read_stationxyfile(textfile)
    
    
    """
    stationdict = {}
    for line in open(textfile).readlines():
        line = line.split()
        for l in range(1,len(line)):
            try:
                line[l] = float(line[l])
            except:
                pass
        sname = line[0]
        if sname not in list(stationdict.keys()):
            stationdict[sname] = line[1:]
        else:
            if read_duplicates:
                if len(line) <= 2:
                    value = line[1]
                else:
                    value = line[1:]
                stationdict[sname].append(value)

            
    return stationdict

def make_unique_filename(infn):

    fn = op.abspath(infn)
    outfn = fn
    i = 1
    while op.isfile(outfn):
        filebase = op.splitext(fn)[0]
        outfn = filebase +'_%i'%i+ op.splitext(fn)[1]
        i += 1

    return outfn

def make_unique_folder(wd,basename = 'run'):
    """
    make a folder that doesn't exist already.
    """        
  
    # define savepath. need to choose a name that doesn't already exist
    i = 1
    svpath_str = basename
    svpath = svpath_str+'_%02i'%i
    while os.path.exists(op.join(wd,svpath)):
        i += 1
        svpath = svpath_str+'_%02i'%i
        
    savepath = op.join(wd,svpath)
    os.mkdir(savepath)
        
    return savepath


def validate_save_file(savepath=None,savefile=None,basename=None,prioritise_savefile = False):
    """
    Return savepath, savefile and basename, ensuring they are internally
    consistent and populating missing fields from the others or using defaults.
    
    Prioritises savepath and basename. I.e. if savepath, savefile and basename
    are all valid but inconsistent, savefile will be updated to reflect
    savepath and basename
    
    :param savepath: directory to save to
    :param savefile: full file path to save to
    :param basename: base file name to save to
    
    """
    
    if prioritise_savefile:
        if savefile is not None:
            if os.path.exists(savefile):
                if os.path.isdir(savefile):
                    savepath = None
                else:
                    savepath, basename = None, None
                    
    
    # first, check if savepath is a valid directory
    if savepath is not None:
        if not os.path.isdir(savepath):
            if not os.path.exists(savepath):
                savepath = None
            else:
                savepath = os.path.dirname(savepath)

    # second, if basename is None, get it from savefile or set a default
    if basename is None:
        if savefile is not None:
            basename = os.path.basename(savefile)
        else:
            basename = 'default.dat'

    # third, if savepath is None, get it from savefile or set a default
    if savepath is None:
        if savefile is not None:
            savepath = os.path.dirname(savefile)
            if not os.path.isdir(savepath):
                savepath = os.getcwd()
        else:
            savepath = os.getcwd()
    
    # finally, make savefile so it is consistent with savepath and basename
    savefile = os.path.join(savepath,basename)
            
    return savepath, savefile, basename
    
            
def sort_folder_list(wkdir,order_file,indices=[0,9999],delimiter = ''):
    """
    sort subfolders in wkdir according to order in order_file
    
    wkdir = working directory containing subfolders
    order = full path to text file containing order.
            needs to contain a string to search on that is the same length
            for each item in the list
    indices = indices to search on; default take the whole string
    
    returns a list of directories, in order.
    
    """
    order = read1columntext(order_file)

    plst = []
    flst = [i for i in os.listdir(wkdir) if os.path.exists(os.path.join(wkdir,i))]
#    print flst
    for o in order:
        for f in flst:
            if str.lower(f.strip().split(delimiter)[0][indices[0]:indices[1]]) == str.lower(o)[indices[0]:indices[1]]:
                plst.append(os.path.join(wkdir,f))
    return plst


def get_pathlist(masterdir, search_stringlist = None, search_stringfile = None,
                 start_dict = {},split='_',extension='',folder=False):
    """
    get a list of files or folders by searching on a string contained in 
    search_stringlist or alternatively search_stringfile
    
    returns:
    dictionary containing search strings as keys and file/folder as values
    
    masterdir - directory to search in
    search_stringlist = list containing string identifiers for files or folders, 
                        e.g. k0101 will work for edifile k0101.edi or 
                        folder k0101.
    search_stringfile = alternative to search_stringlist (need to provide one)
                        will get search_stringlist from a file, full path
                        or make sure you are in the correct directory!
    start_dict = starting dictionary to append to, default is an empty dict
    split = if no exact match is found, search string will be split using split
            character, useful when matching up edi's to inversion directories
            that both contain additional characters
    extension = file extension, e.g. '.edi'
    
    
    """
    
    if search_stringfile is not None:
        if (search_stringlist is None) or (len(search_stringlist)) == 0:
            search_stringlist = read1columntext(search_stringfile)

    flist = [i for i in os.listdir(masterdir) if i[len(i)-len(extension):] == \
             extension]
    
    if folder:
        flist = [op.join(masterdir,i) for i in flist if op.isdir(op.join(masterdir,i))]

    for s in search_stringlist:
        s = str.lower(s)#.split(split)[indices[0]:indices[1]]
        for d in flist:
            d = str.lower(d)#.split(split)[indices[0]:indices[1]]
            append = False
            if s in os.path.basename(d):
                append = True
            else:
                slst = s.strip().split(split)
                for ss in slst:
                    if ss in d:
                        append = True
            if append:
                start_dict[s] = op.join(masterdir,d)

    return start_dict
                

def get_sampling_interval_fromdatafile(filename, length = 3600):
    """ 
    Find sampling interval from data file.

    Provide data file (purely numerical content) and total 
    data length in seconds (default 3600). Data are read in by 
    the Numpy 'loadtxt'-function, the lentgh of the data array yields
    the sampling interval.  

    Lines beginning with # are ignored.

    """

    fn = op.abspath(op.realpath(filename))
    dd = np.loadtxt(fn)
    sampling_interval = length/float(len(dd))

    return sampling_interval


def EDL_make_Nhour_files(n_hours,inputdir, sampling , stationname = None, outputdir = None):

    """
    See 'EDL_make_dayfiles' for description and syntax.

    Only difference: output files are blocks of (max) N hours, starting to count 
    at midnight (00:00h) each day.

    Conditions:
    
    1.   24%%N = 0
    2.   input data files start on the hour marks

    Not working yet!!

    """
    #print '\n\tnot working yet - code under development !!\n'
    #return

    try:
        if 24%n_hours != 0:
            raise Exception("problem!!!")
    except:
        sys.exit('ERROR - File block length must be on of: 1,2,3,4,6,8,12 \n')
    
    n_hours = int(n_hours)
    #list of starting hours for the data blocks:
    lo_hours = [int(i) for i in np.arange(int(24/n_hours))*n_hours ]
    no_blocks = len(lo_hours)
    
    #build a list that contains the respective groups of starting hours belonging
    # to the same block
    lo_blockhours = []
    counter = 0 
    dummylist = []
    for i in range(24):
        dummylist.append(i)
        counter += 1
        if counter == n_hours:
            lo_blockhours.append(dummylist)
            dummylist = []
            counter = 0 

    # most of the following code is redundant/taken from the 
    # 'EDL_make_dayfiles' function
    # This can be cleaned up later

    try:
        if type(inputdir)==str:
            raise Exception("problem!!!")
        lo_foldernames = [i for i in inputdir]
    except TypeError:
        lo_foldernames = [inputdir]

    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz']

    lo_allfiles = []
    pattern = '*.[ebEB][xyzXYZ]'
    if stationname is not None:
        pattern = '*{0}*.[ebEB][xyzXYZ]'.format(stationname.lower())
    print('\nSearching for files with pattern: ',pattern)

    for folder in lo_foldernames:
        wd = op.abspath(op.realpath(folder)) 
        if not op.isdir(wd):
            #print 'Directory not existing: %s' % (wd)
            lo_foldernames.remove(wd)
            continue    

        lo_dirfiles = [op.abspath(op.join(wd,i))  for i in os.listdir(wd) 
                        if fnmatch.fnmatch(i.lower(),pattern.lower()) is True]
        lo_allfiles.extend(lo_dirfiles)   

    #check, if list of files is empty
    if len(lo_allfiles) == 0:
        if stationname is not None:
            raise MTex.MTpyError_inputarguments('Directory(ies) do(es) not contain'\
            ' files to combine for station {0}:\n {1}'.format(stationname, inputdir))

        raise MTex.MTpyError_inputarguments('Directory does not contain files'\
                                            ' to combine:\n {0}'.format(inputdir))

    #define subfolder for storing dayfiles
    outpath = op.join(os.curdir,'{0}hourfiles'.format(int(n_hours)))    
    if outputdir is not None:
        try:
            outpath = op.abspath(op.join(os.curdir,outputdir))
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise Exception("problem!!!")
            if not os.access(outpath, os.W_OK):
                raise Exception("problem!!!")
        except:
            print('Cannot generate writable output directory {0} - using'\
                    ' generic location "dayfiles" instead'.format(outpath))
            outpath = op.join(wd,'{0}hourfiles'.format(int(n_hours)))    
            pass

    #generate subfolder, if not existing
    if not op.exists(outpath):
        try:
            os.makedirs(outpath)
        except:
            MTex.MTpyError_inputarguments('Cannot generate output'\
                                ' directory {0} '.format(outpath))

    #outer loop over all components
    for comp in components:

        #make list of files for the current component
        lo_files = np.array([op.join(wd,i) for i in lo_allfiles 
                            if (i.lower()[-2:] == comp)])

        #make list of starting times for the respective files
        lo_starttimes = np.array([EDL_get_starttime_fromfilename(f) 
                                    for f in lo_files])
        
        #sort the files by their starting times
        idx_chronologic = np.argsort(lo_starttimes)
        
        #obtain sorted lists of files and starting times
        lo_sorted_files = list(lo_files[idx_chronologic])
        lo_sorted_starttimes = list(lo_starttimes[idx_chronologic])
        idx = 0
        while idx < len(lo_sorted_starttimes):
            try:
                val = lo_sorted_starttimes[idx]
            except:
                break
            if val is None:
                dummy = lo_sorted_files.pop(idx)
                dummy = lo_sorted_starttimes.pop(idx)
            else:
                idx += 1


        #set stationname, either from arguments or from filename
        if stationname is None:
            stationname = EDL_get_stationname_fromfilename(lo_sorted_files[0]).upper()

        #set counting variables - needed for handling of consecutive files

        #flags, checking, if data has to written to file
        sameday = True
        sameblock = True
        fileopen = False
        complete = False
        #numerical index of files for same combination of date and hour
        fileindex = 0
        #current file's daily block number 
        blockindex = 0 

        #allocate a data array to fill
        # this is more memory efficient than extending lists!!
        #cater for potential rounding errors:
        if sampling < 1:
            max_n_data = 3600*int(n_hours) * (int(1./sampling)+1)
        else:
            max_n_data = int(3600.*int(n_hours)/sampling) + 1

        block_data = np.zeros(max_n_data,'int')


        #loop over all (sorted) files for the current component
        for idx_f,f in enumerate(lo_sorted_files):

            try:

                print('Reading file %s' %(f))
                #starting time of current file
                file_start_time = lo_sorted_starttimes[idx_f]

                #get tuple with the starting time of the current file
                file_start = time.gmtime(file_start_time)
            
                #read in raw data
                data_in = []
                Fin = open(f)
                for line in Fin:#.readlines():
                #    try:
                    data_in.append(int(float(line.strip())))
                 #   except:
                  #      pass
                data_in = np.array(data_in)
                Fin.close()
                #data_in = np.loadtxt(f)
            except:
                print('WARNING - could not read file - skipping...')
                continue
            no_samples = len(data_in)

            #time axis of the file read in :
            tmp_file_time_axis = np.arange(no_samples)*sampling+file_start_time

            #end: time of the last sample + 1x sampling-interval
            file_end_time =  tmp_file_time_axis[-1] + sampling
         
            #set the current file's starting time as starting time for output file, 
            #if no output file is open already

            
            if fileopen is False:
                #starttime of output file
                outfile_starttime =  file_start_time
                #output data sample index:                
                arrayindex = 0

                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    block_data[:len(data_in)] = data_in                    
                    
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    block_data[:len(data_in)] = data_in[:,1]
                
                #jump with index to current point on time axis 
                arrayindex += len(data_in)

                #current (virtual) end time of open file
                outfile_endtime = file_end_time

                #find the date code of the outputfile
                file_date = '{0}{1:02}{2:02}'.format(file_start[0],
                                                 file_start[1], file_start[2])
                
                #...aaaand the hour code as well
                data_hour = file_start[3]

                #determine, which of the daily data blocks we are currently 
                #processing/writing
                blockindex = no_blocks-1
                
                file_hour = lo_hours[-1]
                for t in range(len(lo_hours)-1):
                    if lo_hours[t+1]>data_hour:
                        blockindex = t
                        file_hour = lo_hours[blockindex]
                        break


                #define output filename
                new_fn = '{0}_{5}hours_{1}_{2:02d}_{3}.{4}'.format(stationname,
                                                 file_date, file_hour,fileindex, comp,n_hours)
                #absolute filename:
                new_file = op.abspath(op.join(outpath,new_fn))
                
                #open output file 
                F = open(new_file,'w')
                
                #set flag for further loop steps
                fileopen = True

            
            else:
                #check, if the new file ends earlier than data in buffer.
                #if yes, just skip this file:
                if file_end_time < outfile_endtime:
                    continue 

                #if current file starts earlier than the endtime of data in buffer, but extends the time span
                #then delete ambiguous  parts of the buffer:
                elif (outfile_endtime - file_start_time) > epsilon:

                    #find point on the outfile time axis for the beginning of current file:
                    overlap_idx = arrayindex - int((outfile_endtime - file_start_time)/sampling)

                    #set the array index back to the appropriate value corresponding to the 
                    #start of the new file
                    arrayindex = overlap_idx
                   
  
                #append current data                  
                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    block_data[arrayindex:arrayindex+len(data_in)] = data_in                    
                    #outfile_data.extend(data_in.tolist())
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    block_data[arrayindex:arrayindex+len(data_in)] = data_in[:,1]                    
                    #outfile_data.extend(data_in[:,1].tolist())

                #update position in time
                arrayindex += len(data_in)

                #update (virtual) end of outfile data
                outfile_endtime = (arrayindex+1)*sampling + outfile_starttime


            #-----------
            # current file has been read in, data in buffer have been updated
            # now check, if it has to be written to file....

            #check, if there is a next file at all:
            try:
                next_file_start_time = lo_sorted_starttimes[idx_f + 1]
            except:
                complete = True

            #if there is a next file, 
            # - check, if it's the same day
            # - check, if it' the same block
            # - check, if it continues at the end of the current one:

            if complete is False:
                next_file_start_time = lo_sorted_starttimes[idx_f + 1]
                next_file_start = time.gmtime(next_file_start_time)
                
                nextfile_hour = next_file_start[3]

                nextfile_blockindex  = no_blocks-1
                for t in range(len(lo_hours)-1):
                    if lo_hours[t+1]>nextfile_hour:
                        nextfile_blockindex = t
                        break


                nextfile_day = '{0}{1:02}{2:02}'.format(next_file_start[0],
                                                 next_file_start[1], next_file_start[2])
                if  nextfile_day != file_date: 
                    complete = True
                    sameday = False                    
                if  nextfile_blockindex != blockindex:
                    complete = True
                    sameblock = False

            if complete is False:
                if next_file_start_time - file_end_time > epsilon: 
                    complete = True
                    sameblock = True
                    sameday = True

         
            #ipdb.set_trace()

            #check, if the file has to be closed and written now
            if complete is True :

                #define header info
                if outfile_starttime%1==0:
                    outfile_starttime = int(outfile_starttime)

                    headerline = '# {0} {1} {2:.1f} {3} {4} \n'.format(
                                    stationname, comp.lower(), 1./sampling, 
                                    outfile_starttime, arrayindex)
                else:
                    headerline = '# {0} {1} {2:.1f} {3:f} {4} \n'.format(
                                    stationname, comp.lower(), 1./sampling, 
                                    outfile_starttime, arrayindex)

                F.write(headerline)

                #outfile_array = np.zeros((len(outfile_timeaxis),2))
                #outfile_array[:,0] = outfile_timeaxis
                #outfile_array[:,1] = outfile_data
                for i in range(arrayindex):
                    F.write('{0}\n'.format(int(block_data[i])))
                #outstring = '\n'.join(['{0:d}'.format(i) for i in day_data[:arrayindex]])
                #F.write(outstring)
                #np.savetxt(F,day_data[:arrayindex],fmt='%d')
                #np.savetxt(F, np.array(outfile_data))
                arrayindex = 0
                
                F.close()
                print('\t wrote file %s'%(new_file))

                fileopen = False
                complete = False
                #blockindex = (blockindex+1)%no_blocks
    
                if sameday is True and sameblock is True : 
                    fileindex +=1        
                else:
                    fileindex = 0

 


def EDL_make_dayfiles(inputdir, sampling , stationname = None, outputdir = None):
    """

    Concatenate ascii time series to dayfiles (calendar day, UTC reference).

    Data can be within a single directory or a list of directories. 
    However, the files in the directory(ies) 'inputdir' have to be for 
    one station only, and named with a 2 character suffix, defining the channel! 

    If the time series are interrupted/discontinuous at some point, a new file 
    will be started after that point, where the file index 'idx' is increased by 1.
    If no stationname is given, the leading non-datetime characters in the first 
    filename are used.


    Files are named as 'stationname_samplingrate_date_idx.channel'
    Stationname, channel, and sampling are written to a header line.

    Output data consists of a single column float data array. The data are 
    stored into one directory. If 'outputdir' is not specified, a subdirectory 
    'dayfiles' will be created witihn the current working directory. 

    Note: 
    Midnight cannot be in the middle of a file, because only file starts are 
    checked for a new day!!

    """
    try:
        if type(inputdir)==str:
            raise Exception("problem!!!")
        lo_foldernames = [i for i in inputdir]
    except TypeError:
        lo_foldernames = [inputdir]

    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz']

    lo_allfiles = []
    pattern = '*.[ebEB][xyzXYZ]'
    if stationname is not None:
        pattern = '*{0}*.[ebEB][xyzXYZ]'.format(stationname.lower())
    print('\nSearching for files with pattern: ',pattern)

    for folder in lo_foldernames:
        wd = op.abspath(op.realpath(folder)) 
        if not op.isdir(wd):
            #print 'Directory not existing: %s' % (wd)
            lo_foldernames.remove(wd)
            continue    

        lo_dirfiles = [op.abspath(op.join(wd,i))  for i in os.listdir(wd) 
                        if fnmatch.fnmatch(i.lower(),pattern.lower()) is True]
        lo_allfiles.extend(lo_dirfiles)   

    #check, if list of files is empty
    if len(lo_allfiles) == 0:
        if stationname is not None:
            raise MTex.MTpyError_inputarguments('Directory(ies) do(es) not contain'\
            ' files to combine for station {0}:\n {1}'.format(stationname, inputdir))

        raise MTex.MTpyError_inputarguments('Directory does not contain files'\
                                            ' to combine:\n {0}'.format(inputdir))

    #define subfolder for storing dayfiles
    outpath = op.join(os.curdir,'dayfiles')    
    if outputdir is not None:
        try:
            outpath = op.abspath(op.join(os.curdir,outputdir))
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise Exception("problem!!!")
            if not os.access(outpath, os.W_OK):
                raise Exception("problem!!!")
        except:
            print('Cannot generate writable output directory {0} - using'\
                    ' generic location "dayfiles" instead'.format(outpath))
            outpath = op.join(wd,'dayfiles')    
            pass

    #generate subfolder, if not existing
    if not op.exists(outpath):
        try:
            os.makedirs(outpath)
        except:
            MTex.MTpyError_inputarguments('Cannot generate output'\
                                ' directory {0} '.format(outpath))

    #outer loop over all components
    for comp in components:

        #make list of files for the current component
        lo_files = np.array([op.join(wd,i) for i in lo_allfiles 
                            if (i.lower()[-2:] == comp)])

        #make list of starting times for the respective files
        lo_starttimes = np.array([EDL_get_starttime_fromfilename(f) 
                                    for f in lo_files])
        
        #sort the files by their starting times
        idx_chronologic = np.argsort(lo_starttimes)
        
        #obtain sorted lists of files and starting times
        lo_sorted_files = list(lo_files[idx_chronologic])
        lo_sorted_starttimes = list(lo_starttimes[idx_chronologic])
        idx = 0
        while idx < len(lo_sorted_starttimes):
            try:
                val = lo_sorted_starttimes[idx]
            except:
                break
            if val is None:
                dummy = lo_sorted_files.pop(idx)
                dummy = lo_sorted_starttimes.pop(idx)
            else:
                idx += 1


        #set stationname, either from arguments or from filename
        if stationname is None:
            stationname = EDL_get_stationname_fromfilename(lo_sorted_files[0]).upper()

        #set counting variables - needed for handling of consecutive files

        sameday = 0
        fileopen = 0
        incomplete = 0
        fileindex = 0

        #allocate a data array to fill
        # this is more memory efficient than extending lists!!
        #cater for potential rounding errors:
        if sampling < 1:
            max_n_data = 86400 * (int(1./sampling)+1)
        else:
            max_n_data = int(86400./sampling) + 1

        day_data = np.zeros(max_n_data,'int')



        #loop over all (sorted) files for the current component
        for idx_f,f in enumerate(lo_sorted_files):

            try:

                print('Reading file %s' %(f))
                #starting time of current file
                file_start_time = lo_sorted_starttimes[idx_f]


                #get tuple with the starting time of the current file
                file_start = time.gmtime(file_start_time)
            
                #read in raw data
                data_in = []
                Fin = open(f)
                for line in Fin:#.readlines():
                #    try:
                    data_in.append(int(float(line.strip())))
                 #   except:
                  #      pass
                data_in = np.array(data_in)
                Fin.close()
                #data_in = np.loadtxt(f)
            except:
                print('WARNING - could not read file - skipping...')
                continue
            no_samples = len(data_in)

            tmp_file_time_axis = np.arange(no_samples)*sampling+file_start_time
            #file_time_axis = (np.arange(no_samples)*sampling +
            #                 file_start_time).tolist()


            #time of the last sample + 1x sampling-interval
            #file_end_time =  file_time_axis[-1] + sampling
            file_end_time =  tmp_file_time_axis[-1] + sampling
         



            #set the time as starting time for output file, if no output file is open already
            if fileopen == 0:
                outfile_starttime =  file_start_time

                #outfile_timeaxis = file_time_axis
                old_time_axis = tmp_file_time_axis[:]
                

                arrayindex = 0

                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in                    
                    #outfile_data = data_in.tolist()
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in[:,1]
                    #outfile_data = data_in[:,1].tolist()
                
                #jump with index to current point on time axis 
                arrayindex += len(data_in)
                outfile_endtime = file_end_time


                file_date = '{0}{1:02}{2:02}'.format(file_start[0],
                                                 file_start[1], file_start[2]) 


                #define output filename
                new_fn = '{0}_1day_{1}_{2}.{3}'.format(stationname,
                                                 file_date, fileindex, comp)
                
                new_file = op.abspath(op.join(outpath,new_fn))
                
                #open output file 
                F = open(new_file,'w')
                
                fileopen = 1


            
            else:
                #check, if the new file ends earlier than data in buffer.
                #if yes, just skip this file:
                if file_end_time < outfile_endtime:
                    continue 

                #if current file starts earlier than the endtime of data in buffer then delete ambiguous  parts of the buffer:
                #elif (outfile_timeaxis[-1] - file_start_time) > epsilon:
                elif (outfile_endtime - file_start_time) > epsilon:

                    #find point on the outfile time axis for the beginning of current file:
                    overlap_idx = arrayindex - int((outfile_endtime - file_start_time)/sampling)

                    #set the array index back
                    arrayindex = overlap_idx
                   
                    #re-define outfile time axis and data
                    # outfile_timeaxis = np.delete(outfile_timeaxis,
                    #                              np.arange(len(outfile_timeaxis) - 
                    #                              overlap_idx) + 
                    #                              overlap_idx).tolist()


                    # outfile_data = np.delete(outfile_data, 
                    #                             np.arange(len(outfile_data) - 
                    #                             overlap_idx) + 
                    #                             overlap_idx).tolist()
                

                #old_time_axis = tmp_file_time_axis[:]
                #append current file's time axis
                #outfile_timeaxis.extend(file_time_axis)
                    
                #append current data                  
                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in                    
                    #outfile_data.extend(data_in.tolist())
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in[:,1]                    
                    #outfile_data.extend(data_in[:,1].tolist())

                arrayindex += len(data_in)
                print(len(data_in),arrayindex)

                arrayindex += len(data_in)
                outfile_endtime = (arrayindex+1)*sampling + outfile_starttime



            #-----------

            #check, if there is a next file:
            try:
                next_file_start_time = lo_sorted_starttimes[idx_f + 1]
            except:
                incomplete = 1

            #if there is a next file, 
            # - check, if it's the same day
            # - check, if it continues at the end of the current one:
            if incomplete == 0:
                next_file_start_time = lo_sorted_starttimes[idx_f + 1]
                next_file_start = time.gmtime(next_file_start_time)
                
                if next_file_start[2] == file_start[2] :
                    #print 'sameday',file_start[:]
                    sameday = 1
                else:
                    incomplete = 1
                    sameday = 0
                    fileindex = 0
                    #print '\t NOT sameday', fileindex


                if next_file_start_time - file_end_time > epsilon: 
                    incomplete = 1

            if incomplete == 1 and sameday == 1 : 
                fileindex +=1        

           

            #check, if the file has to be closed and written now
            if incomplete == 1 :

                #define header info
                if outfile_starttime%1==0:
                    outfile_starttime = int(outfile_starttime)

                    headerline = '# {0} {1} {2:.1f} {3} {4} \n'.format(
                                    stationname, comp.lower(), 1./sampling, 
                                    outfile_starttime, arrayindex)
                else:
                    headerline = '# {0} {1} {2:.1f} {3:f} {4} \n'.format(
                                    stationname, comp.lower(), 1./sampling, 
                                    outfile_starttime, arrayindex)


                F.write(headerline)

                #outfile_array = np.zeros((len(outfile_timeaxis),2))
                #outfile_array[:,0] = outfile_timeaxis
                #outfile_array[:,1] = outfile_data
                for i in range(arrayindex):
                    F.write('{0}\n'.format(int(day_data[i])))
                #outstring = '\n'.join(['{0:d}'.format(i) for i in day_data[:arrayindex]])
                #F.write(outstring)
                #np.savetxt(F,day_data[:arrayindex],fmt='%d')
                #np.savetxt(F, np.array(outfile_data))
                arrayindex = 0
                

                F.close()
                print('\t wrote file %s'%(new_file))

                fileopen = 0
                incomplete = 0
    




def EDL_get_starttime_fromfilename(filename): 
    """ 
    Return starttime of data file in epoch seconds.

    Starting time is determined by the filename. This has to be of the form
    'somthing/*.stationname.ddmmyyHHMMSS.??'


    """     
    #clip parent paths and structure
    bn = op.basename(filename)
    parts_of_bn = bn.split('.')
    timestamp = parts_of_bn[-2]
    
    try:
        secs = int(float(timestamp[-2:]))
        mins = int(float(timestamp[-4:-2]))
        hours =int(float(timestamp[-6:-4]))
        day =  int(float(timestamp[-8:-6]))
        month =int( float(timestamp[-10:-8]))
        year = int(float(timestamp[-12:-10]))
        if year < 50:
            year += 2000
        else:
            year += 1900

        timetuple = (year, month, day, hours, mins, secs)

        epochtime = calendar.timegm(timetuple)

    except:
        epochtime = None

    return epochtime


def EDL_get_stationname_fromfilename(filename):

    bn = op.basename(filename)
    parts_of_bn = bn.split('.')
    stationtime = parts_of_bn[-2]

    stationname = stationtime[:-12].upper()

    if len(stationname) == 0:
        stationname = 'DUMMYSTATION'


    return stationname


def read_data_header(fn_raw):
    """
    Deprecated!!!
    USE 
              read_ts_header

    INSTEAD


        Read the header line of MTpy TS data files.

    input
    -----
    MTpy TS data file name

    output
    -------
    list of header elements:
    stationname, channel, sampling rate, starttime first sample, 
    starttime last sample, unit, lat, lon, elevation

    """


    fn = op.abspath(op.realpath(fn_raw))

    if not op.isfile(fn):
        raise MTex.MTpyError_inputarguments('Not a file:%s'%fn)
    try:
        F = open(fn, 'r')
    except:
        raise MTex.MTpyError_inputarguments('File not readable:%s'%fn)

    firstline = F.readline().strip().split()
    if not firstline[0][0] == '#':
        raise MTex.MTpyError_ts_data('Time series data file does '
            'not have a proper header:%s'%fn)

    F.close()

    header_list = []

    idx_header = 0

    if len(firstline[0]) > 1:
        header_list.append(firstline[0][1:].upper())
    else:
        header_list.append(firstline[1].upper())
        idx_header += 1

    header_list.append( firstline[idx_header+1].lower() )
    header_list.append( float(firstline[idx_header+2]) )
    header_list.append( float(firstline[idx_header+3]) )
    header_list.append( int(float(firstline[idx_header+4])) )
    header_list.append( firstline[idx_header+5].lower() )
    header_list.append( float(firstline[idx_header+6]) )
    header_list.append( float(firstline[idx_header+7]) )
    header_list.append( float(firstline[idx_header+8]) )


    return header_list


def read_2c2_file(filename):
    """
    Read in BIRRP 2c2 coherence files and return 4 lists 
    containing [period],[freq],[coh],[zcoh]. Note if any of the coherences are 
    negative a value of 0 will be given to them.

    """

    period = []
    freq = []
    coh1 = []
    zcoh1 = []

    F_in = open(filename,'r')
    data_raw = F_in.readlines()
    
    for ii in range(len(data_raw)):

        coh_row = data_raw[ii].strip().split()
        
        try:
            period.append(float(coh_row[0]))
        except:
            period.append(0.)
        try:
            freq.append(  float(coh_row[1]))
        except:
            freq.append(0.)
        try:
            coh1.append(  float(coh_row[2]))
        except:
            coh1.append(0.)
        try:
            zcoh1.append( float(coh_row[3]))
        except:
            zcoh1.append(0.)

    indexorder = np.array(period).argsort()

    period = np.array(period)[indexorder]
    freq = np.array(freq)[indexorder]
    coh1 = np.array(coh1)[indexorder]
    zcoh1 = np.array(zcoh1)[indexorder]

    return period, freq, coh1, zcoh1

def validate_ts_file(tsfile):
    """
        Validate MTpy timeseries (TS) data file
        Return Boolean value True/False .

    """ 
    tsfile = op.abspath(tsfile)

    try:
        header = read_ts_header(tsfile)

        if header['station'] is None:
            #print 'header'
            raise Exception("header has no station")
        if header['channel'] is None:
            #print 'channel'
            raise Exception("header has no channel")
        
        sr = float(header['samplingrate'])
        t0 = float(header['t_min'])
        ns = int(float(header['nsamples']))
        
        data = np.loadtxt(tsfile)
        
        if len(data) != ns:
            #print 'data length'
            raise Exception("data length wrong")
        if data.dtype not in [int, float]:
            #print 'data type'
            raise Exception("data type wrong")

    except:
        #print 'number'
        return False


    return True



def read_ts_header(tsfile):
    """ Read in the header line from MTpy timeseries data files.
        
        Return header as dictionary. Return empty dict, 
        if no header line was found.
    """

    header_dict = {}

    tsfile = op.abspath(tsfile)
    
    if not op.isfile(tsfile):
        raise MTex.MTpyError_inputarguments('Error - '
            'input file not existing: {0}'.format(tsfile))

    try:
        with open(tsfile,'r') as F:
            firstline =''
            #ignoring empty lines or lines with just the '#' character in it
            while len(firstline) == 0:
                firstline = F.readline().strip()
                if firstline == '#':
                    firstline = ''
        if firstline[0] != '#':
            raise Exception("First line does not begin with #")
    except:
        raise MTex.MTpyError_ts_data('No header line found -'
            ' check file: {0}'.format(tsfile))
        

    firstline = firstline.replace('#','')
    headerlist = firstline.split()


    for i in range(len(headerlist)):
        header_dict[lo_headerelements[i]] = headerlist[i]
        #old header had tmax instead of n_samples:
        if ((i == 4) and float(headerlist[4])%1 != 0 
            and float(headerlist[i]) > float(headerlist[i-1])):
            header_dict[lo_headerelements[i]] = int(
                    (float(headerlist[i]) - float(headerlist[i-1])
                    )*float(headerlist[i-2]) )+1

    headerlements = ['samplingrate','t_min','nsamples','lat','lon','elev']

    for h in headerlements:                   
        try: 
            header_dict[h] = float(header_dict[h])
        except:
            pass
        try:
            if header_dict[h]%1==0:
                header_dict[h] = int(header_dict[h])
        except:
            pass

    return header_dict


def get_ts_header_string(header_dictionary):
    """
        Return a MTpy time series data file header string from a dictionary.

    """
    
    header_string = '# '
    for headerelement in lo_headerelements:
        if headerelement in header_dictionary:
            header_string += '{0} '.format(str(header_dictionary[headerelement]))
        else:
            header_string += '\t '   

    header_string += '\n'

    return header_string



def write_ts_file_from_tuple(outfile,ts_tuple, fmt='%.8e'):
    """
        Write an MTpy TS data file, where the content is provided as tuple:

        (station, channel,samplingrate,t_min,nsamples,unit,lat,lon,elev, data)

        todo:
        needs tuple-validation

    """

    
    header_dict = {}
    for i in range(len(ts_tuple) -1):
        if ts_tuple[i] is not None:
            header_dict[lo_headerelements[i]] = ts_tuple[i]

    header_string = get_ts_header_string(header_dict)
    data = ts_tuple[-1]

    outfilename = make_unique_filename(outfile)


    try:
        outF = open(outfilename,'w')
        outF.write(header_string)
        np.savetxt(outF, data, fmt=fmt)
        outF.close()
    except ValueError:
        raise MTex.MTpyError_inputarguments('ERROR - could not write content'
                            ' of TS tuple to file : {0}'.format(outfilename))

    return outfilename


def read_ts_file(mtdatafile):
    """
        Read an MTpy TS data file and provide the content as tuple:

        (station, channel,samplingrate,t_min,nsamples,unit,lat,lon,elev, data)
        If header information is incomplete, the tuple is filled up with 'None'

    """

    infile = op.abspath(mtdatafile)
    if not op.isfile(infile):
        raise MTex.MTpyError_inputarguments('ERROR - Data file not '
                                                'existing: {0}'.format(infile))

    header = read_ts_header(infile)
    if len(header) == 0 :
        raise MTex.MTpyError_inputarguments('ERROR - Data file not valid - '
                                        'header is missing : {0}'.format(infile))

    data = np.loadtxt(infile)
    if len(data) != int(float(header['nsamples'])):
        raise MTex.MTpyError_inputarguments('ERROR - Data file not valid '
                                    '- wrong number of samples in data ({1} '
                                    'instead of {2}): {0}'.format(
                                    infile,len(data) , int(float(
                                        header['nsamples']))) )

    lo_header_contents = []

    for i in lo_headerelements:
        if i in header:
            lo_header_contents.append(header[i])
        else:
            lo_header_contents.append(None)
 
    lo_header_contents.append(data)

    return tuple(lo_header_contents)


def reorient_files(lo_files, configfile, lo_stations = None, outdir = None):

    #read config file
    try:
        config_dict = MTcf.read_survey_configfile(configfile)
    except:
        raise MTex.MTpyError_config_file( 'Config file cannot be read:'
                                                    ' {0}'.format(configfile) )

    if lo_stations is not None:
        try:
            if type(lo_stations) == str:
                raise Exception("problem!!!")
            #check, if it's iterable:
            dummy = [i for i in lo_stations]
        except:
            raise MTex.MTpyError_inputarguments('ERROR - "lo_stations"'
                                                ' argument must be iterable!')
    print('\t re-orienting data for collection of stations:\n{0}'.format(lo_stations))
    #Do not require list of headers as input, as this function can be called directly rather than from a 'calibratefiles.py'-like script - so the list not necessarily exists in beforehand - 
    #collect header lines of files in list
    lo_headers = []
    lo_stationnames = []
    for file_idx, filename in enumerate(lo_files):
        header = read_ts_header(filename)
        station = header['station']
        if station.upper() not in [i.upper() for i in lo_stations]:
            #TODO: check, if this causes problems with the indices for the current loop:
            lo_files.remove(filename)
            continue
        lo_headers.append(header)
        lo_stationnames.append(station.upper())



    if len(lo_headers) == 0 :
        if lo_stations is not None:
            print('ERROR - No files with header lines found for station(s)'\
                                                    ' {0}'.format(lo_stations))
        else:
            print('ERROR - No files with header lines found')
        return 1

    lo_stationnames = list(set(lo_stationnames))

    # set up output directory 
    ori_outdir = op.abspath(op.join(os.curdir,'reoriented'))

    if outdir is not None:
        try:
            ori_outdir = op.abspath(op.join(os.curdir,outdir))
            if not op.isdir(ori_outdir):
                os.makedirs(ori_outdir)
        except:
            print('Output directory cannot be generated: {0} - using generic'\
                                                ' location'.format(ori_outdir))
            ori_outdir = op.abspath(op.join(os.curdir,'reoriented'))
    try:
        if not op.isdir(ori_outdir):
            os.makedirs(ori_outdir)
    except:
        #this only comes up, if the generic location cannot be generated
        raise MTex.MTpyError_inputarguments('Generic directory cannot be'
                                        ' generated: {0}'.format(ori_outdir))

    #----------------------
    #start re-orientation
    #present: list of all files, list of all headers, list of all stations
    
    for sta_idx, sta in enumerate(lo_stationnames):
        #print sta
        try:
            stationconfig = config_dict[sta]
        except:
            print('Warning - No config file entry for station {0} -'\
                                        ' no processing possible'.format(sta))
            continue
        
        declination = float(stationconfig.get('declination',0.))


        for sensor in ['e','b']:
            #TODO:
            # reduce this function to the re-orientation of files that have the same length for X and Y. 
            #Do the puzzlling for varying lengths later!!

            for idx_h_x, header_x in enumerate(lo_headers):
                #looking for one specific station
                if not header_x['station'].upper() == sta.upper():
                    continue
                #looking for the specific sensor type
                if not header_x['channel'].lower()[0] == sensor:
                    continue
                #looking for the X channel (the to-be-North)
                if not header_x['channel'].lower()[1] == 'x':
                    continue

                x_file = lo_files[idx_h_x]
                x_header_string = get_ts_header_string(header_x)

                t0 = float(header_x['t_min'])
                #print t0 
                #now look for the respective y-file and possible z-file - unfortunately by another loop over all headers:
                y_file = None
                z_file = None
                for idx_h_y, header_y in enumerate(lo_headers):
                    if (header_y['station'].upper() == sta.upper()) and \
                        (header_y['channel'].lower()[0] == sensor) and \
                        (float(header_y['t_min']) == float(header_x['t_min'] ) ):
                        if (header_y['channel'].lower()[1] == 'y') :
                            y_file = lo_files[idx_h_y]
                            y_header_string = get_ts_header_string(header_y)

                        elif   (header_y['channel'].lower()[1] == 'z') :
                            z_file = lo_files[idx_h_y]

                    else:
                        continue
                if y_file == None:
                    continue

                x_outfn = op.abspath(op.join(ori_outdir,op.basename(x_file)))
                y_outfn = op.abspath(op.join(ori_outdir,op.basename(y_file)))
                if z_file is not None:
                    z_outfn = op.abspath(op.join(ori_outdir,op.basename(z_file)))
                

                xdata = np.loadtxt(x_file)
                ydata = np.loadtxt(y_file)

                #declination is positive, if magnetic North is east of true North.
                # the measured angles are w.r.t. magnetic North, so the given 
                # azimuths do not include the declination 
                #-> thus the declination value is added to azimuths
                if sensor == 'e':
                    xangle = float(stationconfig.get(
                                    'e_xaxis_azimuth', 0.)) + declination
                    yangle = float(stationconfig.get(
                                    'e_yaxis_azimuth',90.)) + declination
                else:
                    xangle = float(stationconfig.get(
                                    'b_xaxis_azimuth', 0.)) + declination
                    yangle = float(stationconfig.get(
                                    'b_yaxis_azimuth',90.)) + declination                


                newx, newy =  MTcc.reorient_data2D(xdata, ydata, 
                            x_sensor_angle = xangle , y_sensor_angle = yangle)
                #print xdata.shape, ydata.shape, newx.shape, newy.shape 

                #continue
                outFx = open(x_outfn,'w')
                outFx.write(x_header_string)
                np.savetxt(outFx,newx)
                outFx.close()
                outFy = open(y_outfn,'w')
                outFy.write(y_header_string)
                np.savetxt(outFy,newy)
                outFy.close()
                written_files = [x_outfn,y_outfn]
                if z_file is not None:
                    shutil.copyfile(z_file, z_outfn)
                    written_files.append(z_outfn)
                print('\tSuccessfullly written files {0}'.format(written_files))


            
    return 0
