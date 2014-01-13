#!/usr/bin/env python

import os
#I know, it's not nice, but I am lazy:
from numpy import *
import sys
import os.path as op
import fnmatch

def main():

    if len(sys.argv) < 2:
        print """\n\tusage: \n
            python plot_stationmap_from_edis.py <edi folder> [optional_arguments] 

        list of optional arguments:

            - <symbolsize (int,24)> 
            - <mapstretchfactor (float,1.)> 
            - <labelsize (int,14)> 
            - <showlabel (bool,True)> 


        Output:

        2 files containing a geographical map (Mercator) with all station 
        locations, which could be extracted from the EDI files in the folder 
        provided:

        station_locations_map.png/svg
        
        (no overwriting of existing files with these names)


        """
        return

    edifolder = sys.argv[1]
    edilist = []
    try:
        if not op.isdir(edifolder):
            raise
        edilist = fnmatch.filter(os.listdir(edifolder),'*.[Ee][Dd][Ii]')
        edilist = [op.abspath(op.join(edifolder,i)) for i in edilist]
        if len(edilist) == 0:
            raise
    except:
        print 'No EDI files in folder {0}'.format(edifolder)
        return



    try:
        symbolsize = int(float(sys.argv[2]))
    except:
        print 'cannot read symbolsize value - using default'
        symbolsize = 24

    try:
        mapstretchfactor = float(sys.argv[3])
    except:
        print 'cannot read mapstretchfactor value - using default'
        mapstretchfactor = 1
        

    try:
        labelsize = int(float(sys.argv[4]))
    except:
        print 'cannot read labelsize value - using default'
        labelsize = 14
    

    try:
        showlabel = sys.argv[5]
        if not bool(showlabel):
            raise
    except:
        print 'cannot read showlabel value - using default (=True)'
        showlabel = True


    f1,f2 = makemap(edilist,mapstretchfactor,symbolsize,labelsize,showlabel)
 
    print 'wrote map files\n\t {0}\n and\n\t {1}\n'.format(f1,f2)


def makemap(edilist,mapstretchfactor,symbolsize,labelsize,showlabel):

    #import of modules here due to warnings from Matplotlib packages
    #these warnings distract the 'usage' information
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import matplotlib as mpl  
    import mtpy.utils.filehandling as MTfh
    import mtpy.core.edi as EDI
    import mtpy.utils.filehandling as MTfh


    lats = []
    lons =[]
    names=[]

    for i in edilist:
        e = EDI.Edi()
        e.readfile(i)
        lats.append(e.lat)
        lons.append(e.lon)
        names.append(e.head['dataid'].lower())

    coords = zeros((len(edilist),2))
    coords[:,0] = lats
    coords[:,1] = lons

    latrange = max(lats) - min(lats)
    lonrange = max(lons) - min(lons)

    #center point for projection:
    c = [mean(lats),mean(lons)]



    #-----------------------
    #Matplotlib options
    mpl.rcParams['font.size'] = 10.
    mpl.rcParams['axes.labelsize'] = 8.
    mpl.rcParams['xtick.labelsize'] = 6.
    mpl.rcParams['ytick.labelsize'] = 6.

    
    plt.close('all')

    fig = plt.figure()#figsize=(5.7,4.3))
    plt.subplots_adjust(left=0.2,right=0.9,top=0.90,bottom=0.1,wspace=0.15,hspace=0.05)
    ax = plt.subplot(111)


    stretch = float(mapstretchfactor)
    total_latmin = max( (min(lats)-stretch*latrange ), -90)
    total_lonmin = max( (min(lons)-stretch*lonrange ), -180)
    total_latmax = min( (max(lats)+stretch*latrange ), 90)
    total_lonmax = min( (max(lons)+stretch*lonrange ), 180)

    total_latrange = total_latmax - total_latmin
    total_lonrange = total_lonmax - total_lonmin

    #determine number of axes labels:
    maximumlabels = 5
    latnum=maximumlabels
    lonnum=maximumlabels
    lonlat_stretch = total_lonrange/total_latrange
    if int(lonlat_stretch) > 2:
        #significantly more long than lat
        factor = int(int(lonlat_stretch)/2.)
        latnum = int(maximumlabels/factor) + 1
        lonnum = maximumlabels
    elif int(lonlat_stretch) <0.5 :
        #significantly more long than lat
        factor = int(int(1./lonlat_stretch)/2.)
        lonnum = int(maximumlabels/factor) + 1
        latnum = maximumlabels

     
    m = Basemap(
        projection='merc',
        lon_0=c[1],lat_0=c[0],lat_ts=c[0],
        llcrnrlat=total_latmin, urcrnrlat=total_latmax,
        llcrnrlon=total_lonmin, urcrnrlon=total_lonmax,
        rsphere=6371200.,resolution='h')#,ax=ax)

    lons.append(total_lonmin)
    lons.append(total_lonmax)
    lats.append(total_latmin)
    lats.append(total_latmax)

    xgrid, ygrid = m(*meshgrid(lons,lats))
    xdiff = xgrid.max()-xgrid.min()
    ydiff = ygrid.max()-ygrid.min()

    largest_extent = max(ydiff,xdiff)


    m.drawcoastlines(linewidth=0.25,ax=ax)
    m.drawcountries(linewidth=0.25,ax=ax)
    m.fillcontinents(color='coral',lake_color='aqua',ax=ax)
    m.drawmapboundary(fill_color='aqua',ax=ax)

    m.drawparallels(
                    [round(i,3) for i in linspace(total_latmin,total_latmax,latnum+1,False)][1:],
                    labels=[1,0,0,0],
                    fmt='%.1f'
                    ) 
    m.drawmeridians(
                    [round(i,3) for i in linspace(total_lonmin,total_lonmax,lonnum+1,False)][1:],
                    labels=[0,0,0,1],
                    fmt='%.1f'
                    ) 
    m.drawrivers()
    m.etopo()

    m.drawmapscale( total_lonmax-0.15*total_lonrange, total_latmax-0.2*total_latrange,
                    c[1],c[0],
                    2*10**(int(log10(largest_extent/1000.))-1),
                    barstyle='simple',
                    labelstyle='simple',fontsize=12,
                    fontcolor='k', fillcolor1='r', fillcolor2='g',
                    ax=None, format='%d',
                    )


    for x,y,name in zip(xgrid[0,:], ygrid[:,0],names):
        plt.plot(x,y,'v',ms=symbolsize,color='k',label=name)
        if showlabel is True:
            plt.text(x,y,name,fontsize=labelsize,ha='center',va='bottom',color='k',backgroundcolor='r')#, labelfontsize=5)


    plt.title('locations of {0} MT stations'.format(len(names)))


    f1 = 'station_locations_map.png'
    f2 = 'station_locations_map.svg'

    f1 = MTfh.make_unique_filename(f1)
    f2 = MTfh.make_unique_filename(f2)

    plt.savefig(f1,format='png',dpi=200)
    plt.savefig(f2,format='svg',transparent=True)

    return op.abspath(f1),op.abspath(f2)


if __name__=='__main__':
    main()

