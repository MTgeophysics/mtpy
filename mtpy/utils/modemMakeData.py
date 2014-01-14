#!/usr/bin/env python
# Generate data file for ModEM
# by Paul Soeffky 2013
# revised by LK 2014


import os,sys
import mtpy.core.edi as E
import numpy as np
from math import cos, sin, asin, sqrt, radians
import mtpy.utils.conversions as conv
import mtpy.utils.merge_periods as mp
import re
import glob



edipath = 'edis_selected'

if not os.path.isdir(edipath):
    print '\n\tERROR - data path does not exist'
    sys.exit()


#flag for merging closely neighbouring periods: 
merge_periods = False

#merge periods, which do not differ more than this threshold (in percent)
merge_threshold = 5



use_tipper = False

errorfloor = 5

#-----------------------------------------------------------------

header_string = """#Data file \n# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
> Full_Impedance \n> exp(-i\omega t)\n> [mV/km]/[nT]
> 0.00
"""

        
edilist = glob.glob(os.path.join(edipath,'*.[Ee][Dd][Ii]'))
edilist = [os.path.abspath(os.path.join(os.curdir,i)) for i in edilist]

lo_ediobjs = []
coord_list = []
utm_list =[]


for edi in edilist:
        
        e = E.Edi()
        e.readfile(edi)
        
        lo_ediobjs.append(e)

        utm_list.append(conv.LLtoUTM(23, e.lat, e.lon))
        coord_list.append([e.lat,e.lon])

lat0 = np.mean(np.array(coord_list)[:,0])
lon0 = np.mean(np.array(coord_list)[:,1])
coord_list = []

#check f\or utm zones...if different assign zone!
#for the moment assume all have the same zone

for utm in utm_list:
    coord_list.append((utm[1],utm[2]))

#coordinates given in Cartesian system, North, East
coords = np.array(coord_list)

East0 = np.mean(coords[:,0])
North0 = np.mean(coords[:,1])

rel_coords = np.zeros_like(coords)

rel_coords[:,0] = coords[:,0] - East0
rel_coords[:,1] = coords[:,1] - North0

#start Impedance tensor part ---------------------------------------------

header_string += '> {0}  {1}\n'.format(lat0,lon0)

impstring = ''
periodlist = []

components = ['XX','XY','YX','YY']

#loop for reading in periods
# in case merging is requested, updating period
period_dict = {}

for idx_edi, edi in enumerate(lo_ediobjs):
    freq2 = edi.freq
    periods=1/freq2
    periods = [np.round(i,5) for i in periods]
    periodlist.extend(periods)

periodlist = sorted(list(set(periodlist)),reverse=False)
    
if merge_periods == True:
    #mp.plot_merging(periodlist,merge_threshold)
    new_periods = mp.merge_periods(periodlist,merge_threshold)
else:
    new_periods = periodlist[:]
#setting up a dictionary for old and new period
for idx,per in enumerate(periodlist):
    period_dict[str(per)] = new_periods[idx]


periodlist = []


for idx_edi, edi in enumerate(lo_ediobjs):

    freq2 = edi.freq
    periods=1/freq2

    zerr=edi.Z.zerr
    zval=edi.Z.z

    northing = rel_coords[idx_edi,1]
    easting = rel_coords[idx_edi,0]
    
    
    #Generate Impedance Array
    for i in range(len(periods)):

        raw_period = periods[i]
        raw_period = np.round(raw_period,5)
        period = float(period_dict[str(raw_period)])
        periodlist.append(period)

        Z = zval[i]
        Zerr = zerr[i]

        period_impstring = ''

        for i in range(2):
            for j in range(2):
                try:
                    rel_err = Zerr[i,j]/np.abs(Z[i,j])
                    if rel_err < errorfloor/100.:
                        raise
                except:
                    Zerr[i,j] = errorfloor/100. * np.abs(Z[i,j])

                comp = components[2*i+j]
                period_impstring += '{0:.5f}  {1}  '.format(period,edi.station)
                period_impstring += '{0:.3f}  {1:.3f}  '.format(edi.lat,edi.lon)
                period_impstring += '{0:.3f}  {1:.3f}  {2}  '.format(northing, easting,0.)
                period_impstring += 'Z{0}  {1:.5E}  {2:.5E}  {3:.5E}  '.format(comp,float(np.real(Z[i,j])),
                                                float(np.imag(Z[i,j])), Zerr[i,j] )
                period_impstring += '\n'

        impstring += period_impstring


n_periods = len(set(periodlist))


print 'Z periods: ',n_periods ,  'files:', len(lo_ediobjs)

header_string += '> {0} {1}\n'.format(n_periods,len(lo_ediobjs))

#print outstring
data=open(r'ModEMdata.dat', 'w')
data.write(header_string)
data.write(impstring)
data.close()


#start Tipper part ---------------------------------------------

errorfloor *= 2.

#Tipper part        
header_string = ''
header_string += """# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
> Full_Vertical_Components \n> exp(-i\omega t)\n> []
> 0.00
"""

header_string += '> {0}  {1}\n'.format(lat0,lon0)

tipperstring = ''
periodlist = []
n_periods = 0
components = ['X','Y']

stationlist = []


for idx_edi, edi in enumerate(lo_ediobjs):

    freq2 = edi.freq
    periods=1/freq2

    tippererr=edi.Tipper.tippererr
    tipperval=edi.Tipper.tipper

    northing = rel_coords[idx_edi,1]
    easting = rel_coords[idx_edi,0]
    
    #Generate Tipper Array
    for i in range(len(periods)):

        period = periods[i]
        period = np.round(period,5)
        try:
            T = tipperval[i][0]
        except:
            continue
        try:
            Terr = tippererr[i][0]
        except:
            Terr = np.zeros_like(T,'float')

        if np.sum(np.abs(T)) == 0:
            continue

        stationlist.append(e.station)

        periodlist.append(period)


        period_tipperstring = ''


        for i in range(2):
        
            try:
                rel_err = Terr[i]/np.abs(T[i])
                if rel_err < errorfloor/100.:
                    raise
            except:
                Terr[i] = errorfloor/100. * np.abs(T[i])

            comp = components[i]
            period_tipperstring += '{0:.5f}  {1}  '.format(period,edi.station)
            period_tipperstring += '{0:.3f}  {1:.3f}  '.format(edi.lat,edi.lon)
            period_tipperstring += '{0:.3f}  {1:.3f}  {2}  '.format(northing, easting,0.)
            period_tipperstring += 'T{0}  {1:.5E}  {2:.5E}  {3:.5E}  '.format(comp,float(np.real(T[i])),
                                            float(np.imag(T[i])), Terr[i] )
            period_tipperstring += '\n'

        tipperstring += period_tipperstring


n_periods = len(set(periodlist))
n_stations = len(set(stationlist))
if use_tipper is True:
    print 'Tipper periods: ',n_periods, 'stations:', n_stations
else:
    print 'no Tipper information in data file'

header_string += '> {0} {1}\n'.format(n_periods,len(lo_ediobjs))




if (len(tipperstring)>0 ) and (use_tipper is True):
    data = open(r'ModEMdata.dat', 'a')
    data.write(header_string)
    data.write(tipperstring.expandtabs(4))
    data.close()


print "end"