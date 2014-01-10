#!/usr/bin/env python
# Generate data file for ModEM
# by Paul Soeffky 2013
# revised by LK 2014


import os
import mtpy.core.edi as E
import numpy as np
from math import cos, sin, asin, sqrt, radians
import mtpy.utils.conversions as conv
import re


edip = 'edis'
edipath = edip


#find average latitude
filepath = '.'

errorfloor = 5



header_string = """#Data file 
\n# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
> Full_Impedance \n> exp(-i\omega t)\n> [V/m]/[T]
> 0.00
"""
latlist=[]
fedilist=[]
tipstring2=[]
latlist2=[]
lonlist2=[]

coord_list = []

periodlistlen=[]
periodlistlen2=[]
latlon=''
for edi in os.listdir(edipath):
	if edi.endswith('.edi'):
		fullepath=os.path.join(edipath,edi)
		fedilist.append(fullepath)
		
		e=E.Edi()
		e.readfile(fullepath)
		coord_list.append((e.lat,e.lon))
		freq2=e.freq

		nperiod=len(freq2)
		periodlistlen=nperiod
		periodlistlen2.append(periodlistlen)
		periodlist=1/freq2
		periodlist2=str(periodlist)
		periodno=str(nperiod)
			

coord_list = np.array(list(set(coord_list)))
latlist3 = list(coord_list[:,0])
lonlist3 = list(coord_list[:,1])

a=np.mean(latlist3)

Minlat=min(latlist3,key=lambda x:abs(x-a))

Lat=str(Minlat)

maxperno=max(periodlistlen2)
maxperno2=str(maxperno)


UTMElist=[]
UTMNlist=[]

for edi in os.listdir(edipath):
	if edi.endswith('.edi'):
		fullepath=os.path.join(edipath,edi)
		fedilist.append(fullepath)
		
		e=E.Edi()
		e.readfile(fullepath)
		
		lat6=(e.lat)
		lon6=(e.lon)

		
		U=conv.LLtoUTM(23, lat6, lon6)

		
		UTMN=U[2]
		UTME=U[1]
		UTMElist.append(UTME)
		UTMNlist.append(UTMN)
		
		
		



#print UTMNlist

UTMNlist3=[float(i) for i in UTMNlist]
UTMElist3=[float(i) for i in UTMElist]

aveN=np.mean(UTMNlist3)
aveE=np.mean(UTMElist3)




for edi in os.listdir(edipath):
	if edi.endswith('.edi'):
		fullepath=os.path.join(edipath,edi)
		fedilist.append(fullepath)
		
		e=E.Edi()
		e.readfile(fullepath)
		
		lat6=str(e.lat)
		lat7=(lat6[:7])
		lon6=str(e.lon)
		lon7=(lon6[:8])		
		
		if lat7 == Lat:
			latlon +=lat7.ljust(10)
			lat9=lat7
			latlon +=lon7.ljust(10)
			lon9=lon7


header_string += '> {0}  {1}\n'.format(aveN,aveE)

#find Edi files


edilist = [ edifile for edifile in os.listdir(edipath) if edifile.endswith('.edi')]		
#print stationlst

nstat=len(edilist)
#print nstat
edino=str(nstat)
#print number of edi files
header_string += '> {0} {1}\n'.format(maxperno2,edino)

fedilist = []

impstring = ''

for edi in edilist:
	fullepath=os.path.join(edipath,edi)
	fedilist.append(fullepath)

	e=E.Edi()
	e.readfile(fullepath)
	
	lat6=(e.lat)
	lon6=(e.lon)
	
	
	lat2=str(lat6)
	lon2=str(lon6)
	
	U=conv.LLtoUTM(23, lat6, lon6)

	
	UTMN=U[2]
	UTME=U[1]
	
	freq2=e.freq
	elev2=e.elev

	nperiod=len(freq2)
	periodlist=1/freq2

	periodno=str(nperiod)

	zerr=e.Z.zerr
	zval=e.Z.z


	#Generate Impadence Array
	for i in range(nperiod):

		period = periodlist[i]
		Z = zval[i]
		Zerr = zerr[i]
		components = ['XX','XY','YX','YY']

		period_impstring = ''

		easting = UTME - aveE
		northing = UTMN - aveN	

		for i in range(2):
			for j in range(2):
				try:
					rel_err = Zerr[i,j]/np.abs(Z[i,j])
					if rel_err < errorfloor/100.:
						raise
				except:
					#pass
					Zerr[i,j] = errorfloor/100. * np.abs(Z[i,j])

			comp = components[2*i+j]
			period_impstring += '{0:.5E}  {1}  '.format(period,e.station)
			period_impstring += '{0:.3f}  {1:.3f}  '.format(e.lat,e.lon)
			period_impstring += '{0:.3f}  {1:.3f}  {2}  '.format(northing, easting,0.)
			period_impstring += 'Z{0}  {1:E}  {2:.5E}  {3:.5E}  '.format(comp,float(np.real(Z[i,j])),
												float(np.imag(Z[i,j])), Zerr[i,j] )
			period_impstring += '\n'

		impstring += period_impstring

#print outstring
data=open(r'ModEMdata.dat', 'w')
data.write(header_string)
data.write(impstring)
data.close()


#Tipper part		
tiphead = ''
tiphead += """# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
> Full_Vertical_Components \n> exp(-i\omega t)\n> []
> 0.00
>"""
tiphead += Lat+" "
tiphead += lon2
tiphead += "\n> "
tiphead += periodno+" "
tiphead += edino+"\n"



tipper=[]
#Generate Tipper Array
tipstring = ''
fedilist=[]
for edi in edilist:
	fullepath=os.path.join(edipath,edi)
	fedilist.append(fullepath)

	e=E.Edi()
	e.readfile(fullepath)
	
	lat6=(e.lat)
	lon6=(e.lon)


	
	U=conv.LLtoUTM(23, lat6, lon6)

	
	UTMN=U[2]
	UTME=U[1]
	
	tip = e.Tipper.tipper
	tiperr = e.Tipper.tippererr		
	
	freq2=e.freq

	nperiod=len(freq2)
	periodlist=1/freq2


	for i in range(nperiod):
		try:
			period_tipstring = ''

			period = periodlist[i]
			t = tip[i]
			terr = tiperr[i]
			components = ['X','Y']

			easting = UTME - aveE
			northing = UTMN - aveN	
			
			#do not use Tipper entry, if value is zero
			if np.abs(tip[i]) == 0.:
				continue

			for i in range(2):
				try:
					rel_err = tiperr[i]/np.abs(tip[i])
					if rel_err < errorfloor/100.:
						raise
				except:
						tiperr[i] = errorfloor/100. * np.abs(tip[i])

				comp = components[i]
				period_tipstring += '{0:.5E}  {1}  '.format(period,e.station)
				period_tipstring += '\t {0:.3f}  {1:.3f}  '.format(e.lat,e.lon)
				period_tipstring += '\t {0:.3f}  {1:.3f}  {2}  '.format(northing, easting,0.)
				period_tipstring += '\t T{0}  {1:.5E}  {2:.5E}  {3:.5E} '.format(comp,np.real(tip[i]),
													np.imag(tip[i]),tiperr[i] )
				period_tipstring += '\n'

		except:
			continue

		tipstring += period_tipstring

if len(tipstring)>0:
	data = open(r'ModEMdata.dat', 'a')
	data.write(tiphead)
	data.write(tipstring.expandtabs(4))
	data.close()



	



print "end"