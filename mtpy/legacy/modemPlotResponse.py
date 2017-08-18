#!/usr/bin/env python

"""
mtpy/utils/modemPlotResponse.py

Module for plotting response curves and original data points incl. error bars.

Data are plotted directly, i.e. as impedance tensor components Zij.
Each plot contains real and imaginary part for one component.
4 plots are used for visualising one station.

Each plot window contains max 4 stations (i.e. 16 plots).

@UofA/LK Jan 2014

"""
#--------------------------------------------------------------------

import numpy as np
import os,sys
import os.path as op
from pylab import *

#precision in decimal digits (used for comparing periods):
precision = 3

#--------------------------------------------------------------------

components = ['ZXX','ZXY','ZYX','ZYY']
#components = ['ZXY','ZYX']


def checkdatafile(fn):

	state = True

	return state

def checkresponsefile(fn):

	state = True

	return state

def checkdata2response(datafn,responsefn): 

	match = True


	return match



def build_data_dict(datafile, responsefile):

	"""
	Reading files and generate dictionary

	input:
	- original ModEM data file
	- ModEM modelling response file

    output:
    - dictionary
	    [station]: periods, data, error, response
    """

	if checkdatafile(datafile) is False:
		print 'ERROR - invalid data file'
		sys.exit()

	if checkresponsefile(responsefile) is False:
		print 'ERROR - invalid response file'
		sys.exit()

	if checkdata2response(datafile,responsefile) is False:
		print 'ERROR - data and response file do not belong together'
		sys.exit()



	stationlist = []
	indata = []

	data_dict = {}
	
	periods = []
	data = []
	error = []
	response = []

	Fin = open(datafile)
	for line in Fin:
		line = line.strip()
		if (len(line) == 0) or (line[0] in ['#','>']) :
			continue
		line = line.split()	

		sta = line[1].upper()

		if sta in data_dict:
			station_data = data_dict[sta]
		else:
			print 'new station: ',sta
			station_data = [[],[],[],[]]


		periodlist = station_data[0]
		datalist = station_data[1]
		errorlist = station_data[2]
		responselist = station_data[3]


		period = np.round(float(line[0]),precision )
		if period in periodlist:
			idx_period = periodlist.index(period)
			
		else:
			#print 'new period: ',period
			periodlist.append(period)
			data = [None,None,None,None]
			error = [None,None,None,None]
			response = [None,None,None,None]
			datalist.append(data)
			errorlist.append(error)
			responselist.append(response)
			idx_period = len(periodlist) - 1
		
		data = datalist[idx_period]
		error = errorlist[idx_period]
		comp = line[7].upper()
		try:
			idx_comp = components.index(comp)
		except:
			continue

		data[idx_comp] = np.complex(float(line[8]),float(line[9]))
		error[idx_comp] = float(line[10])
		datalist[idx_period] = data
		errorlist[idx_period] = error

		station_data[0] = periodlist
		station_data[1] = datalist 
		station_data[2] = errorlist 
		station_data[3] = responselist 
		data_dict[sta] = [periodlist,datalist,errorlist,responselist]

	Fin.close()
	
	print len(data_dict)
	

	Fin2 = open(responsefile)
	for line in Fin2:
		line = line.strip()
		if line[0] in ['#','>']:
			continue
		line = line.split()	
		sta = line[1].upper()
		if sta in data_dict:
			station_data = data_dict[sta]
		else:
			station_data = [[],[],[],[]]


		periodlist = station_data[0]
		datalist = station_data[1]
		errorlist = station_data[2]
		responselist = station_data[3]

		period = np.round(float(line[0]),precision )

		if period in periodlist:
			idx_period = periodlist.index(period)		
		else:
			#print 'new period: ',period
			periodlist.append(period)
			data = [None,None,None,None]
			error = [None,None,None,None]
			response = [None,None,None,None]
			datalist.append(data)
			errorlist.append(error)
			responselist.append(response)
			idx_period = len(periodlist) - 1

		response = responselist[idx_period]
		comp = line[7].upper()
		try:
			idx_comp = components.index(comp)
		except:
			continue

		response[idx_comp] = np.complex(float(line[8]),float(line[9]))
		responselist[idx_period] = response

		#print 'response:',idx_period, idx_comp,len(periodlist),len(datalist),len(errorlist),len(responselist)

		station_data[0] = periodlist
		station_data[1] = datalist 
		station_data[2] = errorlist 
		station_data[3] = responselist 
		data_dict[sta] = [periodlist,datalist,errorlist,responselist]
	Fin2.close()

	#sorting for ascending periods 
	for sta in data_dict:
		data = data_dict[sta]
		sorted_list_tuple = sorted(zip(*data))
		four_tuples = zip(*sorted_list_tuple) 
		data_dict[sta] = [list(i) for i in four_tuples]
		#print sta
		#print four_tuples


	return data_dict


def plotZ(data_dictionary, no_comps = 4,step=1):


	close('all')
	ion()
	#maximum station data in one figure:
	max_stations = 4

	fignum = 1
	station_counter = 0

	plotcounter = 0

	for sta in sorted(data_dictionary):

		plotcounter += 1
		if plotcounter%step != 0:
			continue
		print sta

		station_counter += 1
		#if station_counter >4 :
		#	break


		figure(fignum)

		if station_counter%max_stations == 0:
			fignum += 1
	
		subfig_index_vertical = (station_counter-1)%max_stations

		data = data_dictionary[sta]

		for idx_c, comp in enumerate(components):
			if no_comps == 2:
				if idx_c in [0,3]:
					continue


			subfig_index_horizontal = idx_c
			
			data_periods = []
			in_data_real = []
			in_data_imag = []
			in_data_error = []

			resp_periods = []
			resp_real = []
			resp_imag = []

			for idx_p, period in enumerate(data[0]): 
				try:
					dr = float(np.real(data[1][idx_p][idx_c]))
					di = float(np.imag(data[1][idx_p][idx_c]))
					if dr is None or di is None:
						raise
					try:
						derr = float(data[2][idx_p][idx_c])
					except:
						derr = 0.
					data_periods.append(period)
					in_data_real.append(dr) 
					in_data_imag.append(di) 
					in_data_error.append(derr)
				#skip, if data cannot be read
				except:
					continue

				try:
					rr = float(np.real(data[3][idx_p][idx_c]))
					ri = float(np.imag(data[3][idx_p][idx_c]))

					if rr is None or ri is None:
						raise

					resp_periods.append(period)
					resp_real.append(rr) 
					resp_imag.append(ri) 
				#skip, if data cannot be read
				except:
					continue

			subfigure_index = subfig_index_vertical *4+subfig_index_horizontal+1
			ax = subplot(max_stations,4, subfigure_index)

			plot(resp_periods[::-1],resp_real,c='b')
			plot(resp_periods[::-1],resp_imag,c='r')

			#scatter(data_periods,in_data_real,c='b',marker='o')
			orig_real = errorbar(data_periods[::-1],in_data_real,yerr=in_data_error,c='b',marker='x',ls='none',)
			#scatter(data_periods,in_data_imag,c='r',marker='x')
			orig_imag = errorbar(data_periods[::-1],in_data_imag,yerr=in_data_error,c='r',marker='x',ls='none',)
			
			ax.set_xscale('log')

			if subfig_index_vertical == 0:
				ax.set_title(comp)
			
			if subfig_index_vertical == (max_stations-1):
				ax.set_xlabel('Period (in s)')


			#print subfig_index_vertical,subfig_index_horizontal,subfigure_index

			if subfig_index_horizontal == 0:
				ax.set_ylabel(sta)
			if no_comps == 2:
				if subfig_index_horizontal == 1:
					ax.set_ylabel(sta)


			if (subfigure_index == 1)  or (no_comps == 2 and subfigure_index == 2):

				ax.legend([orig_real,orig_imag],['RE','IM'],ncol=1,
					numpoints=1,markerscale=0.8,frameon=True,labelspacing=0.3, 
					prop={'size':8},fancybox=True,shadow=False)
			tick_params(axis='both', which='major', labelsize=8)
			tick_params(axis='both', which='minor', labelsize=6)

		tight_layout()

	show()


