#!/usr/bin/env python 


import os,sys
import os.path as op

import mtpy.core.edi as MTedi

def main():

	fn = sys.argv[1]

	if not op.isfile(fn):
		print '\n\tFile does not exist: {0}\n'.format(fn)
		sys.exit()
	saveplot = False
	if len(sys.argv)>2:
		arg2 = sys.argv[2]
		if 's' in arg2.lower():
			saveplot = True

	fn = plotedi(fn,saveplot)

def plotedi(fn, saveplot=False):


	edi = MTedi.Edi()
	try:
		edi.readfile(fn)
	except:
		print '\n\tERROR - not a valid EDI file: {0}\n'.format(fn)
		sys.exit()

	if saveplot is True:
		import matplotlib 
		matplotlib.use('Agg')

	from pylab import *


	res_te=[]
	res_tm=[]
	phi_te=[]
	phi_tm=[]
	reserr_te=[]
	reserr_tm=[]
	phierr_te=[]
	phierr_tm=[]

	for r in edi.Z.resistivity:
		res_te.append(r[0,1])
		res_tm.append(r[1,0])

	for p in edi.Z.phase:
	    phi_te.append(p[0,1]%90)
	    phi_tm.append(p[1,0]%90)

	if np.mean(phi_te)>90 and np.mean(phi_tm)>90:
		phi_te = [i%90 for i in phi_te]
		phi_tm = [i%90 for i in phi_tm]

	for r in edi.Z.resistivity_err:
		reserr_te.append(r[0,1])
		reserr_tm.append(r[1,0])
	for p in edi.Z.phase_err:
	    phierr_te.append(p[0,1])
	    phierr_tm.append(p[1,0])

	periods = 1./edi.freq


	ax1 = subplot(211)
	errorbar(periods,res_te,reserr_te, marker='x',c='b',fmt='x')
	errorbar(periods,res_tm,reserr_tm, marker='x',c='r',fmt='x')
	xscale('log')
	yscale('log')
	minval=min( min(res_te,res_tm))
	maxval=max(max(res_te,res_tm))
	xlim(0.5*min(periods),2*max(periods))

	ylim([0.01,1000])
	ylim([minval/10,maxval*10])


	autoscale(False)

	ylabel('app.res. in Ohm m')
	setp( ax1.get_xticklabels(), visible=False)
	## share x only
	ax2 = subplot(212, sharex=ax1)
	autoscale(False)

	#ylim(-45,135)
	errorbar(periods,phi_te,phierr_te,marker='x',c='b',fmt='x')
	errorbar(periods,phi_tm,phierr_tm,marker='x',c='r',fmt='x')
	ylabel('phase')
	xlabel('period (in s)')
	plot([xlim()[0],xlim()[1]],[45,45],'-.',c='0.7')
	ylim([-50,100])

	tight_layout()
	if saveplot is True:

		ioff()
		outfn = op.splitext(fn)[0]+'.png'
		savefig(outfn, bbox_inches='tight')
		close('all')
		ion()
		return outfn



	else:
		ion()
		show(block=True)

		return None



if __name__=='__main__':
    main()
