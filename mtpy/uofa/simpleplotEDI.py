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

def plotedi(fn, saveplot=False, component=None):


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
		lo_comps = ['n','e']

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

	resplotelement_xy = None
	resplotelement_yx = None

	axes = figure('EDI '+fn)
	ax1 = subplot(211)

	if 'n' in lo_comps:
		resplotelement_xy = errorbar(periods,res_te,reserr_te, marker='x',c='b',fmt='x')
	if 'e' in lo_comps:
		resplotelement_yx =errorbar(periods,res_tm,reserr_tm, marker='x',c='r',fmt='x')
	xscale('log')
	yscale('log')
	minval=min( min(res_te,res_tm))
	maxval=max(max(res_te,res_tm))
	xlim(0.5*min(periods),2*max(periods))

	#ylim([0.1,100])
	ylim([minval/10,maxval*10])


	autoscale(False)

	ylabel(r' $\rho$ (in $\Omega m$)')
	setp( ax1.get_xticklabels(), visible=False)
	## share x only
	ax2 = subplot(212, sharex=ax1)
	autoscale(False)

	#ylim(-45,135)
	if 'n' in lo_comps:
		errorbar(periods,phi_te,phierr_te,marker='x',c='b',fmt='x')
	if 'e' in lo_comps:
		errorbar(periods,phi_tm,phierr_tm,marker='x',c='r',fmt='x')
	ylabel('Phase angle ($\degree$)')
	xlabel('Period (in s)')
	plot([xlim()[0],xlim()[1]],[45,45],'-.',c='0.7')
	ylim([-0,90])

	ax1.legend([resplotelement_xy,resplotelement_yx],['$E_{X}/B_Y$','$E_Y/B_X$'],loc=2,ncol=1,
					numpoints=1,markerscale=0.8,frameon=True,labelspacing=0.3, 
					prop={'size':8},fancybox=True,shadow=False)


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
