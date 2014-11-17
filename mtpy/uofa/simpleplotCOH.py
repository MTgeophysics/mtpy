#!/usr/bin/env python

"""
mtpy/utils/plotcoherence.py


Script for plotting coherence file data (*.coh)

So far supported: converted BIRRP  1 stage processing output data

Input:
Coherence data file name

File must contain column data

periods: col 1
coh Ex : col 3
coh Ey : col 5


@UofA, 2014
(LK)

"""

#=================================================================


import numpy as np
import os.path as op
import sys

def main():
	
	saveflag = False

	fn = sys.argv[1]
	if len(sys.argv)>2:
		arg2 = sys.argv[2]
		if 's' in arg2.lower():
			saveflag = True
	plotcoh(fn,saveflag)



def plotcoh(fn, saveplot=False):

	data = np.loadtxt(fn)

	if saveplot is True:
		import matplotlib 
		matplotlib.use('Agg')

	from pylab import *

	periods = data[:,0]
	if data.shape[1] != 99:

		#ipdb.set_trace()
		
		coh1 = data[:,2]
		coh2 = data[:,4]


		ax1 = figure('coherence')
		
		cohplotelement1 = None
		cohplotelement2 = None
		cohplotelement3 = None
		cohplotelement4 = None


		plot(periods,coh1)
		cohplotelement1=scatter(periods,coh1,marker='x',c='b')
		

		plot(periods,coh2)
		cohplotelement2=scatter(periods,coh2,marker='x',c='r')
		try:
			plot(periods,data[:,6],'g:')
			cohplotelement3 = scatter(periods,data[:,6],marker='d',c='g')
		except:
			pass
		try:
			plot(periods,data[:,8],'y:')
			cohplotelement4 = scatter(periods,data[:,8],marker='d',c='y')
		except:
			pass

		xscale('log')
		ylim([-.1,1.1])
		xlim(0.5*min(periods),2*max(periods))
		autoscale(False)
		xlabel('Period (in s)')
		ylabel('Coherence$^2$')

		eps = 0.05
		if (cohplotelement3 is not None) and (cohplotelement4 is not None):
			ax1.legend([cohplotelement1,cohplotelement2,cohplotelement3,cohplotelement4],
					['$E_{X}$','$E_{Y}$','$B_{X}$','$B_{Y}$'],loc=1,ncol=2,bbox_to_anchor=[1-eps, 1-eps],
					numpoints=1,markerscale=0.8,frameon=True,labelspacing=0.3, 
					prop={'size':9},fancybox=True,shadow=False)

		else:
			ax1.legend([cohplotelement1,cohplotelement2],
					['$E_{X}$','$E_{Y}$'],loc=1,ncol=2,bbox_to_anchor=[1-eps, 1-eps],
					numpoints=1,markerscale=0.8,frameon=True,labelspacing=0.3, 
					prop={'size':9},fancybox=True,shadow=False)



		if saveplot is True:
			ioff()
			outfn = op.splitext(fn)[0]+'.coh.png'
			savefig(outfn, bbox_inches='tight')
			close('all')
			ion()
			return outfn

		else:
			tight_layout()
			ion()
			show(block=True)

			return None



if __name__=='__main__':
    main()




