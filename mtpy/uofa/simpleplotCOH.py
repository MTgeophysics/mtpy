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

import mtpy.core.z as MTz
import mtpy.core.edi as MTedi
import mtpy.analysis.pt as MTpt

import mtpy.processing.coherence as MTcoh
import mtpy.utils.exceptions as MTex 

from pylab import *
import ipdb

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

	periods = data[:,0]
	if data.shape[1] == 6:

		#ipdb.set_trace()
		
		coh1 = data[:,2]
		coh2 = data[:,4]

		ax1 = figure(11)
		
		plot(periods,coh1)
		scatter(periods,coh1,marker='x',c='b')
		

		plot(periods,coh2)
		scatter(periods,coh2,marker='x',c='r')
		
		xscale('log')
		ylim([-.1,1.1])
		xlim(0.5*min(periods),2*max(periods))
		autoscale(False)
		xlabel('period (in s)')
		ylabel('squared Coherence')

		if saveplot is True:
			ioff()
			outfn = op.splitext(fn)[0]+'.png'
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




