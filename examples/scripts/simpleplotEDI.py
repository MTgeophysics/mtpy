#!/usr/bin/env python 

"""

simpleplotEDI.py

A quick and dirty script to plot the MT response in resistivity/phase form from an EDI file (off-diagonal elements only).


Copyright (c) 2014, Lars Krieger/the MTpy team
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""



import os,sys
import os.path as op

import mtpy.core.edi as MTedi

from pylab import *



fn = sys.argv[1]

edi = MTedi.Edi()
try:
	edi.readfile(fn)
except:
	print '\n\tFile does not exist: {0}\n'.format(fn)
	sys.exit()


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
    phi_te.append(p[0,1])
    phi_tm.append(p[1,0])

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
errorbar(periods,res_te,reserr_te, marker='s',c='b',fmt='x')
errorbar(periods,res_tm,reserr_tm, marker='s',c='r',fmt='x')
xscale('log')
yscale('log')
minval=min( min(res_te,res_tm))
maxval=max(max(res_te,res_tm))
ylim([minval/10,maxval*10])
xlim(0.5*min(periods),2*max(periods))
autoscale(False)

ylabel('resistivity')
setp( ax1.get_xticklabels(), visible=False)
## share x only
ax2 = subplot(212, sharex=ax1)
autoscale(False)

ylim(-90,270)
errorbar(periods,phi_te,phierr_te,marker='s',c='b',fmt='x')
errorbar(periods,phi_tm,phierr_tm,marker='s',c='r',fmt='x')
ylabel('phase')
xlabel('period (in s)')

tight_layout()
show(block=True)
