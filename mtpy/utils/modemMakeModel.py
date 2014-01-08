#!/usr/bin/env python

#build ModEM input Model from ModEM data file

import numpy as np

#==============================================================================

# paramters:

n_xpadding = 5
n_ypadding = 6
#factor with which the padding stretches outside the central rectangle grid
padding_stretch = 1.4

n_layers = 11

#determine minimum block sizes
#used in the inner rectangle - constant widths
dx = 500
dy = 400
#region around stations discretised with these sizes
#outside, the grid steps will be extended exponentially
#the size of padding is determined by  the numbers of cells as defined above

#depth of first layer
z0 = 1000

#total model depth
model_depth = 100000

#stretching factor for the whole model extension
model_extension_factor = 1.

#starting resistivity value for homog. halfspace setup
rho0 = 100.

#name of datafile (to be handled as argument later on)
datafile = 'de.dat'

#name of output model file
modelfile = 'THE_modelfile.rho'

#==============================================================================
#==============================================================================

outstring = ''

outstring += '# ModEM model generated with MTpy - layout read from datafile: {0}\n'.format(datafile)

Fin = open(datafile,'r')
data = Fin.readlines()
Fin.close()

coords = []
line_8 = data[7].strip().split()


n_datapoints = float(line_8[1]) * float(line_8[2]) * 4

#read station coordinates
#start in line after header info, determined by starting character '>'
for dataline in data:
	if dataline.strip()[0] in ['#','>']:
		continue
	try:
		line = data[int(l)]
		line = line.strip().split()
		co = (float(line[4]),float(line[5]),float(line[6]))
		coords.append(co)
	except:
		continue

# local, Cartesian coordinates:
coords = np.array(list(set(coords)))

xmin = min(coords[:,0])
xmax = max(coords[:,0])
ymin = min(coords[:,1])
ymax = max(coords[:,1])

x_range = xmax - xmin
y_range = ymax - ymin

n_center_xblocks = int(x_range/dx) + 1
n_center_yblocks = int(y_range/dy) + 1

...

x0 = np.mean(coords[:,0])
y0 = np.mean(coords[:,1])

xmin_eff = xmin - x_range
xmax_eff = xmax + x_range
ymin_eff = ymin - y_range
ymax_eff = ymax + y_range

#--------------------------------------------
# calculate size of x blocks:
# at the moment just use even numbers of blocks:
half_n_xblocks = int(n_xblocks/2) - 1
#leave one block out for side padding
length2cover_x =  (xmax_eff - xmin_eff )/2.
minblock = 2*length2cover_x/(half_n_xblocks * (half_n_xblocks+1))

blockwidths_x = list(minblock * (np.arange(half_n_xblocks)+1))

#add the padding cell:
blockwidths_x.append(3*max(blockwidths_x))

xmin_padded = xmin_eff - blockwidths_x[-1]

xblocks = blockwidths_x[::-1]
xblocks.extend(blockwidths_x)

#--------------------------------------------
# calculate size of y blocks:
# at the moment just use even numbers of blocks:
half_n_yblocks = int(n_yblocks/2) - 1
#leave one block out for side padding
length2cover_y =  (ymax_eff - ymin_eff )/2.
minblock = 2*length2cover_y/(half_n_yblocks * (half_n_yblocks+1))

blockwidths_y = list(minblock * (np.arange(half_n_yblocks)+1))

#add the padding cell:
blockwidths_y.append(3*max(blockwidths_y))

ymin_padded = ymin_eff - blockwidths_y[-1]

yblocks = blockwidths_y[::-1]
yblocks.extend(blockwidths_y)

#--------------------------------------------

#build block depths:

n_layers_eff = n_layers - 2
one padding and one spilt uppermost layer

log_part_thickness = model_depth - (n_layers_eff-1) * z0
depths = np.logspace( np.log10(z0), np.log10(log_part_thickness), n_layers_eff ) + \
					 np.arange(n_layers_eff) * z0

depths = list(depths)
depths.append(3*max(depths))
thicknesses = [z0/2.]
for i, layer in enumerate(depths):
	if i == 0 :
		t = layer/2.
	else:
		t = layer - depths[i-1]
	thicknesses.append(t)
thicknesses.append(3*max(thicknesses))


outstring += '{0}    {1}    {2}    {3}    {4}\n'.format(len(xblocks),len(yblocks),
													len(thicknesses), 0,'LOGE')

xstring = ''
for block in xblocks:
	xstring += '{0:.3f}  '.format(block)
xstring += '\n'

outstring += xstring

ystring = ''
for block in yblocks:
	ystring += '{0:.3f}  '.format(block)
ystring += '\n'

outstring += ystring


zstring = ''
for block in thicknesses:
	zstring += '{0:.3f}  '.format(block)
zstring += '\n'

outstring += zstring


for idx_z in range(len(thicknesses)):
	z_string = ''
	#empty line before each layer:
	z_string += '\n'

	for idx_y in range(len(yblocks)):
		y_string = ''
		for idx_x in range(len(xblocks)):
			x_string = '{0:.5E}  '.format(np.log(rho0))
			y_string += x_string
		y_string += '\n'
		z_string += y_string
	outstring += z_string


co_reference = '{0}  {1}  {2} \n'.format(xmin_padded,ymin_padded,0)

outstring += co_reference

outstring += '0 \n'

Fout= open(modelfile,'w')
Fout.write(outstring)
Fout.close()



