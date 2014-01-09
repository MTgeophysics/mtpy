#!/usr/bin/env python

#build ModEM input Model from ModEM data file

import numpy as np
import sys
#==============================================================================

# paramters:

n_xpadding = 5
n_ypadding = 6

#number of vertical padding layers is set to 3 !
#factor with which the padding stretches outside the central rectangle grid
padding_stretch = 1.5

n_layers = 15

#determine minimum block sizes
#used in the inner rectangle - constant widths
dx = 470
dy = 500
#region around stations discretised with these sizes
#outside, the grid steps will be extended exponentially
#the size of padding is determined by  the numbers of cells as defined above

#number of trys to shift the grid for getting own cells for each station
n_maximum_gridshifts = 70

#depth of first layer
z0 = 500

#total model depth in meters
model_depth = 100000

#stretching factor for the whole model extension
model_extension_factor = 1.

#starting resistivity value for homog. halfspace setup
rho0 = 100.

#name of datafile (to be handled as argument later on)
datafile = 'Modular_NLCG_002.data'

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
#line_8 = data[7].strip().split()


#n_datapoints = int(float(line_8[1]) * float(line_8[2]) * 4)

#read station coordinates
#start in line after header info, determined by starting character '>'
for dataline in data:
	line = dataline.strip().split()
	if (len(line) == 0) or line[0].strip()[0] in ['#','>']:
		continue
	try:
		line = dataline.strip().split()
		co = (float(line[4]),float(line[5]),float(line[6]))
		coords.append(co)
	except:
		continue

# local, Cartesian coordinates:
coords = np.array(list(set(coords)))

#reduce grid to 2D - assuming all stations are at the surface
xmin = min(coords[:,0])
xmax = max(coords[:,0])
ymin = min(coords[:,1])
ymax = max(coords[:,1])

x_range = xmax - xmin
y_range = ymax - ymin

n_center_xblocks = int(x_range/dx) + 3
n_center_yblocks = int(y_range/dy) + 3

center_widthX = n_center_xblocks * dx
center_widthY = n_center_yblocks * dy

surplusX = center_widthX - x_range
surplusY = center_widthY - y_range

all_points_in_single_cell = False
n_shifts = 0 
x_shifts = 0 
y_shifts = 0 


while all_points_in_single_cell is False:
	if n_shifts > n_maximum_gridshifts:
		break

	shifting_fraction = np.sqrt(n_maximum_gridshifts) + 1

	offset_x = x_shifts * dx/shifting_fraction
	offset_y = y_shifts * dy/shifting_fraction

	if n_shifts >0:
		print '{0} shift(s): x-offset {1} m - y-offset {2} m'.format(n_shifts,offset_x,offset_y)


	center_x0 = xmin - surplusX/2. + offset_x
	center_y0 = ymin - surplusY/2. + offset_y

	grid_x_points = (np.arange(n_center_xblocks+1) * dx) + center_x0
	grid_y_points = (np.arange(n_center_yblocks+1) * dy) + center_y0


	station_cells = []

	for idx_sta,co in enumerate(coords):
		idx_x = np.argmin(np.abs(grid_x_points-co[0]))
		if (grid_x_points-co[0])[idx_x] == 0:
			# coordinate lies on a node line => need to shift
			print 'station coordinates lie on cell nodes'
			break
		#otherwise, shift the index to correspond with the row of blocks, if necessary:
		if grid_x_points[idx_x] > co[0] : 
			idx_x -= 1
		idx_y = np.argmin(np.abs(grid_y_points-co[1]))
		if (grid_y_points-co[1])[idx_y] == 0:
			# coordinate lies on a node line => need to shift
			break
		#otherwise, shift the index to correspond with the row of blocks, if necessary:
		if grid_y_points[idx_y] > co[1] : 
			idx_y -= 1

		#cells enumerated West->East first, then northwards
		cell_index = idx_x * n_center_xblocks + idx_y
		station_cells.append(cell_index)

	if len(set(station_cells)) == len(coords):
		all_points_in_single_cell = True


	#shift the grid 
	x_shifts += 1
	if x_shifts >= (shifting_fraction - 1):
		x_shifts = 0
		y_shifts += 1

	n_shifts += 1
	


if all_points_in_single_cell < 1:
	print 'ERROR - cannot build grid having each station in a single cell!\n'\
	'change the values for dx,dy or remove stations'
	sys.exit()


#Now the inner grid is  well distributed over the stations
#add padding to the sides:

grid_x_points = list(grid_x_points)
x_padding_widths = [dx]
for idx_pad in range(n_xpadding):
	pad = x_padding_widths[-1] * padding_stretch
	grid_x_points.insert(0,grid_x_points[0]-pad)
	grid_x_points.append(grid_x_points[-1]+pad)
	x_padding_widths.append(pad)

grid_y_points = list(grid_y_points)
y_padding_widths = [dy]
for idy_pad in range(n_ypadding):
	pad = y_padding_widths[-1] * padding_stretch
	grid_y_points.insert(0,grid_y_points[0]-pad)
	grid_y_points.append(grid_y_points[-1]+pad)
	y_padding_widths.append(pad)


xmin_padded = grid_x_points[0]
ymin_padded = grid_y_points[0]

# transfer the block coordinates into block widths
xblocks = []
for idx_x in range(len(grid_x_points)-1):
	xblocks.append(grid_x_points[idx_x+1] - grid_x_points[idx_x])
yblocks = []
for idy_y in range(len(grid_y_points)-1):
	yblocks.append(grid_y_points[idy_y+1] - grid_y_points[idy_y])

#---------------------------------------------------------------------
n_zpadding = 3

#build block depths:

n_layers_eff = n_layers - 1
#splitted uppermost layer

log_part_thickness = model_depth - (n_layers_eff-1) * z0
depths = np.logspace( np.log10(z0), np.log10(log_part_thickness), n_layers_eff ) + \
					 np.arange(n_layers_eff) * z0


depths = list(depths)

thicknesses = [z0/2.]
for i, layer in enumerate(depths):
	if i == 0 :
		t = layer/2.
	else:
		t = layer - depths[i-1]
	thicknesses.append(t)
print thicknesses

padding = [thicknesses[-1]*padding_stretch]
for idx_pad in range(n_zpadding-1):
	padding.append(padding[-1]*padding_stretch)
total_padding = np.sum(padding)
pad_ratio = total_padding/model_depth

if pad_ratio < 1:
	padding = list(np.array(padding)/pad_ratio)
if pad_ratio >2 :
	padding = list(np.array(padding)/pad_ratio*2)
print pad_ratio ,total_padding,model_depth
print padding, np.sum(padding)


thicknesses.extend(padding)


grid_z_points = [0]
for t in thicknesses:
	grid_z_points.append(grid_z_points[-1]+t)


#some information for the user:
print '\n\t Model set up - dimensions: {0:.1f}x{1:.1f}x{2:.1f} km^3 ({3}x{4}x{5} cells)\n'.format(
	(grid_x_points[-1]-grid_x_points[0])/1000.,(grid_y_points[-1]-grid_y_points[0])/1000.,
	depths[-1]/1000.,len(grid_x_points)-1,len(grid_y_points)-1,len(grid_z_points)-1)


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


def plotgrid(stations,grid_x,grid_y,grid_z=None, n_xpadding = None, n_y_padding=None, n_zpadding_layers = None):
	from pylab import *
	ion()
	close('all')
	

	grid_x = [i/1000. for i in grid_x]
	grid_y = [i/1000. for i in grid_y]

	# Note: X and Y are swapped - mathematical definition used in the plotting functions!!!
	#fig = figure(1)
	#ax = fig.gca()
	fig = figure(figsize=(8, 6))
	if 	grid_z is not None:
		ax = subplot2grid((1, 3), (0, 0), colspan=2,aspect='equal')
	else:
		ax = subplot2grid((1, 3), (0, 0), colspan=3,aspect='equal')
	

	#ax = subplot(1,2,1)
	ax.scatter(stations[:,1]/1000.,stations[:,0]/1000.,c='r')
	ax.scatter([ymin_padded/1000.],[xmin_padded/1000.],c='b',marker='x',s=40)
	outline_x = [min(grid_x),min(grid_x),max(grid_x),max(grid_x),min(grid_x)]
	outline_y = [min(grid_y),max(grid_y),max(grid_y),min(grid_y),min(grid_y)]
	ax.plot(outline_y,outline_x,c='r')
	

	if n_xpadding is not None and n_ypadding is not None: 
		regular_x = [grid_x[n_xpadding],grid_x[n_xpadding],
					grid_x[-n_xpadding],grid_x[-n_xpadding],grid_x[n_xpadding]]
		regular_y = [grid_y[n_ypadding],grid_y[-n_ypadding],
					grid_y[-n_ypadding],grid_y[n_ypadding],grid_y[n_ypadding]]
		ax.plot(regular_y,regular_x,c='b')



	extension_factor = 0.1
	x_extent = max(grid_x) - min(grid_x)
	x_extension = extension_factor * x_extent
	ax.set_ylim([min(grid_x) - x_extension,max(grid_x) + x_extension])
	
	y_extent = max(grid_y) - min(grid_y)
	y_extension = extension_factor * y_extent
	ax.set_xlim([min(grid_y) - y_extension,max(grid_y) + y_extension])
	

	ax.set_yticks(grid_x, minor=True)
	ax.yaxis.grid(False, which='major')
	ax.yaxis.grid(True, which='minor',c='g')
	ax.set_xticks(grid_y, minor=True)
	ax.xaxis.grid(False, which='major')
	ax.xaxis.grid(True, which='minor',c='g')
	ax.set_xlabel('Easting (Y-coordinate) in km')
	ax.set_ylabel('Northing (X-coordinate) in km')
	ax.set_title('Model geometry (origin at {0:.1f},{1:.1f})'.format(xmin_padded,ymin_padded))
	ax.set_aspect('equal',adjustable='box')
	draw()

	if grid_z is not None:

		grid_z = [-i/1000. for i in grid_z]
		bottom_index = len(grid_z) - n_zpadding_layers -1
		
		ax2 = subplot2grid((1, 3), (0, 2),aspect='equal')
		#fig2 = figure(2)
		#ax2 = fig2.gca()
		#ax2 = subplot(1,2,2)
		outline_z = [min(grid_z),min(grid_z),max(grid_z),max(grid_z),min(grid_z)]
		outline_y = [min(grid_y),max(grid_y),max(grid_y),min(grid_y),min(grid_y)]
		plot(outline_y,outline_z,c='r')

		plot([min(grid_y),max(grid_y)],[grid_z[bottom_index],grid_z[bottom_index]],c='b')

		ax2.axhline(linewidth=2, color='k')

		extension_factor = 0.1

		z_extent = max(grid_z) - min(grid_z)
		z_extension = extension_factor * z_extent
		ax2.set_ylim([min(grid_z) - z_extension,max(grid_z) + z_extension])
		
		y_extent = max(grid_y) - min(grid_y)
		y_extension = extension_factor * y_extent
		ax2.set_xlim([min(grid_y) - y_extension,max(grid_y) + y_extension])
		#ax2.set_aspect('equal','datalim')
		ax2.set_yticks(grid_z, minor=True)
		ax2.yaxis.grid(False, which='major')
		ax2.yaxis.grid(True, which='minor',c='k')
		ax2.set_xlabel('Easting (Y-coordinate) in km')
		ax2.set_ylabel('Depth in km')
		ax2.set_title('Model layers')

		ax2.set_aspect('equal',adjustable='box')

	#tight_layout()
	show()
	raw_input()

plotgrid(coords,grid_x_points,grid_y_points,grid_z_points,n_xpadding,n_ypadding, n_zpadding)
