#!/bin/env python
"""
Description:
    Class for plotting geological data (Poly and line data)
    in shapefiles.
References:
 
CreationDate:   11/24/17
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     11/24/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os

import numpy as np
import logging, traceback

from mtpy.utils import gis_tools
import matplotlib.pyplot as plt
from matplotlib import colors

from shapely.geometry import LineString, Polygon, MultiPolygon, shape
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch
from collections import defaultdict
import fiona

class Geology:
    def __init__(self, sfn,
                 symbolkey='SYMBOL',
                 
                 minLon=None, maxLon=None,
                 minLat=None, maxLat=None):
        '''
        Class for plotting geological (Poly/Line) data in shapefiles

        :param sfn: shape file name
        :param symbolkey: key in shapefile for map symbol
        :param minLon: minimum longitude of bounding box for clipping shapefile contents; this
                       is necessary to ensure that the map-legend contains only elements
                       that pertain to the region of interest
        :param maxLon: minimum longitude of bounding box
        :param minLat: minimum latitude of bounding box
        :param maxLat: maximum latitude of bounding box
        '''

        self._sfn = sfn
        self._properties = []
        self._geometries = []
        self._symbolkey = symbolkey
        self._hasLUT = False

        self._boundingPoly = None
        if(minLon != None and maxLon != None and minLat != None and maxLat != None):
            self._boundingPoly = Polygon([(minLon,minLat), (maxLon,minLat),
                                          (maxLon,maxLat), (minLon, maxLat)])

        # Load file
        sf = None
        try:
            sf = fiona.open(sfn)
        except Exception as err:
            print('Failed to read %s' % (sfn))
            logging.error(traceback.format_exc())
            exit(-1)

        for feature in sf:
            if feature['geometry'] is not None:
                g = shape(feature['geometry'])
    
                # filter geometry based on intersection with region
                # of interest
                if(self._boundingPoly != None):
                    if (self._boundingPoly.intersects(g)):
                        self._properties.append(feature['properties'])
                        self._geometries.append(g)
                else:
                    self._properties.append(feature['properties'])
                    self._geometries.append(g)
            # end for
        # end for
        sf.close()
    # end func

    def processLUT(self, lutfn, lut_delimiter=' '):
        '''
        Loads a colour lookup table of the following format:

        Header Line
        <unitname> r,g,b
        <unitname> r,g,b
        .
        .
        .

        The first line in the file is a header and subsequent lines list the geological
        unit name, followed by a rgb triplet.

        :param lutfn: Look-up-table file name
        '''

        self._lutfn = None
        self._lutDict = None
        self._hasLUT = False
        # Process look-up-table
        if (lutfn is not None):
            self._lutfn = lutfn
            self._lutDict = defaultdict(list)

            f = None
            try:
                f = open(lutfn)
            except Exception as err:
                print('Failed to read %s' % (lutfn))
                logging.error(traceback.format_exc())
                exit(-1)

            lines = f.readlines()
            for i, line in enumerate(lines):
                if (i == 0): continue  # skip header
                key = line.split(lut_delimiter)[0]
                color = np.array(line.strip().split(lut_delimiter)[1].split(','))
                try:
                    color = color.astype(np.float) / 256.
                except ValueError:
                    # if string value is provided (named matplotlib color), convert to a rgb tuple
                    color = np.array(colors.to_rgba(color[0])[:3])
                        

                self._lutDict[key] = color
            # end for
            self._hasLUT = True
        # end if
    # end func

    def plot(self, ax, m, lutfn=None, lut_delimiter=' ',
             default_polygon_color='grey', **kwargs):
        '''
        Plots a shapefile. This function assumes that a shapefile containing polygonal
        data will have an attribute column named 'SYMBOL', which is used to pick corresponding
        color values from the colour lookup table, described in function processLUT.

        :param ax: plot axis
        :param m: basemap instance
        :param lutfn: colour look-up-table file name, which if not provided, polygons
                      are coloured by the default colour (default_polygon_color). This
                      parameter is ignored for shapefiles that contain only line data
        :param default_polygon_color: default color for polygons; overridden by colors
                                      provided in look-up-table, if given
        :param kwargs: list of relevant matplotlib arguments, e.g. alpha, zorder, color, etc.
        :return:
            legend_handles: legend handles for polygonal data; empty list for line data
            legend_labels: symbol names for polygonal data; empty list for line data
        '''

        # Populate lookup table
        self.processLUT(lutfn,lut_delimiter=lut_delimiter)
        print("plotting geology")

        patches = []
        legend_handles = []
        legend_labels = []
        handles = set()
        

        ecolor_is_fcolor = False
        if ('edgecolor' in list(kwargs.keys()) and kwargs['edgecolor'] == 'face'):
            ecolor_is_fcolor = True
        # Process geometry
        
        for i, feature in enumerate(self._geometries):
            fcolor = None
            symbol = ''
            if (self._hasLUT):
                symbol = self._properties[i][self._symbolkey]
                fcolor = self._lutDict[symbol]
            if ((fcolor is None) or (np.iterable(fcolor) and len(fcolor) == 0)): 
                fcolor = default_polygon_color
            
            if (isinstance(feature, Polygon)):
                polygon = feature
                x, y = polygon.exterior.coords.xy
                if m is None:
                    px, py = x, y
                else:
                    px, py = m(x, y)
                # end if

                holes = []
                for interior in polygon.interiors:
                    x, y = interior.coords.xy
                    if m is None:
                        ipx, ipy = x, y
                    else:
                        ipx, ipy = m(x, y)
                    # end if

                    holes.append(list(zip(ipx, ipy)))
                # end for

                ppolygon = Polygon(shell=list(zip(px, py)), holes=holes)
                
                if (fcolor is not None): 
                    kwargs['facecolor'] = fcolor
                if ('edgecolor' not in list(kwargs.keys()) and not ecolor_is_fcolor):
                    kwargs['edgecolor'] = 'none'
                elif ecolor_is_fcolor:
                    kwargs['edgecolor'] = fcolor

                if ('fill') not in list(kwargs.keys()): kwargs['fill'] = True

                pp = PolygonPatch(ppolygon, **kwargs)
                patches.append(pp)

                # filter duplicates
                if (symbol not in handles):
                    handles.add(symbol)
                    legend_handles.append(pp)
                    legend_labels.append(symbol)

            elif (isinstance(feature, MultiPolygon)):
                multiPolygon = feature

                for polygon in multiPolygon:
                    x, y = polygon.exterior.coords.xy
                    if m is None:
                        px, py = x, y
                    else:
                        px, py = m(x, y)
                    # end if

                    holes = []
                    for interior in polygon.interiors:
                        x, y = interior.coords.xy
                        if m is None:
                            ipx, ipy = x, y
                        else:
                            ipx, ipy = m(x, y)
                        # end if

                        holes.append(list(zip(ipx, ipy)))
                    # end for

                    ppolygon = Polygon(shell=list(zip(px, py)), holes=holes)

                    if (fcolor is not None): kwargs['facecolor'] = fcolor
                    if ('edgecolor' not in list(kwargs.keys()) and not ecolor_is_fcolor):
                        kwargs['edgecolor'] = 'none'
                    elif ecolor_is_fcolor:
                        kwargs['edgecolor'] = fcolor
                    if ('fill') not in list(kwargs.keys()): kwargs['fill'] = True
    
                    pp = PolygonPatch(ppolygon, **kwargs)
                    patches.append(pp)

                # filter duplicates
                if (symbol not in handles):
                    handles.add(symbol)
                    legend_handles.append(pp)
                    legend_labels.append(symbol)
                # end for
            elif (isinstance(feature, LineString)):
                line = feature
                x, y = line.coords.xy
                if m is None:
                    px,py = x,y
                else:
                    px, py = m(x, y)
                ax.plot(px, py, **kwargs)
            # end if
        # end for
        if (len(patches)):
            ax.add_collection(PatchCollection(patches, match_original=True))

        return legend_handles, legend_labels
    # end func

    def _xy_to_local(self,x,y,epsg_from,epsg_to,centre_shift,scale_factor):
        '''
        
        '''
        xl,yl = gis_tools.epsg_project(x,y,epsg_from,epsg_to)
        xl = (np.array(xl) + centre_shift[0])/scale_factor
        yl = (np.array(yl) + centre_shift[1])/scale_factor
        
        return xl,yl
    # end func

    def plotlocal(self, epsg_from, epsg_to, centre_shift=[0.,0.], ax=None,
                  map_scale='m', default_polygon_color='grey', **kwargs):
        '''
        Plots a shapefile as lines in local coordinates (for overlaying on
        existing depth slice plotting functions, for example)
        :epsg_from: source epsg
        :epsg_to: target epsg
        :centre_shift: option to shift by [x,y]
        :ax: axes instance to plot on
        :map_scale: 'km' or 'm'
        :default_polygon_color: default color
        :kwargs: key word arguments to the matplotlib plot function
        local coordinates.
        :return:
        '''
        # set default line colour to black
        if 'color' not in list(kwargs.keys()):
            kwargs['color'] = 'k'
        
        if ax is None:
            ax = plt.subplot(111)
            
        if map_scale == 'km':
            scale_factor = 1000.
        else:
            scale_factor = 1.

        patches = []
        legend_handles = []
        legend_labels = []
        handles = set()

        ecolor_is_fcolor = False
        if ('edgecolor' in list(kwargs.keys()) and kwargs['edgecolor'] == 'face'):
            ecolor_is_fcolor = True
        # Process geometry
        for i, feature in enumerate(self._geometries):
            fcolor = None
            symbol = ''
            if (self._hasLUT):
                symbol = self._properties[i][self._symbolkey]
                fcolor = self._lutDict[symbol]
            if (fcolor == []): 
                fcolor = default_polygon_color

            if (isinstance(feature, Polygon)):
                polygon = feature
                x, y = polygon.exterior.coords.xy
                
                px, py = self._xy_to_local(x, y, epsg_from, epsg_to,
                                           centre_shift, scale_factor)
                ppolygon = Polygon(list(zip(px, py)))

                if (fcolor is not None): kwargs['facecolor'] = fcolor
                if ('edgecolor' not in list(kwargs.keys()) and not ecolor_is_fcolor):
                    kwargs['edgecolor'] = 'none'
                else:
                    kwargs['edgecolor'] = fcolor
                if ('fill') not in list(kwargs.keys()): kwargs['fill'] = True

                pp = PolygonPatch(ppolygon, **kwargs)
                patches.append(pp)

                # filter duplicates
                if (symbol not in handles):
                    handles.add(symbol)
                    legend_handles.append(pp)
                    legend_labels.append(symbol)

            elif (isinstance(feature, MultiPolygon)):
                multiPolygon = feature

                for polygon in multiPolygon:
                    x, y = polygon.exterior.coords.xy
                    px, py = self._xy_to_local(x, y, epsg_from, epsg_to,
                                               centre_shift, scale_factor)
                    ppolygon = Polygon(list(zip(px, py)))
                    
                    if (fcolor is not None): kwargs['facecolor'] = fcolor
                    if ('edgecolor' not in list(kwargs.keys()) and not ecolor_is_fcolor):
                        kwargs['edgecolor'] = 'none'
                    else:
                        kwargs['edgecolor'] = fcolor
                    if ('fill') not in list(kwargs.keys()): kwargs['fill'] = True

                    pp = PolygonPatch(ppolygon, **kwargs)
                    patches.append(pp)

                    # filter duplicates
                    if (symbol not in handles):
                        handles.add(symbol)
                        legend_handles.append(pp)
                        legend_labels.append(symbol)
                        # end for
            elif (isinstance(feature, LineString)):
                line = feature
                x, y = line.coords.xy
                px, py = self._xy_to_local(x, y, epsg_from, epsg_to,
                                           centre_shift, scale_factor)
                ax.plot(px, py, **kwargs)
                # end if
        # end for
        if (len(patches)):
            ax.add_collection(PatchCollection(patches, match_original=True))

        return ax, legend_handles, legend_labels
    # end func
# end class

def main():
    """
    define main function
    :return:
    """
    imaging = os.path.dirname(os.path.abspath(__file__))
    mtpy = os.path.dirname(imaging)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'geology')
    polyDataShapefile = os.path.join(ModEM_files, 'NT_LithInterp_2500K_region.shp')
    lineDataShapefile = os.path.join(ModEM_files, 'NT_Fault_2500K_polyline.shp')

    pg1 = Geology(polyDataShapefile) # Polygon data
    pg2 = Geology(lineDataShapefile) # Line data

    from mpl_toolkits.basemap import Basemap
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(16,12)
    m = Basemap(ax=ax, llcrnrlon=110,llcrnrlat=-50,urcrnrlon=155.,urcrnrlat=-10,
                 resolution='l', projection='merc', lat_0 = 39.5, lon_0 = 1)

    m.drawcoastlines()
    pg1.plot(ax, m)
    pg2.plot(ax, m)
    plt.savefig('/tmp/a.pdf')

    return
    # end


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()
