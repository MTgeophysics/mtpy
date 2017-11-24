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

from shapely.geometry import LineString, Polygon, MultiPolygon, shape
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch
from collections import defaultdict
import fiona

class plot_geology:
    def __init__(self, sfn):
        '''
        Class for plotting geological (Poly/Line) data in shapefiles

        :param sfn: shape file name
        '''
        self._sfn = sfn
        self._properties = []
        self._geometries = []
        self._hasLUT = False

        # Load file
        sf = fiona.open(sfn)
        for feature in sf:
            self._properties.append(feature['properties'])
            self._geometries.append(shape(feature['geometry']))
        # end for
        sf.close()
    # end func

    def processLUT(self, lutfn):
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

            f = open(lutfn)
            lines = f.readlines()
            for i, line in enumerate(lines):
                if (i == 0): continue  # skip header
                key = line.split(' ')[0]
                color = np.array(line.split(' ')[1].split(','))
                color = color.astype(np.float) / 256.

                self._lutDict[key] = color
            # end for
            self._hasLUT = True
        # end if
    # end func

    def plot(self, ax, m, lutfn=None,
             default_polygon_color='grey', **kwargs):
        '''
        Plot shapefile

        :param ax: plot axis
        :param m: basemap instance
        :param lutfn: colour look-up-table file name, which if not provided, polygons
                      are coloured by the default colour (default_polygon_color). This
                      parameter is ignored for shapefiles that contain only line data
        :param default_polygon_color: default color or polygons; overridden by colors
                                      provided in look-up-table, if given
        :param kwargs: ist of matplotlib relevant keywords, e.g. alpha, zorder, color, etc.
        :return:
        '''

        # Populate lookup table
        self.processLUT(lutfn)

        patches = []
        # Process geometry
        for i, feature in enumerate(self._geometries):

            fcolor = None
            if (self._hasLUT): fcolor = self._lutDict[self._properties[i]['SYMBOL']]
            if (fcolor == []): fcolor = default_polygon_color

            if (isinstance(feature, Polygon)):
                polygon = feature
                x, y = polygon.exterior.coords.xy
                px, py = m(x, y)
                ppolygon = Polygon(zip(px, py))

                if (fcolor is not None): kwargs['facecolor'] = fcolor
                if ('edgecolor' not in kwargs.keys()): kwargs['edgecolor'] = 'none'
                if ('fill') not in kwargs.keys(): kwargs['fill'] = True
                patches.append(PolygonPatch(ppolygon,
                                            **kwargs))
            elif (isinstance(feature, MultiPolygon)):
                multiPolygon = feature

                for polygon in multiPolygon:
                    x, y = polygon.exterior.coords.xy
                    px, py = m(x, y)
                    ppolygon = Polygon(zip(px, py))

                    if (fcolor is not None): kwargs['facecolor'] = fcolor
                    if ('edgecolor' not in kwargs.keys()): kwargs['edgecolor'] = 'none'
                    if ('fill') not in kwargs.keys(): kwargs['fill'] = True
                    patches.append(PolygonPatch(ppolygon,
                                                **kwargs))
                    # end for
            elif (isinstance(feature, LineString)):
                line = feature
                x, y = line.coords.xy
                px, py = m(x, y)
                ax.plot(px, py, **kwargs)
                # end if
        # end for
        if (len(patches)):
            ax.add_collection(PatchCollection(patches, match_original=True))
    # end func
# end class

def main():
    """
    define main function
    :return:
    """
    imaging2 = os.path.dirname(__file__)
    mtpy = os.path.dirname(imaging2)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'geology')
    polyDataShapefile = os.path.join(ModEM_files, 'NT_LithInterp_2500K_region.shp')
    lineDataShapefile = os.path.join(ModEM_files, 'NT_Fault_2500K_polyline.shp')

    pg1 = plot_geology(polyDataShapefile)
    pg2 = plot_geology(lineDataShapefile)

    return
    # end


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()