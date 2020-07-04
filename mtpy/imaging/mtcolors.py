# -*- coding: utf-8 -*-
"""
Created on Tue May 14 18:05:59 2013

@author: jpeacock-pr
"""

import matplotlib.colors as colors
from matplotlib import cm
import numpy as np

#==============================================================================
# Make some color maps for plotting
#==============================================================================
#yellow to red
ptcmapdict = {'red':((0.0, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),

              'green':((0.0, 0.0, 1.0),
                       (1.0, 0.0, 1.0)),

              'blue':((0.0, 0.0, 0.0),
                      (1.0, 0.0, 0.0))}

mt_yl2rd=colors.LinearSegmentedColormap('mt_yl2rd', ptcmapdict, 256)

#blue to yellow to red
skcmapdict = {'red':((0.0, 0.0, 0.0),
                     (.5, 1.0, 1.0),
                     (0.5, 0.0, 1.0),
                     (1.0, 1.0, 1.0)),
              'green':((0.0, 1.0, 0.0),
                       (.5, 1.0, 0.0),
                       (.5, 0.0, 1.0),
                       (1.0, 0.0, 1.0)),
              'blue':((0.0, 0.0, 1.0),
                      (.5, 0.0, 1.0),
                      (0.5, 0.1, 0.1),
                      (1.0, 0.1, 0.1))}

mt_bl2yl2rd=colors.LinearSegmentedColormap('mt_bl2yl2rd', skcmapdict, 256)

#blue to white to red
skcmapdict2 = {'red':  ((0.0, 0.0, 0.0),
                        (0.25,0.0, 0.0),
                        (0.5, 0.8, 1.0),
                        (0.75,1.0, 1.0),
                        (1.0, 0.4, 1.0)),

               'green': ((0.0, 0.0, 0.0),
                         (0.25,0.0, 0.0),
                         (0.5, 0.9, 0.9),
                         (0.75,0.0, 0.0),
                         (1.0, 0.0, 0.0)),

               'blue':  ((0.0, 0.0, 0.4),
                         (0.25,1.0, 1.0),
                         (0.5, 1.0, 0.8),
                         (0.75,0.0, 0.0),
                         (1.0, 0.0, 0.0))}

mt_bl2wh2rd=colors.LinearSegmentedColormap('mt_bl2wh2rd', skcmapdict2, 256)


#blue to white to red in segmented colors
mt_seg_bl2wh2rd = colors.ListedColormap(((0, 0, 1), (.5, .5, 1), (.75, .75, 1),
                                         (.9, .9, 1), (1, 1, 1), (1.0, .9, .9),
                                         (1, .75, .75), (1, .5, .5),(1, 0, 0)))

#white to blue
ptcmapdict3 = {'red':((0.0, 1.0, 1.0),
                      (1.0, 0.0, 0.0)),

               'green':((0.0, 1.0, 1.0),
                        (1.0, 0.0, 0.0)),

               'blue':((0.0, 1.0, 1.0),
                       (1.0, 1.0, 1.0))}
mt_wh2bl = colors.LinearSegmentedColormap('mt_wh2bl', ptcmapdict3, 256)

#white to orange
cmapdict_wh2or = {'red':((0.0, 1.0, 1.0),
                         (1.0, .95, 0.0)),

                  'green':((0.0, 1.0, 1.0),
                           (1.0, .45, .95)),

                  'blue':((0.0, 1.0, 1.0),
                          (1.0, 0, 0))}
mt_wh2or = colors.LinearSegmentedColormap('mt_wh2or', cmapdict_wh2or, 256)

#red to blue
rtcmapdict = {'red':((0.0, 0.0, 1.0),
                     (1.0, 0.0, 1.0)),

              'green':((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

              'blue':((0.0, 1.0, 0.0),
                      (1.0, 1.0, 0.0))}
mt_rd2bl = colors.LinearSegmentedColormap('mt_rd2bl', rtcmapdict, 256)

#blue to green to red
ptcmapdict4 = {'red':  ((0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.5, 0.9, 1.0),
                        (0.75,1.0, 1.0),
                        (1.0, 0.45, 1.0)),

               'green': ((0.0, 0.0, 0.0),
                         (0.25, 0.5, 0.5),
                         (0.5, 1.0, 1.0),
                         (0.75,0.5, 0.5),
                         (1.0, 0.0, 0.0)),

              'blue':  ((0.0, 0.0, 0.45),
                        (0.25, 1.0, 1.0),
                        (0.5, 1.0, 0.9),
                        (0.75,0.0, 0.0),
                        (1.0, 0.0, 0.0))}
mt_bl2gr2rd = colors.LinearSegmentedColormap('mt_bl2gr2rd', ptcmapdict4, 256)

#red to green to blue
ptcmapdict4 = {'red':  ((0.0, 0.0, 0.45),
                        (0.25, 1.0, 1.0),
                        (0.5, 1.0, 0.9),
                        (0.75,0.0, 0.0),
                        (1.0, 0.0, 0.0)),

               'green': ((0.0, 0.0, 0.0),
                         (0.25, 0.5, 0.5),
                         (0.5, 1.0, 1.0),
                         (0.75,0.5, 0.5),
                         (1.0, 0.0, 0.0)),

              'blue':  ((0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.5, 0.9, 1.0),
                        (0.75,1.0, 1.0),
                        (1.0, 0.45, 1.0))}
mt_rd2gr2bl = colors.LinearSegmentedColormap('mt_rd2gr2bl', ptcmapdict4, 256)

#mtcmapdict2 = {'red':  ((0.0, 0.5, 0.2),
#                        (0.25, 1.0, 1.0),
#                        (0.50, 1.0, 1.0),
#                        (0.55, 0.9, 0.9),
#                        (1.0, 0.0, 0.0),
#                        (1.0, 0.2, 0.2)),
#
#               'green': ((0.0, 0.0, 0.0),
#                         (0.25, 0.3, 0.3),
#                         (0.50, 1.0, 1.0),
#                         (0.55, 0.9, 0.9),
#                         (1.0, 0.0, 0.0),
#                         (1.0, 0.3, 0.3)),
#
#               'blue':  ((0.0, 0.0, 0.0),
#                         (0.25, 0.0, 0.0),
#                         (0.50, 0.9, 0.9),
#                         (0.55, 1.0, 1.0),
#                         (1.0, 0.50, 0.50),
#                         (1.0, 0.4, 0.75))}
mtcmapdict2 = {'red':  ((0.0, 0.4, 0.2),
                        (0.20, 0.992, 0.992),
                        (0.50, 0.953, 0.953),
                        (0.62, 0.384, 0.384),
                        (1.0, 0.12, 0.12)),

               'green': ((0.0, 0.0392, 0.0392),
                         (0.20, 0.3098, 0.3098),
                         (0.50, 0.953, 0.953),
                         (0.62, 0.529, 0.529),
                         (1.0, 0.094, 0.094)),

               'blue':  ((0.0, 0.00, 0.00),
                         (0.20, 0.0, 0.0),
                         (0.50, 0.953, 0.953),
                         (0.62, 1.0, 1.0),
                         (1.0, 0.45, 0.45))}
mt_rd2wh2bl = colors.LinearSegmentedColormap('mt_rd2wh2bl', mtcmapdict2, 256)

mtcmapdict3 = {'red':  ((0.0, 0.2, 0.2),
                        (0.25, 0.2, 0.2),
                        (0.5, 1.0, 1.0),
                        (0.75, 1.0, 1.0),
                        (1.0, 0.2, 0.5)),

               'green': ((0.0, 0.3, 0.3),
                         (0.25, 0.3, 0.3),
                         (0.5, 1.0, 1.0),
                         (0.75, 0.3, 0.3),
                         (1.0, 0.0, 0.0)),

               'blue':  ((0.0, 0.75, 0.45),
                         (0.25, 0.85, 0.85),
                         (0.5, 1.0, 1.0),
                         (0.75, 0.0, 0.0),
                         (1.0, 0.0, 0.0))}

mt_rd2wh2bl_r = colors.LinearSegmentedColormap('mt_rd2wh2bl_r', mtcmapdict3, 256)

rdylbu_data = {'blue': [[0.0, 0.000, 0.000],
                        [0.1, 0.000, 0.000],
                        [0.2, 0.000, 0.000],
                        [0.3, 0.000, 0.000],
                        [0.4, 0.000, 0.000],
                        [0.5, 1.000, 1.000],
                        [0.6, 0.990, 0.990],
                        [0.7, 0.950, 0.950],
                        [0.8, 0.900, 0.900],
                        [0.9, 0.550, 0.550],
                        [1.0, 0.250, 0.250]],
               'green': [[0.0, 0.000, 0.000],
                         [0.1, 0.100, 0.100],
                         [0.2, 0.400, 0.400],
                         [0.3, 0.800, 0.800],
                         [0.4, 0.900, 0.900],
                         [0.5, 1.000, 1.000],
                         [0.6, 0.900, 0.900],
                         [0.7, 0.800, 0.800],
                         [0.8, 0.400, 0.400],
                         [0.9, 0.100, 0.100],
                         [1.0, 0.000, 0.000]],
                'red':[[0.0, 0.250, 0.250],
                       [0.1, 0.550, 0.550],
                       [0.2, 0.900, 0.900],
                       [0.3, 0.950, 0.950],
                       [0.4, 0.990, 0.990],
                       [0.5, 1.000, 1.000],
                       [0.6, 0.000, 0.000],
                       [0.7, 0.000, 0.000],
                       [0.8, 0.000, 0.000],
                       [0.9, 0.000, 0.000],
                       [1.0, 0.000, 0.000]]}

mt_rdylbu = colors.LinearSegmentedColormap('mt_rdylbu', rdylbu_data, 256)


cmapdict = {'mt_yl2rd' : mt_yl2rd,
            'mt_bl2yl2rd' : mt_bl2yl2rd,
            'mt_wh2bl' : mt_wh2bl,
            'mt_rd2bl' : mt_rd2bl,
            'mt_bl2wh2rd' : mt_bl2wh2rd,
            'mt_seg_bl2wh2rd' : mt_seg_bl2wh2rd,
            'mt_bl2gr2rd' : mt_bl2gr2rd,
            'mt_rd2gr2bl' : mt_rd2gr2bl,
            'mt_wh2or' : mt_wh2or,
            'mt_rd2wh2bl': mt_rd2wh2bl,
            'mt_rd2wh2bl_r': mt_rd2wh2bl_r,
            'mt_rdylbu': mt_rdylbu}
# add matplotlib built-in colormaps
cmapdict.update(cm.cmap_d)

#make functions for getting the color from each map according to the variable
#cvar

def get_color(cvar,cmap):
    """
    gets the color to plot for the given color map
    
    """
    if cmap == 'mt_yl2rd':
        plot_color = get_mt_yl2rd(cvar)
        return plot_color

    elif cmap == 'mt_wh2bl':
        plot_color = get_mt_wh2bl(cvar)
        return plot_color

    elif cmap == 'mt_bl2wh2rd' or cmap=='mt_seg_bl2wh2rd':
        plot_color = get_mt_bl2wh2rd(cvar)
        return plot_color

    elif cmap == 'mt_bl2yl2rd':
        plot_color = get_mt_bl2yl2rd(cvar)
        return plot_color

    elif cmap == 'mt_bl2gr2rd':
        plot_color = get_mt_bl2gr2rd(cvar)
        return plot_color

    elif cmap == 'mt_rd2gr2bl':
        plot_color = get_mt_rd2gr2bl(cvar)
        return plot_color

    elif cmap == 'mt_wh2or':
        plot_color = get_mt_wh2or(cvar)
        return plot_color

    elif cmap == 'mt_rd2wh2bl':
        plot_color = get_mt_rd2wh2bl(cvar)
        return plot_color

    elif cmap == 'mt_rd2wh2bl_r':
        plot_color = get_mt_rd2wh2bl_r(cvar)
        return plot_color
    
    else:
        try:
            return get_matplotlib_cval(cmap,cvar)

        except:
            print('Color map: {0} is not supported yet.'.format(cmap))


def get_matplotlib_cval(cmap,cvar):
    """
    gets the color for any matplotlib colormaps
    
    """
    return cm.get_cmap(cmap)(cvar)
    

def get_mt_yl2rd(cvar):
    """
    gets color for the color map that goes from yellow to red
    
    """

    if cvar >= 1:
        plot_color = (1, 0, 0)
    elif cvar <= 0:
        plot_color = (1, 1, 0)
    else:
        plot_color = (1, 1-abs(cvar), 0.1)

    return plot_color

def get_mt_wh2bl(cvar):
    """
    gets color for the color map that goes from white to blue
    
    """

    if cvar >= 1:
        plot_color = (0, 0, 1)
    elif cvar <= 0:
        plot_color = (1, 1, 1)
    else:
        plot_color = (1-abs(cvar), 1-abs(cvar), 1)

    return plot_color

def get_mt_wh2or(cvar):
    """
    gets color for the color map that goes from white to orange
    
    """

    if cvar >= 1:
        plot_color = (1, .5, 0)
    elif cvar <= 0:
        plot_color = (1, 1, 1)
    else:
        plot_color = (1, abs(cvar)*.5+.5, abs(cvar))

    return plot_color

def get_mt_bl2wh2rd(cvar):
    """
    gets color for the color map that goes from blue to white to red
    
    """

    if cvar < 0 and cvar > -1:
        plot_color = (1+cvar, 1+cvar, 1)
    elif cvar <= -1:
        plot_color = (0, 0, 1)
    elif cvar >= 0 and cvar < 1:
        plot_color = (1, 1-cvar, 1-cvar)
    elif cvar >= 1:
        plot_color = (1, 0, 0)

    return plot_color

def get_mt_bl2yl2rd(cvar):
    """
    gets color for the color map that goes from blue to yellow to red
    
    """

    if cvar < 0 and cvar > -1:
        plot_color = (1+cvar, 1+cvar, -cvar)
    elif cvar <= -1:
        plot_color = (0, 0, 1)
    elif cvar >= 0 and cvar < 1:
        plot_color = (1, 1-cvar, .01)
    elif cvar >= 1:
        plot_color = (1, 0, 0)

    return plot_color

def get_mt_bl2gr2rd(cvar):
    """
    gets color for the color map that goes from blue to greenish to red
    
    """

    if cvar < 0 and cvar > -1:
        plot_color = (1+cvar, 1+cvar/2, 1)
    elif cvar <= -1:
        plot_color = (0, 0, 1)
    elif cvar >= 0 and cvar < 1:
        plot_color = (1, 1-cvar/2, 1-cvar)
    elif cvar >= 1:
        plot_color = (1, 0, 0)

    return plot_color

def get_mt_rd2gr2bl(cvar):
    """
    gets color for the color map that goes red to greenish to blue
    
    """

    if cvar < 0 and cvar > -1:
        plot_color = (1, 1+cvar/2, 1+cvar)
    elif cvar <= -1:
        plot_color = (1, 0, 0)
    elif cvar >= 0 and cvar < 1:
        plot_color = (1-cvar, 1-cvar/2, 1)
    elif cvar >= 1:
        plot_color = (0, 0, 1)

    return plot_color

def get_mt_rd2wh2bl(cvar):
    """
    gets color for the color map that goes red to white to blue
    
    """

    if cvar < 0 and cvar > -1:
        plot_color = (1, 1+cvar/3, 1+cvar)
    elif cvar <= -1:
        plot_color = (1, 0, 0)
    elif cvar >= 0 and cvar < 1:
        plot_color = (1-cvar, 1-cvar/3, 1)
    elif cvar >= 1:
        plot_color = (0, 0, 1)

    return plot_color

def get_mt_rd2wh2bl_r(cvar):
    """
    gets color for the color map that goes red to white to blue
    
    """
    # blue
    if cvar < -.5 and cvar > -1:
        plot_color = (.2, .3, 1.5+cvar)
    if cvar < 0 and cvar > -.5:
        plot_color = (1.6*cvar+1, 1.4*cvar+1, .5*cvar+1)
    elif cvar <= -1:
        plot_color = (.2, .3, .5)
    elif cvar >= 0 and cvar < .5:
        plot_color = (1, -1.2*cvar+1, -2*cvar+1)
    # red
    elif cvar >= .5 and cvar < 1:
        plot_color = (-cvar+1.5, -.6*cvar+.6, 0)
    elif cvar >= 1:
        plot_color = (.5, 0, 0)

    return plot_color



def get_plot_color(colorx, comp, cmap, ckmin=None, ckmax=None, bounds=None):
    """
    gets the color for the given compnent, color array and cmap

    Note: we now use the linearSegmentedColorMap objects, instead of the get_color function
    """

    #get face color info
    if comp in ['phimin', 'phimax', 'phidet', 'ellipticity', 'geometric_mean',
                'azimuth', 'strike']:
        if ckmin is None or ckmax is None:
            raise IOError('Need to input min and max values for plotting')

        '''
        cvar = (colorx-ckmin)/(ckmax-ckmin)
        if cmap == 'mt_bl2wh2rd' or cmap == 'mt_bl2yl2rd' or \
           cmap == 'mt_bl2gr2rd' or cmap == 'mt_rd2gr2bl' or \
           cmap == 'mt_rd2wh2bl' or cmap == 'mt_rd2wh2bl_r':
            cvar = 2*cvar-1

        return get_color(cvar, cmap)
        '''
        norm = colors.Normalize(ckmin, ckmax)
        if(cmap in list(cmapdict.keys())):
            return cmapdict[cmap](norm(colorx))
        else:
            return cm.get_cmap(cmap)(norm(colorx))
    elif comp == 'skew' or comp == 'normalized_skew':
        '''
        cvar = 2*colorx/(ckmax-ckmin)
        return get_color(cvar, cmap)
        '''

        norm = colors.Normalize(ckmin, ckmax)
        if (cmap in list(cmapdict.keys())):
            return cmapdict[cmap](norm(colorx))
        else:
            return cm.get_cmap(cmap)(norm(colorx))
    
    elif comp == 'skew_seg' or comp == 'normalized_skew_seg':
        if bounds is None:
            raise IOError('Need to input bounds for segmented colormap')

        '''
        for bb in range(bounds.shape[0]):
            if colorx >= bounds[bb] and colorx < bounds[bb+1]:
                cvar = float(bounds[bb])/bounds.max()
                return get_color(cvar, cmap)

            #if the skew is extremely negative make it blue
            elif colorx < bounds[0]:
                cvar = -1.0
                return get_color(cvar, cmap)

            #if skew is extremely positive make it red
            elif colorx > bounds[-1]:
                cvar = 1.0
                return get_color(cvar, cmap)
        '''
        norm = colors.Normalize(bounds[0], bounds[-1])
        step = abs(bounds[1] - bounds[0])
        ### need to get the color into a bin so as to not smear the colors.
        
        if colorx > max(bounds):
            colorx = max(bounds)
        elif colorx < min(bounds):
            colorx = min(bounds)
        elif abs(colorx) <= step:
            colorx = 0
        else:
            colorx = int(step * round(float(colorx - np.sign(colorx) * (abs(colorx) % step))/ step))

        if (cmap in list(cmapdict.keys())):
            return cmapdict[cmap](norm(colorx))
        else:
            return cm.get_cmap(cmap)(norm(colorx))
    else:
        raise NameError('color key '+comp+' not supported')

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
      
         cmap: colormap instance, eg. cm.jet. 
         N: number of colors.
     
     Example
         x = resize(arange(100), (5,100))
         djet = cmap_discretize(cm.jet, 5)
         imshow(x, cmap=djet)
    """

    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in range(N+1)]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)




