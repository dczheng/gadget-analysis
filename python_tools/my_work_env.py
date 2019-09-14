#!/usr/bin/env python3

import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib import cm
import matplotlib.colors as mplc
import cos_and_const as mycc
import matplotlib.ticker as tik

from matplotlib.font_manager import FontProperties
#https://matplotlib.org/2.0.0/examples/pylab_examples/fonts_demo.html

plt.style.use( 'ggplot' )

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

#figs_dir = "./"
#data_dir = "./data/"
figs_dir = "/mnt/ddnfs/data_users/dczheng/simulation/plots/figs/"
data_dir = "/mnt/ddnfs/data_users/dczheng/simulation/plots/data/"

ds9a = {'red': lambda v : np.interp(v, [0, 0.25, 0.5, 1],
                                        [0, 0, 1, 1]),
         'green': lambda v : np.interp(v, [0, 0.25, 0.5, 0.77, 1],
                                          [0, 1, 0, 0, 1]),
         'blue': lambda v : np.interp(v, [0, 0.125, 0.5, 0.64, 0.77, 1],
                                         [0, 0, 1, 0.5, 0, 0])}
ds9b = {'red': lambda v : 4 * v - 1,
        'green': lambda v : 4 * v - 2,
        'blue': lambda v : np.select([v < 0.25, v < 0.5, v < 0.75, v <= 1],
                                      [4 * v, -4 * v + 2, 0, 4 * v - 3])}

# Note that this definition slightly differs from ds9cool, but make more sense to me...
ds9cool = {'red': lambda v : 2 * v - 1,
           'green': lambda v : 2 * v - 0.5,
           'blue': lambda v : 2 * v}


ds9i8 = {'red': lambda v : np.where(v < 0.5, 0, 1),
        'green': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                       v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                      [0, 1, 0, 1, 0, 1, 0, 1]),
        'blue': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                      v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                      [0, 0, 1, 1, 0, 0, 1, 1])}

ds9aips0 = {'red': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0.475, 0, 0.373, 0, 0, 1, 1, 1]),
            'green': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0, 0, 0.655, 0.596, 0.965, 1, 0.694, 0]),
            'blue': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0.608, 0.785, 0.925, 0, 0, 0, 0, 0])}

ds9rainbow = {'red': lambda v : np.interp(v, [0, 0.2, 0.6, 0.8, 1], [1, 0, 0, 1, 1]),
              'green': lambda v : np.interp(v, [0, 0.2, 0.4, 0.8, 1], [0, 0, 1, 1, 0]),
              'blue': lambda v : np.interp(v, [0, 0.4, 0.6, 1], [1, 1, 0, 0])}

# This definition seems a bit strange...
ds9he = {'red': lambda v : np.interp(v, [0, 0.015, 0.25, 0.5, 1],
                                        [0, 0.5, 0.5, 0.75, 1]),
         'green': lambda v : np.interp(v, [0, 0.065, 0.125, 0.25, 0.5, 1],
                                          [0, 0, 0.5, 0.75, 0.81, 1]),
         'blue': lambda v : np.interp(v, [0, 0.015, 0.03, 0.065, 0.25, 1],
                                         [0, 0.125, 0.375, 0.625, 0.25, 1])}

ds9heat = {'red': lambda v : np.interp(v, [0, 0.34, 1], [0, 1, 1]),
           'green': lambda v : np.interp(v, [0, 1], [0, 1]),
           'blue': lambda v : np.interp(v, [0, 0.65, 0.98, 1], [0, 0, 1, 1])}




# Set aliases, where colormap exists in matplotlib
cm.cmap_d['ds9bb'] = cm.cmap_d['afmhot']
cm.cmap_d['ds9grey'] = cm.cmap_d['gray']

# Register all other colormaps
cm.register_cmap('ds9b', data=ds9b)
cm.register_cmap('ds9cool', data=ds9cool)
cm.register_cmap('ds9a', data=ds9a)
cm.register_cmap('ds9i8', data=ds9i8)
cm.register_cmap('ds9aips0', data=ds9aips0)
cm.register_cmap('ds9rainbow', data=ds9rainbow)
cm.register_cmap('ds9he', data=ds9he)
cm.register_cmap('ds9heat', data=ds9heat)

def real2tex( x, b=0, n=0 ):
    if x == '':
        return x
    x = ('%e'%x).split( 'e' )
    if float(x[0]) == 1:
        x = r'10^{%i}'%( int( x[1] )+b )
    else:
        fmt = '%%.%if'%(n)
        fmt = r'%s \times'%(fmt) + '\, 10^{%i}'
        x = (fmt)%( float(x[0]), int( x[1] )+b )
    return x

def remove_tick_labels( ax, xy, bar=0 ):
    if 'x' in xy:
        t = [ l.get_text() for l in ax.get_xticklabels() ]
        print( "remove xtick labels:" )
        print( t )
        ll = [''] * len(t)
        ax.set_xticklabels( ll )
        if bar == 0:
            ax.tick_params( axis='x', width=0 )
        t = [ l.get_text() for l in ax.get_xticklabels() ]
        print( t )
    if 'y' in xy:
        t = [ l.get_text() for l in ax.get_yticklabels() ]
        print( "remove ytick labels:" )
        print( t )
        ll = [''] * len(t)
        ax.set_yticklabels( ll )
        if bar == 0:
            ax.tick_params( axis='y', width=0 )

def fmt_tick_labels( ax, xy ):
    if 'x' in xy:
        t = [ l.get_text() for l in ax.get_xticklabels() ]
        tt = []
        for x in t:
            xx = x.split( '$' )
            if len(xx) < 2:
                tt.append( xx )
            else:
                tt.append( r'$%s$'%real2tex( float(xx[1]) ) )
        ax.set_xticklabels( tt )
    if 'y' in xy:
        t = [ l.get_text() for l in ax.get_yticklabels() ]
        tt = []
        for x in t:
            xx = x.split( '$' )
            if len(xx) < 2:
                tt.append( xx )
            else:
                tt.append( r'$%s$'%real2tex( float(xx[1]) ) )
        ax.set_yticklabels( tt )

def make_log_ticks( xmin, xmax, n, a=0, axis=None ):
    xmin = np.log10( xmin )
    xmax = np.log10( xmax )
    if a == 0:
        a = 1
        while (xmax-xmin)/a > 10:
            a += 1
    loc = np.arange( int(xmin)-1, int(xmax)+2, a )
    fmt = [ r'$10^{%.0f}$'%i for i in loc ] 
    loc = ( loc - xmin ) / ( xmax-xmin ) * ( n-1 )
    t1 = []
    t2 = []
    for i in range( len(loc) ):
        if loc[i]>0 and loc[i]<n:
            t1.append( loc[i] )
            t2.append( fmt[i] )
    loc = tik.FixedLocator( t1 )
    fmt = tik.FixedFormatter( t2 )
    if axis:
        axis.set_major_locator( loc )
        axis.set_major_formatter( fmt )
    return ( loc, fmt )

def set_tick_params( ax, ls, axis='both', pad=5 ):
    ax.tick_params( axis=axis, direction='in', pad=pad, labelsize=ls )
    ax.minorticks_off()
