#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker as tik
#import matplotlib.colors as mplc
#from matplotlib import cm
import sys

if ( len( sys.argv ) != 2 ):
    print( "Please give data file!" )
    exit()

data = np.loadtxt( sys.argv[1] )
m, n = data.shape
print( "( %i, %i )"%( m, n ) )

v = data[ 0, 2:-1 ]
for i in range( 1, m ):
    m = data[ i, 1 ]
    if m < 1e14:
        break
    p = data[ i, 2:-1 ]
    if i % 2:
        ss = '.-'
    else:
        ss = '-'
    plt.plot( v, p, ss, label="%.2e"%( m ) )

plt.xscale( 'log' )
plt.yscale( 'log' )
plt.xlabel( r'$\nu \; [MHz]$' )
plt.ylabel( r'$ I_{\nu}\; [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1} ]$' )
plt.legend()
plt.title( 'z = ' + sys.argv[1][-8:-4] )
png_fn = sys.argv[1][:-4] + '.png'
print( 'save fig to ' + png_fn )
plt.savefig( png_fn )
