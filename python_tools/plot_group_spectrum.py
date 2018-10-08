#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser()

parser.add_argument( '-f', type=str, help='file name', required=True )
parser.add_argument( '-t', type=str, help='data type', required=True )

args = parser.parse_args()

fn = args.f
data_type = args.t
data_type_all = [ 'ele', 'rad' ]

if not(data_type in data_type_all):
    print( "ERROR" )
    exit()

data = np.loadtxt( fn )
m, n = data.shape
print( "shape: ( %i, %i )"%( m, n ) )

v = data[ 0, 2:-1 ]
for i in range( 1, m ):
    m = data[ i, 1 ]
    index = data[ i, 0 ]
    p = data[ i, 2:-1 ]
    if i % 2:
        ss = '.-'
    else:
        ss = '-'
    t = "%.2e"%m
    t = t.split( 'e' )
    if t[1][0] == '+':
        t[1] = t[1][1:]
    plt.plot( v, p, ss, label=r"$(%04i)%s \times 10^{%s}$"%(index, t[0], t[1] ) )

plt.xscale( 'log' )
plt.yscale( 'log' )

if ( data_type == 'rad' ):
    plt.xlabel( r'$\nu \; [MHz]$' )
    plt.ylabel( r'$ I_{\nu}\; [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1} ]$' )

if ( data_type == 'ele' ):
    plt.xlabel( r'$q$' )
    plt.ylabel( r'$f \; [cm^{-3}]$' )

plt.legend()
plt.grid()

plt.title( 'z = ' + fn[-8:-4] )

png_fn = fn[:-4] + '.png'
print( 'save fig to ' + png_fn )
plt.savefig( png_fn )
