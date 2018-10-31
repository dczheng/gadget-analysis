#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import argparse

parser = argparse.ArgumentParser()

parser.add_argument( 'f' )
parser.add_argument( '-t', type=str, help='data type', required=True )
parser.add_argument( '-fit', help='fitting', action='store_true' )
parser.add_argument( '-pdf', help='save pdf', action='store_true' )
parser.add_argument( '-nosr', help='no sr', action='store_true' )

args = parser.parse_args()

fn = args.f
data_type = args.t
data_type_all = [ 'ele', 'rad' ]

if not(data_type in data_type_all):
    print( "ERROR" )
    exit()

data = np.loadtxt( fn )

if len(data) == 0:
    print( "No Group!" )
    exit()

m, n = data.shape

#v = data[ 0, 2:-1 ]
v = data[ 0, 2: ]
logv = np.log(v)

def f( x, a, b ):
    return x*a+b

for i in range( 1, m ):
    m = data[ i, 1 ]
    index = data[ i, 0 ]
    #p = data[ i, 2:-1 ]
    p = data[ i, 2: ]

    p = p * 1e23

    if i % 2:
        ss = '.-'
    else:
        ss = '-'
    t = "%.2e"%m
    t = t.split( 'e' )
    if t[1][0] == '+':
        t[1] = t[1][1:]

    if args.fit:
        print( '[%i] fitting:'%i )
        logp = np.log(p)
        x = logv[ logp != -np.inf ]
        y = logp[ logp != -np.inf ]
        r = curve_fit( f, x, y )
        print( r )
        plt.plot( v, p, ss, label=r"$[%i]\ %s\times 10^{%s} \; , \alpha: %.2f$"%(index, t[0], t[1], np.abs( r[0][0] ) ) )
    else:
        plt.plot( v, p, ss, label=r"$[%i]\ %s\times 10^{%s}\,M_\odot$"%(index, t[0], t[1] ) )


plt.xscale( 'log' )
plt.yscale( 'log' )

if ( data_type == 'rad' ):
    plt.xlabel( r'$\nu \; [MHz]$' )

    if args.nosr:
        plt.ylabel( r'$ I_{\nu}\; [Jy]$' )
    else:
        plt.ylabel( r'$ I_{\nu}\; [Jy/sr]$' )

if ( data_type == 'ele' ):
    plt.xlabel( r'$q$' )
    plt.ylabel( r'$f \; [cm^{-3}]$' )

plt.legend()
plt.grid()

plt.title( 'z = ' + fn[-8:-4] )

if args.pdf:
    fn = fn[:-4] + '.pdf'
else:
    fn = fn[:-4] + '.png'
print( 'save fig to ' + fn )

plt.savefig( fn )
