#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )

import numpy as np
import argparse
import matplotlib.pyplot as plt
#import matplotlib.ticker as tik
#import matplotlib.colors as mplc

parser = argparse.ArgumentParser()

parser.add_argument( 'f' )
parser.add_argument( '-g', help='plot gadget data', action='store_true' )
parser.add_argument( '-eps', help='save eps', action='store_true' )
parser.add_argument( '-nolog', help='disable log axis', action='store_true' )
parser.add_argument( '-K0', help='minimum k',type=float )
parser.add_argument( '-K1', help='maximu k',type=float )
parser.add_argument( '-P0', help='minimum P',type=float )
parser.add_argument( '-P1', help='maximu P',type=float )
parser.add_argument( '-Mpc', help='Mpc Input', action='store_true' )

args = parser.parse_args()

if args.g:
    print( 'load gadget data ...' )
    data = np.loadtxt( args.f, skiprows=2008 )
    #print( data[0:2,:] )
else:
    print( 'load gadget-tools data ...' )
    data = np.loadtxt( args.f )

k = data[:,0]
p = data[:,1]

fig = plt.figure()
ax = fig.add_subplot( 111 )

if not( args.Mpc ):
    k = k * 1000

k = k[ p>0 ]
p = p[ p>0 ]

ax.plot( k, p, '.' )
ax.plot( k, p, '-.' )

if not(args.nolog):
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )

ax.set_xlabel( r'$k\,[h/Mpc]$' )
ax.set_ylabel( r'$\Delta^2(k)\,[\frac{4 \pi k^3}{(2 \pi L)^3} P(k)]$' )

if args.K0 or args.K1:
    klim = list(plt.xlim())

    if args.K0:
        klim[0] = args.K0

    if args.K1:
        klim[1] = args.K1

    plt.xlim( klim )

if args.P0 or args.P1:
    plim = list(plt.ylim())

    if args.P0:
        plim[0] = args.P0

    if args.P1:
        plim[1] = args.P1

    plt.ylim( plim )

ax.set_title( 'Power Spectrum' )


if args.eps:
    file_suffix = '.eps'
else:
    file_suffix = '.png'

fn = args.f[0:-4] + file_suffix

plt.savefig( fn )
