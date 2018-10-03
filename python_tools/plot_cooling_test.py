#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.loadtxt( sys.argv[1] )
q = data[ :,0 ]
df = data[ :,1 ]
f = data[ :,2 ]
ff = data[ :,3 ]
df = np.abs( df )

qmin = q[-2]
qmax = q[-1]
fmin = f[-2]
fmax = f[-1]

print( qmin, qmax, fmin, fmax )

q = q[:-2]
f = f[:-2]
df = df[:-2]
ff = ff[:-2]

t = sys.argv[1][:-4].split( '_' )
rho = int( t[2] )
a   = float( t[3] )
c   = float( t[4] )
B   = float( t[5] )
print( "rho: %i, a: %g, c: %g, B: %g"%( rho, a, c, B ) )

plt.plot( q, df, '-.', label="df" )
plt.plot( q, f, label="f" )
plt.plot( q, f-df, '--', label="f-df" )
plt.plot( q, ff, '.', label="ff" )
plt.plot( qmin, fmin, '*', markersize=8 )
plt.plot( qmax, fmax, '*', markersize=8 )

plt.title( '%i'%(rho) )

plt.xscale( 'log' )
plt.yscale( 'log' )
plt.legend()
plt.grid()

plt.savefig( sys.argv[1][:-4] + '.png' )


