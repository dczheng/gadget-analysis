#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def k1( q ):
    if 0 <= q <= 0.5:
        return 8.0 / np.pi * ( 1 - 6*q*q + 6*q*q*q )
    elif q <= 1:
        return 8.0 / np.pi * ( 2.0 * np.power( 1-q, 3.0 ) )
    else:
        return 0

fig = plt.figure()
ax_k1 = fig.add_subplot( 221 )
ax_k1_2d = fig.add_subplot( 222 )
N = 100
x = np.linspace( 0, 1.2, N )
y = np.zeros( N )
for i in range( N ):
    y[i] = k1( x[i] )
ax_k1.plot( x, y )
ax_k1.set_title( 'k1' )


N = 256
img = np.zeros( (N,N) )
for i in range( N ):
    for j in range( N ):
        x = (j-N/2) / (N/2)
        y = (i-N/2) / (N/2)
        r = np.sqrt( x*x + y*y )
        img[ i, j ] = k1(r)
k_2d=ax_k1_2d.imshow( img )
ax_k1_2d.set_title( 'k1_2d' )
plt.colorbar( k_2d )

plt.show()
