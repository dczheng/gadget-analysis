#!/usr/bin/env python3
import sys
import numpy as np

f = open( sys.argv[1], 'r' )
data = f.readlines()
f.close

fa = open( 'time_list.txt', 'w' )
fz = open( 'z_list.txt', 'w' )
a = np.array( [] )
z = np.array( [] )
for i in data:
    t = i.split()
    if ( len( t )<3 ): continue
    z1 = float( t[0] )
    z2 = float( t[1] )
    zn = int( t[2] )
    if ( z2 == z1 ): continue
    if ( z2 < z1 ):
        z = z2
        z2 = z1
        z1 = z
    print( "z1 = %f, z2 = %g, zn = %i"%(z1, z2, zn) )
    zl = np.linspace( z1, z2, zn+1 )
    zl = zl[:-1]
    al = 1 / ( zl+1 )
    a = np.hstack( [a, al] )
    z = np.hstack( [z, zl] )
a = np.sort( a )
z = -np.sort( -z )
print( a )
print( z )
for i in a:
    fa.write( "%f\n"%(i) )
for i in z:
    fz.write( "%f\n"%(i) )
fa.close()
fz.close()
