#!/usr/bin/env python3
import h5py
import numpy as np
import sys

f = h5py.File( sys.argv[1], 'r' )
z = f['Header'].attrs[ "Redshift" ]
h = f['/PartType0/SmoothingLength'][()]
m = f[ 'PartType0/Masses' ][()]
d = f[ 'PartType0/Density' ][()]
divb = f[ 'PartType0/DivergenceOfMagneticField' ][()]
b = f[ 'PartType0/MagneticField' ][()]
b = np.sqrt(b[:,0]**2 + b[:,1]**2 + b[:,2]**2)
f.close()


idx = (b>0) * (divb>0) 
h = h[idx]
m = m[idx]
b = b[idx]
d = d[idx]
divb = np.fabs(divb[idx])
bmu = b * 1e6   # to muG

v = m / d
err = divb * h / b 

print( '-' * 30 + '%.2f'%z + '-'*30 )

def f( v1, v2, v ):
    idx = (v>v1) * (v<v2)
    vv = v[idx]
    print( "[%8.2f~%-8.2f [muG]]: %i"%( v1, v2, len(vv) ) )

for i in range(-2, 4):
    f( 10**(i), 10**(i+1), bmu )

print( "[divb error] min: %g, max: %g, mean: %g"%(\
    err.min(), err.max(), (err*v).sum() / v.sum()) )

print( '-' * 30 + '%.2f'%z + '-'*30 )
