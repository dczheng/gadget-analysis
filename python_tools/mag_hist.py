#!/usr/bin/env python3

import h5py
from my_work_env import *

fn = sys.argv[1]
f = h5py.File( fn, 'r' )
bb = f[ '/PartType0/MagneticField' ].value
#b = [ np.sqrt(x[0]**2+x[1]**2+x[2]**2) for x in bb ]
m, n = bb.shape
b = bb.reshape(m*n, 1)
b = b[b>0] * 1e6 # to muG

print( b.max() )
print( "> 100 muG" )
t = b[b>100] # 100 muG
print( t )
print( len(t) )
print( "> 1000 muG" )
t = b[b>1000] # 1000 muG
print( t )
print( len(t) )

b = np.log10(b)
bins = np.linspace( -20, 3, 1000 )
plt.hist( b, bins=bins )
plt.yscale('log')
plt.savefig( fn+'_mag_hist.png' )

