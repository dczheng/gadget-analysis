#!/usr/bin/env python3

import h5py
from my_work_env import *

fn = sys.argv[1]
f = h5py.File( fn, 'r' )
b = f[ '/PartType0/MagneticField' ].value
#print( b )
b = b[b>0]
b = np.log10(b)
bins = np.linspace( -26, -2, 50 )
plt.hist( b, bins=bins )
plt.yscale('log')
plt.savefig( fn+'_mag_hist.png' )

