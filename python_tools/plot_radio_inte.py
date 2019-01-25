#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )
d = np.loadtxt( './radio_inte.dat' )

p = d[:,0]
f = d[:,1]

plt.loglog( p, f )
plt.legend()

plt.savefig( './radio_inte.pdf', figsize=(5,5) )

