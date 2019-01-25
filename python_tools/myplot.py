#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import sys

#data_0001 = np.loadtxt( './hge256_data/hge_256_100_0.001_0.01_Spec_Tot.dat' )
data = np.loadtxt( sys.argv[1] )

#plt.plot( data_0001[:,0], data_0001[:,1]*1e23, label="CRE100_005" )
plt.plot( data[:,0], data[:,1]*1e23 )

plt.xscale( 'log' )
plt.yscale( 'log' )

plt.xlabel( r'$\nu \; [MHz]$' )
plt.ylabel( r'$I_{\nu} \; [Jy/sr]$' )

plt.title( "Radio Spectrum" )

plt.legend()
plt.grid()

plt.savefig( sys.argv[1][:-4] + '.png' )

