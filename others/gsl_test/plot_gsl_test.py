#!/usr/bin/env python3

import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt( './gsl_bessel_Knu.dat' )
x = d[:,0]
K = d[:,1]
IK = d[:,2]
fig, axs = plt.subplots( 1,2 )
axs[0].loglog( x, K )
axs[1].loglog( x, IK )
plt.savefig( './gsl_bessel_Knu.png' )
plt.close()
