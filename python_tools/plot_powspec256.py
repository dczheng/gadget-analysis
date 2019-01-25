#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

BoxSize = 100
Nmesh  = 256

def ps_interp( k0, k, p ):

    if ( k0 < k[0] or k0 > k[-1] ):
        return 0

    if ( k0 == k[0] ):
        return p[0]

    i = 0
    while( k[i] < k0 ):
        i = i+1

    i = i-1

    dk = ( k0 - k[i] ) / ( k[i+1] - k[i] )

    return p[i] * ( 1-dk ) + p[i+1] * dk


gs = gridspec.GridSpec( 2, 1, height_ratios=[3,1] )
ax1 = plt.subplot( gs[0] )
ax2 = plt.subplot( gs[1] )

kmin = 2 * np.pi / BoxSize
kmax = 2 * np.pi / ( BoxSize / Nmesh )



ps = np.loadtxt( './hge_256_100_PowSpec.dat' )
ps_cre = np.loadtxt( './hge_256_100_0.01_0.01_PowSpec.dat' )

k = ps[:,0] * 1000
Delta = ps[:,1]
k = k[ Delta > 0]
Delta = Delta[ Delta > 0 ]

if k.min() > kmin:
    kmin = k.min()
if k.max() < kmax:
    kmax = k.max()

k_cre = ps_cre[:,0] * 1000
Delta_cre = ps_cre[:,1]
k_cre = k_cre[ Delta_cre > 0 ]
Delta_cre = Delta_cre[ Delta_cre > 0 ]

if k_cre.min() > kmin:
    kmin = k_cre.min()
if k_cre.max() < kmax:
    kmax = k_cre.max()

ax1.plot( k, Delta_cre, '-.', label="CRE256" )
ax1.plot( k_cre, Delta_cre, '.', label="CRE256_01_01", markersize=2 )

#ax1.set_xlabel( r'$k \; [hMpc^{-1}]$' )
ax1.set_ylabel( r'$\Delta^2(k)$' )
ax1.set_xscale( 'log' )
ax1.set_yscale( 'log' )
ylim = ax2.get_ylim()

#ax1.vlines( 2*np.pi/BoxSize, ylim[0], ylim[1], linestyle='-.', label=r'$\frac{2 \pi}{L}$' )
#ax1.vlines( 2*np.pi/BoxSize*Nmesh, ylim[0], ylim[1], linestyle='-.', label=r'$\frac{2 \pi N}{L}$' )

ax1.legend()
ax1.set_title( "Power Spectrum ( z=0 )" )

kk = np.linspace( np.log10( kmin ), np.log10( kmax ), 100 )
kk = np.power( 10, kk )
dd = []

for i in kk:
    d1 = ps_interp( i, k, Delta )
    d2 = ps_interp( i, k_cre, Delta_cre )
    if d1 == d2:
        dd.append( 0 )
    else:
        dd.append( ( d2-d1 ) / d1 )

ax2.plot( kk, dd )
ax2.set_xlabel( r'$k \; [hMpc^{-1}]$' )
ax2.set_ylabel( r'$\|\Delta^2-\Delta^2_{NoCRE}| \; / \; \Delta^2_{NoCRE}$' )
#ax2.legend(fontsize=4)
ax2.set_xscale( 'log' )

plt.savefig( "hge_256_100_powspec.pdf" )
