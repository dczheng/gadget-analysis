#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

BoxSize = 50
Nmesh  = 128
my_markersize = 2
kmin = 2 * np.pi / BoxSize
kmax = 2 * np.pi / ( BoxSize / Nmesh )

gs = gridspec.GridSpec( 2, 1, height_ratios=[3,1] )
ax1 = plt.subplot( gs[0] )
ax2 = plt.subplot( gs[1] )

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



'''
sigma8_camb = 0.8371
sigma8 = 0.817
halo_fit = np.loadtxt( './decay_matterpower.dat' )
k = halo_fit[:, 0]
Delta = halo_fit[:, 1] * k**3 / ( 2 * np.pi ** 2 )
Delta = Delta / ( sigma8_camb**2 ) * sigma8 ** 2
Dmin = Delta.min()
Dmax = Delta.max()
ax1.plot( k, Delta, label='HaloFit' )
'''

zeta = 0.01
chi_list = [ 0.001, 0.01, 0.02, 0.03 ]

ps0 = np.loadtxt( './hge_128_50_PowSpec.dat' )
k0 = ps0[:,0] * 1000
D0 = ps0[:,1]

k0 = k0[ D0 > 0 ]
D0 = D0[ D0 > 0 ]
ax1.plot( k0, D0, '-.', label="NoCRE", markersize=my_markersize )

if k0.min() > kmin:
    kmin = k0.min()
if k0.max() < kmax:
    kmax = k0.max()

kk = np.linspace( np.log10( kmin ), np.log10( kmax ), 100 )
kk = np.power( 10, kk )

for chi in chi_list:
    fn = "hge_128_50_%g_%g_PowSpec.dat"%(zeta, chi)

    ps = np.loadtxt( fn )
    k = ps[:,0] * 1000
    D = ps[:,1]

    k = k[ D>0 ]
    D = D[ D>0 ]

    ax1.plot( k, D, '-.', label=r'$\zeta=%g, \; \chi=%g$'%(zeta, chi), markersize=my_markersize )

    dd = []
    for i in kk:
        d1 = ps_interp( i, k0, D0 )
        d2 = ps_interp( i, k, D )

        if d1 == d2:
            dd.append( 0 )
        else:
            dd.append( ( d2-d1 ) / d1 )

    ax2.plot( kk, dd, label=r'$\zeta=%g, \; \chi=%g$'%(zeta, chi) )


#ax1.set_xlabel( r'$k \; [hMpc^{-1}]$' )
ax1.set_ylabel( r'$\Delta^2(k)$' )
ax1.set_xscale( 'log' )
ax1.set_yscale( 'log' )
ax1.legend()
ax1.set_title( "Power Spectrum" )

ax2.set_xlabel( r'$k \; [hMpc^{-1}]$' )
ax2.set_ylabel( r'$\frac{|\Delta^2-\Delta^2_{NoCRE}|}{\Delta^2_{NoCRE}}$' )
ax2.legend(fontsize=4)

ax2.set_xscale( 'log' )

plt.savefig( "hge_128_50_powspec.pdf" )
