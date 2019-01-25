#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

BoxSize = 100
Nmesh  = 256

sigma8_camb = 0.8371
sigma8 = 0.817

my_markersize = 2

gs = gridspec.GridSpec( 2, 1, height_ratios=[4,1] )

ax1 = plt.subplot( gs[0] )
ax2 = plt.subplot( gs[1] )

kmin = 2 * np.pi / BoxSize
kmax = 2 * np.pi / ( BoxSize / Nmesh )

halo_fit = np.loadtxt( './decay_matterpower.dat' )

k = halo_fit[:, 0]
Delta = halo_fit[:, 1] * k**3 / ( 2 * np.pi ** 2 )

Delta = Delta / ( sigma8_camb**2 ) * sigma8 ** 2

Dmin = Delta.min()
Dmax = Delta.max()

ax1.plot( k, Delta, label='HaloFit' )

zeta = 0.01
chi_list = [ 0.001, 0.01, 0.02, 0.03 ]

fn_list = []
pic_label_list = []
for chi in chi_list:
    fn = "hge_128_50_%g_%g_PowSpec.dat"%(zeta, chi)
    fn_list.append( fn )
    pic_label_list.append( r'$\zeta=%g, \; \chi=%g$'%(zeta, chi) )

fn_list.append( "hge_128_50_PowSpec.dat" )
pic_label_list.append( "NoCRE" )


for i in range( len(fn_list) ):

    fn = fn_list[i]
    pic_label = pic_label_list[i]

    data = np.loadtxt( fn )

    k = data[:, 0]
    k = k * 1000 #convert to Mpc
    Delta = data[:, 1]

    k = k[ Delta>0 ]
    Delta = Delta[ Delta> 0 ]

    ax1.plot( k, Delta, '-.', label=pic_label, markersize=my_markersize )

    print( "K: %.3f, %.3f, %.3f, %.3f"%(kmin, kmax, k.min(), k.max()) )
    print( "Delta: %e, %e, %e, %e"%( Dmin, Dmax, Delta.min(), Delta.max() ) )

    if k.min() > kmin:
        kmin = k.min()
    if k.max() < kmax:
        kmax = k.max()

    if Delta.min() > Dmin:
        Dmin = Delta.min()
    if Delta.max() < Dmax:
        Dmax = Delta.max()


kmax = 6
#ax1.set_xlabel( r'$k \; [hMpc^{-1}]$' )
ax1.set_ylabel( r'$\Delta^2(k)$' )

ax1.set_xscale( 'log' )
ax1.set_yscale( 'log' )

print( "kmin: %g, kmax: %g"%(kmin, kmax) )
print( "Dmin: %g, Dmax: %g"%(Dmin, Dmax) )

ax1.set_xlim( [kmin, kmax] )
ax1.set_ylim( [Dmin, Dmax] )

ax1.legend()

ax1.set_title( "Power Spectrum" )


data = np.loadtxt( fn_list[-1] )
k_no_cre = data[:,0]
k_no_cre = k_no_cre * 1000 #convert to Mpc
Delta_no_cre = data[:,1]

for i in range( len(fn_list)-1 ):

    fn = fn_list[i]
    pic_label = pic_label_list[i]

    data = np.loadtxt( fn )

    k = data[:, 0]
    k = k * 1000 #convert to Mpc
    Delta = data[:, 1]

    dDelta = np.abs( Delta-Delta_no_cre )

    index = dDelta > 0
    k = k[ index ]
    dDelta = dDelta[index] / Delta_no_cre[index]

    ax2.plot( k, dDelta, label=pic_label )

ax2.set_xlabel( r'$k \; [hMpc^{-1}]$' )
ax2.set_ylabel( r'$\frac{|\Delta^2-\Delta^2_{NoCRE}|}{\Delta^2_{NoCRE}}$' )
ax2.legend(fontsize=4)

ax2.set_xscale( 'log' )

ax2.set_xlim( [kmin, kmax] )

plt.savefig( "hge_128_50_powspec.pdf" )
