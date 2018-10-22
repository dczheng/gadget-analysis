#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt

my_markersize = 2

flag256 = 0

if flag256:
    BoxSize = 50
    Nmesh  = 128
else:
    BoxSize = 100
    Nmesh  = 256

kmin = 2 * np.pi / BoxSize
kmax = 2 * np.pi / ( BoxSize / Nmesh )

halo_fit = np.loadtxt( './decay_matterpower.dat' )

sigma8_camb = 0.8371
sigma8 = 0.817

k = halo_fit[:, 0]
Delta = halo_fit[:, 1] * k**3 / ( 2 * np.pi ** 2 )

#plt.plot( k, Delta, label='HaloFit_bef' )

Delta = Delta / (sigma8_camb**2) * sigma8**2

Dmin = Delta.min()
Dmax = Delta.max()

plt.plot( k, Delta, label='HaloFit' )


if flag256:
    fn = "hge_256_100_PowSpec.dat"
else:
    fn = "hge_128_50_PowSpec.dat"

data = np.loadtxt( fn )

k = data[:, 0]
k = k * 1000 #convert to Mpc
Delta = data[:, 1]
k = k[ Delta>0 ]
Delta = Delta[ Delta> 0 ]

print( Delta )

plt.plot( k, Delta, '.-', label='gadget-tools', markersize=my_markersize )
#plt.plot( k, Delta, '.' )

print( "K: ", kmin, kmax, k.min(), k.max() )
print( "Delta: ", Dmin, Dmax, Delta.min(), Delta.max() )

if k.min() > kmin:
    kmin = k.min()
if k.max() < kmax:
    kmax = k.max()

if Delta.min() > Dmin:
    Dmin = Delta.min()
if Delta.max() < Dmax:
    Dmax = Delta.max()

if flag256:
    gdata = np.loadtxt( './hge_256_100_powerspec_057.txt', skiprows=2008 )
else:
    gdata = np.loadtxt( './hge_128_50_powerspec_057.txt', skiprows=2008 )

k = gdata[:, 0]
k = k * 1000
Delta = gdata[:, 1]

k = k[Delta>0]
Delta = Delta[Delta>0]


print( Delta )
plt.plot( k, Delta, '.', label='gadget', markersize=my_markersize )

plt.xlabel( r'$k \; [hMpc^{-1}]$' )
plt.ylabel( r'$\Delta^2(k)$' )

plt.xscale( 'log' )
plt.yscale( 'log' )

print( "kmin: %g, kmax: %g"%(kmin, kmax) )
print( "Dmin: %g, Dmax: %g"%(Dmin, Dmax) )

plt.xlim( [kmin, kmax] )
plt.ylim( [Dmin, Dmax] )

plt.legend()
plt.grid()

plt.title( "Power Spectrum", family='monospace', style='normal' )

if flag256:
    fn = "hge_256_100_compare_powspec.pdf"
    print( 'save to ' + fn )
else:
    fn = "hge_128_50_compare_powspec.pdf"
    print( 'save to ' + fn )
plt.savefig( fn )
