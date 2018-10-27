#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

my_markersize = 2

flag256 = 1

def top_hat_filter( k, R ):
    kR = k * R;

    if ( kR < 1e-4 ):
        w = 1
    else:
        w = 3.0 * ( np.sin(kR)/kR**3 - np.cos(kR)/kR**2 )

    return w

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

def plot_top_hat_filter():

    R = 8
    h = 0.7

    k = np.linspace( -4, np.log10( 350/(R/h) ), 100 )
    k = np.power( 10, k )

    w = []
    for i in k:
        w.append( top_hat_filter( i, R/h ) )

    plt.plot( k, w, label="R: %g, h:%g"%(R,h) )

    plt.xscale( 'log' )
    plt.yscale( 'log' )
    plt.xlabel( 'k' )
    plt.ylabel( 'w[k]' )

    plt.legend()
    plt.savefig( "top_hat_filter.pdf" )
    plt.close()

def test_ps_interp(k, p):

    kk = np.linspace( np.log( k.min() ), np.log( k.max() ), 100 )
    kk = np.exp( kk )

    pp = []
    for i in kk:
        pp.append(ps_interp( i, k, p ))

    plt.plot( k, p,  '.', label='origin' )
    plt.plot( kk, pp, '*', label='interp' )

    plt.xscale( 'log' )
    plt.yscale( 'log' )

    plt.legend()
    plt.savefig( 'test_ps_interp.pdf' )
    plt.close()


def calc_sigma8( k, p ):

    print( "calc_sigma8 ..." )

    R = 8

    plot_top_hat_filter()
    test_ps_interp( k, p )

    f = lambda x: ps_interp( x, k, p ) * x**2 * top_hat_filter( x, R )**2

    k0 = 1e-99 / R
    k1 = 350 / R

    r = quad( f, k0, k1 )

    print( r )

    sigma8 = np.sqrt( r[0] / ( 2 * np.pi**2 ) )

    print( "sigma8: ", sigma8 )

    return sigma8



#sigma8_camb = 0.8371
#sigma8_camb = 0.80
#sigma8 = 0.817
#sigma8 = 0.837


sigma8 = 0.817

if flag256:
    print( "*256*" )
    BoxSize = 100
    Nmesh  = 256
else:
    print( "*128*" )
    BoxSize = 50
    Nmesh  = 128

kmin = 2 * np.pi / BoxSize
kmax = 2 * np.pi / ( BoxSize / Nmesh )

print( 'load halofit... ' )
halo_fit = np.loadtxt( './decay_matterpower.dat' )

sigma8_camb = calc_sigma8( halo_fit[:, 0], halo_fit[:,1] )

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

#########################################################
print( 'load gadtet-tools...' )
data = np.loadtxt( fn )

k = data[:, 0]
k = k * 1000 #convert to Mpc
Delta = data[:, 1]
k = k[ Delta>0 ]
Delta = Delta[ Delta> 0 ]

#print( Delta )

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

#########################################################
print( 'load gadtet...' )
if flag256:
    gdata = np.loadtxt( './hge_256_100_powerspec_057.txt', skiprows=2008 )
else:
    gdata = np.loadtxt( './hge_128_50_powerspec_057.txt', skiprows=2008 )

k = gdata[:, 0]
k = k * 1000
Delta = gdata[:, 1]

k = k[Delta>0]
Delta = Delta[Delta>0]


#print( Delta )
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
