#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from astropy.cosmology import LambdaCDM
m_e = 9.10953e-28
e_e = 4.8032e-10
c   = 2.9979e10
erg_2_Gev = 624.15
sig_t = 6.56246e-25
pc = 2.0857e18
Jansky = 10e-23

B = [ -5, -6, -7, -8, -9 ]
p = np.linspace( 100, 10000, 1000 )
g = np.sqrt( 1+p*p )
E = m_e * g * c * c
#print( m_e * c * c * erg_2_Gev * 1e3 )
#print( E * erg_2_Gev )
for b in B:
    bb = 10 ** b
    #v_c = 3 / np.pi / 4 * e_e * bb / m_e * ( E / m_e ) ** 2
    #print( bb )
    v_c = e_e * bb / 2 / np.pi / m_e / c * g * g
    plt.plot( E * erg_2_Gev, v_c / 1e6, label=r'$B=e^{%i}\mu G$'%(b) )
plt.xlabel( 'Gev' )
plt.ylabel( 'MHz' )
#plt.xscale( 'log' )
plt.yscale( 'log' )
plt.legend()
#plt.show()
plt.savefig( '1.png' )
plt.close()

alpha = 3.1
v_c = 1.4e9
n0 = 1e-9
g0 = 1
L = 1000 * 1000 * pc
v = np.linspace( 1, 5000, 100 ) #MHz

cos = LambdaCDM( H0=70, Om0=0.3, Ode0=0.7 )
dis = cos.luminosity_distance( 0.2 ).value * 1e6 * pc
print( dis )
for b in B:
    bb = 10 ** b
    u_b = bb * bb / ( 8 * np.pi )
    E_tot = 2 / 3 * sig_t * u_b * n0 / v_c * np.power( v / v_c, -(alpha-1)/2 ) * np.power( L, 3 )
    #print( E_tot.min() )
    flux = E_tot / ( 4*np.pi*dis*dis ) / Jansky * 1000
    #print( 4 * np.pi * dis * dis )
    #plt.plot( v, E_tot, label=r'$B=e^{%i}\mu g$'%(b) )
    plt.plot( v, flux, label=r'$B=e^{%i}\mu g$'%(b) )
plt.ylabel( r'$(mJy)$' )
plt.xlabel( 'MHz' )
plt.yscale( 'log' )
#plt.xscale( 'log' )
plt.legend()
#plt.show()
plt.savefig( '2.png' )

