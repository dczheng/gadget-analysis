#!/usr/bin/env python3
import numpy as np
import scipy.constants as sc
from scipy.integrate import quad
import matplotlib.pyplot as plt


Omega_m = 0.3175#0.3
Omega_Lambda = 0.6825#0.7
Omega_k = 0.0
Omega_r = 0.0
H0 = 70#67.11#70

def E( a ):
    return np.sqrt( Omega_m / np.power( a, 3.0 ) + \
                    Omega_r / np.power( a, 4.0 ) + \
                    Omega_Lambda +\
                    Omega_k / np.power( a, 2.0 ) )

def dt_da( a ):
    return 1.0 / E(a) / a

H0_in_SI = H0 / ( 1.0e6 * sc.parsec ) * 1.0e3

print( "H0 in SI: ", H0_in_SI )

N = 1000
a_min_i = 1e-30
a_min = 1.0e-9
a = np.linspace( np.log10(a_min), 0, N )
a = np.power(10,a)
t = np.zeros( N )
t_err = np.zeros( N )
#tmp = quad( dt_da, a_min_i, 1e-5 )
#print( tmp )
for i in range( N ):
    tmp = quad( dt_da, a_min_i, a[i] )
    t[i] = tmp[0]
    t_err[i] = tmp[1]
t = t / H0_in_SI / sc.year / 1.0e9
dtda = dt_da( a )
#plt.plot( a, t )
#plt.plot( a, t_err )
fig = plt.figure()
ax_t = fig.add_subplot( 221 )
ax_dt_da = fig.add_subplot( 222 )
ax_H_a = fig.add_subplot( 223 )
ax_p = fig.add_subplot( 224 )

#ax_t.plot( a, t  )
z = 1/a - 1
ax_t.plot(a, t)
#ax_t.loglog( a, t )
ax_t.set_title( r'$age\ (Gyr)$' )
#ax_t.set_xlabel( r'$a$' )
ax_t.grid()


ax_dt_da.plot( a, dtda )
ax_dt_da.set_title( r'$\frac{dt}{da}\ (Gyr)$' )
#ax_dt_da.set_xlabel( r'$a$' )
ax_dt_da.grid()

a = np.linspace( 0.5, 1 )
H_a = H0 * E( a )
ax_H_a.plot( a, H_a )
ax_H_a.set_ylabel( r'$H(a)\ ( Km/s/Mpc )$' )
ax_H_a.set_xlabel( r'$a$' )
ax_H_a.grid()

ax_p.text( 0.3, 0.85, r'$\Omega_m = %.3f$'%(Omega_m) )
ax_p.text( 0.3, 0.70, r'$\Omega_\Lambda = %.3f$'%(Omega_Lambda) )
ax_p.text( 0.3, 0.55, r'$\Omega_k = %.3f$'%(Omega_k) )
ax_p.text( 0.3, 0.40, r'$\Omega_r = %.3f$'%(Omega_r) )
ax_p.text( 0.3, 0.25, r'$H_0 = %.3f$'%(H0) )


#plt.savefig( './cosm_t_%.3f_%.3f_%.3f_%.3f_%.3f.png'%( Omega_m, Omega_Lambda, Omega_k, Omega_r, H0  ) )
plt.savefig( 'time.png' )

