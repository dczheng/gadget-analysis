#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import scipy.special as ss
from scipy.integrate import quad

import tools_and_constants as tc
from my_work_env import *

def _K( x ):
    return ss.kv( 5.0/3.0, x )

def _F( x ):

    if x < 1e-5:
        return 0
    r = quad( _K, x, np.inf, epsabs=0, epsrel=1e-3 )
    #print( x, r )
    return x * r[0]

def _f( p, c, a, pmin, pmax ):
    if p<pmin or p>pmax:
        return 0
    return c * np.power( p, -a )

def dp2dVdv( B, nu, c, a, pmin, pmax ):

    fac = np.sqrt(3) * np.pi / 4 * tc.e3 * B / tc.m_ec2
    #print( fac )

    ff = lambda p: _f( p, c, a, pmin, pmax ) \
            * _F( nu / ( 3 * ( 1+p*p ) * tc.e * B / ( 16.0 * tc.m_ec) ) )

    r = quad( ff, pmin, pmax, epsabs=0, epsrel=1e-3 )
    #print( r )

    return r[0]*fac

def plot_F():

    x = np.power( 10, np.linspace( -5, 1, 1000 ) )
    y = [ _F(i) for i in x ]
    #print( y[0] )
    plt.plot( x, y )
    plt.savefig( output_dir + 'F_x_py.pdf' )


print( "e: %e"%tc.e )
print( "e3: %e"%tc.e3 )
print( "m: %e"%tc.m_e )
print( "c: %e"%tc.c )
print( "mc: %e"%tc.m_ec )
print( "mc2: %e"%tc.m_ec2 )

#plot_F()

def calc_test():
    B = 1e-8
    z = 0.09
    a = 3.1
    pmin = 3.125
    pmax = 191200
    C = 1.21e-10

    #L = 0.1 * tc.Mpc
    L = 5 * tc.Kpc
    lc = tc.D_c( z )
    la = tc.D_a( z )
    ll = tc.D_l( z )
    nu = 1e9
    flux =  dp2dVdv( B, nu, C, a, pmin, pmax )
    flux = flux * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / ( L*L/(lc*lc) ) # tc.mJy

    flux = flux * 4500 + flux * 10 * 2258 + flux*100 * 500
    print( "flux: %g"%flux )


def plot_flux():
    B = 1e-6
    z = 0.1
    a = 1 / ( 1+z )
    pmin = 0.1
    alpha = 3

    #nu = 1.4 * 1e9
    nua = np.power( 10, np.linspace( 0, 3, 50 ) ) # MHz
    pmaxa = np.power( 10, np.linspace(4,5, 8) )
    L = 0.1 * tc.Mpc

    lc = tc.D_c( z )
    la = tc.D_a( z )
    ll = tc.D_l( z )

    print( "z: %g, a: %g"%(z, a) )
    print( "lc: %e, la: %e, ll: %e"%( lc, la, ll ) )
    print( "la/lc: %g, lc/ll: %g"%(la/lc, lc/ll) )

    fig, ax = plt.subplots( 1, 1 )

    C = 1
    for pmax in pmaxa:
        print( pmin, pmax )
        r = np.array([ dp2dVdv( B, nu*1e6, C, alpha, pmin, pmax ) for nu in nua ])
        r = r * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / tc.mJy
        t = "%e"%pmax
        t = t.split('e')
        ax.plot( nua, r, label=r'$p_{\rm max} = %.1f \times 10^{%i}$'%(float(t[0]), int(t[1])) )

    ax.set_xlabel( r'$\nu \, [\rm MHz]$', fontsize=20 )
    ax.set_ylabel( r'$\frac{1}{C}\frac{d^{2}P(\nu)}{dVd\nu}\,[\rm erg\,s^{-1}\,Hz^{-1}]$', fontsize=20)
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.minorticks_off()

    ax.tick_params( axis='both', direction='in', labelsize=20, pad=5 )
    ax.grid()
    ax.legend()
    fig.tight_layout()

    fig.savefig( output_dir + 'Flux.pdf', figsize=(5,5) )

calc_test()
