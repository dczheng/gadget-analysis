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

def dp2dVdv2( B, nu, c, a, pmin, pmax ):

    fac1 = 2 * np.pi * c * tc.e2 * nu / (np.sqrt(3) * tc.c)
    fac2 = 16 * tc.m_ec * nu / ( 3 * tc.e * B  )
    print( 'fac1: %g, fac2: %g'%(fac1, fac2) )

    ff = lambda t: np.power( fac2 / t - 1, -(1+a)/2 ) * _F(t) / (t*t)

    '''
    tt = np.power( 10, np.linspace( -5, 5, 20 ) )
    fff = []
    for t in tt:
        fff.append( ff(t) )
    fig, ax = plt.subplots( 1,1 )
    ax.loglog( tt, fff )
    fig.savefig( 'x.pdf' )
    '''

    tmin = fac2 / ( 1 + pmax**2 )
    tmax = fac2 / ( 1 + pmin**2 )
    print( 'tmin: %g, tmax: %g'%(tmin, tmax) )

    r = quad( ff, tmin, 1000, epsabs=0, epsrel=1e-3 )
    #print( r )

    return r[0]*fac1


def test_B():
    nu = 1e9
    c = 1
    a = 3.1
    pmin = 0.1
    pmax = 1e8
    Ba = np.power( 10, np.linspace( -6, -9, 10 ) )

    ff = lambda p: _f( p, c, a, pmin, pmax ) \
            * _F( nu / ( 3 * ( 1+p*p ) * tc.e * B / ( 16.0 * tc.m_ec) ) )

    fff = []
    for B in Ba:
        r = quad( ff, pmin, pmax, epsabs=0, epsrel=1e-3 )
        print( r )
        fff.append( r[0] )

    fff = np.array( fff )

    fig, axs = plt.subplots( 1, 2 )
    ax = axs[0]
    ax.plot( Ba, fff )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xlabel( 'B' )
    ax.set_ylabel( 'P' )

    ax = axs[1]
    ax.plot( Ba, fff*Ba )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xlabel( 'B' )
    ax.set_ylabel( 'PB' )

    fig.tight_layout()
    fig.savefig( 'B_fff.png' )

    print( fff )
    print( fff*Ba )


def plot_F():

    x = np.power( 10, np.linspace( -5, 1, 1000 ) )
    y = [ _F(i) for i in x ]
    #print( y[0] )
    plt.plot( x, y )
    plt.savefig( output_dir + 'F_x_py.pdf' )


#plot_F()


def plot_B_C():

    d = np.loadtxt( sys.argv[1] )
    index = d[:,-1] > 0
    dd = d[index,:]
    B = dd[:,3]
    C = dd[:,4]
    a = dd[:,5]
    pmin = dd[:,6]
    pmax = dd[:,7]
    n = dd[:,8]
    e = dd[:,9]
    z = 0.09

    print( "mean, C: %g, B: %g"%(C.mean(), B.mean()) )

    L = 5 * tc.Kpc
    lc = tc.D_c( z )
    la = tc.D_a( z )
    ll = tc.D_l( z )

    fac1 = np.sqrt(3) * np.pi / 4 * tc.e3 / tc.m_ec2
    fac1 = fac1 * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / tc.mJy
    fac2 = 3 * tc.e / ( 16.0 * tc.m_ec)
    print( len(e) )
    print( 'fac1: %g, fac2: %g'%(fac1, fac2) )

    fig,axs = plt.subplots( 1, 2 )
    #ax.plot( B, C, '.', ms=1 )

    ax = axs[0]
    ax.hist2d( np.log10(B), np.log10(C), bins=30 )
    #plt.colorbar( img )
    #ax.set_xscale( 'log' )
    #ax.set_yscale( 'log' )
    ax.set_xlabel( 'B' )
    ax.set_ylabel( 'C' )

    ax = axs[1]
    CB = C * np.power( B, (a+1)/2 )
    print( "CB sum: %g"%CB.sum() )

    ax.hist( np.log10( CB ), bins=30 )
    ax.set_xlabel( r'$C B^{\frac{\alpha+1}{2}}$' )
    ax.set_ylabel( 'N' )
    ax.set_yscale( 'log' )

    fig.tight_layout()
    fig.savefig( "B_C.png" )

def calc_test():
    B = 1e-10
    z = 0.09
    a = 2.1
    pmin = 0.1
    pmax = 1e7
    C = 1e-9

    #L = 0.1 * tc.Mpc
    L = 1 * tc.Kpc
    lc = tc.D_c( z )
    la = tc.D_a( z )
    ll = tc.D_l( z )
    nu = 1e9
    print( "nu: %g, B: %g, z: %g, a: %g, pmin: %g, " \
            "pmax: %g, c: %g, L^2: %g, ll^2: %g"
            %(nu, B, z, a, pmin, pmax, C, L**3, ll**2 ) )
    flux =  dp2dVdv( B, nu, C, a, pmin, pmax )
    print( "dp2dVdv: %g"%flux )
    flux = flux * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / tc.mJy #/  ( L*L/(lc*lc) ) # tc.mJy
    print( "flux: %g"%flux )

    '''
    flux =  dp2dVdv2( B, nu, C, a, pmin, pmax )
    flux = flux * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / tc.mJy #/  ( L*L/(lc*lc) ) # tc.mJy
    print( "flux2: %g"%flux )
    '''

def calc_group_radio():

    #from mpi4py import MPI

    '''
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    '''

    d = np.loadtxt( sys.argv[1] )
    index = d[:,-1] > 0
    dd = d[index,:]
    B = dd[:,3]
    C = dd[:,4]
    a = dd[:,5]
    pmin = dd[:,6]
    pmax = dd[:,7]
    n = dd[:,8]
    e = dd[:,9]

    print( np.histogram( a ) )

    #print( dd[:,3:] )
    #print( len(C[C>1e-9]), C.max() )
    #exit()

    z = 0.09
    L = 5 * tc.Kpc / (1+z) / 0.68
    lc = tc.D_c( z )
    la = tc.D_a( z )
    ll = tc.D_l( z )
    N = len(e)

    '''
    if rank == 0:
        print( "N: %i"%N )
    '''

    nu1 = 1e6
    nu2 = 1e9

    flux1 = 0
    flux2 = 0
    for i in range( N ):
        if ( a[i] > 3 ):
            continue
        pmax[i] = 1e7
        #a[i] = 3.1
        #C[i] = 1e-9
        '''
        if rank == 0:
            if i % 10 == 0:
                print( '[%8i], flux[%g]: %g, flux[%g]: %g'
                        %(i, nu1, flux1, nu2, flux2))
        if i % size:
            continue
        '''
        '''
        if i % 10 == 0:
            print( '[%8i, %8i], flux[%g]: %g, flux[%g]: %g'
                 %(i, N, nu1, flux1, nu2, flux2))
        '''

        print( '[%8i, %8i]B: %e, C: %e, a: %g, pmin: %g, pmax: %e' \
                    %(i, N, B[i], C[i], a[i], pmin[i], pmax[i] ) )

        r = dp2dVdv( B[i], nu1, C[i], a[i], pmin[i], pmax[i] )
        r = r * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / tc.mJy
        print( "r: %g, flux1: %g"%(r, flux1) )
        flux1 += r

        r = dp2dVdv( B[i], nu2, C[i], a[i], pmin[i], pmax[i] )
        r = r * ( 4/3 * np.pi * L**3 ) / ( np.pi * 4 * ll * ll ) / tc.mJy
        print( "r: %g, flux2: %g"%(r, flux2) )
        flux2 += r

    print( "nu: %g, f: %g"%(nu1, flux1) )
    print( "nu: %g, f: %g"%(nu2, flux2) )
    '''
    flux1 = np.array(comm.gather( flux1, root=0 ))
    flux2 = np.array(comm.gather( flux2, root=0 ))

    if rank == 0:
        print( "nu: %g, f: %g"%(nu1, flux1.sum()) )
        print( "nu: %g, f: %g"%(nu2, flux2.sum()) )
    '''


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

def part_angle():
    L = 0.1 * tc.Mpc
    z = 0.09
    lc = tc.D_c( z )
    print( "angle: %g"%(L**2 / lc**2) )

#print( "e: %e"%tc.e )
#print( "e3: %e"%tc.e3 )
#print( "m: %e"%tc.m_e )
#print( "c: %e"%tc.c )
#print( "mc: %e"%tc.m_ec )
#print( "mc2: %e"%tc.m_ec2 )

print( 16 / 3.0 / 0.3 * tc.m_ec / tc.e )
#test_B()
#plot_B_C()


#part_angle()
#calc_test()
#calc_group_radio()
