#!/usr/bin/env python3

from my_work_env import *

from scipy.special import beta as sbeta
from scipy.special import betainc as sbetainc
from scipy.integrate import quad


gamma = lambda p: np.sqrt( 1 + p*p )
beta = lambda p: p / gamma(p)

coul_fac1 = lambda b, n_e: 3 * mycc.sigma_t * mycc.c * n_e / ( 2 * (b**2 ) )

omega_plasma = lambda n_e: np.sqrt( 4 * np.pi * mycc.e2 * n_e / mycc.m_e  )

coul_fac2 = lambda g, b, n_e: mycc.m_ec2 * b * np.sqrt( g - 1 ) / ( mycc.hbar * omega_plasma( n_e ) )

def beta_func( x, a, b ):

    if x<1e-30:
        return 0

    #print( x )

    if b>0:
        r = sbeta(a,b) * sbetainc(a,b,x)
        return r
    else:
        f = lambda x: x**(a-1) * (1-x)**(b-1)
        r = quad( f, 1e-20*x, x )
        return r[0]

def cre_n( c, a, q ):
    return c * np.power( q, 1-a ) / ( a-1 )

def cre_n2( c, a, q1, q2 ):
    return cre_n( c, a, q1 ) - cre_n( c, a, q2 )

def cre_Tbar( a, q ):

    aa = ( a-2 )/2
    bb = ( 3-a ) / 2
    x = 1 / ( 1+q*q )

    return ( np.power( q, a-1 ) / 2 * beta_func( x, aa, bb ) + \
             np.sqrt( 1+q*q ) - 1 ) * mycc.m_ec2

def cre_Tbar_gadget( a, q ):
    Tbar = cre_Tbar( a, q )
    return Tbar / mycc.gadget_energy_in_erg

def cre_e( c, a, q ):
    return cre_n( c, a, q ) * cre_Tbar( a, q )

def cre_e_gadget( c, a, q ):
    e = cre_n( c, a, q ) * cre_Tbar( a, q ) / mycc.m_e
    return e / ( mycc.gadget_energy_in_erg / mycc.gadget_mass_in_g )

def cre_e_gadget_2( c, a, q1, q2 ):
    return cre_e_gadget( c, a, q1 ) - cre_e_gadget( c, a, q2 )

def cre_T( q ):
    return mycc.m_ec2 * ( np.sqrt( 1+q*q ) - 1 )

def cre_f( c, a, q, qmin ):

    if q<qmin:
        return 0
    return c * np.power( q, -a )

def cre_dedt( c, a, qmin, dpdt ):

    return -cre_f( c, a, qmin, qmin ) * \
            cre_T( qmin ) * \
            dpdt( qmin )

def cre_dndt( c, a, qmin, dpdt ):

    return -cre_f( c, a, qmin, qmin ) * \
            dpdt( qmin )

def plot_cre_e():
    q = np.power( 10, np.linspace( -1, 7, 100 ) )
    a = 3.01
    c = 1
    e = [ cre_e(c, a, qq) for qq in q ]
    plt.plot( q, e, label='C: %g, a: %g'%( c, a ) )
    plt.xlabel( r'$q$' )
    plt.xscale( 'log' )
    plt.ylabel( r'$\epsilon$' )
    plt.yscale( 'log' )
    plt.savefig( 'e.pdf' )



def dp_dt_coul ( p, n_e ):

    g = gamma(p)
    b = beta( p )

    fac1 = coul_fac1( b, n_e )
    fac2 = coul_fac2( g, b, n_e )

    return fac1* ( \
                np.log( fac2 ) \
                - np.log(2) * ( 0.5 * b**2 + 1/g ) \
                + 0.5 \
                + ( 0.25 * ( 1 - 1/g ) )**2 \
                )

def dp2_dpdt_coul( p, n_e ):

    g = gamma(p)
    g2 = g**2
    g3 = g**3
    b = beta( p )
    fac1 = coul_fac1( b, n_e )

    return -dp_dt_coul( p, n_e ) * 2 / ( b * g3  ) \
                                + fac1 * ( \
                                            0.5 * b / (g - 1) + 1 / ( b * g3 ) \
                                            - np.log(2) * b / g2 * ( 1/g -1 ) \
                                            + 0.125 *  b / g3 * ( g-1 )\
                                            )
def dp_dt_syn( p, B ):

    g = gamma( p )

    Ccool = B**2 / ( 8 * np.pi ) * 4 * mycc.sigma_t / ( 3 * mycc.m_ec )

    return Ccool * g * p

def dp2_dpdt_syn( p, B ):

    return dp_dt_syn( p, B ) * ( p / (1+p**2) + 1 / p )

def dp_dt_bre( p, n ):

    return 6.622e-22 * n * 1e6 * ( gamma(p) - 1 ) / beta(p)

def dp2_dpdt_bre( p, n ):

    g = gamma( p )
    b = beta( p )

    return 3.66e-22 * n * 1e6 * ( 1 - ( g-1 ) / ( p**2 * g ) )

def comp_T_with_c():

    c_dat = np.loadtxt( output_dir + 'CRE_T.dat' )
    p = c_dat[:,0]
    a = 3
    T = [ cre_Tbar( a, i ) for i in p ]
    T = np.array( T )
    T = T / mycc.m_ec2

    del_t = ( c_dat[:,1] - T ) / T

    fig, axs = plt.subplots( 2,1 )

    axs[0].plot( p, T )
    axs[1].plot( p, del_t )

    for i in range( 2 ):
        ax = axs[i]
        ax.set_xscale( 'log' )
        ax.set_xlabel( r'$p$' )

    axs[0].set_ylabel( r'$T$' )
    axs[0].set_yscale( 'log' )
    axs[1].set_ylabel( r'$\Delta T$' )

    fig.savefig( output_dir + 'T_py_c.pdf' )

    fd = open( output_dir + 'T_py.dat', 'w' )

    for i in range( len(p) ):
        fd.write( "%g %g\n"%( p[i], T[i] ) )

    fd.close()

def plot_cooling_timescale():

    c = 1
    a = 2.1
    B = 1e-6
    n = 1000 * mycc.rho_bar_crit0 / mycc.m_p
    title = r'$\rho=10^3\,\rho_{\rm bar}, B=1\mu \rm G$'
    ne = n
    pmin = 1e-2
    pmax = 1e9
    p = np.power( 10, np.linspace( np.log10(pmin), np.log10(pmax), 1000 ) )

    dpdt_f = [\
            lambda p: dp_dt_coul( p, ne ),\
            lambda p: dp_dt_coul( p, ne ),\
            lambda p: dp_dt_syn( p, B ),\
            lambda p: dp_dt_syn( p, B )
            ]
    labels = [ \
            r'$\tau_{\epsilon, \rm coul}$',\
            r'$\tau_{n, \rm coul}$',\
            r'$\tau_{\epsilon, \rm syn}$',\
            r'$\tau_{n, \rm syn}$'\
            ]
    ls = [\
            'k-',\
            'b-',\
            'k--',\
            'b--'\
            ]

    dx_dt_f = [\
            cre_dedt, \
            cre_dndt, \
            cre_dedt, \
            cre_dndt
            ]

    dx_f = [\
            cre_e, \
            cre_n, \
            cre_e, \
            cre_n
            ]

    fig, ax = plt.subplots( 1,1 )

    for i in range( 4 ):
        tau = np.array([ dx_f[i]( c, a, pp) / np.abs(dx_dt_f[i]( c, a, pp,  dpdt_f[i] )) for pp in p ])
        tau = tau / mycc.Myr
        ax.plot( p, tau, ls[i], label=labels[i] )
        #print( tau )

    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_ylabel( r'$\tau [\rm Myr]$' )
    ax.set_xlabel( r'$p$' )
    ax.axhline( y = 1, linestyle='-.', color='y', linewidth=1, label=r'$\tau=1 \, \rm Myr$' )
    ax.legend( title = title, \
            prop={'size':6}, framealpha=0.1 )

    fig.savefig( output_dir + 'tau.pdf', figsize=(5,5) )



def calc_dpdt_dp2dpdt( p, n, ne, B ):

    dpdt_coul = dp_dt_coul( p, ne )
    dpdt_syn = dp_dt_syn( p, B )
    dpdt_bre = dp_dt_bre( p, n )

    dp2dpdt_coul = dp2_dpdt_coul( p, ne )
    dp2dpdt_syn  = dp2_dpdt_syn( p, B )
    dp2dpdt_bre  = dp2_dpdt_bre( p, n )

    dpdt = dpdt_coul + dpdt_syn + dpdt_bre
    dp2dpdt = dp2dpdt_coul + dp2dpdt_syn + dp2dpdt_bre

    return  ( dpdt_coul, dpdt_syn, dpdt_bre, dp2dpdt_coul, dp2dpdt_syn, dp2dpdt_bre )

def plot_dp_dt( dpdt_dp2dpdt, p ):

    a = dpdt_dp2dpdt
    dpdt_coul = a[0]
    dpdt_syn  = a[1]
    dpdt_bre  = a[2]
    dp2dpdt_coul = a[3]
    dp2dpdt_syn  = a[4]
    dp2dpdt_bre  = a[5]

    axN = 2
    fig, axs = plt.subplots( 1, axN )

    ls = [ '-k', '--r', '-.y' ]
    ax = axs[0]
    ax.plot( p, dpdt_coul, ls[0],  label='coul' )
    ax.plot( p, dpdt_syn, ls[1], label='rad' )
    #ax.plot( p, dpdt_bre, ls[2], label='bre' )

    ax.set_ylabel( r'$\frac{dp}{dt}~[s^{-1}]$', fontsize=20 )

    ax = axs[1]

    ax.plot( p, dp2dpdt_coul, ls[0], label='coul' )
    ax.plot( p, dp2dpdt_syn, ls[1], label='rad' )
    #ax.plot( p, dp2dpdt_bre, ls[2], label='bre' )

    '''
    t0 = p[ dp2dpdt_coul < 0 ]
    t1 = dp2dpdt_coul[ dp2dpdt_coul < 0 ]
    t1 = np.abs( t1 )
    ax.plot( t0, t1, label='- coul' )
    '''

    ax.set_ylabel( r'$\frac{d^2p}{dpdt}~[s^{-1}]$', fontsize=20 )


    for i in range( axN ):
        axs[i].set_xscale( 'log' )
        axs[i].set_yscale( 'log' )
        axs[i].set_xlabel( r'$p$', fontsize=20)
        #axs[i].grid()
        axs[i].legend()
        axs[i].tick_params( axis='both', direction='in', labelsize=15, pad=5, length=0 )


    #plt.show()
    fig.tight_layout()
    fig.savefig( output_dir + 'dp_dt.pdf' )
    #fig.savefig( output_dir + "/elec_coul_cool.pdf" )
    plt.close()

def plot_dn( p, dpdt_dp2dpdt, dt, alpha, ax, title ):

    a = dpdt_dp2dpdt
    dpdt_coul = a[0]
    dpdt_syn  = a[1]
    dpdt_bre  = a[2]

    #1: bre
    #2: syn
    #3: coul
    #4: coul + syn + bre
    #5: coul + bre
    #6: coul + syn
    flag = 0b100110
    ls = ['b-', 'r--', 'y-.' ]

    dpdt = dpdt_coul + dpdt_syn
    r = -(alpha-1) / p * dpdt
    dn = np.exp( r * dt )
    ax.plot( p, dn, ls[0], label='coul+syn+ic')

    dpdt = dpdt_coul
    r = -(alpha-1) / p * dpdt
    dn = np.exp( r * dt )
    ax.plot( p, dn, ls[1], label='coul')

    dpdt = dpdt_syn
    r = -(alpha-1) / p * dpdt
    dn = np.exp( r * dt )
    ax.plot( p, dn, ls[2], label='syn+ic')

    ax.set_xscale( 'log' )
    #ax.set_yscale( 'log' )

    fs = 20
    ax.set_xlabel( r'$p$', fontsize=fs )
    ax.set_ylabel( r'$\exp{(\sum_{i}{-\frac{\alpha-1}{p}\frac{dp}{dt}|_i)}}$', fontsize=fs )

    ax.set_title( title, fontsize=fs )
    ax.tick_params( axis='both', labelsize=20, direction='in', pad=5 )
    ax.legend()


def plot_df( p, dpdt_dp2dpdt, dt, alpha, ax, title ):

    a = dpdt_dp2dpdt
    dpdt_coul = a[0]
    dpdt_syn  = a[1]
    dpdt_bre  = a[2]
    dp2dpdt_coul = a[3]
    dp2dpdt_syn  = a[4]
    dp2dpdt_bre  = a[5]

    #1: bre
    #2: syn
    #3: coul
    #4: coul + syn + bre
    #5: coul + bre
    #6: coul + syn
    flag = 0b100110
    ls = ['b-', 'k.', 'k.', 'r--', 'y-.', 'y-.']
    if flag & 0b100000:
        dpdt = dpdt_coul + dpdt_syn
        dp2dpdt = dp2dpdt_coul + dp2dpdt_syn
        r = -alpha / p * dpdt + dp2dpdt
        df = np.exp( r * dt )
        ax.plot( p, df, ls[0], label='coul+rad' )

    if flag & 0b10000:
        dpdt = dpdt_coul + dpdt_bre
        dp2dpdt = dp2dpdt_coul + dp2dpdt_bre
        r = -alpha / p * dpdt + dp2dpdt
        df = np.exp( r * dt )
        ax.plot( p, df, ls[1], label='syn+bre')

    if flag & 0b01000:
        dpdt = dpdt_coul + dpdt_syn + dpdt_bre
        dp2dpdt = dp2dpdt_coul + dp2dpdt_syn + dp2dpdt_bre
        r = -alpha / p * dpdt + dp2dpdt
        df = np.exp( r * dt )
        ax.plot( p, df, ls[2], label='all'  )

    if flag & 0b00100:
        dpdt = dpdt_coul
        dp2dpdt = dp2dpdt_coul
        r = -alpha / p * dpdt + dp2dpdt
        df = np.exp( r * dt )
        ax.plot( p, df, ls[3], label='coul' )

    if flag & 0b00010:
        dpdt = dpdt_syn
        dp2dpdt = dp2dpdt_syn
        r = -alpha / p * dpdt + dp2dpdt
        df = np.exp( r * dt )
        ax.plot( p, df, ls[4], label='rad' )

    if flag & 0b00001:
        dpdt = dpdt_bre
        dp2dpdt = dp2dpdt_bre
        r = -alpha / p * dpdt + dp2dpdt
        df = np.exp( r * dt )
        ax.plot( p, df, ls[5], label='bre' )

    ax.set_xscale( 'log' )
    #ax.set_yscale( 'log' )

    fs = 20
    ax.set_xlabel( r'$p$', fontsize=fs )
    ax.set_ylabel( r'$\exp{\left[-\Delta t \sum_{i}\left({\frac{\alpha}{p}\frac{dp}{dt}|_i - \frac{d2p}{dpdt}|_i}\right)\right]}$', fontsize=fs )

    ax.set_title( title, fontsize=fs )
    ax.tick_params( axis='both', labelsize=20, direction='in', pad=5 )
    ax.legend( prop={'size':15}, framealpha=0.1 )
    ax.grid()

    #plt.show()
    #fig.savefig( 'tt.pdf' )

def plot_dfs_dns():
    Bs = [ 1e-4, 1e-5, 1e-6 ]
    ni = [ 1, 100, 10000 ]

    tles1 = [ \
            r'$\rho = 1 \, \rho_{\rm bar}, B = 1\mu \rm G$',  \
            r'$\rho = 100 \, \rho_{\rm bar}, B = 1\mu \rm G$',  \
            r'$\rho = 10000 \, \rho_{\rm bar}, B = 1\mu \rm G$',  \
            ]

    tles2 = [ \
            r'$\rho = \rho_{\rm bar}, B = 100 \, \mu \rm G$',  \
            r'$\rho = \rho_{\rm bar}, B = 10 \, \mu \rm G$',  \
            r'$\rho = \rho_{\rm bar}, B = 1 \, \mu \rm G$',  \
            ]


    ns = np.array( ni ) * mycc.rho_bar_crit0 / mycc.m_p
    dt = 0.1 * mycc.Myr
    alpha = 3

    pmin = 1e-1
    pmax = 1e8
    p = np.power( 10, np.linspace( np.log10(pmin), np.log10(pmax), 1000 ) )

    N = 3

    fig, axs = plt.subplots( 2, N, figsize=(15, 10) )
    fign, axsn = plt.subplots( 2, N, figsize=(15, 10) )

    n = ns[1]
    B = 1e-4
    ne = n
    dpdt_dp2dpdt = calc_dpdt_dp2dpdt( p, n, ne, B )
    plot_dp_dt( dpdt_dp2dpdt, p )

    for i in range( N ):

        n = ns[i]
        ne = n
        B = 1e-6

        dpdt_dp2dpdt = calc_dpdt_dp2dpdt( p, n, ne, B )
        tle = tles1[i]

        plot_df( p, dpdt_dp2dpdt, dt, alpha, axs[0,i], tle )
        plot_dn( p, dpdt_dp2dpdt, dt, alpha, axsn[0,i], tle )

        B = Bs[i]
        n  = ns[1]
        ne = n

        dpdt_dp2dpdt = calc_dpdt_dp2dpdt( p, n, ne, B )
        tle = tles2[i]
        plot_df( p, dpdt_dp2dpdt, dt, alpha, axs[1,i], tle )
        plot_dn( p, dpdt_dp2dpdt, dt, alpha, axsn[1,i], tle )


    fig.tight_layout()
    fign.tight_layout()

    fig.savefig( output_dir + 'dfs.pdf' )
    fign.savefig( output_dir + 'dns.pdf' )

def plot_for_paper():

    dt = 0.1 * mycc.Myr
    alpha = 2.1
    B = 1e-6
    n = 100 * mycc.rho_bar_crit0 / mycc.m_p
    ne = n
    pmin = 1e-1
    pmax = 1e9
    p = np.power( 10, np.linspace( np.log10(pmin), np.log10(pmax), 1000 ) )

    dpdt_dp2dpdt = calc_dpdt_dp2dpdt( p, n, ne, B )
    plot_dp_dt( dpdt_dp2dpdt, p )

    fig, ax = plt.subplots( 1,1 )
    plot_df( p, dpdt_dp2dpdt, dt, alpha, ax, '' )

    fig.tight_layout()
    fig.savefig( output_dir + 'df.pdf', figsize=(3,3) )


def compare_dpdt_dp2dpdt_with_c(  ):

    B = 1e-5
    #n1 = mycc.rho_bar_crit0 / mycc.m_p
    n = 1 / mycc.m_p
    ne = n * mycc.Xh * ( 1+mycc.Xh ) / ( 2*mycc.Xh )

    c_dat = np.loadtxt( output_dir + 'c_dpdt_dp2dpdt.dat' )

    p = c_dat[:,0]
    py_dat = calc_dpdt_dp2dpdt( p, n, ne, B )

    fig, axs = plt.subplots( 2,2 )

    t = (c_dat[:,1] - py_dat[0]) / py_dat[0]
    axs[0,0].plot( p, t, '--', label='coul dpdt' )

    t = (c_dat[:,2] - py_dat[1]) / py_dat[1]
    axs[0,1].plot( p, t, '-', label='syn dpdt' )

    t = (c_dat[:,3] - py_dat[3]) / py_dat[3]
    axs[1,0].plot( p, t, '--', label='coul dp2dpdt' )

    t = (c_dat[:,4] - py_dat[4]) / py_dat[4]
    axs[1,1].plot( p, t, '-', label='syn dp2dpdt' )

    for i in range(2):
        for j in range(2):
            ax = axs[i,j]
            ax.set_xscale( 'log' )
            ax.set_xlabel( r'$p$' )
            ax.set_ylabel( r'$\delta$' )
            ax.legend()

    fig.tight_layout()
    fig.savefig( output_dir + 'com_dpdt_dp2dpdt_py_c.pdf' )

    fd = open( output_dir + 'py_dpdt_dp2dpdt.dat', 'w' )

    for i in range( len(p) ):
        fd.write( "%g %g %g %g %g\n"%( \
                p[i],
                py_dat[0][i],
                py_dat[1][i],
                py_dat[3][i],
                py_dat[4][i],
                ) )

    fd.close()

def rad_cooling():

    B = 1e-6
    n = 100 * mycc.rho_bar_crit0 / mycc.m_p
    ne = n
    pmin = 1e-1
    pmax = 1e9
    c = 1
    a = 2.5

    ps = np.power( 10, np.linspace( np.log10(pmin), np.log10(pmax), 1000 ) )
    fs = np.array([ cre_f( c, a, p, pmin ) for p in ps ])
    plt.loglog( ps, fs, label='0' )

    dt = 0.1 * mycc.Myr
    print( "dt: ", dt )

    Ccool = B**2 / ( 8 * np.pi ) * 4 * mycc.sigma_t / ( 3 * mycc.m_ec )
    print( "Ccool: ", Ccool )

    pmax = 1 / ( dt * Ccool * ( a-1 ) )
    print( "pmax: ", pmax )
    fs[ ps>pmax ] = 0
    plt.loglog( ps, fs, label='1' )

    dpdt = lambda p: dp_dt_syn(p, B)
    n = cre_n2( c, a, pmin, pmax )
    print( n )

    index = 1
    for i in range( 100 ):
        dndt = cre_dndt( c, a, pmax, dpdt ) - cre_dndt( c, a, pmin, dpdt )
        n = n - dndt * dt
        pmax = (cre_n( c, a, pmin) - n) / c * ( a-1 )
        pmax = np.power( pmax, 1/(1-a) )
        index += 1
        if i % 10 == 0:
            fs[ ps>pmax ] = 0
            plt.loglog( ps, fs, label='%i'%(index) )
        print( pmax )
        if pmax < pmin:
            exit()

    #print( cre_dndt( c, a, pmin, dpdt ) )
    #print( cre_dndt( c, a, pmax, dpdt ) )
    #for i in range(100):

    plt.legend()
    plt.savefig( figs_dir + 'rad_cooling.pdf' )

    exit()

def cooling_func_output():

    B = 1e-6
    pmin = 1e-4
    pmax = 1e8
    pN = 100
    n_e = 1e-3

    p = np.logspace( np.log10(pmin), np.log10(pmax), pN )
    rad = dp_dt_syn( p, B )
    rad2 = dp2_dpdt_syn( p, B )
    coul = dp_dt_coul( p, n_e )
    coul2 = dp2_dpdt_coul( p, n_e )

    fd = open( 'func_output.dat', 'w' )
    for i in range( len(p) ):
        fd.write( "%g %g %g %g %g\n"%(p[i], rad[i], rad2[i], coul[i], coul2[i]) )
    fd.close()



def main():
    #plot_dfs_dns()
    #compare_dpdt_dp2dpdt_with_c()

    #plot_for_paper()
    #comp_T_with_c()
    #plot_cooling_timescale()
    plot_cre_e();
    #rad_cooling()
    #cooling_func_output()



if __name__ == '__main__':
    main()
