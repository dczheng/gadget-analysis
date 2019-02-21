#!/usr/bin/env python3

from my_work_env import *
import tools_and_constants as tc
from scipy.optimize import curve_fit

def my_plot1():
    rad_spec_file1 = './cre_256_no_turb_Spec_0.09_nosr.dat'
    rad_spec_file2 = './cre_256_Spec_0.09_nosr.dat'

    ele_spec_file1 = './cre_256_no_turb_Elec_Spec_0.09.dat'
    ele_spec_file2 = './cre_256_Elec_Spec_0.09.dat'

    dat_r1 = np.loadtxt( rad_spec_file1 )
    dat_r2 = np.loadtxt( rad_spec_file2 )

    dat_e1 = np.loadtxt( ele_spec_file1 )
    dat_e2 = np.loadtxt( ele_spec_file2 )

    fig_r, axs_r= plt.subplots( 1, 2 )
    fig_e, axs_e= plt.subplots( 1, 2 )

    m, n = dat_r1.shape

    v = dat_r1[ 0, 2: ]
    p = dat_e1[ 0, 2: ]

    for i in range( 1, m ):
        M = dat_r1[ i, 1 ] * 1e10
        index = dat_r1[ i, 0 ]
        P1 = dat_r1[ i, 2: ] / tc.mJy
        P2 = dat_r2[ i, 2: ] / tc.mJy

        F1 = dat_e1[i, 2:]
        F2 = dat_e2[i, 2:]

        t = "%.2e"%M
        t = t.split( 'e' )
        if t[1][0] == '+':
            t[1] = t[1][1:]

        if i % 2:
            ss = '--'
        else:
            ss = '-'

        #label=r"$[%i]\ %s\times 10^{%s}\,M_\odot$"%(index, t[0], t[1] )
        if ( i == 1 ):
            ss = 'k-'
            label=r"$%s\times 10^{%s}\,M_\odot$"%(t[0], t[1] )
        else:
            label = ''

        axs_r[0].loglog( v, P1, ss, label=label )
        axs_r[1].loglog( v, P2, ss, label=label )

        axs_e[0].loglog( p, F1, ss, label=label )
        axs_e[1].loglog( p, F2, ss, label=label )


    for i in range( 2 ):
        axs_r[i].legend( framealpha=0.1 )
        axs_r[i].set_xlabel( r'$\nu \, [\rm MHz]$' )
        axs_r[i].set_ylabel( r'$S \, [\rm mJy]$' )


        #ylim = list(axs_r[i].get_ylim())
        #ylim[0] = 1e-5
        #axs_r[i].set_ylim( ylim )
        axs_r[i].grid()
        axs_r[i].minorticks_off()
        axs_r[i].locator_params( numticks=10 )

        axs_e[i].legend( framealpha=0.1 )
        axs_e[i].set_xlabel( r'$p$' )
        axs_e[i].set_ylabel( r'$F$' )

        axs_e[i].grid()
        axs_e[i].minorticks_off()
        axs_e[i].locator_params( numticks=10 )


    axs_r[0].set_title( r'$S2$' )
    axs_r[1].set_title( r'$S3$' )

    axs_e[0].set_title( r'$S2$' )
    axs_e[1].set_title( r'$S3$' )

    if ( len(sys.argv) < 2 ):
        file_pre = ''
    else:
        file_pre = sys.argv[1] + '-'

    fig_r.tight_layout()
    fig_e.tight_layout()

    fig_r.savefig( output_dir + file_pre + 'rad_spec.pdf' )
    fig_e.savefig( output_dir + file_pre + 'ele_spec.pdf' )

def fit_f( x, a, b ):
    return x*a + b

def my_plot0():

    ele_spec_file = './cre_256_no_turb_Elec_Spec_0.00.dat'

    dat_e = np.loadtxt( ele_spec_file )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )

    m, n = dat_e.shape

    print( m, n )

    p = dat_e[ 0, 2: ]
    t = np.abs(p-1e3)
    ii3 = np.where( t == t.min() )[0][0]
    t = np.abs(p-1e5)
    ii5 = np.where( t == t.min() )[0][0]

    print( "ii3: %i [%g], ii5: %i [%g]"%(ii3, p[ii3], ii5, p[ii5] ) )

    for i in range( 1, m ):

        M = dat_e[ i, 1 ] * 1e10
        index = dat_e[ i, 0 ]

        F = dat_e[i, 2:]

        t = "%.2e"%M
        t = t.split( 'e' )
        if t[1][0] == '+':
            t[1] = t[1][1:]

        if i % 2:
            ss = '--'
        else:
            ss = '-'

        #label=r"$[%i]\ %s\times 10^{%s}\,M_\odot$"%(index, t[0], t[1] )
        if ( i == 1 ):
            ss = 'k-'

        y = np.log10( F[ii3:ii5] )
        x = np.log10( p[ii3:ii5] )
        r = curve_fit( fit_f, x, y )
        print( 'ele: ', r[0], r[1] )
        alpha_e = r[0][0]

        label_e = r'$%.2f$'%(-alpha_e)

        if ( i > 10 ):
            continue

        ax.loglog( p, F, ss, label=label_e )


    ax.grid()
    ax.minorticks_off()

    ax.set_xlabel( r'$p$', fontsize=20 )
    ax.set_ylabel( r'$f(p) \, [\rm cm^{-3}]$', fontsize=20 )
    ax.legend( framealpha=0.1, title=r'$\alpha$')
    ylim = list(ax.get_ylim())
    ylim[0] = 1e-25
    ax.set_ylim( ylim )
    ax.tick_params( axis='both', pad=5, direction='in', labelsize=15 )

    fig.tight_layout()

    fig.savefig( output_dir + 'ele0_spec.pdf', figsize=(5,5) )


def my_plot2():

    rad_spec_file = './cre_256_no_turb_Spec_0.09_nosr.dat'
    ele_spec_file = './cre_256_no_turb_Elec_Spec_0.09.dat'

    dat_r = np.loadtxt( rad_spec_file )
    dat_e = np.loadtxt( ele_spec_file )

    fig, axs= plt.subplots( 1, 2 )

    m, n = dat_r.shape

    print( m, n )

    v = dat_r[ 0, 2: ]
    t = np.abs(v-300)
    i300 = np.where( t == t.min() )[0][0]
    t = np.abs(v-1400)
    i1400 = np.where( t == t.min() )[0][0]

    print( "i300: %i [%g]), i1400: %i [%g])"%(i300, v[i300], i1400, v[i1400] ) )
    p = dat_e[ 0, 2: ]

    t = np.abs(p-1e3)
    ii3 = np.where( t == t.min() )[0][0]
    t = np.abs(p-1e5)
    ii5 = np.where( t == t.min() )[0][0]

    print( "ii3: %i [%g], ii5: %i [%g]"%(ii3, p[ii3], ii5, p[ii5] ) )

    ele_alphas = []
    rad_alphas = []

    for i in range( 1, m ):

        M = dat_r[ i, 1 ]
        index = dat_r[ i, 0 ]

        P = dat_r[ i, 2: ] / tc.mJy
        F = dat_e[i, 2:]

        t = "%.2e"%M
        t = t.split( 'e' )
        if t[1][0] == '+':
            t[1] = t[1][1:]

        if i % 2:
            ss = '--'
        else:
            ss = '-'

        #label=r"$[%i]\ %s\times 10^{%s}\,M_\odot$"%(index, t[0], t[1] )
        if ( i == 1 ):
            ss = 'k-'

        y = np.log10( F[ii3:ii5] )
        x = np.log10( p[ii3:ii5] )
        r = curve_fit( fit_f, x, y )
        print( 'ele: ', r[0], r[1] )
        alpha_e = r[0][0]
        ele_alphas.append( -alpha_e )

        y = np.log10( P[i300:i1400] )
        x = np.log10( v[i300:i1400] )
        r = curve_fit( fit_f, x, y )
        print( 'rad: ', r[0], r[1] )
        alpha_r = r[0][0]
        rad_alphas.append( -alpha_r )

        label_e = r'$%.2f$'%(-alpha_e)
        label_r = r'$%.2f$'%(-alpha_r)

        if ( i > 6 ):
            continue

        axs[0].loglog( p, F, ss, label=label_e )
        axs[1].loglog( v, P, ss, label=label_r )


    for i in range( 2 ):
        axs[i].grid()
        axs[i].minorticks_off()
        #axs[i].locator_params( numticks=10 )

        if i == 0:
            axs[i].set_xlabel( r'$p$', fontsize=20 )
            axs[i].set_ylabel( r'$f(p) \, \rm [cm^{-3}]$', fontsize=20 )
            axs[i].legend( framealpha=0.1, title=r'$\alpha$')
            ylim = list(axs[i].get_ylim())
            ylim[0] = 1e-25
            axs[i].set_ylim( ylim )
        else:
            axs[1].set_xlabel( r'$\nu \, [\rm MHz]$', fontsize=20 )
            axs[1].set_ylabel( r'$S \, [\rm mJy]$', fontsize=20 )
            axs[i].legend( framealpha=0.1, title=r'$\alpha_{\rm rad}$' )
        axs[i].tick_params( axis='both', pad=5, direction='in', labelsize=15 )

    if ( len(sys.argv) < 2 ):
        file_pre = ''
    else:
        file_pre = sys.argv[1] + '-'

    fig.tight_layout()

    fig.savefig( output_dir + file_pre + 'ele_rad_spec.pdf', figsize=(5,5) )

    #return
    fig, axs = plt.subplots( 1, 2 )

    axs[0].hist( ele_alphas )
    axs[1].hist( rad_alphas )

    axs[0].set_xlabel( r'$\alpha$' )
    axs[0].set_ylabel( r'$Number$' )

    axs[1].set_xlabel( r'$\alpha_{\rm rad}$' )
    axs[1].set_ylabel( r'$Number$' )

    fig.tight_layout()
    fig.savefig( output_dir + file_pre + 'ele_rad_spec_index.pdf', figsize=(8,8) )


my_plot2()
my_plot0()
