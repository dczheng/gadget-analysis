#!/usr/bin/env python3

from my_work_env import *
from scipy.optimize import curve_fit

def fit_f( x, a, b ):
    return x*a + b

sss = [ '-', '--', '-.'  ]
def my_plot0():

    ele_spec_file =  sys.argv[1]
    fn_out = sys.argv[2]

    dat_e = np.loadtxt( ele_spec_file )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )

    m, n = dat_e.shape

    print( m, n )

    p = dat_e[ 0, 2: ]
    t = np.abs(p-1e3)
    ii3 = np.where( t == t.min() )[0][0]
    t = np.abs(p-1e5)
    #ii5 = np.where( t == t.min() )[0][0]
    ii5 = -1

    print( "ii3: %i [%g], ii5: %i [%g]"%(ii3, p[ii3], ii5, p[ii5] ) )

    for i in range( 1, 50 ):

        M = dat_e[ i, 1 ] * 1e10
        index = dat_e[ i, 0 ]

        F = dat_e[i, 2:]

        t = "%.2e"%M
        t = t.split( 'e' )
        if t[1][0] == '+':
            t[1] = t[1][1:]

        ss = sss[ i % len(sss) ]

        #label=r"$[%i]\ %s\times 10^{%s}\,M_\odot$"%(index, t[0], t[1] )
        if ( i == 1 ):
            ss = 'k-'

        y = np.log10( F[ii3:ii5] )
        x = np.log10( p[ii3:ii5] )
        r = curve_fit( fit_f, x, y )
        print( x.min(), y.min() )
        print( 'ele: ', r[0], r[1] )
        alpha_e = r[0][0]

        label_e = r'$%.2f$'%(-alpha_e)

        if ( i > 20 ):
            continue

        ax.loglog( p, F, ss, label=label_e )


    ax.minorticks_off()

    ax.set_xlabel( r'$p$', fontsize=20 )
    ax.set_ylabel( r'$f(p) \, [\rm cm^{-3}]$', fontsize=20 )
    ax.legend( framealpha=0.1, title=r'$\alpha$')
    ylim = list(ax.get_ylim())
    ylim[0] = 1e-25
    ax.set_ylim( ylim )
    ax.tick_params( axis='both', pad=5, direction='in', labelsize=15 )

    fig.tight_layout()

    fig.savefig( fn_out, figsize=(5,5) )


def my_plot2():

    rad_spec_file = sys.argv[1]
    ele_spec_file = sys.argv[2]
    fn_out = sys.argv[3]

    dat_r = np.loadtxt( rad_spec_file )
    dat_e = np.loadtxt( ele_spec_file )

    fig, axs= plt.subplots( 1, 2, figsize=(2*4,4) )

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
    ii5 = -1
    #ii5 = np.where( t == t.min() )[0][0]

    print( "ii3: %i [%g], ii5: %i [%g]"%(ii3, p[ii3], ii5, p[ii5] ) )

    ele_alphas = []
    rad_alphas = []

    for i in range( 1, m ):

        M = dat_r[ i, 1 ]
        index = dat_r[ i, 0 ]

        P = dat_r[ i, 2: ] / mycc.mJy
        F = dat_e[i, 2:]
#
        t = "%.2e"%M
        t = t.split( 'e' )
        if t[1][0] == '+':
            t[1] = t[1][1:]

        ss = sss[ i % len(sss) ]

        #label=r"$[%i]\ %s\times 10^{%s}\,M_\odot$"%(index, t[0], t[1] )
        if ( i == 1 ):
            ss = 'k-'

        #print( F[ii3:ii5] )
        #print( p[ii3:ii5] )
        y =  F[ii3:ii5]
        x =  p[ii3:ii5]
        x = x[ y>0 ]
        y = y[ y>0 ]
        y = np.log10( y )
        x = np.log10( x )
        r = curve_fit( fit_f, x, y )
        #print( 'ele: ', r[0], r[1] )
        alpha_e = r[0][0]
        ele_alphas.append( -alpha_e )

        y = np.log10( P[i300:i1400] )
        x = np.log10( v[i300:i1400] )
        r = curve_fit( fit_f, x, y )
        #print( 'rad: ', i, r[0], r[1] )
        alpha_r = r[0][0]
        rad_alphas.append( -alpha_r )

        print( "[alpha] ele: %.2f, rad: %.2f"%(-alpha_e, -alpha_r) )

        if ( i > 10 ):
            continue

        label_e = r'$G_{%i}:\,%.2f$'%(i, -alpha_e)
        label_r = r'$G_{%i}:\,%.2f$'%(i, -alpha_r)
        axs[0].loglog( p, F, ss, label=label_e )
        xlim = list(axs[0].get_xlim())
        xlim[0] = 0.3 
        axs[0].set_xlim(xlim)
        axs[1].loglog( v, P, ss, label=label_r )
        #print( axs[1].get_xlim() )


    for i in range( 2 ):
        axs[i].minorticks_off()
        #axs[0].set_xscale( 'log' )
        #axs[i].locator_params( numticks=10 )

        if i == 0:
            axs[i].set_xlabel( r'$p$', fontsize=20 )
            axs[i].set_ylabel( r'$f(p) \, \rm [cm^{-3}]$', fontsize=20 )
            axs[i].legend( framealpha=0.1, title=r'$\alpha$')
            #ylim = list(axs[i].get_ylim())
            #ylim[0] = 1e-25
            #axs[i].set_ylim( ylim )
            #axs[i].set_xlim( [1, 1e7] )
        else:
            axs[i].set_xlabel( r'$\nu \, [\rm MHz]$', fontsize=20 )
            axs[i].set_ylabel( r'$S \, [\rm mJy]$', fontsize=20 )
            axs[i].legend( framealpha=0.1, title=r'$\alpha_{\rm rad}$' )

        axs[i].tick_params( axis='both', pad=5, direction='in', labelsize=10 )

    if ( len(sys.argv) < 2 ):
        file_pre = ''
    else:
        file_pre = sys.argv[1] + '-'

    fig.tight_layout()

    fig.savefig( fn_out, figsize=(5,5) )

    #return
    #fig, axs = plt.subplots( 1, 2 )

    #axs[0].hist( ele_alphas )
    #axs[1].hist( rad_alphas )

    #axs[0].set_xlabel( r'$\alpha$' )
    #axs[0].set_ylabel( r'$Number$' )

    #axs[1].set_xlabel( r'$\alpha_{\rm rad}$' )
    #axs[1].set_ylabel( r'$Number$' )

    #fig.tight_layout()
    #fig.savefig( figs_dir + file_pre + 'ele_rad_spec_index.pdf', figsize=(8,8) )


my_plot2()
#my_plot0()
