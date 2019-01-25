#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

#chi_list = [ 0.001, 0.01, 0.02, 0.03 ]
#chi_list = [ 0.01, 0.02, 0.03, 0.1 ]
#chi_list = [ 0.03, 0.1, 0.2, 0.3 ]
#chi_list = [ 0.3, 0.4, 0.5 ]
chi_list = [ 0.02, 0.03, 0.04, 0.05 ]
#chi_list = [ 0.02, 0.05, 0.1, 0.2 ]
#kchi_list = [ 0.03, 0.05, 0.1, 0.2 ]

zeta = 0.01
z = 0.09

#fig = plt.figure( figsize=[8,2] )

def f( x, a, b ):
    return x*a+b

fig = plt.figure( figsize=(8,8) )

for ci in range( len(chi_list) ):

    fn = "./hge128_data/hge128_%g_%g_spec_%g.dat"%(zeta,  chi_list[ci], z)
    data = np.loadtxt( fn )

    v = data[ 0, 2: ]

    m, n = data.shape

    #ax = fig.add_subplot( 140 + ci + 1 )
    ax = fig.add_subplot( 220 + ci + 1 )

    for i in range( 1, m ):

        m = data[ i, 1 ]

        index = data[ i, 0 ]
        p = data[ i, 2: ]

        if i % 2:
            ss = '--'
        else:
            ss = '-'

        t = "%.2e"%m
        t = t.split( 'e' )

        if t[1][0] == '+':
            t[1] = t[1][1:]

        pmax = p.max()
        p[ p<pmax*1e-8 ] = 0

        logv = np.log10( v )
        logp = np.log10( p )
        fit_params = curve_fit( f, logv, logp )
        #print( fit_params )

        a = ( logp[0] - logp[-1] ) / (logv[0]-logv[-1])
        print( a )

        #ax.loglog( v, p, ss, label=r"$[%i]\ %s\times 10^{%s}\,M\_\odot, \quad \alpha: %.2f$"%(index, t[0], t[1], fit_params[0][0] ), linewidth=0.6 )
        ax.loglog( v, p, ss, label=r"$\rm [%i]\ %s\times 10^{%s}\,M_\odot,\; \alpha=%.2f$"%(index, t[0], t[1], a ) )

    if ci == 0 or ci == 2:
        ax.set_ylabel( r'$I_{\nu}\;\rm [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1}]$', fontsize=10 )

    if ci > 1:
        ax.set_xlabel( r'$\nu \; \rm [MHz]$', fontsize=10 )

    tik = ax.xaxis.get_major_ticks()
    for t in tik:
        t.label.set_fontsize(5)

    tik = ax.yaxis.get_major_ticks()
    for t in tik:
        t.label.set_fontsize(5)

    ax.minorticks_off()
    ax.tick_params( axis='both', direction='in', pad=5, labelsize=15 )

    ax.grid()

    str_chi = str( chi_list[ci] )
    str_chi = str_chi.split( '.' )

    ax.legend( prop={'size':10}, framealpha=0.1 )
    ax.set_title( r'$\rm CRE50\_%s \;$'%(str_chi[1]), fontsize=10)

fig.tight_layout()
plt.savefig( "./hge128_out/chi_test.pdf" )
