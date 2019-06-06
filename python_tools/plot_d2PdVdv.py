#!/usr/bin/env python3

from my_work_env import *

d = np.loadtxt( output_dir + './d2PdVdv_qmax.dat' )

m,n = d.shape
nu = d[ 0, 1: ]

lss = [ '-', '-.', '--' ]
lssN = len( lss )

for i in range(1,m):

    t = d[i,:]

    qmax = t[0]
    t = t[1:]

    qmax = "%.1e"%qmax
    qmax = qmax.split( 'e' )
    qmax[1] = qmax[1][1:]
    while( qmax[1][0] == '0' ):
        qmax[1] = qmax[1][1:]
    print( qmax )

    #plt.loglog( nu/1e6, t, '.', label=r"$ p_{max}: %s \times 10^{%s}$"%(qmax[0], qmax[1]), markersize=1 );
    #plt.loglog( nu/1e6, t, ls=lss[i%lssN], label=r"$ p_{max}: %s \times 10^{%s}$"%(qmax[0], qmax[1]) );
    #plt.loglog( nu/1e6, t, ls=lss[i%lssN], label=r"$ p_{\rm max}: %s \times 10^{%s}$"%(qmax[0], qmax[1]) );
    plt.loglog( nu/1e6, t, label=r"$ p_{\rm max}= %s \times 10^{%s}$"%(qmax[0], qmax[1]) );

plt.xlabel( r'$\nu \; [\rm MHz]$', fontsize=20 )
plt.ylabel( r'$\frac{d^2P(\nu)}{dVd\nu}[\rm erg\, s^{-1}\, cm^{-3}\, Hz^{-1}]$', fontsize=20 )
#plt.grid()

ax = plt.gca()
#ax.margins( 0.5, 0.5 )

ax.tick_params( axis='both', direction='in', labelsize=20, pad=5 )
ax.minorticks_off()

plt.legend(prop={'size':13}, framealpha=0.1)
plt.tight_layout()
plt.savefig( output_dir + './d2PdVdv_qmax.pdf', figsize=(5,5) )

