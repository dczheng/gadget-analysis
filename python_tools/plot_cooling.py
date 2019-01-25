#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

cool = np.loadtxt( './cooling.dat' )
qmax_B = np.loadtxt( './qmax_B.dat' )
qmin_rho = np.loadtxt( './qmin_rho.dat' )

N = 1
fig, axs = plt.subplots( 1, N, figsize=(N*4,4) )
#ax = axs[0]
ax = axs


q = cool[ 0, 4: ]
f = cool[ 1, 4: ]
df = cool[ 2, 4: ]
ff = cool[ 3, 4: ]
qmin = cool[3,0]
fqmin = cool[3,1]
qmax = cool[3,2]
fqmax = cool[3,3]

ax.plot( q, f, label=r'$f_t$', linewidth=1 )
ax.plot( q, -df, '-.', label=r'$-df_t$', linewidth=1 )
ax.plot( q, f+df, '--', label=r'$f_t+df_t$', linewidth=1 )
ax.plot( q, ff, '.', label=r'$f_{\rm{mod}, t+\Delta t}$', markersize=1 )
ax.plot( [qmin], [fqmin], '*', markersize=5, label=r'$ p_{\rm{min}}$' )
ax.plot( [qmax], [fqmax], '*', markersize=5, label=r'$ p_{\rm{max}}$' )
ax.set_xscale( 'log' )
ax.set_yscale( 'log' )
ax.set_xlabel( r'$p$', fontsize=20 )
#ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=3))
ax.yaxis.set_major_locator(plt.LogLocator(numticks=5))
ax.xaxis.set_major_locator(plt.LogLocator(numticks=7))
ax.tick_params( axis='both', labelsize=13, direction='in', pad=5 )
ax.legend(prop={'size':13}, framealpha=0.1)
ax.minorticks_off()
ax.grid()

'''
ax = axs[1]
B = qmax_B[:,0]*1e6
qmax = qmax_B[:,1]
ax.plot( B, qmax )
ax.set_xscale( 'log' )
ax.set_yscale( 'log' )
ax.set_xlabel( r'$B [\mu \rm G]$', fontsize=20 )
ax.set_ylabel( r'$p_{\rm max}$', fontsize=20 )
ax.tick_params( axis='both', labelsize=20, direction='in', pad=5 )
ax.minorticks_off()
ax.grid()

ax = axs[2]
rho= qmin_rho[:,0]
qmin = qmin_rho[:,1]
ax.loglog( rho, qmin )
#ax.set_xlabel( r'$\frac{\rho}{\rho_{\rm bar}}$', fontsize=20 )
ax.set_xlabel( r'${\rho}/{\rho_{\rm bar}}$', fontsize=20 )
ax.set_ylabel( r'$p_{\rm min}$', fontsize=20 )
ax.xaxis.set_major_locator(plt.LogLocator(numticks=6))
ax.tick_params( axis='both', labelsize=20, direction='in', pad=5 )
ax.minorticks_off()
ax.grid()
'''

plt.tight_layout()
plt.savefig( './cooling.pdf' )

exit()

data = np.loadtxt( './cooling_test.dat' )
q = data[0, 6:]
f = data[1, 6:]

data = data[2:, :]

m, n = data.shape

mm = m // 2

fig = plt.figure()

for i in range( 4 ):

    df_i = i * 2
    ff_i = i * 2 + 1

    rho   = data[ff_i, 0]
    B     = data[ff_i, 1]
    qmin  = data[ff_i, 2]
    fqmin = data[ff_i, 3]
    qmax  = data[ff_i, 4]
    fqmax = data[ff_i, 5]

    df    = data[df_i, 6:]
    ff    = data[ff_i, 6:]

    ax = fig.add_subplot( 220 + i + 1 )

    ax.plot( q, f, label=r'$f_t$', linewidth=1 )
    ax.plot( q, -df, '-.', label=r'$-df_t$', linewidth=1 )
    ax.plot( q, f+df, '--', label=r'$f_t+df_t$', linewidth=1 )
    ax.plot( q, ff, '.', label=r'$f_{t+\Delta t}$', markersize=1 )
    ax.plot( [qmin], [fqmin], '*', markersize=5, label=r'$ p_{min}$' )
    ax.plot( [qmax], [fqmax], '*', markersize=5, label=r'$ p_{max}$' )

    rho = "%.2e"%rho
    rho = rho.split( 'e' )
    #print( rho )

    if rho[1][0] == '-':
        t = '-'
    else:
        t = ''

    a = 1
    while( rho[1][a] == '0' ):
        a += 1
    t = t + rho[1][a:]
    rho[1] = t


    '''
    B = "%.2e"%B
    B = B.split( 'e' )
    if ( B[1][0] == '-' ):
        t = '-'
    else:
        t = ''

    a = 1
    while( B[1][a] == '0' ):
        a += 1
    t = t + B[1][a:]
    B[1] = t

    ax.set_title( r'$\rho: 10^{%s}\rho_{bar}\quad B: %s \times 10^{%s}G$'%( rho[1], B[0], B[1] ), fontsize=5 )
    '''
    ax.set_title( r'$\rho: 10^{%s}\rho_{bar}\quad B: %.0f \mu G$'%( rho[1], B*1e6 ), fontsize=5 )

    if ( i == 3 or i == 2 ):
        ax.set_xlabel( 'q', fontsize=6 )


    '''
    tik = ax.xaxis.get_major_ticks()
    for t in tik:
        t.label.set_fontsize(4)

    tik = ax.yaxis.get_major_ticks()
    for t in tik:
        t.label.set_fontsize(4)
    '''

    ax.legend( fontsize=5 )
    ax.grid()


