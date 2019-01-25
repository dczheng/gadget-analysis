#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

qmin_chi = np.loadtxt( './qmin_chi.dat' )
qmax_chi= np.loadtxt( './qmax_chi.dat' )

chi = qmin_chi[0, 1:]
qmin_chi = qmin_chi[1:, :]
m,n = qmin_chi.shape
for i in range( m ):
    rho = qmin_chi[i,0]
    print( "rho:", rho )
    qmin = qmin_chi[i,1:]
    if rho < 1:
        plt.plot( chi, qmin, label='%.1f'%rho )
    else:
        plt.plot( chi, qmin, label='%.0f'%rho )
plt.xscale( 'log' )

plt.xlabel( r'$\chi_{\rm coul}, \chi_{\rm syn+ic}$', fontsize=20 )
plt.ylabel( r'$p_{\rm min}$', fontsize=20 )
plt.minorticks_off()

ax = plt.gca()
ax.tick_params( axis='both', labelsize=20, direction='in', pad=5  )

ax.legend(title=r'$p_{\rm min}, \rho/\rho_{\rm bar}$', loc="upper left", prop={'size':15}, framealpha=0.1)
ax2 = plt.twinx( ax )

chi = qmax_chi[0, 1:]
qmax_chi = qmax_chi[1:, :]
m,n = qmax_chi.shape
for i in range( m ):
    B = qmax_chi[i,0]
    print( "B:", B )
    qmax = qmax_chi[i,1:]
    ax2.plot( chi, qmax, '--', label='%.0f'%B )
ax2.legend(title=r'$p_{\rm max}, B/B_{\rm cmb.0}$', loc="upper right", prop={'size':15}, framealpha=0.1)
ax2.set_yscale( 'log' )
ax2.set_ylabel( r'$p_{\rm max}$', fontsize=20 )
ax2.minorticks_off()
ax2.tick_params( axis='both', labelsize=20, direction='in', pad=5  )

plt.tight_layout()
plt.savefig( 'q_chi.pdf', figsize=(4,4) )
