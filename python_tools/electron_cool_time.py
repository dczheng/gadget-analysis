#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

def t_cool( g_l, z ):
    return 2.32e3 / g_l / np.power( z+1, 4 )


z = [ 0.0, 0.2, 0.4, 0.5, 1.0, 2.0 ]
gamma = np.linspace( 20, 60, 1000 )

fig = plt.figure()
ax_syn = fig.add_subplot( 221 )
ax_ion = fig.add_subplot( 223 )

for i in z:
    f = lambda x: t_cool( x, i )
    t_syn = f( gamma )
    ax_syn.plot( gamma, t_syn, label='z=%.1f'%(i) )


#ax_syn.set_xlabel( r'$\gamma_L$')
ax_syn.set_ylabel( 'Gyr' )
ax_syn.set_title( 'Synchroton' )
ax_syn.legend()

n = 1e6
t_ion = 0.511 * 1e6 * gamma / ( 7.64e-15 * n * ( 3 * np.log( gamma ) + 19.8 ) ) \
        / sc.year / 1e6
ax_ion.plot( gamma, t_ion, label='Ionisation' )
ax_ion.set_ylabel( 'Myr' )
ax_ion.set_xlabel( r'$\gamma$' )
ax_ion.legend()



#plt.show()
plt.savefig( './electron_cool_time.png' )

