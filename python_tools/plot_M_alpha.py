#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import matplotlib.pyplot as plt
import numpy as np

m = 1.5
a = ( 4 + m*m*5 ) / ( 2 * ( m*m -1 ) )
print( a )

output_dir = "/mnt/ddnfs/data_users/dczheng/simulation/paper_plot/"

data = np.loadtxt( output_dir + 'M_alpha.dat' )
M = data[:,0]
a = data[:,1]
egy = data[:,2]

index = M < 10
M = M[index]
a = a[index]
egy = egy[index]

print( "Mmin: %g, Mmax: %g"%( M.min(), M.max() ) )
print( "amin: %g, amax: %g"%( a.min(), a.max() ) )
print( "Emin: %g, Emax: %g"%( egy.min(), egy.max() ) )

print( "aN: %i"%(len(a)) )
b = a[ a>10 ]
print( "bN: %i"%(len(b)) )

fig, axs = plt.subplots( 3,1 )
ax = axs[0]
ax.scatter( M, a, s=1 )
#ax.set_yscale( 'log' )
#ax.set_xlim( [0, 10] )
#ax.set_ylim( [0, 10] )
ax.set_xlabel( r'$M$' )
ax.set_ylabel( r'$\alpha$' )

ax = axs[1]
ax.scatter( M, egy, s=1 )
#ax.set_xlim( [0, 10] )
#ax.set_ylim( [0, 10] )
ax.set_yscale( 'log' )

ax.set_xlabel( r'$M$' )
ax.set_ylabel( r'$E$' )

m = np.linspace( 2, 10 )
a = ( 4 + m*m*5 ) / ( 2 * ( m*m -1 ) )
ax = axs[2]
ax.scatter( m, a, s=1 )
ax.set_xlabel( r'$M$' )
ax.set_ylabel( r'$\alpha$' )


fig.savefig( output_dir + 'M_alpha.pdf' )

