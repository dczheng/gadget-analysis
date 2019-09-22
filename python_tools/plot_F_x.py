#!/usr/bin/env python3

from my_work_env import *

d = np.loadtxt(  sys.argv[1] )

#plt.loglog( d[:,0], d[:,1] )
x = d[:,0]
F = d[:,1]
plt.plot( d[:,0], d[:,1], 'b' )
#plt.loglog( d[:,0], d[:,1], '.' )

#plt.xscale( 'log' )
#plt.yscale( 'log' )

plt.xlabel( 'x', fontsize=20 )
plt.xlim( [0, 4] )
ylim = plt.ylim()
ylim = list( ylim )
ylim[0] = 0
plt.ylim( ylim )

plt.ylabel( r'$F(x)$', fontsize=20 )


x0 = x[ F == F.max() ][0]
print( x0 )

plt.axvline( x = x0, linestyle='-.', color = 'r', label="x=%.2f"%x0 )
ax = plt.gca()
ax.xaxis.set_major_locator( plt.LinearLocator(numticks=6) )
ax.tick_params( axis='both', direction='in', labelsize=20, pad=5, length=0 )
#plt.axis( 'equal' )

plt.legend(prop={'size':20}, framealpha=0.1)
#plt.text( x0-0.15,0-0.055, "%.2f"%x0 )

#plt.grid()
plt.tight_layout()

plt.savefig( 'F_x.pdf', figsize=(4,4) )
