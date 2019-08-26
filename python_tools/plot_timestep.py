#!/usr/bin/env python3

from my_work_env import *

lines = os.popen( 'grep Sync-Point %s'%sys.argv[1] ).readlines()

z = []
ts = []

for l in lines:
    if "Sync-Point" in l:
        ll = l.split( ',' )
        z.append( float(ll[2].split()[1]) )
        ts.append( float(ll[3].split()[1]) )

z = np.array( z )
ts = np.array( ts )

plt.plot( z, ts )
ax = plt.gca()
ax.set_xlabel( 'z' )
ax.set_ylabel( 'timestep' )

#ax2 = ax.twiny()
#ax2.plot( 1/(1+z), ts, '.' )
#ax2.set_xlabel( 'a' )

plt.yscale( 'log' )
#plt.xscale( 'log' )
plt.xlim( [0,10] )
plt.savefig( 'timestep.png' )
