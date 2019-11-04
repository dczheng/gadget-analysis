#!/usr/bin/env python3

from my_work_env import *

fn = sys.argv[1]
fn_o = sys.argv[2]

d = np.loadtxt( fn )
x = d[1:,0]
y = d[1:,1]

#u = 1e5
#y = y / u
#plt.step( x, y )
plt.plot( x, y )
plt.xscale( 'log' )
#plt.yscale( 'log' )
fs = 20
plt.xlabel( r'$\epsilon/\epsilon_{\rm bar}$', fontsize=fs )
plt.ylabel( r'$\rm PDF$', fontsize=fs )
plt.xlim( [1e-8, 1e-1] )
ax = plt.gca()
set_tick_params( ax, 15 )
plt.savefig( fn_o )
