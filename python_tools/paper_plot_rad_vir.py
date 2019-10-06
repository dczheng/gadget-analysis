#!/usr/bin/env python3

from my_work_env import *

d = np.loadtxt( sys.argv[1] )
fn_out = sys.argv[2]

ep = d[:,0]
ek = d[:,1]
disp = d[:,2]
I = d[:,3] * 1e-7

plt.plot( np.abs(ek/ep), I, '*' )
plt.xlabel( r'$\rm T/U$' )
plt.ylabel( r'$\rm I_{1.4G} \, [W Hz^{-1}]$' )
plt.yscale( 'log' )
plt.savefig( fn_out )
