#!/usr/bin/env python3

from my_work_env import *

gr = np.loadtxt( sys.argv[1] )
fn_out = sys.argv[3]
cs = [ 'k', 'r', 'b', 'y' ]
ls = [ '-', '-.' ]
ns = [ 'diffuse_cool', 'diffuse_warm', 'diffuse_hot', 'dense' ]

z = gr[:,0]
for j in range(4):
    plt.plot( z, gr[:,j+1], cs[j], label=ls[j] )

plt.xlim( [0, 4] )
plt.yscale( 'log' )

plt.savefig( fn_out )

