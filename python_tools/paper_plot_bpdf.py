#!/usr/bin/env python3

from my_work_env import *

d = np.loadtxt( sys.argv[1] )
fn_out = sys.argv[2]

ep = d[:,0]
ek = d[:,1]
disp = d[:,2]
I = d[:,3] * 1e-7

plt.plot( ek/ep, I, '*' )
for i in range( len(ep) ):
    a = 1.01
    b = 1.01
    if i == 25:
        a = 0.99
        b = 0.99
    plt.text( ek[i]/ep[i]*a, I[i]*b, r'${%i}$'%i )

plt.xlabel( r'$E_k/E_p$' )
plt.ylabel( r'$P_{\rm 1.4GHz} \, [W Hz^{-1}]$' )
plt.yscale( 'log' )
plt.savefig( fn_out )
