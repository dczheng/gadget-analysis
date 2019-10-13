#!/usr/bin/env python3

from my_work_env import *
import pandas as pd

df = pd.read_csv( sys.argv[1] )
fn_out = sys.argv[2]

M = df['mtot'].to_numpy()
vr = df['vr'].to_numpy()
I =  df['Lum'].to_numpy() * 1e-7
index = (np.abs(vr)<1)
index = (np.abs(vr)<1) * ( M>5e3 )

#vr = np.abs( vr+0.5 )
M  = M[index] 
vr = vr[index]
I = I[index]

#I = np.log10( I )
#plt.plot( M, I, '*' )
plt.plot( vr, I, '*' )
#plt.hist2d( vr, I, bins=15 )
#plt.xscale( 'log' )
plt.yscale( 'log' )
for i in range( len(vr) ):
    a = 1.01
    b = 1.01
    if i == 25:
        a = 0.99
        b = 0.99
    #plt.text( vr[i]*a, I[i]*b, r'${%i}$'%i )

plt.xlabel( r'$E_k/E_p$' )
plt.ylabel( r'$P_{\rm 1.4GHz} \, [W Hz^{-1}]$' )
plt.savefig( fn_out )
