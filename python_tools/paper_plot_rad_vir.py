#!/usr/bin/env python3

from my_work_env import *
import pandas as pd

df = pd.read_csv( sys.argv[1] )
fn_out = sys.argv[2]

M = df['mtot'].to_numpy()
vr = df['vr'].to_numpy()
ep = df['ep'].to_numpy()
ek = df['ek'].to_numpy()
ee = df['ee'].to_numpy()
B = df['B'].to_numpy()
I =  df['Lum'].to_numpy() * 1e-7
index = (np.abs(vr)<1)
index = (np.abs(vr)<1) * ( M>1e4 )

#vr = np.abs( vr+0.5 )
M  = M[index] 
vr = vr[index]
I = I[index]
ee = ee[index]
ep = ep[index]
ek = ek[index]
B = B[index]

print( np.abs(ep)/ek )
print( ek/np.abs(ep) )
fig, axs = plt.subplots( 2, 1, figsize=(5,2*5) )
ax_e = axs[0]
#ax_B = axs[1]
ax_I = axs[1]
#I = np.log10( I )
#plt.plot( M, I, '*' )
#vr = np.abs(ep) / ek
#idx = (I>1e20) 
#I =  I[idx]
#vr = vr[idx]
#ee = ee[idx]
vr = np.abs( vr+0.5 )
ax_e.plot( vr, ee, '*' )
#ax_B.plot( vr, B, '*' )
ax_I.plot( vr, I, '*' )

for i in range( len(vr) ):
    a = 1.01
    b = 1.01
    if i == 25:
        a = 0.99
        b = 0.99
    ax_I.text( vr[i]*a, I[i]*b, r'${%i}$'%i )

xlabels = [\
r'$E_k/E_p$',\
r'$E_k/E_p$',\
]
ylabels = [\
r'$\epsilon_{\rm CRE}/\epsilon$',\
r'$P_{\rm 1.4GHz} \, [W Hz^{-1}]$',\
]

ax_I.set_yscale( 'log' )
for i in range(2):
    axs[i].set_xlabel( xlabels[i] )
    axs[i].set_ylabel( ylabels[i] )

fig.savefig( fn_out )
