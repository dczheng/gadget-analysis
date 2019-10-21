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

M = M[index]
vr = vr[index]
ep = ep[index]
ek = ek[index]
ee = ee[index]
I = I[index]

x = np.abs((vr+0.5))
#x = np.abs( ep / ek + 2 )
#x = vr
y = I
x = x[:25]
y = y[:25]
print( len(x) )

plt.plot( x, y, '*' )

xlabel = r'$|E_k/E_p+0.5|$'
ylabel = r'$P_{\rm 1.4GHz} \, [W Hz^{-1}]$'

plt.xlabel( xlabel )
plt.ylabel( ylabel )

plt.yscale( 'log' )

plt.savefig( fn_out, figsize=(5,5) )
