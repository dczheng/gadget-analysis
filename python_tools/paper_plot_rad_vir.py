#!/usr/bin/env python3

from my_work_env import *
import pandas as pd

df = pd.read_csv( sys.argv[1] )
fn_out = sys.argv[2]

M = df['mtot'].to_numpy()
vr = df['vr'].to_numpy().astype( 'float' )
ep = df['ep'].to_numpy()
ek = df['ek'].to_numpy()
ee = df['ee'].to_numpy()
B = df['B'].to_numpy()
I =  df['Lum'].to_numpy() * 1e-7
idxs = np.arange(len(M))
print( vr )
#index = (np.abs(vr)<1)
#index = (np.abs(vr)<1) * ( M>1e4 )
#index = (M/0.68>1e4) * (I>5e23)
index = (M>1e4) * (I>5e23)
#index = (I>1e24)
#index = (M>1e3)

M = M[index]
vr = vr[index]
ep = ep[index]
ek = ek[index]
ee = ee[index]
I = I[index]
idxs = idxs[index]

x = np.abs((vr+0.5))
#x = vr
#x = np.abs( ep / ek + 2 )
#x = vr
y = I
#x = x[:25]
#y = y[:25]
#index = [ True ] * len(x)
#index[20] = False
#x = x[index]
#y = y[index]

print( len(x) )
for i in range(len(x)):
    plt.text( x[i]*1.01, y[i]*1.02, '%.2f'%(M[i]/1e4) )

plt.plot( x, y, '*' )

xlabel = r'$|E_k/E_p+0.5|$'
ylabel = r'$P_{\rm 1.4GHz} \, [\rm W Hz^{-1}]$'

plt.xlabel( xlabel )
plt.ylabel( ylabel )


plt.yscale( 'log' )
#plt.xscale( 'log' )
plt.xlim( [0,0.18] )


ax = plt.gca()
set_tick_params( ax, 15 )
#ax.axis( 'equal' )

plt.savefig( fn_out, figsize=(5,5) )
