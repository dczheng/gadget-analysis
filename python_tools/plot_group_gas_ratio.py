#!/usr/bin/env python3

from my_work_env import *

fns = sys.argv[1:5]
fn_out = sys.argv[5]

ds = []
for f in fns:
    print( 'load %s ...'%f )
    ds.append( np.loadtxt( f ) )

fs = 5
x_pad = 0.05
y_pad = 0.1
fig = plt.figure( figsize=( 2*fs/(1-3*x_pad), fs/(1-2*y_pad) ) )
dx = ( 1-3*x_pad ) / 2
dy = 1-2*y_pad
axs = [ fig.add_axes( [ x_pad + i*(dx+x_pad), y_pad, dx, dy ] ) for i in range(2) ]

labels1 = [ '1', '2' ]
labels2 = [ 'r_gas', 'r_hot', 'r_warmhot', 'r_condensed', 'r_diffuse' ]
bins = np.linspace( 0.0, 0.12, 10 )
for i in range(2,4):
    ii = i // 2
    ax = axs[ii]
    d = ds[i]
    for j in range(1):
        n = d[:,j+1]
        ax.hist( n, histtype='step', bins=bins, label=labels1[i%2] )

for i in range(2):
    ax.legend()

fig.savefig( fn_out )

