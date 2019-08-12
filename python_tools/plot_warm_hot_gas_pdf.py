#!/usr/bin/env python3

from my_work_env import *

fns = sys.argv[1:5]
ds = [ np.loadtxt( f )[3:,:] for f in fns ]
fn_out = sys.argv[5]

t_pad = 0.2
t_diff = 0.3
fs = 6
fig = plt.figure( figsize=( fs, fs ) )
dx = 1-t_pad
dy = 1-t_pad-t_diff
ax      = fig.add_axes( [t_pad/2, t_pad/2+t_diff, dx, dy ] )
ax_diff = fig.add_axes( [t_pad/2, t_pad/2, dx, t_diff ] )

labels = [\
r'$\rm z: 2$', r'$\rm z:2 \, CRE$',\
r'$\rm z: 0$', r'$\rm z:0 \, CRE$'\
]

ls = [ 'b-', 'b--', 'k-', 'k--' ]
ls2 = [ '.', '*', '.', '*' ]

Nmin = 10000
Yunit = 1e5
Yunit_str =  r'\times 10^5'

for i in range(4):
    d = ds[i]
    x = d[:,0]
    y = d[:,2]
    index = y > Nmin
    x = x[index]
    y = y[index] / Yunit
    ax.plot( x, y, ls[i], label=labels[i] )
    #ax.plot( x, y, ls2[i] )

ax.set_xscale( 'log' )
ax.minorticks_off()
fig.canvas.draw_idle()
remove_tick_labels( ax, 'x' )
ax.set_ylabel( r'$\rm dN/dlog_{10} \rho \, [%s]$'%Yunit_str )

labels = [\
r'$\rm z: 2$',\
r'$\rm z: 0$'
]
ls = [ 'b-.', 'k-.' ]

for i in range(2):
    x = ds[i*2][:,0]
    y = ds[i*2][:,2]
    y_cre = ds[i*2+1][:,2]
    index1 = y > Nmin 
    index2 = y_cre > Nmin 
    index = index1 * index2
    x = x[index]
    y = y[index]
    y_cre = y_cre[index]
    y_diff = ( y_cre-y ) / y
    ax_diff.plot( x, y_diff, ls[i], label=labels[i] )
#ax_diff.set_yscale( 'log' )
ax_diff.set_xlabel( r'$\rho/\bar{\rho}$' )
ax_diff.set_xscale( 'log' )

axs = [ ax, ax_diff ]

for a in axs:
    a.legend()
    set_tick_params( a, 10 )

fig.savefig( fn_out )
