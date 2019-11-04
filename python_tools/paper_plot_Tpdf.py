#!/usr/bin/env python3

from my_work_env import *

fn      = sys.argv[1]
fn_cre  = sys.argv[2] 
fn_out  = sys.argv[3]

d      = np.loadtxt( fn )
d_cre  = np.loadtxt( fn_cre )

T       = d[1:,0]
Nfof    = d[1:,1]
N       = d[1:,2]
Nfof4   = d[1:,3]
Nfof3   = d[1:,4]

T_cre       = d_cre[1:,0]
Nfof_cre    = d_cre[1:,1]
N_cre       = d_cre[1:,2]
Nfof4_cre   = d_cre[1:,3]
Nfof3_cre   = d_cre[1:,4]


if (T-T_cre).sum() != 0:
    print( "error" )
    exit()

fs = 5
t_pad = 0.1
t_err = 0.3
dx = 1 - 2 * t_pad
dy = 1 - 2 * t_pad - t_err
fig = plt.figure( figsize=(fs, fs) )
ax      = fig.add_axes( [ 1.5*t_pad, 1.2*t_pad+t_err, dx, dy ] )
ax_err  = fig.add_axes( [ 1.5*t_pad, 1.2*t_pad, dx, t_err ] )

labels = [ '', 'no', 'fof',\
                #'> 10^{14} \, h^{-1}  M_\odot',\
                #'10^{12} \sim 10^{13} \, h^{-1} \, M_\odot',\
                '> 10^{14}',\
                '10^{12} \sim 10^{13}',\
                ]
cs     = [ '', 'k', 'y', 'b', 'r' ]
#cols = [ 1, 2, 3, 4 ]
cols = [ 3, 4 ]
for i in cols:
    ax.plot( d[1:,0], d[1:,i], cs[i]+'-',\
        label=r'$%s\, [\rm SIM]$'%labels[i] )
    ax.plot( d_cre[1:,0], d_cre[1:,i], cs[i]+'--',\
        label=r'$%s \, [\rm SIM-CRE]$'%labels[i] )

for i in cols:
    idx = d[1:,i] > 0
    x = T[idx]
    y1 = d[1:,i][idx]
    y2 = d_cre[1:,i][idx]
    y  = (y2-y1) / y2
    ax_err.plot( x, y, cs[i],\
        label=r'$%s$'%labels[i] )

ax.set_xscale( 'log' )
ax.set_yscale( 'log' )
ax.legend( framealpha=0.1 )
xlim = ax.get_xlim()
remove_tick_labels(ax, 'x', )
set_tick_params( ax, 10 )
ax.set_ylabel( r'$\rm PDF$' )

ax_err.set_xscale( 'log' )
ax_err.set_xlim( xlim )
ax_err.legend( framealpha=0.1 )
set_tick_params( ax_err, 10 )

ax_err.set_xlabel( r'$T \, \rm [K]$' )
ax_err.set_ylabel( r'$\rm Difference$' )
ax_err.set_ylim( [-0.4, 0.3] )


fig.savefig( fn_out )
