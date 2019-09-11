#!/usr/bin/env python3

from my_work_env import *

fn1= sys.argv[1]
fn1_cre =sys.argv[2]
fn2= sys.argv[3]
fn2_cre =sys.argv[4]
fn_out = sys.argv[5]
fn_all = [ fn1, fn1_cre, fn2, fn2_cre ]

fs = 6 
t_pad = 0.15
t_err = 0.3
fig =  plt.figure( figsize=(fs*(2+3*t_pad), fs*(1+1*t_err+2*t_pad)) )
dx = 1 / ( 2 + 3*t_pad )
dy = 1 / ( 1 + 1*t_err + 2*t_pad )
axs = [ fig.add_axes([ dx*t_pad+i*dx*(1+t_pad), dy*(t_pad + t_err),\
    dx, dy]) for i in range(2) ]
axs_err = [ fig.add_axes([ dx*t_pad+i*dx*(1+t_pad), dy*t_pad,\
    dx, dy*t_err]) for i in range(2) ]
axs_all = axs + axs_err

ds = []
for f in fn_all:
    print( 'load %s ...'%f )
    ds.append( np.loadtxt(f) )

cs = [ 'k', 'b', 'y', 'r' ]
names = [ 'Diffuse-Condense', 'Warm-hot', 'Hot' ]

Nmin = 1e4
for i in range(2):
    if ds[i*2][0,0] != ds[i*2+1][0,0]:
        print( 'error' )
        exit()
    d = ds[i*2][1:,2:]
    d_cre =ds[i*2+1][1:,2:]
    z = ds[i*2][0,0]
    ax = axs[i]
    ax_err = axs_err[i]
    xmin = ds[i*2][1:,0].min()
    xmax = ds[i*2][1:,0].max()
    for j in range(3):
        x = ds[i*2][1:,0]
        y = d[:,j]
        y_cre = d_cre[:,j]

        if not(y.max() > 0 and y_cre.max()>0):
            continue

        index1 = y > 0 
        index2 = y_cre > 0 
        index = index1 * index2
        x = x[index]
        y = y[index]
        y_cre = y_cre[index]

        ax.loglog( x, y, '-', color=cs[j], label=names[j] )
        ax.loglog( x, y, '.', color=cs[j] )
        ax.loglog( x, y_cre, '--', color=cs[j], label='CRE-' + names[j] )
        ax.loglog( x, y_cre, '.', color=cs[j] )

        y_err = (y_cre-y) / y
        idx1 = y > Nmin 
        idx2 = y_cre > Nmin 
        idx = idx1 * idx2
        x = x[idx]
        y_err = y_err[idx]
        print( y_err )

        ax_err.plot( x, y_err, '-', color=cs[j], label= names[j] )
        ax_err.plot( x, y_err, '.', color=cs[j])
        ax_err.set_xscale( 'log' )
        #ax_err.set_yscale( 'symlog', linthreshy=0.001)

    remove_tick_labels( ax, 'x' )
    set_tick_params( ax, 15 )
    set_tick_params( ax_err, 15 )
    ax.legend( prop={'size':15}, framealpha=0.1 )
    ax.axhline( Nmin, color='r' )
    ax_err.legend( prop={'size':13}, framealpha=0.1 )
    ax_err.set_xlabel( r'$\rm \rho/\bar{\rho}$', fontsize=20 )
    ax.set_ylabel( r'$\rm dN/dlog(\rho)$', fontsize=20 )
    ax.set_xlim( [xmin, xmax] )
    ax_err.set_xlim( [xmin, xmax] )
    ax_err.set_ylabel( 'Error' )

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    print( xlim )
    print( ylim )
    fx = xlim[0] * np.exp(np.log(xlim[1]/xlim[0]) * 0.01)
    fy = ylim[0] * np.exp(np.log(ylim[1]/ylim[0]) * 0.85)
    ax.text( fx, fy, r"$\rm z: %.0f$"%z,fontsize=30, color='g' )

fig.savefig( fn_out )
