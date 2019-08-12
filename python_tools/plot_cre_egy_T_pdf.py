#!/usr/bin/env python3

from my_work_env import *

ds = [ np.loadtxt( f ) for f in sys.argv[1:3] ]
hs = []
for i in range(2):
    hs.append( ds[i][0,:] )
    ds[i] = ds[i][1:,:]
fn_out = sys.argv[3]

fs = 5 
x_pad = 0.05
y_pad = 0.1
t_cbar = 0.03

fig = plt.figure( figsize=( fs/(0.5*(1-2*x_pad-t_cbar)), fs/(1-2*y_pad) ) )
dx = (1 - 2*x_pad-t_cbar) / 2
dy = 1 - 2*y_pad
ax = [ fig.add_axes( [ x_pad+i*dx, y_pad, dx, dy ] ) for i in range(2) ]
ax_cbar = fig.add_axes( [ x_pad+2*dx, y_pad, t_cbar, dy ] )

for i in range(2):
    ax[i].set_xticks( [] )
    ax[i].set_yticks( [] )

vmin = np.min( [ ds[0][ds[0]>0].min(), ds[1][ds[1]>0].min() ] )
vmax = np.max( [ ds[0].max(), ds[1].max() ] )
m, n = ds[0].shape
for i in range(2):
    ds[i][m-1,0] = vmin
    ds[i][m-1,1] = vmax

Nmin = 0 
for i in range(2):
    ds[i][ ds[i]<Nmin ] = 0
    img = ax[i].imshow( ds[i], norm=mplc.LogNorm(), cmap=cm.jet )
    ax[i].invert_yaxis()
    ax[i].grid()
    xloc, xfmt = make_log_ticks( hs[i][1], hs[i][2], n, a=2 )
    yloc, yfmt = make_log_ticks( hs[i][3], hs[i][4], m, a=1 )
    #xloc = tik.FixedLocator(xloc)
    #yloc = tik.FixedLocator(yloc)
    #xfmt = tik.FixedFormatter(xfmt)
    #yfmt = tik.FixedFormatter(yfmt)

    ax[i].xaxis.set_major_locator( xloc )
    ax[i].xaxis.set_major_formatter( xfmt )

    ax[i].yaxis.set_major_locator( yloc )

    if i == 0:
        plt.colorbar( img, cax=ax_cbar )
        ax[i].yaxis.set_major_formatter( yfmt )
    else:
        remove_tick_labels( ax[i], 'y', bar=1 )

for i in range(2):
    ax[i].tick_params( axis='both', direction='in', labelsize=15, pad=5 )
ax_cbar.tick_params( axis='both', direction='in', labelsize=15, pad=5 )

fig.savefig( fn_out )


