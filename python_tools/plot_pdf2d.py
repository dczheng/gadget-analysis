#!/usr/bin/env python3

from my_work_env import *

fn = sys.argv[1]
fn_out = sys.argv[2]

if len( sys.argv ) < 4:
    mynorm = mplc.LogNorm()
else:
    mynorm = None

d = np.loadtxt( fn )
h = d[ 0, : ]
d = d[ 1:, : ]
m, n = d.shape

z = h[0]
xmin = h[1]
xmax = h[2]
if ( h[3] ):
    xmin = np.power( 10, xmin )
    xmax = np.power( 10, xmax )
ymin = h[4]
ymax = h[5]
if ( h[6] ):
    ymin = np.power( 10, ymin )
    ymax = np.power( 10, ymax )

print( "xmin: %g, xmax: %g, ymin: %g, ymax: %g"%\
       ( xmin, xmax, ymin, ymax ) )


t_cbar = 0.05
t_pad = 0.1
fs = 5
fig = plt.figure( figsize=( fs/(1-t_cbar-2*t_pad), fs/( 1-2*t_pad ) ) )

ax = fig.add_axes( [ t_pad, t_pad, 1-t_cbar-2*t_pad, 1-2*t_pad ] )
ax_cbar = fig.add_axes( [ 1-t_cbar-t_pad, t_pad, t_cbar, 1-2*t_pad ] )

img = ax.imshow( d, norm=mynorm, cmap=cm.jet )
plt.colorbar( img, cax=ax_cbar )

ax.grid()
ax.invert_yaxis()
make_log_ticks( xmin, xmax, n, axis=ax.xaxis )
make_log_ticks( ymin, ymax, n, axis=ax.yaxis )
set_tick_params( ax, 15 )
set_tick_params( ax_cbar, 15 )

fig.savefig( fn_out )
