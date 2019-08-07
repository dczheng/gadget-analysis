#!/usr/bin/env python3

from my_work_env import *
import matplotlib.ticker as tik

#matplotlib.style.use( 'default' )

ss =  [ '', '_condensed', '_diffuse', '_warm-hot', '_hot' ]
ds = []
for s in ss:
    fn = "%s%s_%03i.dat"%( sys.argv[1], s, int(sys.argv[2] ) )
    print( 'load %s ...'%fn )
    ds.append( np.loadtxt(fn)[1:,:] )

fn_out = "%s_%03i.png"%( sys.argv[1], int(sys.argv[2]) )

fs = 5 
fig = plt.figure( figsize=( 3*fs, 2*fs) )
dx = 1 / 3.0
dy = 1 / 2.0
rbar = 0.05

axs = [ fig.add_axes( [0, (1-dy-rbar)/2+rbar, dx, dy] ) ]
ax_cbar = fig.add_axes( [0, (1-dy-rbar)/2, dx, rbar] )

for i in range(4):
    ii = i // 2 + 1
    jj = i %  2
    axs.append( fig.add_axes( [ ii*dx, jj*dy, dx, dy ] ) )

for a in axs:
    a.set_xticks( [] )
    a.set_yticks( [] )

vmin = 1e10
vmax = -vmin 
for d in ds:
    if not( d.max() > 0 ):
        continue
    vmin = np.min( [d[d>0].min(), vmin] )
    vmax = np.max( [d.max(), vmax] )
print( vmin, vmax )

m, n = ds[0].shape
for i in range( 5 ):
    d = ds[i]
    ax = axs[i]
    d[0,0] = vmin
    d[0,1] = vmax
    d[ d<vmax * 1e-7 ] = vmax * 1e-8
    img = ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.plasma )
    if i == 0:
        plt.colorbar( img, cax=ax_cbar, orientation='horizontal' )
        ax_cbar.minorticks_off()
        ax_cbar.tick_params( direction='in', labelsize=15 )
    else:
        ax.text( m*0.05, n*0.05, ss[i][1:], color='y', size=30 )

    ax.invert_yaxis()

fig.savefig( fn_out, dpi=300 )
