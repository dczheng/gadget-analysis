#!/usr/bin/env python3

from my_work_env import *

fn_dens_pre = sys.argv[1]
ss = ['x', 'y', 'z']
fns = [ "%s_%s.dat"%(fn_dens_pre, ss[i]) for i in range(3) ]
fns.append( sys.argv[2] )
fns.append( sys.argv[3] )
fn_out = sys.argv[4]

ds = []
for f in fns:
    print( "load %s ..."%f )
    ds.append( np.loadtxt( f ) )

for i in range(3):
    ds[i] = ds[i][1:,:]

s = 4
t_pad = 0.1
t_diff = 0.3
fig = plt.figure( figsize=(2*s, 2*s) )
dx = 1/2
dy = 1/2

axs = []
for i in range(2):
    for j in range(2):
        if i == 1 and j == 0:
            continue
        axs.append( fig.add_axes( [i*dx, j*dy, dx, dy] ) )

axs.append( fig.add_axes( [ dx+t_pad/2, t_pad/2+(dy-t_pad)*t_diff, dx-t_pad, (dy-t_pad)*(1-t_diff) ] ) )
axs.append( fig.add_axes( [ dx+t_pad/2, t_pad/2, dx-t_pad, (dy-t_pad)*t_diff ] ) )

for i in range(3):
    ax = axs[i]
    ax.imshow( ds[i], norm=mplc.LogNorm(), cmap=cm.plasma )
    #ax.imshow( ds[i], cmap=cm.plasma )
    ax.grid()
    ax.set_xticks( [] )
    ax.set_yticks( [] )

ax = axs[3]
ax_diff = axs[4]
labels = [ 'no', 'CRE' ]
Tunit = 1e7
Nmin = 1000
for i in range(2):
    ds[3+i][:,0] = ds[3+i][:,0] / 1000
    ds[3+i][:,1] = ds[3+i][:,1] / Tunit
for i in range(2):
    x = ds[3+i][:,0]
    y = ds[3+i][:,1]
    n = ds[3+i][:,2]
    index = n > Nmin
    x = x[index]
    y = y[index]
    ax.plot( x, y, label=labels[i] )
rmin = ds[3][:,0].min()
rmax = ds[3][:,0].max()
ax.legend()
ax.set_xlim( [rmin, rmax] )
ax.set_xscale( 'log' )

index1 = ds[3][:,1] > 0
index2 = ds[4][:,1] > 0
index3 = ds[3][:,2] > Nmin
index4 = ds[4][:,2] > Nmin
index = index1 * index2 * index3 * index4
x = ds[3][:,0][index]
y = ds[3][:,1][index]
y_cre = ds[4][:,1][index]
y_err = ( y_cre - y ) / y
ax_diff.plot( x, y_err )
ax_diff.set_xscale( 'log' )
ax_diff.set_xlim( [rmin, rmax] )


fig.savefig( fn_out )
    
