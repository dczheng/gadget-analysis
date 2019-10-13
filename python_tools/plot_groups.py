#!/usr/bin/env python3

from my_work_env import *

data_dir = sys.argv[1]
ds_name  = sys.argv[2]
snapidx  = int( sys.argv[3] )
NG       = int( sys.argv[4] )
Proj     = sys.argv[5]

if len( sys.argv ) > 6:
    cols = sys.argv[6]
else:
    cols = 5

rows = NG // cols
if NG % cols:
    rows += 1

fs = 2
fig = plt.figure( figsize=(fs*cols, fs*rows) ) 
dx = 1 / cols
dy = 1 / rows
axs = [ fig.add_axes([ j*dx, (rows-1-i)*dy, dx, dy ])\
        for i in range(rows) for j in range(cols) ]

for i in range(NG):
    ax = axs[i]
    fn = '%s/%s/%s_%03i_%04i_%s.dat'%(\
            data_dir, ds_name, ds_name, snapidx, i, Proj  )
    print( 'plot %s ...'%fn )
    
    d = np.loadtxt( fn )
    m = d[0, 8:8+6].sum()
    d = d[1:,:]
    nm, nn = d.shape

    mx = real2tex( m, 10, 2 )
    #print( mx )
    ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.hot )
    #ax.text( 0.02*nm, 0.02*nn, r'$%s\, h^{-1}M_\odot$'%mx  )
    ax.text( 0.02*nm, 0.02*nn, r'$%s$'%mx  )


for a in axs:
    a.set_xticks( [] )
    a.set_yticks( [] )
    a.invert_yaxis()
    a.grid()

fn_out = './%s_%03i_%s.png'%( ds_name, snapidx, Proj )
fig.savefig( fn_out )
