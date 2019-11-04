#!/usr/bin/env python3
 
from my_work_env import *


g_dir = sys.argv[1]
si    = int(sys.argv[2])
N = int( sys.argv[3] )

pp = [ 'x', 'y', 'z' ]
cols = 5 
rows =  N // cols
if N % cols:
    rows += 1
fs = 2
dx = 1.0/cols
dy = 1.0/rows

for p in range(3):
    fig = plt.figure( figsize=(fs/dx, fs/dy) )
    axs = [ fig.add_axes( [j*dx, (rows-1-i)*dy, dx, dy] )\
        for i in range(rows)\
            for j in range(cols) ]
    for a in axs:
        a.set_xticks( [] )
        a.set_yticks( [] )

    for gi in range(N):

        ax = axs[gi]
        fn = "%s/Density/Density_%03i_%04i_%s.dat"%( g_dir, si, gi, pp[p] )
        print( fn )
        d = np.loadtxt( fn )[ 1:, : ]
        #print( d[d>0].min(), d.max() )
        #ax.imshow( d )
        #ax.set_title( '%i'%gi )
        m, n = d.shape
        ax.text( 0.1*n, 0.3*m, '%i'%gi, fontsize=30 )
        ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.jet )
        ax.grid()
    fig.savefig( "g_%s.png"%(pp[p]) )
    plt.close()


