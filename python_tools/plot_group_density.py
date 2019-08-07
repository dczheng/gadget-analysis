#!/usr/bin/env python3

from my_work_env import *

fn_pre = sys.argv[1]
N = int(sys.argv[2])
fn_suf = sys.argv[3]
dir_out = sys.argv[4]

M = 5
dx = 1/M
dy = 1/M


index = 0
ss = [ 'tot', 'dm', 'gas',  '',  '', 'star', 'bh' ]

while index < N:
    fig = plt.figure( figsize=(10,10) )
    axs = [ fig.add_axes( [ii//M*dx, ii%M*dy, dx, dy] ) for ii in range( M*M ) ]
    for a in axs:
        a.set_xticks( [] )
        a.set_yticks( [] )
    for i in range( M*M ):
        ii = index + i
        if ii > N-1:
            continue
        fn = fn_pre + '%04i'%ii + fn_suf
        print( 'plot %s ...'%(fn) )
        d = np.loadtxt( fn )
        h = d[0,:]
        d = d[1:,:]
        m, n = d.shape

        axs[i].imshow( d, norm=mplc.LogNorm() )

        ms = h[15:15+6]
        r200 = h[15+6]
        ms = np.array([ ms.sum() ] + list(ms))

        dd = 1
        for j in range(1):
            if not( ms[j] > 0 ):
                continue
            t = r'$M_{%s}: \, %s M_{\odot}$'%( ss[j], real2tex( ms[j], n=2 ) )
            axs[i].text( 0.1*n, 0.1*dd*m, t )
            dd += 1
        axs[i].text( 0.1*n, 0.1*dd*m, 'r200: %.1f Mpc'%( r200/1000 ) )


    fn_out = dir_out + '/gs_%03i.png'%index 
    fig.savefig( fn_out )
    plt.close()

    index += M*M

