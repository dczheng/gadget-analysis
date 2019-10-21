#!/usr/bin/env python3

from my_work_env import *
from matplotlib.patches import Circle
from scipy.signal import convolve2d

DataDir     = sys.argv[1]
snap_idx    = int( sys.argv[2] )
g_idx       = int( sys.argv[3] )
ds_name = [  "Density", "MagneticField", "Mach", "Cre_e",\
             "Radio1",
             "RadioIndex",
             "Radio",
            ]
fig_name = [  r"$ {\rho}/{\bar{\rho}}$", \
              r"$ B \,[\rm \mu G]$",\
              r"$\rm Mach$",  \
              r"$ \epsilon / \epsilon_{\rm bar}$", \
              r"$ I_{\rm 1.4G}\, [\rm Jy\,arcmin^{-2}]$",\
              r'$ \alpha_{\rm rad} \; [\rm 350MHz-1.4GHz]$',
              r"$ I_{\rm 350M}\, [\rm Jy\,arcmin^{-2}]$",\
              ]
norms   = [  mplc.LogNorm,\
             mplc.LogNorm,\
             None,\
             mplc.LogNorm,\
             mplc.LogNorm,
             None,
             mplc.LogNorm,
             ]
cmaps   = [ \
        cm.hot,\
        #cm.viridis,\
        #cm.gist_ncar,\
        #cm.gnuplot2,\
        #cm.spectral,\
        #plt.get_cmap( 'ds9b' ),\
        cm.magma,\
        plt.get_cmap( 'ds9a' ),\
        #cm.plasma,\
        cm.spectral,\
        cm.magma,\
        #cm.magma,\
        #cm.jet,\
        #cm.ocean,\
        #cm.magma,\
        #cm.summer,\
        #cm.winter,\
        #cm.inferno,\
        #cm.gnuplot,\
        #cm.nipy_spectral,\
        cm.spectral,\
        cm.magma,\
        cm.gist_heat,\
        ]

ds = []
nds = len(ds_name)
for i in range(nds):
    t = []
    for c in ['x', 'y', 'z']:
        fn = "%s/%s/%s_%03i_%04i_%s.dat"%( \
              DataDir, \
              ds_name[i], ds_name[i], \
              snap_idx, g_idx, c )
        print( 'load `%s`...'%fn )
        t.append( np.loadtxt(fn) )
    ds.append(t)

fs = 5 
t_cbar = 0.2 
n = 3
m = nds-1
fig = plt.figure( figsize=(fs*m, fs*(n+t_cbar)) )
dx = 1 / m
dy = 1 / (n+t_cbar) 
axs = []
for j in range(m):
    aa = [ fig.add_axes( [ j*dx, i*dy+dy*t_cbar, dx, dy ] ) for i in range(n) ]
    axs.append( aa )
r = 2/3.0 
cbar_axs = [ fig.add_axes( [j*dx, t_cbar*dy*r, dx, t_cbar*dy*(1-r)] )\
                    for j in range(m) ]
for aa in axs:
    for a in aa:
        a.set_xticks( [] )
        a.set_yticks( [] )

L = []
for i in range(n):
    L.append(np.abs(ds[0][i][0,1]) * 2 / 1000)

print( np.abs(L) )
for j in range(m+1):
    for i in range(n):
        ds[j][i] = ds[j][i][1:,:]

r_cre_e = 1e-4
r_rad = 1e-12

idx_rad = []
for j in range(m+1):
    if "Radio" in ds_name[j] and "RadioIndex" not in ds_name[j]:
        idx_rad.append(j)
print( idx_rad )

xxx = [ 'Y-X', 'Z-X', 'X-Y' ]
for j in range(m):
    print( 'plot %s'%ds_name[j] )
    vmin = np.min( [ dd[dd>0].min() for dd in ds[j] ] )
    vmax = np.max( [ dd.max() for dd in ds[j]] )
    print( "vmin: %g, vmax: %g"%(vmin, vmax) )
    cmap = cmaps[j]
    norm = norms[j]
    cax = cbar_axs[j]

    if "Cre_e"  in ds_name[j]:
        for i in range(n):
            d = ds[j][i]
            d[d<d.max()*r_cre_e] = 0
            ds[j][i] = d
        #vmin = vmax * 1e-4

    if "Radio" in ds_name[j] and "RadioIndex" not in ds_name[j]:
        #print( ds_name[j] )
        for i in range(n):
            d = ds[j][i]
            d[d<d.max()*r_rad] = 0
            d = d / mycc.Jy
            ds[j][i] = d
        #vmin = vmax * 1e-10

    if "RadioIndex" in ds_name[j]:
        for i in range(n):
            d = ds[idx_rad[0]][i]
            idx = d<d.max()*r_rad
            for kk in idx_rad[1:]:
                d = ds[kk][i]
                idx1 = d<d.max()*r_rad
                idx = idx + idx1

            d = ds[j][i]
            d[idx] = 0
            print( len(d[d>3]), len(d[d>0]), d[d>3] )
            d[d>3] = 0

            ds[j][i] = d



    vmin = np.min( [ dd[dd>0].min() for dd in ds[j] ] )
    vmax = np.max( [ dd.max() for dd in ds[j]] )
    print( "[new] vmin: %g, vmax: %g"%(vmin, vmax) )

    for i in range(n):
        d = ds[j][i]
        ax = axs[j][i]
        if not norm:
            d[ d==0 ] = np.nan

        if norm:
            img = ax.imshow( d, norm=norm(vmin=vmin, vmax=vmax), cmap=cmap)
        else:
            img = ax.imshow( d, vmin=vmin, vmax=vmax, cmap=cmap)

        if i == 0:
            cbar = plt.colorbar( img, cax=cax, orientation='horizontal' )
            #cbar.set_ticks( [] )
            set_tick_params( cax, 35 )
            #cax.set_xlabel( fig_name[j], fontsize=30 )

            mm, nn = d.shape
            rx = 0.4
            ry = 0.05
            fs = 40

            if "MagneticField" in ds_name[j]:
                make_log_ticks( vmin, vmax, 2, a = 2, axis=cax.xaxis )

            if "Density" in ds_name[j]:
                make_log_ticks( vmin, vmax, 2, a = 1, axis=cax.xaxis )

            if "Cre_e" in ds_name[j]:
                make_log_ticks( vmin, vmax, 2, a = 2, axis=cax.xaxis )

            if "Radio" in ds_name[j] and "RadioIndex" not in ds_name[j]:
                make_log_ticks( vmin, vmax, 2, a = 3, axis=cax.xaxis )
                rx = 0.1

            if "RadioIndex" in ds_name[j]:
                rx = 0.1
                fs = 30

            ax.text(rx*nn, ry*mm, fig_name[j], fontsize=fs )
        ax.invert_yaxis()
        if j == 0:
            mm, nn = d.shape
            fs = 30
            ax.text(0.01*nn, 0.6*mm, r'$%.0f \, h^{-1}\, {\rm Mpc} $'%L[i],\
                        rotation=90, fontsize=fs )

            ax.text(0.1*nn, 0.8*mm, '%s'%xxx[i], fontsize=fs )

fig.savefig( 'g_%03i_%04i.pdf'%(snap_idx, g_idx) )
