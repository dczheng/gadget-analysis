#!/usr/bin/env python3

from my_work_env import *
from matplotlib.patches import Circle
from scipy.signal import convolve2d
from scipy import interpolate

DataDir     = sys.argv[1]
snap_idx    = int( sys.argv[2] )
g_idx       = int( sys.argv[3] )
ds_name = [  "Density", "MagneticField", "Mach", "Cre_e",\
            "Cre_qmax",
            "Cre_qmax"
             #"Radio1",
             #"RadioIndex",
             #"Radio",
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
             mplc.LogNorm,\
             #None,
             mplc.LogNorm,\
             mplc.LogNorm,
             None,
             mplc.LogNorm,
             ]
cmaps   = [ \
        #cm.viridis,\
        #cm.brg,\
        #cm.PuRd_r,\
        #cm.gist_earth,\
        #cm.ocean,\
        #cm.winter,\
        cm.jet,\
        cm.jet,\
        #cm.magma,\
        cm.jet,\
        #plt.get_cmap( 'ds9a' ),\
        #cm.jet,\
        #cm.gnuplot,\
        cm.spectral,\
        cm.magma,\
        cm.spectral,\
        cm.gnuplot,\
        #cm.hot,\
        cm.winter,\
        cm.plasma,\
        cm.jet,\
        cm.gist_heat,\
        cm.magma,\
        #cm.ocean,\
        #cm.magma,\
        #cm.summer,\
        #cm.inferno,\
        #cm.nipy_spectral,\
        #cm.viridis,\
        #cm.gist_ncar,\
        #cm.gnuplot2,\
        #cm.spectral,\
        #plt.get_cmap( 'ds9b' ),\
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

r_cre_e = 1e-5
r_rad = 1e-10

idx_rad = []
for j in range(m+1):
    if "Radio" in ds_name[j] and "RadioIndex" not in ds_name[j]:
        idx_rad.append(j)
print( idx_rad )

xxx = [ 'Y-X', 'Z-X', 'X-Y' ]

for j in range(m):
    if "Mach" in ds_name[j]:
        mmax = 5 
        mmin = 1.1
        mach_index = j
        for i in range(n):
            d = ds[j][i]
            mm, nn = d.shape
            print( len(d[d>mmax]), len(d[d>mmax]) / ( mm*nn ), d[d>mmax] )
            #d[d>mmax] = mmax
            d[d<mmin] = 0

for j in range(m):
    print( '-' * 30 )
    print( 'plot %s'%ds_name[j] )
    vmin = np.min( [ dd[dd>0].min() for dd in ds[j] ] )
    vmax = np.max( [ dd.max() for dd in ds[j]] )
    print( "vmin: %g, vmax: %g"%(vmin, vmax) )
    cmap = cmaps[j]
    norm = norms[j]
    cax = cbar_axs[j]

    if "Cre_e"  in ds_name[j]:
        emax = 1e-3
        for i in range(n):
            d = ds[j][i]
            #print( len(d[d>emax]), d[d>emax] )
            #d[d>emax] = emax
            d[d<d.max()*r_cre_e] = 0
            ds[j][i] = d
        #vmin = vmax * 1e-4

    if "MagneticField"  in ds_name[j]:
        r_B = 1e-8
        for i in range(n):
            d = ds[j][i]
            #print( len(d[d>emax]), d[d>emax] )
            #d[d>emax] = emax
            #d[d<d.max()*r_B] = 0
            ds[j][i] = d
        #vmin = vmax * 1e-4

    if "Radio" in ds_name[j] and "RadioIndex" not in ds_name[j]:
        print( ds_name[j] )
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
            amax = 3.1
            print( len(d[d>amax]), len(d[d>0]), d[d>amax] )
            d[d>amax] = 0

            ds[j][i] = d

    #if "Mach" in ds_name[j]:
    #    mmax = 10 
    #    mmin = 1.1
    #    mach_index = j
    #    for i in range(n):
    #        d = ds[j][i]
    #        print( len(d[d>mmax]), d[d>mmax] )
    #        d[d>mmax] = mmax
    #        d[d<mmin] = 0

    #for i in range(n):
    #    ds[j][i][ np.where( ds[mach_index][i] == 0 ) ] =  0

    vmin = np.min( [ dd[dd>0].min() for dd in ds[j] ] )
    vmax = np.max( [ dd.max() for dd in ds[j]] )
    print( "[new] vmin: %g, vmax: %g"%(vmin, vmax) )

    #for dd in ds[j]:
    #    dd[dd==0] = dd[dd>0].min() / 10 


    for i in range(n):
        d = ds[j][i]
        ax = axs[j][i]
        if "Mach" in ds_name[j]:
            mm, nn = ds[j][i].shape
            Y, X = np.meshgrid( range(mm), range(nn) )
            #f = interpolate.interp2d(X, Y, ds[j][i], kind='cubic')
            #f = interpolate.RegularGridInterpolator( range(mm), range(nn), ds[j][i] )
            y = np.linspace( 0, mm, mm*4 ) 
            x = np.linspace( 0, nn, nn*4 ) 
            #d = f(x, y)
        d = ds[j][i]


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
                make_log_ticks( vmin, vmax, 2, a = 2, axis=cax.xaxis )

            if "Cre_e" in ds_name[j]:
                make_log_ticks( vmin, vmax, 2, a = 1, axis=cax.xaxis )

            #if "Mach" in ds_name[j]:
            #    make_ticks( vmin, vmax, nn=1, a=1, r=0.9, axis=cax.xaxis )

            if "Radio" in ds_name[j] and "RadioIndex" not in ds_name[j]:
                make_log_ticks( vmin, vmax, 2, a = 2, axis=cax.xaxis )
                rx = 0.1

            if "RadioIndex" in ds_name[j]:
                rx = 0.1
                fs = 30

            ax.text(rx*nn, ry*mm, fig_name[j], fontsize=fs, color='black' )
        ax.invert_yaxis()
        if j == m-1:
            mm, nn = d.shape
            fs = 40
            ax.text(nn-0.2*nn, 0.6*mm, r'$%.0f \, h^{-1}\, {\rm Mpc} $'%L[i],\
                        rotation=90, fontsize=fs, color='black' )

            ax.text(0.1*nn, 0.8*mm, '%s'%xxx[i], fontsize=fs )

fig.savefig( 'g_%03i_%04i.pdf'%(snap_idx, g_idx) )
