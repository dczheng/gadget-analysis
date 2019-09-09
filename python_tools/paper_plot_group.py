#!/usr/bin/env python3

from my_work_env import *

DataDir     = sys.argv[1]
snap_idx    = int( sys.argv[2] )
Projs       = sys.argv[3:]
#tProjs       = [ 'z', 'x', 'x', 'x', 'x', 'x' ]
NGroup  = len( Projs )
ds_name = [  "Density", "MagneticField", "Mach", "Cre_e",\
             "Radio"
            ]
fig_name = [  r"$\rm {\rho}/{\bar{\rho}}$", \
              r"$\rm B \,[\mu G]$",\
              r"$\rm Mach$",  \
              r"$\rm \epsilon / \epsilon_{\rm bar}$", \
              r"$\rm I_{1.4G}\, [mJy\,sr^{-1}]$"\
              ]
norms   = [  mplc.LogNorm(),\
             mplc.LogNorm(),\
             mplc.LogNorm(),\
             mplc.LogNorm(),\
             mplc.LogNorm()
             ]
cmaps   = [ \
        #cm.jet,\
        #cm.viridis,\
        #cm.ocean,\
        #cm.magma,\
        #cm.summer,\
        #cm.winter,\
        cm.hot,\
        #plt.get_cmap( 'ds9b' ),\
        #cm.inferno,\
        #cm.gnuplot,\
        cm.nipy_spectral,\
        #cm.gist_ncar,\
        plt.get_cmap( 'ds9a' ),\
        cm.plasma,\
        cm.magma,\
        cm.gist_heat,\
        ]
m = len(ds_name)
n = NGroup

ds = []
for i in range(m):
    t = []
    for j in range(n):
        fn = "%s/%s/%s_%03i_%04i_%s.dat"%( \
              DataDir, \
              ds_name[i], ds_name[i], \
              snap_idx, j, Projs[j] )
        print( "load `%s` ..."%fn )
        t.append( np.loadtxt(fn)[1:,:] )
    ds.append(t)

for i in range(m):
    if "Cre_e" in ds_name[i]:
        for j in range(n):
            d = ds[i][j]
            #d[ d<d.max()*1e-5 ] = 0
    if "Radio" in ds_name[i]:
        for j in range(n):
            ds[i][j] = ds[i][j] / mycc.mJy
    if "Mach" in ds_name[i]:
        for j in range(n):
            ds[i][j][0,0] = 1

for i in range(m):
    vmin = 1e100
    vmax = -vmin

    for j in range(n):
        v = ds[i][j].max()
        if v > vmax:
            vmax = v
        v = ds[i][j][ ds[i][j]>0 ].min()
        if v < vmin:
            vmin = v
    for j in range(n):
        ds[i][j][0,0] = vmin
        ds[i][j][0,1] = vmax
    print( "%-15s min:%g, max:%g"%( '['+ds_name[i]+']', vmin, vmax ) )

for i in range(m):
    if not(norms[i]):
        for j in range(n):
            d = ds[i][j]
            d[ d==0 ] = np.nan


fs = 10 
t_cbar = 0.4 
fig = plt.figure( figsize=(fs*(n+t_cbar), fs*m) )
dx = 1 / (n+t_cbar)
dy = 1 / m
axs = []
for i in range(m):
    aa = [ fig.add_axes( [ j*dx, i*dy, dx, dy ] ) for j in range(n) ]
    axs.append( aa )
for i in range(m*n):
    a = axs[i//n][i%n]
    a.set_xticks( [] )
    a.set_yticks( [] )

cbar_axs = [ fig.add_axes( [n*dx, i*dy, t_cbar*dx/4, dy] ) for i in range(m) ]

for i in range(m):
    cmap = cmaps[i]
    norm = norms[i]
    cax = cbar_axs[m-1-i]
    for j in range(n):
        fn = "%s/%s/%s_%03i_%04i_%s.dat"%( \
              DataDir, \
              ds_name[i], ds_name[i], \
              snap_idx, j, Projs[j] )
        print( "plot `%s` ..."%fn )
        ax = axs[m-1-i][j]
        img = ax.imshow( ds[i][j], norm=norm, cmap=cmap )
        if j == 0:
            cbar = plt.colorbar( img, cax=cax )
            #cbar.set_ticks( [] )
            set_tick_params( cax, 60 )
            cax.set_ylabel( fig_name[i], fontsize=60 )
            #if "Mach" in ds_name[i]:
            #    mmin = int( np.log10( ds[i][j][ds[i][j]>0].min() ) )
            #    mmax = int( np.log10( ds[i][j].max() ) )
            #    tik = [ 1, 2.5, 5, 10, 100 ]
            #    cbar.set_ticks([])
        ax.invert_yaxis()

mm, nn = ds[0][0].shape
#font = FontProperties()
#font.set_size( 'xx-large' )
#font.set_weight( 'medium' )
for i in range(n):
    ax = axs[m-1][i]
    ax.text( 0.05*nn, 0.05*mm, r"$G_{%i}$"%(i+1), fontsize=90 )

#fig.savefig( str(snap_idx) + '_' + ''.join(Projs) + '.png' )
fig.savefig( str(snap_idx) + '_' + ''.join(Projs) + '.pdf' )

