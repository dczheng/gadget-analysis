#!/usr/bin/env python3

from my_work_env import *

DataDir     = sys.argv[1]
snap_idx    = int( sys.argv[2] )
Projs       = sys.argv[3:]
#tProjs       = [ 'z', 'x', 'x', 'x', 'x', 'x' ]
NGroup  = len( Projs )
ds_name = [  "Density", "MagneticField", "Mach", "Cre_e"]
fig_name = [  r"$\rm Density$", r"$\rm MagneticField$", r"$Mach$", r"$\epsilon_{\rm CRE}$"]
norms   = [  mplc.LogNorm(), mplc.LogNorm(), mplc.LogNorm(), mplc.LogNorm() ]
cmaps   = [ \
        #cm.jet,\
        plt.get_cmap( 'ds9b' ),\
        plt.get_cmap( 'ds9a' ),\
        plt.get_cmap( 'ds9a' ),\
        cm.jet,\
        cm.jet,\
        cm.gist_heat,\
        cm.nipy_spectral,\
        cm.jet,\
        cm.jet \
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


fs = 5
t_cbar = 0.2 
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

cbar_axs = [ fig.add_axes( [n*dx, i*dy, t_cbar*dx/3, dy] ) for i in range(m) ]

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
            cbar_axs[m-1-i].set_ylabel( fig_name[i], fontsize=25 )

fig.savefig( str(snap_idx) + '_' + ''.join(Projs) + '.png' )

