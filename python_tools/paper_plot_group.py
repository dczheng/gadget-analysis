#!/usr/bin/env python3

from my_work_env import *
from matplotlib.patches import Circle
from scipy.signal import convolve2d

DataDir     = sys.argv[1]
snap_idx    = int( sys.argv[2] )
Projs       = sys.argv[3:]
#tProjs       = [ 'z', 'x', 'x', 'x', 'x', 'x' ]
NGroup  = len( Projs )
ds_name = [  "Density", "MagneticField", "Mach", "Cre_e",\
             "Radio",
             #"Radio1",
             #"Radio2"
            ]
fig_name = [  r"$\rm {\rho}/{\bar{\rho}}$", \
              r"$\rm B \,[\mu G]$",\
              r"$\rm Mach$",  \
              r"$\rm \epsilon / \epsilon_{\rm bar}$", \
              r"$\rm I_{350M}\, [mJy\,arcsec^{-2}]$",\
              r"$\rm I_{1.4G}\, [mJy\,arcsec^{-2}]$",\
              r'$\rm \alpha_{rad} \; [350MHz-1.4GHz]$',
              ]
norms   = [  mplc.LogNorm,\
             mplc.LogNorm,\
             None,\
             mplc.LogNorm,\
             mplc.LogNorm,
             mplc.LogNorm,
             None
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
        cm.gist_heat,\
        ]
m = len(ds_name)
n = NGroup

ds = []
r200 = []
r200_L = []
r200_c = [ 'w', 'w', 'w', 'k', 'k', 'k', 'k' ]
res = []
Ls  = []

for i in range(m):
    t = []
    for j in range(n):
        fn = "%s/%s/%s_%03i_%04i_%s.dat"%( \
              DataDir, \
              ds_name[i], ds_name[i], \
              snap_idx, j, Projs[j] )
        print( "load `%s` ..."%fn )
        d = np.loadtxt(fn)
        t.append( d[1:,:] )
        
        xmin = d[0,1]
        xmax = d[0,2]
        ymin = d[0,4]
        ymax = d[0,5]
        
        N = d[1:, :].shape[0]
        dL = (xmax - xmin) / N 
        Ls.append( (xmin, xmax, ymin, ymax) )
        res.append( dL )
        #print( dL )
        r200_L.append( d[0, 15] )
        r200.append( d[0, 15] / dL )

    ds.append(t)

print( Ls )
print( res )
print( r200 )

print( Ls[:NGroup] )
print( res[:NGroup] )
print( r200[:NGroup] )

name2index = {}
for i in range(m):
    name2index[ ds_name[i] ] = i
print( name2index )

i = name2index[ "Cre_e" ]
for j in range(n):
    d = ds[i][j]
    d[ d<d.max()*1e-4 ] = 0

i = name2index[ "Radio" ]
for j in range(n):
    ds[i][j] = ds[i][j] / mycc.mJy
    idx1 = ds[name2index['Density']][j] > 1e4
    idx2 = ds[name2index['Mach']][j] > 2 
    idx = idx1 * idx2
    #ds[i][j][ idx ] = 0 
    #idx_mag = ds[name2index['MagneticField']][j] > 10 
    #ds[i][j][ idx_mag ] = 0
    ds[i][j][ds[i][j]<ds[i][j].max()*1e-13] = 0

'''
i = name2index[ "Radio1" ]
for j in range(n):
    ds[i][j] = ds[i][j] / mycc.mJy
    idx1 = ds[name2index['Density']][j] > 1e4
    idx2 = ds[name2index['Mach']][j] > 2 
    idx = idx1 * idx2
    #ds[i][j][ idx ] = 0 
    #idx_mag = ds[name2index['MagneticField']][j] > 10 
    #ds[i][j][ idx_mag ] = 0
    ds[i][j][ds[i][j]<ds[i][j].max()*1e-13] = 0

#i = name2index[ "RadioIndex" ]
#for j in range(n):
#    ds[i][j] *= -1
#    #ds[i][j][ds[name2index['Radio']][j]==0] = 0
#    ds[i][j][ds[i][j] > 3] = 0
'''


i = name2index[ "Mach" ]
for j in range(n):
    ds[i][j][0,0] = 1
    #ds[i][j][ds[i][j]>4] = 1
    #ds[i][j][ idx_mag ] = ds[i][j].min()
    #idx1 = ds[name2index['Density']][j] > 1e5
    #idx2 = ds[i][j] > 2 
    #ds[i][j][ idx1*idx2 ] = 1

i = name2index[ "MagneticField" ]
for j in range(n):
    ds[i][j][ds[i][j]<1e-5] = 0
    idx_mag = ds[name2index['MagneticField']][j] > 10 
    print( len(ds[i][j][ idx_mag ]), ds[i][j][idx_mag] )
    #ds[i][j][ idx_mag ] = 10

vmin = []
vmax = []
for i in range(m):
    vmin_local = 1e100
    vmax_local = -vmin_local

    for j in range(n):
        v = ds[i][j].max()
        if v > vmax_local:
            vmax_local = v
        v = ds[i][j][ ds[i][j]>0 ].min()
        if v < vmin_local:
            vmin_local = v
    vmin.append( vmin_local )
    vmax.append( vmax_local )
print( 'vmin:', vmin )
print( 'vmax:', vmax )

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

def f( x, sigma ):
    return 1/(np.sqrt( 2 * np.pi ) * sigma) * np.exp( -x**2/ ( 2*sigma**2 ) )

def gen_kernel( N, kres, sigma ):
    ker = np.zeros( [N,N] )
    xo, yo = N//2, N//2
    for i in range(N):
        for j in range(N):
            r = np.sqrt( (j-xo)**2 + (i-yo)**2 ) * kres
            #print( f(r, sigma) )
            ker[i,j] = f( r, sigma )
            #if r == 0:
            #    print( i, j, r, f(r,sigma) )
    return ker 

#print( gen_kernel( 4, 10, 10 ) )
#print( gen_kernel( 8, 10, 10 ) )
#exit()

conv_flag = 0
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

        if ds_name[i] == 'Radio' and conv_flag:
            if conv_flag:
                sigma = res[j] / ( 2 * np.sqrt( np.log(2) ) )
                kernel = gen_kernel( 4, res[j], sigma )
                print( kernel )
                t = convolve2d( ds[i][j], kernel ) 
                if norm:
                    img = ax.imshow( t, norm=norm(vmin=vmin[i], vmax=vmax[i]),\
                    cmap=cmap)
                else:
                    img = ax.imshow( t, norm=norm, cmap=cmap,\
                    vmin=vmin[i], vmax=vmax[i] )
        else:
            if norm:
                img = ax.imshow( ds[i][j],\
                norm=norm(vmin=vmin[i], vmax=vmax[i]),\
                cmap=cmap)
            else:
                img = ax.imshow( ds[i][j], norm=norm, cmap=cmap,\
                vmin=vmin[i], vmax=vmax[i] )

        mm, nn = ds[i][j].shape

        cir = Circle( xy=(nn/2, mm/2), radius=r200[j], fill=False, color=r200_c[i] );
        ax.add_patch( cir )

        #cir = Circle( xy=(nn/2, mm/2), radius=r200[j]/10, fill=False );
        #ax.add_patch( cir )
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
    ax = axs[0][i]
    ax.text( 0.05*nn, 0.05*mm, r"$R_{200}=%.2f \, h^{-1}{\rm Mpc}$"%( r200_L[i]/1000), fontsize=60 )

#fig.savefig( str(snap_idx) + '_' + ''.join(Projs) + '.png' )
fig.savefig( str(snap_idx) + '_' + ''.join(Projs) + '.pdf' )

