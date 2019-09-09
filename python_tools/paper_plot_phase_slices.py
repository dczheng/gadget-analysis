#!/usr/bin/env python3

from my_work_env import *

#matplotlib.style.use( 'default' )

data_dir = sys.argv[1]
fn_out = sys.argv[2]


#d = np.loadtxt( "./Density_diffuse.dat"  )[1:,:]
#plt.imshow( d, norm=mplc.LogNorm(), cmap=cm.hot  )
#ax = plt.gca()
#ax.set_xticks( []  )
#ax.set_yticks( []  )
#ax.invert_yaxis()
#plt.savefig( fn_out, figsize=(5,5)  )
#exit()

N = 4
ss =  [\
'Diffuse',\
'Warm-hot',\
'Hot',\
'Condensed',\
]
ds = []
for s in ss:
    fn = "%s/Density_%s.dat"%( data_dir, s )
    print( 'load %s ...'%fn )
    ds.append( np.loadtxt(fn)[1:,:] )


NNN = 2
vmin = 1e10
vmax = -vmin 
for d in ds:
    if not( d.max() > 0 ):
        continue
    vmin = np.min( [d[d>0].min(), vmin] )
    vmax = np.max( [d.max(), vmax] )
print( vmin, vmax )

for i in range(N):
    if i < NNN:
        ds[i][ds[i]<vmin] = vmin
    else:
        ds[i][0,1] = vmin
    ds[i][0,0] = vmax

for d in ds:
    print( d.min(), d.max(), d.shape )

fs = 5 
t_cbar = 0.4
fig = plt.figure( figsize=( fs*(N+t_cbar), fs) )
dx = 1 / (N+t_cbar) 
dy = 1
axs = [ fig.add_axes([i*dx, 0, dx, dy]) for i in range(N) ]
ax_cbar = fig.add_axes([N*dx+dx*t_cbar*0.05, 0, dx*t_cbar/5, dy])
for a in axs:
    a.set_xticks( [] )
    a.set_yticks( [] )

#ds[0] = ds[0] + ds[1] + ds[2] + ds[3]
m, n = ds[0].shape
cmap = cm.hot
#cmap = cm.jet
tx = 0.1 * n
ty = 0.8 * m
for i in range(N):
    d = ds[i]
    ax = axs[i]
    img = ax.imshow( d, norm=mplc.LogNorm(), cmap=cmap )
    if i == 1:
        cbar = plt.colorbar( img, cax=ax_cbar )
        ax_cbar.minorticks_off()
        set_tick_params( ax_cbar, 30 )
        vmin = int( np.log10(vmin) )
        vmax = int( np.log10(vmax) )
        tik = [ 10**v for v in range(vmin, vmax+1, 2) ]
        cbar.set_ticks( tik )
        cbar.ax.set_ylabel( r'$\rho/\bar{\rho}$', fontsize=30 )
    ax.invert_yaxis()
    ax.spines['top'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    if i < NNN:
        ax.text( tx, ty, ss[i], color='white', fontsize=40 )
    else:
        ax.text( tx, ty, ss[i], color='black', fontsize=40 )
        ax.set_facecolor( 'white' )

fig.savefig( fn_out, dpi=300 )
