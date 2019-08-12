#!/usr/bin/env python3

from my_work_env import *
import phys_and_const as tc


matplotlib.style.use( 'default' )
plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )
d100 = np.loadtxt( data_dir + 'radio-slice-0.2-100.dat' )[1:, :]
d500 = np.loadtxt( data_dir + 'radio-slice-0.2-500.dat' )[1:, :]
d1400 = np.loadtxt( data_dir + 'radio-slice-0.2-1400.dat' )[1:, :]
ds = [ d100, d500, d1400 ]
Names = [
        r'$\rm 100 \, MHz\,[mJy sr^{-1}]$',
        r'$\rm 500 \, MHz\,[mJy sr^{-1}]$',
        r'$\rm 1400 \, MHz\,[mJy sr^{-1}]$'
        ]

N = len( ds )
bar_r = 0.05
tit_r = 0.1
figl = 5
fig = plt.figure( figsize=(N*figl/(1-bar_r), figl/(1-tit_r)) )

dx = (1-bar_r) / N
axs = [ fig.add_axes( [i*dx, tit_r, dx, 1-tit_r] ) for i in range(N) ]

t = 1/3
bar_ax = fig.add_axes( [N*dx, tit_r, bar_r*t, 1-tit_r] )

vmin = ds[0].min()
vmax = ds[0].max()

for d in ds:
    print( d.min(), d.max() )
    if d.min() < vmin:
        vmin = d.min()
    if d.max() > vmax:
        vmax = d.max()
for d in ds:
    d[0,0] = vmin
    d[0,1] = vmax

cmap = cm.jet

axs_all = np.hstack( (axs, bar_ax) )
for ax in axs_all:
    ax.set_xticks( [] )
    ax.set_yticks( [] )

for i in range(N):
    ds[i][ ds[i]<ds[i].max() * 1e-20 ] = 0
    ds[i] = ds[i] / tc.mJy

for i in range(N):
    img = axs[i].imshow( ds[i], norm=mplc.LogNorm(), cmap=cmap )
    axs[i].set_xlabel( Names[i], fontsize=20)
    if i == 0:
        cbar = plt.colorbar( img, cax=bar_ax )

        vmin = int( np.log10(ds[i][ds[i]>0].min()) )
        vmax = int( np.log10(ds[i].max()) )
        t = range( vmin-1, vmax+1, 2 )
        tik = [ 10.0**v for v in t ]
        print( tik )
        #cbar.set_ticks( tik )

        bar_ax.minorticks_off()
        bar_ax.tick_params( axis='x', direction='in', labelsize=15, top=True, bottom=True )

fig.savefig( figs_dir + 'radio-slice.pdf', dpi=300 )
