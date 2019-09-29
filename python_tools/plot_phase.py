#!/usr/bin/env python3

from my_work_env import *
import matplotlib.ticker as tik

#matplotlib.style.use( 'default' )

ds = [ np.loadtxt( sys.argv[i+1] ) for i in range(4) ]
# 0: d1, 1:d1_cre, 2:d2, 3:d2_cre
hs = []
f_out = sys.argv[5]

for i in range(4):
    hs.append( ds[i][0,:] )
    ds[i] = ds[i][1:,:]

zs = [hs[0][0], hs[2][0]]
if zs[0] != hs[1][0] or zs[1] != hs[3][0]:
    print( "error1" )
    exit()

print( zs )

xmin, xmax, xlog, ymin, ymax, ylog = hs[0][1:7] 
m, n = ds[i].shape
for h in hs:
    mm, nn = ds[i].shape
    if xmin != hs[i][1] or \
       xmax != hs[i][2] or \
       ymin != hs[i][4] or \
       ymax != hs[i][5] or\
       m != mm or\
       n != nn:
       print( "error2" )
       exit()

Nmin = 1e5 
contour_levles = [ 1e4, 1e5, 2e5, 1e6 ]
print( "contour levels: ", contour_levles )
Errmin = 1e-2
err_lognorm = 0
zfmt = r'$z:%.0f$'

#ds = [ ds[1], ds[0], ds[3], ds[2] ]
for i in range(2):
    index = np.where( ds[i*2] < Nmin )
    t1 = ds[i*2].copy()
    t2 = ds[i*2+1].copy()
    t1[index] = 1
    t2[index] = 1
    t = (t2 - t1) / t1
    ds[i*2+1] = t

vmin = [0] * 4
vmax = [0] * 4
for i in range(2):
    if i == 0:
        vmin[i] = np.min( [ ds[i][ds[i]>0].min(), ds[i+2][ds[i+2]>0].min() ] )
    else:
        vmin[i] = np.min( [ ds[i].min(), ds[i+2].min() ] )
    vmax[i] = np.max( [ ds[i].max(), ds[i+2].max() ] )

    if vmin[i] * vmax[i] < 0:
        vv = np.max( [ np.abs(vmin[i]), np.abs(vmax[i]) ] )
        vmin[i] = -vv
        vmax[i] = vv

    vmin[i+2] = vmin[i]
    vmax[i+2] = vmax[i]

#print( vmin, vmax )

if not(err_lognorm):
    for i in range(2):
        ds[i*2+1][ds[i*2+1] == 0] = np.nan


d1 = ds[0]
d1_err = ds[1] 

d2 = ds[2]
d2_err = ds[3] 


fs = 8
rc = 0.1

fig = plt.figure( figsize=(fs/(1-3*rc),fs/(1-2*rc)) )
rrc = 0.3 
raw = ( 1-3*rc ) / 2.0
rah = ( 1-2*rc ) / 2.0 
ax1      = fig.add_axes(      [rc,          rah+rc,  raw,     rah] )
ax1_cbar = fig.add_axes(      [rc+raw,      rah+rc,  rc*rrc,  rah] )
ax1_err = fig.add_axes(       [2*rc+raw,    rah+rc,  raw,     rah] )
ax1_err_cbar = fig.add_axes(  [2*rc+2*raw,  rah+rc,  rc*rrc,  rah] )

ax2     = fig.add_axes(       [rc,          rc, raw,    rah] )
ax2_cbar = fig.add_axes(      [rc+raw,      rc, rc*rrc, rah] )
ax2_err = fig.add_axes(       [2*rc+raw,    rc, raw,    rah] )
ax2_err_cbar = fig.add_axes(  [2*rc+2*raw,  rc, rc*rrc, rah] )

axs = [ ax1, ax1_err, ax2, ax2_err ]
axs_cbar = [ ax1_cbar, ax1_err_cbar, ax2_cbar, ax2_err_cbar ]

cmap1 = cm.jet
#cmap2 = cm.jet
#cmap2 = cm.viridis
#cmap2 = cm.tab20
#cmap2 = cm.tab10
cmap2 = cm.PiYG
cmap2 = plt.get_cmap( 'ds9a' )
cmap2 = plt.get_cmap( 'ds9b' )
cmap2 = plt.get_cmap( 'ds9cool' )
cmap2 = plt.get_cmap( 'ds9he' )
cmap2 = cm.PRGn
cmap2 = cm.seismic
cmaps = [ cmap1, cmap2, cmap1, cmap2 ]

norm1 = mplc.LogNorm()
if err_lognorm:
    norm2 = mplc.SymLogNorm( linthresh=Errmin )
else:
    norm2 = None

norms = [ norm1, norm2, norm1, norm2 ] 

print( "xmin: %g, xmax: %g, ymin: %g, ymax: %g"%( xmin, xmax, ymin, ymax ) )
xloc = np.arange( int(xmin-1), int(xmax+2), 2 )
xfmt = [ r"$10^{%.0f}$"%i for i in xloc ]
xloc = ( xloc - xmin ) / (xmax-xmin) * (n-1) 
t1 = []
t2 = []
for i in range( len(xloc) ):
    if xloc[i] > 0 and xloc[i] < n:
        t1.append( xloc[i] )
        t2.append( xfmt[i] )
xloc = t1
xfmt = t2
print( "xticks: ", xfmt )

yloc = np.arange( int(ymin-1), int(ymax+2), 2 )
yfmt = [ r"$10^{%.0f}$"%i for i in yloc ]
yloc = ( yloc - ymin ) / (ymax-ymin) * (n-1) 
t1 = []
t2 = []
for i in range( len(yloc) ):
    if yloc[i] > 0 and yloc[i] < n:
        t1.append( yloc[i] )
        t2.append( yfmt[i] )
yloc = t1
yfmt = t2
print( "yticks: ", yfmt )

xloc = tik.FixedLocator(xloc)
xfmt = tik.FixedFormatter(xfmt)
yloc = tik.FixedLocator(yloc)
yfmt = tik.FixedFormatter(yfmt)

x3 = (3 - xmin) / ( xmax - xmin ) * (n-1)
y5 = (5 - ymin) / ( ymax - ymin ) * (m-1)
y7 = (7 - ymin) / ( ymax - ymin ) * (m-1)

fc = 'black'
font = FontProperties()
#font.set_family('monospace')
font.set_size( 'xx-large' )
font.set_weight('medium')

fx = [ x3*0.4, (n-x3)*0.1+x3, n*0.2, n*0.85 ]
fy = [ y5*0.1, y5*0.3, (y7-y5)*0.6+y5, (m-y7)*0.5+y7 ]
ft = [ 'Diffuse', 'Condensed', 'Warm-hot', 'Hot' ]

for i in range(4):

    if i % 2 == 0:
        con = axs[i].contour( ds[i], levels=contour_levles, linestyles='dashdot', linewidths=1, colors=['black'] )

        for l in con.collections:
            ps = l.get_paths()
            for p in ps:
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
                axs[i+1].plot( x, y, 'b-.', linewidth=1 )

    if i % 2:
        t = ds[i]
        index = ds[i-1] < Nmin 
        t[index] = np.nan

    ds[i][0,0] = vmin[i]
    ds[i][0,1] = vmax[i]

    img = axs[i].imshow( ds[i], norm=norms[i], cmap=cmaps[i] )
    plt.colorbar( img, cax=axs_cbar[i] )
    if i % 2 == 0:
        axs_cbar[i].set_ylabel( r'$\rm {N}/[{dlog(\rho)dlog(T)}]$', fontsize=25 )

    if i > 1:
        axs[i].xaxis.set_major_locator( xloc )
        axs[i].xaxis.set_major_formatter( xfmt )
        axs[i].set_xlabel( r"${\rho}/{\bar{\rho}}$", fontsize=20 )
    else:
        axs[i].set_xticks( [] )

    if i % 2 ==  0:
        axs[i].yaxis.set_major_locator( yloc )
        axs[i].yaxis.set_major_formatter( yfmt )
        axs[i].set_ylabel( r"$T[K]$", fontsize=20 )
    else:
        axs[i].set_yticks( [] )

    set_tick_params( axs_cbar[i], 15 )
    set_tick_params( axs[i], 15 )

    axs[i].invert_yaxis()
    axs[i].grid()

    axs[i].plot( [x3, x3], [0, y5],  'c--' )
    axs[i].plot( [0, n-1], [y5, y5], 'c--' )
    axs[i].plot( [0, n-1], [y7, y7], 'c--' )

    if i == 0:
        axs[i].set_title( "Gas Phase", fontsize=20 )
    if i == 1:
        axs[i].set_title( "Relative Difference", fontsize=20 )


    for kk in range(len(fx)):
        axs[i].text( fx[kk], fy[kk], ft[kk], \
            fontproperties=font )

    axs[i].text( 0, 1, zfmt%zs[i//2], fontproperties=font  )


fig.savefig( f_out )

#my_cmap = cm.Set1
#my_cmap = cm.cubehelix
#my_cmap = cm.ocean
#my_cmap = cm.Paired
#my_cmap = cm.tab10
#my_cmap = cm.tab20b
#my_cmap = cm.tab20c
#my_cmap = cm.Set2
#my_cmap = cm.Set3
#my_cmap = cm.PiYG
#my_cmap = cm.BrBG
#my_cmap = cm.RdYlBu
#my_cmap = cm.RdBu
#my_cmap = cm.Spectral
