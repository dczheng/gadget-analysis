#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tik

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

data_dir = './group/'
data_name = [
        'Density',
        'Hge_n',
        'Hge_qmax',
        'Mach',
        'MagneticField',
        'MagneticField'
        ]

titles = [
        'Density',
        'n',
        r'q_{max}',
        'Mach Number',
        'Magnetic Field',
        'Magnetic Field'
        ]

cmaps = [
        cm.jet,
        #cm.gist_ncar,
        cm.jet,
        #cm.jet,
        cm.gist_ncar,
        cm.gist_ncar,
        cm.jet,
        cm.jet,
        ]

contour_flags = [ 1, 1, 1, 1, 1, 1 ]

dN = len( data_name )

figm = 2
fign = 3

fig, axs = plt.subplots( figm,fign )

redshift = 0.09
group_index = 0
projection = 'y'

fn_suffix = '_%.2f_%04i_%s.dat'%( redshift, group_index, projection )

ds = []

for dn in data_name:
    f = data_dir + dn + '/' + dn + fn_suffix
    ds.append( np.loadtxt( f ) )
    print( 'load %s ...'%dn )

gp = ds[0][0,:]

for i in range( dN ):
    dd = ds[i][1:,:]
    ax = axs[i//fign][i%fign]
    my_cmap = cmaps[i]
    title = titles[i]
    cf = contour_flags[i]

    img = ax.imshow( dd, norm=mplc.LogNorm(), cmap=my_cmap )
    cbar = plt.colorbar( img, ax=ax, shrink=0.75, pad=0 )

    if cf:
        plt.contour( dd )

    ax.invert_yaxis()
    ax.set_title( title )

plt.savefig( 'Group.pdf' )

DensMin = data_info[1]
DensMax = data_info[2]
TempMin = data_info[3]
TempMax = data_info[4]
GlobDensMin = data_info[5]
GlobDensMax = data_info[6]
GlobTempMin = data_info[7]
GlobTempMax = data_info[8]
z = data_info[10]

#plt.savefig( 't.png' )

m,n = ds[0].shape
xmin = int( DensMin ) + 1
xmax = int( DensMax )
xl = DensMax - DensMin
ymin = int( TempMin ) + 1
ymax = int( TempMax )
yl = TempMax - TempMin

print( "xmin: %i, xmax: %i, ymin: %i, ymax: %i"%( xmin, xmax, ymin ,ymax ) )

xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
t = [ r"$10^{%.0f}$"%i for i in xLocList ]
xFmtList = t * dN
xLocList = ( xLocList - DensMin ) / xl * (n-1)
xLocList = np.hstack( [xLocList, n+xLocList, 2*n+xLocList] )
print( "xticks: ", xFmtList )

yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
yFmtList = [ r"$10^{%.0f}$"%i for i in yLocList ]
yLocList = ( yLocList - TempMin ) / yl * (m-1)
print( "yticks: ", yFmtList )

xloc = tik.FixedLocator(xLocList)
xfmt = tik.FixedFormatter(xFmtList)
yloc = tik.FixedLocator(yLocList)
yfmt = tik.FixedFormatter(yFmtList)

#fig, axs = plt.subplots( 1, dN, figsize=(3*dN,3) )
gs = gridspec.GridSpec( 1, 2, width_ratios=(1,dN-1) )

axs = [plt.subplot(gs[0, 0]), \
        plt.subplot(gs[0, 1]) ]

ds2 = [ ds[0], np.hstack( ds[1:] ) ]
fig = plt.gcf()
fig.set_size_inches(3*dN, 3)

contour_levels = np.linspace( 3, 6, 4 )
contour_levles = np.flip( 1 / np.power( 10, contour_levels ), 0)
print( "contour levels: ", contour_levles )
con = []

for i in range( 2 ):
    ax = axs[i]
    dd = ds2[i]
    if i == 0:
        my_cmap = cm.jet
        #my_cmap = cm.tab10
        #my_cmap = cm.tab20
        my_norm = mplc.LogNorm()
    else:
        #my_cmap = cm.jet
        #my_cmap = cm.Set1
        #my_cmap = cm.cubehelix
        #my_cmap = cm.ocean
        #my_cmap = cm.magma
        #my_cmap = cm.Paired
        #my_cmap = cm.tab10
        my_cmap = cm.tab20
        #my_cmap = cm.tab20b
        #my_cmap = cm.tab20c
        #my_cmap = cm.Set2
        #my_cmap = cm.Set3
        #my_cmap = cm.Set1
        my_norm = mplc.SymLogNorm( linthresh=1e-3 )

    img = axs[i].imshow( dd, norm=my_norm, cmap=my_cmap )
    cbar = plt.colorbar( img, ax=ax, pad=0, shrink=0.785 )

    t = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels( t, fontsize=7 )
    #exit()

    ax.invert_yaxis()

    ax.xaxis.set_major_locator( xloc )
    ax.xaxis.set_major_formatter( xfmt )
    ax.xaxis.set_tick_params(labelsize=5)

    axs[i].yaxis.set_major_locator( yloc )
    if i == 0:
        ax.yaxis.set_major_formatter( yfmt )
        ax.set_ylabel( r"$T[K]$" )
        ax.yaxis.set_tick_params(labelsize=5)
        con = ax.contour( dd, levels=contour_levles, linestyles='dashdot', linewidths=1, colors=['black'] )
        #print( con )
    else:
        ax.yaxis.set_major_formatter( tik.FixedFormatter([]) )
        #t = t + ' ' * n
        if ( dN > 2 ):
            ax.axvline( x = i*n, color='k', linewidth=1 )

        for l in con.collections:
            ps = l.get_paths()
            for p in ps:
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
                ax.plot( x, y, 'k-.', linewidth=1 )

    ax.set_title( titlelist[i] )
    ax.set_xlabel( r"$\frac{\rho}{\rho_{Bar}}$" )


    #axs[i].tick_params[ 'fontsize' ] = 5


#plt.savefig( 'test.png' )
#plt.subplots_adjust( vspace=0.001 )
plt.savefig( fn_out )
