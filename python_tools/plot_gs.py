#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tik


#axs = []


#mask_r = 1e-2
#mask_r = 0
plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

data_dir = "./hge256_data/"
data_name = [ 'hge256_GasState_0.00.dat', \
          #      'hge256_0.001_0.05_GasState_0.00.dat', \
                'hge256_0.01_0.05_GasState_0.00.dat' \
                ]
titlelist =[ "CRE100" , \
#"CRE100\_001            CRE100\_01", \
"CRE100\_01" \
]
dN = len( data_name )

fn_out = 'hge256_out/GasPhase.pdf'

ds = []
for dn in data_name:
    print( 'load %s ...'%(data_dir+dn) )
    ds.append( np.loadtxt( data_dir + dn ) )

data_info = ds[0][ 0, : ]

for i in range( dN ):
    ds[i] = ds[i][1:,:]

for i in range( 1, dN ):
    ds[i] = (ds[i] - ds[0] ) / ds[0]


#ds.append( np.abs(ds[1] - ds[0] ) / ds[0] )
#ds.append( np.abs(ds[2] - ds[0] ) / ds[0] )
#ds.append( (ds[1] - ds[0] ) / ds[0] )
#ds.append( (ds[2] - ds[0] ) / ds[0] )

#print( np.isnan(ds[3]) )


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

#xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
xLocList = np.arange( xmin, xmax+1, 2 )
t = [ r"$10^{%.0f}$"%i for i in xLocList ]
xFmtList = t * dN
xLocList = ( xLocList - DensMin ) / xl * (n-1)
xLocList = np.hstack( [xLocList, n+xLocList, 2*n+xLocList] )
print( "xticks: ", xFmtList )

#yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
yLocList = np.arange( ymin, ymax+1, 2)
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

contour_levels = np.linspace( 1, 4, 4 )
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
        #my_cmap = cm.tab20
        #my_cmap = cm.tab20b
        #my_cmap = cm.tab20c
        #my_cmap = cm.Set2
        #my_cmap = cm.Set3
        #my_cmap = cm.Set1
        #my_cmap = cm.viridis
        my_cmap = cm.PiYG
        my_cmap = cm.BrBG
        my_cmap = cm.RdYlBu
        my_cmap = cm.RdBu
        my_cmap = cm.Spectral

        my_norm = mplc.SymLogNorm( linthresh=1e-3 )
        dd[ dd>0.15 ] = np.nan
        dd[ dd<-0.15 ] = np.nan

    img = axs[i].imshow( dd, norm=my_norm, cmap=my_cmap )
    cbar = plt.colorbar( img, ax=ax, pad=0, shrink=0.75 )

    t = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels( t, fontsize=10 )
    cbar.ax.tick_params( direction='in', width=0.3, length=1.5, labelsize=10 )
    #exit()

    ax.invert_yaxis()

    ax.xaxis.set_major_locator( xloc )
    ax.xaxis.set_major_formatter( xfmt )
    #ax.xaxis.set_tick_params(labelsize=5)

    ax.yaxis.set_major_locator( yloc )
    ax.yaxis.set_major_formatter( yfmt )
    ax.yaxis.set_tick_params(labelsize=5)
    if i == 0:
        #cbar.set_label( r'$\frac{dF(<T, <\frac{\rho}{\rho_{bar}})}{dlog_{10}(T)dlog_{10}({\rho}/{\rho_{Bar}})}$',
        #        fontsize=6)
        ax.set_title( r'$f(T, {\rho}/{\rho_{\rm bar}})$',)
        ax.set_ylabel( r"$T[K]$" )
        con = ax.contour( dd, levels=contour_levles, linestyles='dashdot', linewidths=1, colors=['black'] )
        #print( con )
    else:
        #ax.yaxis.set_major_formatter( tik.FixedFormatter([]) )
        ax.set_title( r'$\rm Relative \, difference$',)
        #t = t + ' ' * n
        if ( dN > 2 ):
            ax.axvline( x = i*n, color='k', linewidth=1 )

        for l in con.collections:
            ps = l.get_paths()
            for p in ps:
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
                ax.plot( x, y, 'b-.', linewidth=1 )
        '''
        ddd = dd.copy()
        ddd[ ddd<0 ] = np.nan
        ddd[ ddd>0.1 ] = np.nan
        ddd[ np.isnan( ddd ) ] = 1
        #print( ddd )

        fig2 = plt.figure()
        ax2 = fig2.add_subplot( 111 )

        img2 = ax2.imshow( ddd )
        plt.colorbar( img2, ax=ax2 )

        con2 = ax2.contour( ddd, linewidths=1, colors=['red'] )

        for l in con2.collections:
            ps = l.get_paths()
            for p in ps:
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
                if ( len(x)<50 ):
                    continue
                ax.plot( x, y, 'k-', linewidth=1 )
            break
        #ax2.invert_yaxis()
        #fig2.savefig( 'x.png' )

        ddd = dd.copy()
        ddd[ ddd>0 ] = np.nan
        ddd[ ddd<-0.1 ] = np.nan
        ddd[ np.isnan( ddd ) ] = 1
        #print( ddd )

        fig2 = plt.figure()
        ax2 = fig2.add_subplot( 111 )

        img2 = ax2.imshow( ddd )
        plt.colorbar( img2, ax=ax2 )

        con2 = ax2.contour( ddd, linewidths=1, colors=['red'] )

        for l in con2.collections:
            ps = l.get_paths()
            for p in ps:
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
                if ( len(x)<50 ):
                    continue
                ax.plot( x, y, 'w-', linewidth=1 )
            break
        #ax2.invert_yaxis()
        #fig2.savefig( 'x.png' )
        '''

    #ax.set_title( titlelist[i] )
    #ax.set_xlabel( r"$\frac{\rho}{\rho_{bar}}$" )
    ax.set_xlabel( r"${\rho}/{\rho_{\rm bar}}$" )
    ax.tick_params( axis='both', direction='in', pad=5, labelsize=10 )


    #axs[i].tick_params[ 'fontsize' ] = 5


#plt.savefig( 'test.png' )
#plt.subplots_adjust( vspace=0.001 )
plt.tight_layout()
fig.savefig( fn_out )


'''
        ds2[1] = np.hstack( [ds[0], ds[1], ds[2]] )

for i in range(2):

    if i == 0:
        dd = np.hstack( [ds[0], ds[1], ds[2]] )

        for j in range(3):

    else:
        dd = np.hstack( [ds[3], ds[4]] )
        for j in range(2):
            axs[i].text( n/7 + j*n, n+2, titlelist[j+3],fontsize=7 )

    axs[i].axvline( x = n, color='k', linewidth=1 )

    img = axs[i].imshow( dd, norm=my_norm, cmap=my_cmap )
    #cbar = plt.colorbar( img, ax=axs[i], pad=0, extend='neither' )
    #cbar = plt.colorbar( img, ax=axs[i], pad=0, shrink=0.788 )
    cbar = plt.colorbar( img, ax=axs[i], pad=0 )


    if i == 1:

    #cbar.ax.tick_params( labelsize=8 )

#plt.tight_layout()
#plt.show()
plt.savefig( 'hge256_out/GasPhase.pdf' )
'''

'''
fn_png = sys.argv[1][:-4] + '.png'
print( 'save image to ' + fn_png )
plt.savefig( fn_png )
#plt.show()
'''
