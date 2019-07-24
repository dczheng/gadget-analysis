#!/usr/bin/env python3

from my_work_env import *
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tik


#matplotlib.style.use( 'default' )
data_name = [ 'Phase.dat', \
                'cre_Phase.dat' \
                ]
titlelist =[ "S1" , \
"S2" \
]
dN = len( data_name )

ds = []
for dn in data_name:
    print( 'load %s ...'%(data_dir+dn) )
    ds.append( np.loadtxt( data_dir + dn ) )

data_info = ds[0][ 0, : ]

for i in range( dN ):
    ds[i] = ds[i][1:,:]

for i in range( 1, dN ):
    index = np.where( ds[0] < 0.1 )
    t1 = ds[i].copy()
    t2 = ds[0].copy()
    t1[index] = 0
    t2[index] = 0
    ds[i] = (t1 - t2) / t2


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

x3 = (3 - DensMin) / ( DensMax - DensMin ) * n
y5 = (5 - TempMin) / ( TempMax - TempMin ) * m
y7 = (7 - TempMin) / ( TempMax - TempMin ) * m

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
        #dd[ dd>0.15 ] = np.nan
        #dd[ dd<-0.15 ] = np.nan

    img = axs[i].imshow( dd, norm=my_norm, cmap=my_cmap )
    cbar = plt.colorbar( img, ax=ax, pad=0, shrink=0.7 )

    t = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels( t, fontsize=10 )
    cbar.ax.tick_params( direction='in', width=0.3, length=1.5, labelsize=10 )
    #exit()

    ax.invert_yaxis()
    ax.grid()

    ax.xaxis.set_major_locator( xloc )
    ax.xaxis.set_major_formatter( xfmt )
    #ax.xaxis.set_tick_params(labelsize=5)

    ax.yaxis.set_major_locator( yloc )
    ax.yaxis.set_major_formatter( yfmt )
    ax.yaxis.set_tick_params(labelsize=5)

    ax.plot( [x3, x3], [0, y5], 'w--' )
    ax.plot( [0, n], [y5, y5], 'w--' )
    ax.plot( [0, n], [y7, y7], 'w--' )
    fc = 'black'
    fs = 10
    font = FontProperties()
    #font.set_family('monospace')
    font.set_size( 'large' )
    font.set_weight('medium')

    fx = [ x3*0.4, (n-x3)*0.1+x3, n*0.1, n*0.85 ]
    fy = [ y5*0.1, y5*0.3, (y7-y5)*0.75+y5, (m-y7)*0.3+y7 ]
    ft = [ 'Diffuse', 'Condensed', 'Warm-hot', 'Hot' ]

    for kk in range(len(fx)):
        ax.text( fx[kk], fy[kk], ft[kk], \
            fontproperties=font )

    if i == 0:
        #cbar.set_label( r'$\frac{dF(<T, <\frac{\rho}{\rho_{bar}})}{dlog_{10}(T)dlog_{10}({\rho}/{\rho_{Bar}})}$',
        #        fontsize=6)
        ax.set_title( r'$f(T, {\rho}/{\bar{\rho}})$',)
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
    ax.set_xlabel( r"${\rho}/{\bar{\rho}}$" )
    ax.tick_params( axis='both', direction='in', pad=5, labelsize=10 )


    #axs[i].tick_params[ 'fontsize' ] = 5


plt.tight_layout()
fig.savefig( figs_dir + 'GasPhase.pdf' )
