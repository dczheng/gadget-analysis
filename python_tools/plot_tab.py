#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mplc
from matplotlib import cm
import matplotlib.ticker as tik

data_list = [ np.loadtxt( 'beta.dat' ),
             np.loadtxt( 'tab_beta.dat' ),
             np.loadtxt( 'E.dat' ),
             np.loadtxt( 'tab_E.dat' ),
             np.loadtxt( 'T.dat' ),
             np.loadtxt( 'tab_T.dat' ) ]
data_name = [ 'beta', 'tab_beta', 'E', 'tab_E', 'T', 'tab_T' ]

m,n = data_list[0].shape
sep_str = '-' * 20

print( sep_str )
print( "Shape:" )
for i in range( len( data_list ) ):
    print( "%10s"%data_name[i], data_list[i].shape )
    if ( (m,n) != data_list[i].shape ):
        print( "ERROR" )
        exit()
print( sep_str )

amin = data_list[0][0,0]
amax = data_list[0][0,1]
qmin = data_list[0][0,2]
qmax = data_list[0][0,3]
qmin = np.log10( qmin )
qmax = np.log10( qmax )
print( "amin: %g, amax: %g, qmin: %g, qmax: %g"%( amin, amax, qmin, qmax ) )

for i in range( len( data_list ) ):
    data_list[i] = data_list[i][1:, :]

m, n = data_list[0].shape

for i in range( len(data_list) ):
    if ( (m,n) != data_list[i].shape ):
        print( "ERROR" )
        exit()

print( sep_str )
for i in range( len(data_list) ):
    print( "%10s, min: %e, max: %e"%( data_name[i], data_list[i].min(), data_list[i].max() ) )



xmin = int( qmin ) + 1
xmax = int( qmax )

xLocList = np.linspace( xmin, xmax, xmax-xmin+1 )
xFmtList = []
for i in xLocList:
    xFmtList.append( r'$10^{%.0f}$'%i )
xLocList = ( xLocList - xmin ) / ( xmax-xmin ) * ( n-1 )

yLocList = np.linspace( amin, amax, 5 )
yFmtList = []
for i in yLocList:
    yFmtList.append( '%.1f'%i )
yLocList = ( yLocList - amin ) / ( amax-amin ) * ( m-1 )

print( 'xLocList:', xLocList )
print( 'xFmtList:', xFmtList )
print( 'yLocList:', yLocList )
print( 'yFmtList:', yFmtList )

xloc = tik.FixedLocator( xLocList )
xfmt = tik.FixedFormatter( xFmtList )
yloc = tik.FixedLocator( yLocList )
yfmt = tik.FixedFormatter( yFmtList )

def my_plot( data, norm, colorbar_label, fn ):
    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    img = ax.imshow( data, norm=norm, cmap=cm.jet )

    cbar = fig.colorbar( img )
    cbar.set_label( colorbar_label )

    ax.xaxis.set_major_locator( xloc )
    ax.xaxis.set_major_formatter( xfmt )
    ax.yaxis.set_major_locator( yloc )
    ax.yaxis.set_major_formatter( yfmt )
    ax.invert_yaxis()
    ax.set_xlabel( 'q' )
    ax.set_ylabel( r'$\alpha$' )
    #ax.axis( 'equal' )
    #fig.tight_layout()

    plt.savefig( fn )
    plt.close()

pic_labels = [ r'$\beta_{\frac{1}{1+q^2}}\left(\frac{\alpha-2}{2}, \frac{3-\alpha}{2}\right)$',
              r'$\beta_{\frac{1}{1+q^2}}\left(\frac{\alpha-2}{2}, \frac{3-\alpha}{2}\right)$',
              'E',
              'tab_E',
              'T',
              'tab_T' ]



print( sep_str )

for i in range( len(data_list) ):
    print( 'plot %s ...'%data_name[i] )
    #data_list[i][ data_list[i] > data_list[i].min() * 1e8 ] = 0
    my_plot( data_list[i],
        mplc.LogNorm(),
        pic_labels[i],
        data_name[i]+'.png' )


