#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tik
import matplotlib.colors as mplc
from matplotlib import cm
import sys

if ( len( sys.argv ) != 3 ):
    print( "Please give data file and fig file!" )
    exit()


data = np.loadtxt( sys.argv[1] )
data_info = data[ 0, : ]
data = np.flipud(data[ 1:, : ])
z = data_info[0]
XMin = data_info[1] / 1000.0
XMax = data_info[2] / 1000.0
YMin = data_info[3] / 1000.0
YMax = data_info[4] / 1000.0
#print( data_info )
m,n = data.shape
print( "PicSize: (%i,%i)"%( m, n ) )
print( "XMin: %g, XMax: %g, YMin: %g, YMax: %g"%\
        (XMin, XMax, YMin, YMax ) )
#print( data_info )

NX = 5
NY = 5

if ( XMin - int( XMin ) != 0  ):
    xmin = int( XMin ) + 1
else:
    xmin = XMin

xmax = int( XMax )
xl = XMax - XMin
if ( YMin - int( YMin ) != 0 ):
    ymin = int( YMin ) + 1
else:
    ymin = YMin
ymax = int( YMax )
yl = YMax - YMin

print( "xmin: %i, xmax: %i, ymin: %i, ymax: %i"%( xmin, xmax, ymin ,ymax ) )

if ( xmin > xmax ):
    pass
else:
    xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
    xFmtList = xLocList
    xLocList = ( xLocList - XMin ) / xl * n
    print( "xticks: ", xFmtList )
if ( ymin > ymax ):
    pass
else:
    yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
    yFmtList = yLocList
    yLocList = m - ( yLocList - YMin ) / yl * m
    print( "yticks: ", yFmtList )



fig = plt.figure()
ax = fig.add_subplot( 111 )
norm = mplc.LogNorm()
img = ax.imshow( data, norm=norm, cmap=cm.jet )

ax.set_xlabel( "Mpc" )
ax.set_ylabel( "" )
ax.set_title( "cluster(z=%.2f)"%(z) )

xloc = tik.FixedLocator(xLocList)
xfmt = tik.FixedFormatter(xFmtList)
yloc = tik.FixedLocator(yLocList)
yfmt = tik.FixedFormatter(yFmtList)
ax.xaxis.set_major_locator( xloc )
ax.xaxis.set_major_formatter( xfmt )
ax.yaxis.set_major_locator( yloc )
ax.yaxis.set_major_formatter( yfmt )

'''
ax.set_xticks( np.linspace( 0, n, NX ) )
ax.set_xticklabels( DensList )
ax.set_yticks( np.linspace( 0, m, NY ) )
ax.set_yticklabels( TempList )
'''
cbar = fig.colorbar( img )
cbar.set_label( "Temperature" )

plt.savefig( sys.argv[2] )
#plt.show()
