#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tik
import matplotlib.colors as mplc
from matplotlib import cm
import sys

if ( len( sys.argv ) != 2 ):
    print( "Please give data file!" )
    exit()


data = np.loadtxt( sys.argv[1] )
data_info = data[ 0, : ]
data = np.flipud(data[ 1:, : ])
XMin = data_info[1] / 1000.0
XMax = data_info[2] / 1000.0
YMin = data_info[3] / 1000.0
YMax = data_info[4] / 1000.0
GlobXMin = data_info[5] / 1000.0
GlobXMax = data_info[6] / 1000.0
GlobYMin = data_info[7] / 1000.0
GlobYMax = data_info[8] / 1000.0
proj = int( data_info[9] )
z = data_info[10]
mgas = data_info[ 15+0 ]
mdm = data_info[ 15+1 ]
ms = data_info[ 15+4 ]
mbh = data_info[ 15+5 ]

#data[ data < data.max() * 1e-5 ] = 0

if ( proj == 0 ):
    ProDire = 'x'
if ( proj == 1 ):
    ProDire = 'y'
if ( proj == 2 ):
    ProDire = 'z'
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
    xLocList = []
    xFmtList = []
else:
    xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
    xFmtList = xLocList
    xLocList = ( xLocList - XMin ) / xl * n
    print( "xticks: ", xFmtList )
if ( ymin > ymax ):
    yLocList = []
    yFmtList = []
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

pt1y = 0.1
pt2y = pt1y + pt1y * 0.5
pt3y = pt1y + pt1y * 1
pt4y = pt1y + pt1y * 1.5

pt1x = pt1y
pt2x = pt1x
pt3x = pt1x
pt4x = pt1x

ax.text( m * pt1x, n * pt1y, r"$ M_{gas}: %.2e$"%( mgas) )
ax.text( m * pt2x, n * pt2y, r"$ M_{dm} : %.2e$"%( mdm ) )
ax.text( m * pt3x, n * pt3y, r"$ M_{s}  : %.2e$"%( ms  ) )
ax.text( m * pt4x, n * pt4y, r"$ M_{bh} : %.2e$"%( mbh ) )

ax.grid()

'''
ax.set_xticks( np.linspace( 0, n, NX ) )
ax.set_xticklabels( DensList )
ax.set_yticks( np.linspace( 0, m, NY ) )
ax.set_yticklabels( TempList )
'''

cbar = fig.colorbar( img )
cbar.set_label( r"$Radio \; (%s) $"%(ProDire) )

fn_png = sys.argv[1][:-4] + '.png'
print( "save image to " +  fn_png )
plt.savefig( fn_png )
#plt.show()
