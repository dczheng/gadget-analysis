#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as tik
import matplotlib.colors as mplc
from matplotlib import cm

parser = argparse.ArgumentParser()

parser.add_argument( '-f',   type=str,    help='file name', required=True )
parser.add_argument( '-log', help='log plot', action="store_true" )
parser.add_argument( '-m',   help='mass', action="store_true" )
parser.add_argument( '-eps',  help='save eps picture', action="store_true" )
parser.add_argument( '-cu',  type=float,  help='upper cutoff' )
parser.add_argument( '-cd',  type=float,  help='lower cutoff' )
parser.add_argument( '-cl',  type=str,    help='colorbar label' )

args = parser.parse_args()

fn = args.f

data = np.loadtxt( fn )

data_info = data[ 0, : ]
data = data[ 1:, : ]
XMin =     data_info[1] / 1000.0
XMax =     data_info[2] / 1000.0
YMin =     data_info[3] / 1000.0
YMax =     data_info[4] / 1000.0
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

if ( proj == 0 ):
    ProDire = 'x'
if ( proj == 1 ):
    ProDire = 'y'
if ( proj == 2 ):
    ProDire = 'z'

m,n = data.shape

print( "PicSize: (%i,%i)"%( m, n ) )
print( "XMin: %g, XMax: %g, YMin: %g, YMax: %g"%\
        (XMin, XMax, YMin, YMax ) )

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
    xLocList = ( xLocList - XMin ) / xl * (n-1)
    print( "xticks: ", xFmtList )
if ( ymin > ymax ):
    yLocList = []
    yFmtList = []
else:
    yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
    yFmtList = yLocList
    yLocList = ( yLocList - YMin ) / yl * (m-1)
    print( "yticks: ", yFmtList )

xloc = tik.FixedLocator(xLocList)
xfmt = tik.FixedFormatter(xFmtList)
yloc = tik.FixedLocator(yLocList)
yfmt = tik.FixedFormatter(yFmtList)


fig = plt.figure()
ax = fig.add_subplot( 111 )

if args.log:
    if args.cu != None:
        cu = np.power( 10, args.cu )
        data[data>=cu] = 0

    if args.cd != None:
        cd = np.power( 10, args.cd )
        data[data<=cd] = 0

    norm = mplc.LogNorm()
    img = ax.imshow( data, norm=norm, cmap=cm.jet )
else:
    if args.cu != None:
        data[data>=args.cu] = np.nan

    if args.cd != None:
        data[data<=args.cd] = np.nan

    img = ax.imshow( data, cmap=cm.jet )


ax.set_xlabel( "Mpc" )
ax.set_ylabel( "" )
ax.set_title( "cluster (z=%.2f)"%(z) )
ax.xaxis.set_major_locator( xloc )
ax.xaxis.set_major_formatter( xfmt )
ax.yaxis.set_major_locator( yloc )
ax.yaxis.set_major_formatter( yfmt )
ax.invert_yaxis()
ax.grid()

if args.cl:
    clabel = '$%s\ (%s)$'%( args.cl, ProDire )
    #print( clabel )
    #exit()
else:
    clabel = r'$%s$'%(ProDire)

cbar = fig.colorbar( img )
cbar.set_label( clabel )

if args.m:

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

fn_out = fn[:-4] + '.png'

if args.eps:
    fn_out = fn[:-4] + '.eps'

print( "save image to " +  fn_out )
plt.savefig( fn_out )
