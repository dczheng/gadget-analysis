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

parser.add_argument( 'f' )
parser.add_argument( '-log', help='log plot', action="store_true" )
parser.add_argument( '-m',   help='mass', action="store_true" )
parser.add_argument( '-eps',  help='save eps picture', action="store_true" )
parser.add_argument( '-ng',  help='negtive ', action="store_true" )
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
mtot = mgas + mdm + ms + mbh

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

if args.ng:
    data = -data

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
    clabel = '$%s\ [%s]$'%( args.cl, ProDire )
    #print( clabel )
    #exit()
else:
    clabel = r'$%s$'%(ProDire)

cbar = fig.colorbar( img )
cbar.set_label( clabel )

def put_text( ax, x, y, s, v ):
    fmt = r'$%s: %s \times 10^{%s}$'
    t = "%.2e"%v
    t = t.split( 'e' )
    if t[1][0] == '+':
        t[1] = t[1][1:]
    ax.text( x, y, fmt%(s, t[0], t[1]) )

if args.m:

    p0x = 0.01
    p0y = 0.9
    offset = 0.05

    px = p0x
    py = p0y

    put_text( ax, px*n, py*m, 'Total', mtot )
    py = py - offset

    put_text( ax, px*n, py*m, 'Dark', mdm )
    py = py - offset

    put_text( ax, px*n, py*m, 'Gas', mgas )
    py = py - offset

    put_text( ax, px*n, py*m, 'Star', ms )
    py = py - offset

    put_text( ax, px*n, py*m, 'BH', mbh )




fn_out = fn[:-4] + '.png'

if args.eps:
    fn_out = fn[:-4] + '.eps'

print( "save image to " +  fn_out )
plt.savefig( fn_out )
