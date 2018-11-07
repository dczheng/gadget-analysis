#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tik
import matplotlib.colors as mplc
import argparse
from matplotlib import cm

parser = argparse.ArgumentParser()

parser.add_argument( 'f' )
parser.add_argument( '-pdf', help='save pdf', action='store_true' )

args = parser.parse_args()

data = np.loadtxt( args.f)

data_info = data[ 0, : ]
data = np.flipud(data[ 1:, : ])
DensMin = data_info[1]
DensMax = data_info[2]
BMin = data_info[3]
BMax = data_info[4]
GlobDensMin = data_info[5]
GlobDensMax = data_info[6]
GlobBMin = data_info[7]
GlobBMax = data_info[8]
z = data_info[10]
vmax = data.max()
#data[ data < 1e-3 * vmax ] = 0
#print( data_info )
m,n = data.shape
print( "PicSize: (%i,%i)"%( m, n ) )
print( "DensMin: %g, DensMax: %g, BMin: %g, BMax: %g"%\
        (DensMin, DensMax, BMin, BMax ) )
#print( data_info )

NDens = 5
NB = 5

'''
DensList = []
BList = []
for i in np.linspace( DensMin, DensMax, NDens ):
    DensList.append( "%.2f"%(i) )
for i in np.flip( np.linspace( BMin, BMax, NB ), 0 ):
    BList.append( "%.2f"%(i) )
'''

xmin = int( DensMin )
xmax = int( DensMax )
xl = DensMax - DensMin
ymin = int( BMin ) + 1
ymax = int( BMax )
yl = BMax - BMin

print( "xmin: %i, xmax: %i, ymin: %i, ymax: %i"%( xmin, xmax, ymin ,ymax ) )

if ( xmin > xmax ):
    pass
else:
    xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
    xFmtList = xLocList
    xLocList = ( xLocList - DensMin ) / xl * (n-1)
    print( "xticks: ", xFmtList )
if ( ymin > ymax ):
    pass
else:
    yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
    yFmtList = yLocList
    yLocList = m - ( yLocList - BMin ) / yl * (m-1)
    print( "yticks: ", yFmtList )



fig = plt.figure()
ax = fig.add_subplot( 111 )
norm = mplc.LogNorm()
img = ax.imshow( data, norm=norm, cmap=cm.jet )

ax.set_xlabel( r"$log[\frac{\rho}{\rho_{Bar}}]$" )
ax.set_ylabel( r"$log[B/\mu G]$" )
ax.set_title( "z=%.2f"%(z) )

xloc = tik.FixedLocator(xLocList)
xfmt = tik.FixedFormatter(xFmtList)
yloc = tik.FixedLocator(yLocList)
yfmt = tik.FixedFormatter(yFmtList)
ax.xaxis.set_major_locator( xloc )
ax.xaxis.set_major_formatter( xfmt )
ax.yaxis.set_major_locator( yloc )
ax.yaxis.set_major_formatter( yfmt )

'''
ax.set_xticks( np.linspace( 0, n, NDens ) )
ax.set_xticklabels( DensList )
ax.set_yticks( np.linspace( 0, m, NB ) )
ax.set_yticklabels( BList )
'''
cbar = fig.colorbar( img )
cbar.set_label( "Probability" )

if args.pdf:
    fn = args.f[0:-4] + '.pdf'
else:
    fn = args.f[0:-4] + '.png'

print( 'save image to ' + fn )
plt.savefig( fn )

