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
DensMin = data_info[1]
DensMax = data_info[2]
TempMin = data_info[3]
TempMax = data_info[4]
#print( data_info )
m,n = data.shape
print( "PicSize: (%i,%i)"%( m, n ) )
print( "DensMin: %g, DensMax: %g, TempMin: %g, TempMax: %g"%\
        (DensMin, DensMax, TempMin, TempMax ) )
#print( data_info )

NDens = 5
NTemp = 5

'''
DensList = []
TempList = []
for i in np.linspace( DensMin, DensMax, NDens ):
    DensList.append( "%.2f"%(i) )
for i in np.flip( np.linspace( TempMin, TempMax, NTemp ), 0 ):
    TempList.append( "%.2f"%(i) )
'''

xmin = int( DensMin ) + 1
xmax = int( DensMax )
xl = DensMax - DensMin
ymin = int( TempMin ) + 1
ymax = int( TempMax )
yl = TempMax - TempMin

print( "xmin: %i, xmax: %i, ymin: %i, ymax: %i"%( xmin, xmax, ymin ,ymax ) )

if ( xmin > xmax ):
    pass
else:
    xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
    xFmtList = xLocList
    xLocList = ( xLocList - DensMin ) / xl * n
    print( "xticks: ", xFmtList )
if ( ymin > ymax ):
    pass
else:
    yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
    yFmtList = np.flip( yLocList, 0 )
    yLocList = ( yLocList - TempMin ) / yl * m
    print( "yticks: ", yFmtList )



fig = plt.figure()
ax = fig.add_subplot( 111 )
norm = mplc.LogNorm()
img = ax.imshow( data, norm=norm, cmap=cm.jet )

ax.set_xlabel( r"$\frac{\rho}{\rho_{Bar}}(10^x)$" )
ax.set_ylabel( r"$Temperature (10^x K)$" )
ax.set_title( "Gas State(z=%.2f)"%(z) )

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
ax.set_yticks( np.linspace( 0, m, NTemp ) )
ax.set_yticklabels( TempList )
'''
cbar = fig.colorbar( img )
cbar.set_label( "Probability" )

plt.savefig( sys.argv[2] )
#plt.show()