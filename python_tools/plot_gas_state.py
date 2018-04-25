#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tik
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

DensList = []
TempList = []
for i in np.linspace( DensMin, DensMax, NDens ):
    DensList.append( "%.2f"%(i) )
for i in np.flip( np.linspace( TempMin, TempMax, NTemp ), 0 ):
    TempList.append( "%.2f"%(i) )

fig = plt.figure()
ax = fig.add_subplot( 111 )
img = ax.imshow( data )

ax.set_xlabel( "Density ($10^x gcm^-3$)" )
ax.set_ylabel( "Temperature ($10^x$ K)" )
ax.set_title( "Gas State(z=%.2f)"%(z) )

xloc = tik.LinearLocator(NDens)
xfmt = tik.FixedFormatter(DensList)
yloc = tik.LinearLocator(NTemp)
yfmt = tik.FixedFormatter(TempList)
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
cbar.set_label( "Number $(10^x)$" )

plt.show()
