#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tik
import matplotlib.colors as mplc
import sys

if ( len( sys.argv ) != 2 ):
    print( "Please give data file!" )
    exit()


data = np.loadtxt( sys.argv[1] )
data_info = data[ 0, : ]
data = data[ 1:,: ]
z = data_info[0]
XMin = data_info[1]
XMax = data_info[2]
YMin = data_info[3]
YMax = data_info[4]
#print( data_info )
m,n = data.shape
print( "PicSize: (%i,%i)"%( m, n ) )
print( "XMin: %g, XMax: %g, YMin: %g, YMax: %g"%\
        (XMin, XMax, YMin, YMax ) )
#print( data_info )

NX = 5
NY = 5

XList = []
YList = []
for i in np.linspace( XMin, XMax, NX ):
    XList.append( "%.f"%(i/1000) )
for i in np.linspace( YMin, YMax, NY ):
    YList.append( "%.f"%(i/1000) )

fig = plt.figure()
ax = fig.add_subplot( 111 )

norm = mplc.LogNorm()
img = ax.imshow( data, norm=norm )

ax.set_xlabel( "Mpc" )
ax.set_ylabel( "" )
ax.set_title( "Gas Density(z=%.2f)"%(z) )

xloc = tik.LinearLocator(NX)
xfmt = tik.FixedFormatter(XList)
yloc = tik.LinearLocator(NY)
yfmt = tik.FixedFormatter(YList)
ax.xaxis.set_major_locator( xloc )
ax.xaxis.set_major_formatter( xfmt )
ax.yaxis.set_major_locator( yloc )
ax.yaxis.set_major_formatter( yfmt )

ax.invert_yaxis()

cbar = fig.colorbar( img )
cbar.set_label( r"$gcm^{-2}$" )

fn_png = sys.argv[1][:-4] + '.png'
print( 'save image to ' + fn_png )
plt.savefig( fn_png )
#plt.show()
