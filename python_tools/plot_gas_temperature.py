#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )
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
for i in np.flip( np.linspace( YMin, YMax, NY ), 0 ):
    YList.append( "%.f"%(i/1000) )

fig = plt.figure()
ax = fig.add_subplot( 111 )
img = ax.imshow( data )

ax.set_xlabel( "Mpc" )
ax.set_ylabel( "" )
ax.set_title( "Gas Temperature(z=%.2f)"%(z) )

xloc = tik.LinearLocator(NX)
xfmt = tik.FixedFormatter(XList)
yloc = tik.LinearLocator(NY)
yfmt = tik.FixedFormatter(YList)
ax.xaxis.set_major_locator( xloc )
ax.xaxis.set_major_formatter( xfmt )
ax.yaxis.set_major_locator( yloc )
ax.yaxis.set_major_formatter( yfmt )

'''
ax.set_xticks( np.linspace( 0, n, NX ) )
ax.set_xticklabels( XList )
ax.set_yticks( np.linspace( 0, m, NY ) )
ax.set_yticklabels( YList )
'''
cbar = fig.colorbar( img )
cbar.set_label( r"$10^x Kcm$" )

plt.savefig( sys.argv[2] )
#plt.show()
