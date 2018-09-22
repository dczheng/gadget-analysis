#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mplc
from matplotlib import cm
import matplotlib.ticker as tik

beta = np.loadtxt( './beta.txt' )

amin = beta[0,0]
amax = beta[0,1]
qmin = np.log10( beta[0,2] )
qmax = np.log10( beta[0,3] )

beta = beta[1:, :]

print( "amin: %g, amax: %g, qmin: %g, qmax: %g"%( amin, amax, qmin, qmax ) )

m,n = beta.shape
print( "( %i, %i )"%(m,n) )

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


fig = plt.figure()
ax = fig.add_subplot( 111 )

xloc = tik.FixedLocator( xLocList )
xfmt = tik.FixedFormatter( xFmtList )
yloc = tik.FixedLocator( yLocList )
yfmt = tik.FixedFormatter( yFmtList )

ax.xaxis.set_major_locator( xloc )
ax.xaxis.set_major_formatter( xfmt )
ax.yaxis.set_major_locator( yloc )
ax.yaxis.set_major_formatter( yfmt )


norm = mplc.LogNorm()
img = ax.imshow( beta, norm=norm, cmap=cm.jet )
#img = ax.imshow( beta, cmap=cm.jet )

ax.invert_yaxis()

cbar = fig.colorbar( img )
cbar.set_label( r'$\beta_{\frac{1}{1+q^2}}\left(\frac{\alpha-2}{2}, \frac{3-\alpha}{2}\right)$' )

plt.savefig( './beta.png' )
plt.close()
