#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mplc
from matplotlib import cm
import matplotlib.ticker as tik

aerr = np.loadtxt( './find_root_alpha_err.txt' )
qerr = np.loadtxt( './find_root_q_err.txt' )
Aerr = np.loadtxt( './find_root_A_err.txt' )
Berr = np.loadtxt( './find_root_B_err.txt' )
aa = np.loadtxt( './find_root_aa.txt' )
qq = np.loadtxt( './find_root_qq.txt' )

amin = aerr[0,0]
amax = aerr[0,1]
qmin = np.log10(aerr[0,2])
qmax = np.log10(aerr[0,3])

aerr = aerr[ 1:, : ]
qerr = qerr[ 1:, : ]
Aerr = Aerr[ 1:, : ]
Berr = Berr[ 1:, : ]
aa = aa[ 1:, : ]
qq = qq[ 1:, : ]

am, an = aerr.shape
qm, qn = qerr.shape

print( am, an, qm, qn )

print( "alpha: %f %f"%(amin,amax) )
print( "q: %f %f"%(qmin,qmax) )

print( "aerr max: ", aerr.max() )
print( "qerr max: ", qerr.max() )
print( "Aerr max: ", Aerr.max() )
print( "Berr max: ", Berr.max() )
print( "aa max: ", aa.max() )
print( "qq max: ", qq.max() )


norm = mplc.LogNorm()

NX = 5
NY = 5

xmin = int( qmin ) + 1
xmax = int( qmax )


#--------------------------------------------------
xLocList = np.linspace( xmin, xmax, xmax-xmin+1 )
xFmtList = []
for i in xLocList:
    xFmtList.append( r"$10^{%.0f}$"%i )
xLocList = (xLocList - xmin) / ( xmax-xmin ) * (an-1)

yLocList = np.linspace( amin, amax, NY )
yFmtList = []
for i in yLocList:
    yFmtList.append( "%.1f"%i )
yLocList = (yLocList - amin) / ( amax-amin ) * (am-1)

print( 'xLocList:', xLocList )
print( 'xLocFmt:', xFmtList )
print( 'yLocList:', yLocList )
print( 'yLocFmt:', yFmtList )

def set_global_pic_params( ax ):
    xloc = tik.FixedLocator( xLocList )
    xfmt = tik.FixedFormatter( xFmtList )
    yloc = tik.FixedLocator( yLocList )
    yfmt = tik.FixedFormatter( yFmtList )

    ax.xaxis.set_major_locator( xloc )
    ax.xaxis.set_major_formatter( xfmt )
    ax.yaxis.set_major_locator( yloc )
    ax.yaxis.set_major_formatter( yfmt )

    ax.invert_yaxis()



#-------------------------------------------------

print( 'plot alpha error ...' )
fig = plt.figure()
ax = fig.add_subplot( 111 )

img_aerr = ax.imshow( aerr, norm=norm, cmap=cm.jet )

ax.set_title( r'$\alpha \; Error$' )
cbar = fig.colorbar( img_aerr )
cbar.set_label( r'$\frac{\alpha_{calc}-\alpha_0}{\alpha_0}$' )

set_global_pic_params( ax )
ax.set_xlabel( r'$q_{0}$' )
ax.set_ylabel( r'$\alpha_{0}$' )



plt.savefig( './find_root_alpha_err.png' )
plt.close()
#--------------------------------------------------

print( 'plot q error ...' )
fig = plt.figure()
ax = fig.add_subplot( 111 )

img_qerr = ax.imshow( qerr, norm=norm, cmap=cm.jet )

ax.set_title( r'$q \; Error $' )
cbar = fig.colorbar( img_qerr )
cbar.set_label( r'$\frac {q_{calc}-q_0}{q_0}$' )

set_global_pic_params( ax )
ax.set_xlabel( r'$q_{0}$' )
ax.set_ylabel( r'$\alpha_{0}$' )

ax.invert_yaxis()

plt.savefig( './find_root_q_err.png' )
plt.close()

#--------------------------------------------------

print( 'plot A error ...' )
fig = plt.figure()
ax = fig.add_subplot( 111 )

img_Aerr = ax.imshow( Aerr, norm=norm, cmap=cm.jet )

ax.set_title( r'$A \; Error $' )
cbar = fig.colorbar( img_Aerr )
cbar.set_label( r'$\frac {A_{calc}-A_{inj}}{A_{inj}}$' )

set_global_pic_params( ax )
ax.set_xlabel( r'$q_{inj}$' )
ax.set_ylabel( r'$\alpha_{inj}$' )


plt.savefig( './find_root_A_err.png' )
plt.close()

#--------------------------------------------------

print( 'plot B error ...' )
fig = plt.figure()
ax = fig.add_subplot( 111 )

img_Berr = ax.imshow( Berr, norm=norm, cmap=cm.jet )

ax.set_title( r'$B \; Error $' )
cbar = fig.colorbar( img_Berr )
cbar.set_label( r'$\frac {B_{calc}-B_{inj}}{B_{inj}}$' )
ax.set_xlabel( r'$q_{inj}$' )
ax.set_ylabel( r'$\alpha_{inj}$' )

set_global_pic_params( ax )

plt.savefig( './find_root_B_err.png' )
plt.close()

#--------------------------------------------------

print( 'plot aa ...' )
fig = plt.figure()
ax = fig.add_subplot( 111 )

#img_aa = ax.imshow( aa, norm=norm, cmap=cm.jet )
img_aa = ax.imshow( aa, cmap=cm.jet )

ax.set_title( r'$\alpha$' )
cbar = fig.colorbar( img_aa )
cbar.set_label( r'$\alpha$' )
ax.set_xlabel( r'$q_{inj}$' )
ax.set_ylabel( r'$\alpha_{inj}$' )

set_global_pic_params( ax )

plt.savefig( './find_root_aa.png' )
plt.close()


#--------------------------------------------------

print( 'plot qq ...' )
fig = plt.figure()
ax = fig.add_subplot( 111 )

#img_qq = ax.imshow( qq, norm=norm, cmap=cm.jet )
img_qq = ax.imshow( qq, cmap=cm.jet )

ax.set_title( r'$q$' )
cbar = fig.colorbar( img_qq )
cbar.set_label( r'$q$' )

set_global_pic_params( ax )
ax.set_xlabel( r'$q_{inj}$' )
ax.set_ylabel( r'$\alpha_{inj}$' )

plt.savefig( './find_root_qq.png' )
plt.close()


