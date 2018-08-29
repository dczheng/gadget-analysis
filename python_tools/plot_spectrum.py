#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker as tik
#import matplotlib.colors as mplc
#from matplotlib import cm
import sys
import os

if ( len( sys.argv ) != 2 ):
    print( "Please give data file!" )
    exit()

DataFileList = os.listdir( sys.argv[1] )

for f in DataFileList:

    if f[-4:] != '.dat':
        continue

    print( 'plot %s ...'%(f) )

    d = np.loadtxt( f )
    plt.loglog( d[0,:], d[1,:] )
    plt.loglog( d[0,:], d[1,:], '*' )
    plt.xlabel( r'$\nu \; [MHz]$' )
    plt.ylabel( r'$ I_{\nu}\; [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1} ]$' )
    #plt.legend()

    t = f.split( '_' )[-1][:-4]
    if t == 'Tot':
        plt.title( 'Total' )
    else:
        plt.title( 'Redshift %s'%t )

    png_fn = f[:-4] + '.png'
    print( 'save fig to ' + png_fn )
    plt.savefig( png_fn )
    plt.close()

    print('')
