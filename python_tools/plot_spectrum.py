#!/usr/bin/env python3
import matplotlib as mpl
mpl.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker as tik
#import matplotlib.colors as mplc
#from matplotlib import cm
import sys

if ( len( sys.argv ) != 2 ):
    print( "Please give data file!" )
    exit()

plt.xlabel( r'$\nu \; [MHz]$' )
plt.ylabel( r'$ I_{\nu}\; [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1} ]$' )

png_fn = f[:-4] + '.png'
print( 'save fig to ' + png_fn )
plt.savefig( png_fn )
plt.close()

    print('')
