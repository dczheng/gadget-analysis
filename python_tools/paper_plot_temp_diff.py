#!/usr/bin/env python3

from my_work_env import *

cmaps = [
        plt.get_cmap( 'ds9b' ),\
        plt.get_cmap( 'ds9a' ),\
        cm.gist_heat,\
        cm.nipy_spectral,\
        plt.get_cmap( 'ds9heat' ),\
        plt.get_cmap( 'ds9a' ),\
        plt.get_cmap( 'ds9he' ),\
        plt.get_cmap( 'ds9cool' ),\
        plt.get_cmap( 'ds9rainbow' ),\
        ]

FileNames1 = data_dir + '/temp_141.dat'
FileNames2 = data_dir + '/cre_temp_141.dat'
data1 = np.loadtxt( FileNames1 )
data2 = np.loadtxt( FileNames2 )
N =  len( FileNames1 )

#FileNames1_little = data_dir + '/temp_zoom_141.dat'
#FileNames2_little = data_dir + '/cre_temp_zoom_141.dat'
#data1_little = np.loadtxt( FileNames1_little )
#data2_little = np.loadtxt( FileNames2_little )

norm = mplc.SymLogNorm( linthresh=1e-3 )
cmap = cmaps[0]

N = 3
fig = plt.figure(  )
fig, ax = plt.subplots( 1, 1, figsize=(5,5) )
ax.set_xticks( [] )
ax.set_yticks( [] )

dd = data1 - data2
dd = dd / data1
print( dd.min(), dd.max() )

img = ax.imshow( dd, norm=norm, cmap=cmap )
plt.colorbar( img, ax=ax )

plt.savefig( figs_dir + 'temp-diff.pdf', dpi=300 )
#plt.show()

