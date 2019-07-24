#!/usr/bin/env python3

from my_work_env import *

cmaps = [
        plt.get_cmap( 'ds9b' ),\
        plt.get_cmap( 'ds9a' ),\
        cm.seismic,\
        cm.PiYG,\
        cm.gist_heat,\
        plt.get_cmap( 'ds9cool' ),\
        plt.get_cmap( 'ds9b' ),\
        plt.get_cmap( 'ds9heat' ),\
        cm.nipy_spectral,\
        plt.get_cmap( 'ds9a' ),\
        plt.get_cmap( 'ds9he' ),\
        plt.get_cmap( 'ds9cool' ),\
        plt.get_cmap( 'ds9rainbow' ),\
        ]

FileName0 = data_dir + '/Density_141.dat'
FileName1 = data_dir + '/temp_141.dat'
FileName2 = data_dir + '/cre_temp_141.dat'
data0 = np.loadtxt( FileName0 )[1:,:]
data1 = np.loadtxt( FileName1 )[1:,:]
data2 = np.loadtxt( FileName2 )[1:,:]
data_all = [ data0, data1, data2 ]

norms = [
        mplc.LogNorm(),
        mplc.LogNorm(),
        mplc.SymLogNorm( linthresh=1e-2 )
        ]

N = 3
l = 5
ty = 0.1
fig = plt.figure( figsize=( l*N, l/(1-ty) ) )
dx = 1/N
dy = ( 1-ty )
axs = [ fig.add_axes( [ i*dx, ty, dx, dy ] ) for i in range(N) ]

tt = 1/3
for i in range(N):
    axs.append( fig.add_axes( [i*dx, ty*(1-tt), dx, ty*tt] ) )

for a in axs:
    a.set_xticks( [] )
    a.set_yticks( [] )
    a.minorticks_off()

data_all[2] = (data_all[2] - data_all[1]) / data_all[1]
m, n = data_all[0].shape

a = 10
a_str = '>10: '
print( a_str, len(data_all[2][data_all[2]>a]) / ( m*n )  )
data_all[2][data_all[2]>a] = np.nan


#a = 1e-2
#a_str = '<1e-2: '
#print( a_str, len(data_all[2][data_all[2]<a]) / ( m*n )  )
#data_all[2][data_all[2]<a] = np.nan

#data_all[2][np.isnan(data_all[2])] = \
#data_all[2][ np.logical_not( np.isnan( data_all[2] ) ) ].min() * 1e-2

r = 1e-7
for i in range(N):
    if i == 0:
        v = data_all[i].max() * r
        data_all[i][ data_all[i]<v ] = v
    else:
        data_all[i][ data_all[0] == v ] = np.nan
    img = axs[i].imshow( data_all[i], norm=norms[i], cmap=cmaps[i] )
    plt.colorbar( img, cax=axs[i+N], orientation='horizontal' )
    axs[i+N].tick_params( axis='x', direction='in', labelsize=10, top=True,
            bottom=True )
    axs[i+N].minorticks_off()

plt.savefig( figs_dir + 'temp-diff.pdf', dpi=300 )
exit()

