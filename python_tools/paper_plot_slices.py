#!/usr/bin/env python3

from my_work_env import *


#my_cmap = cm.jet
#my_cmap = cm.Set1
#my_cmap = cm.cubehelix
#my_cmap = cm.ocean
#my_cmap = cm.magma
#my_cmap = cm.Paired
#my_cmap = cm.tab10
#my_cmap = cm.tab20
#my_cmap = cm.tab20b
#my_cmap = cm.tab20c
#my_cmap = cm.Set2
#my_cmap = cm.Set3
#my_cmap = cm.Set1
#my_cmap = cm.viridis
#my_cmap = cm.PiYG
#my_cmap = cm.BrBG
#my_cmap = cm.RdYlBu
#my_cmap = cm.RdBu
#my_cmap = cm.Spectral

N = 5
Names = [
        r'$\rm Density$', \
        r'$\rm Mach \, Number$', \
        r'$\rm Magnetic \, Field$', \
        r'$\rm CRE \,Number \, Density$',\
        r'$\rm Radio$'\
        ]
FileNames = [
        data_dir + '/Density_0.00.dat', \
        data_dir + '/MachNumber_0.00.dat', \
        data_dir + '/256_Density_0.20.dat', \
        data_dir + '/256_Density_0.20.dat', \
        data_dir + '/256_Density_0.20.dat'
                ]
cmaps = [
        cm.viridis, \
        #cm.jet, \
        cm.magma, \
        cm.viridis, \
        cm.viridis, \
        cm.viridis \
        ]

norms = [
        mplc.LogNorm(),\
        #        None, \
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm()
        ]

data = [ np.loadtxt( f ) for f in FileNames ]

data = [ d[1:, :] for d in data ]

t = 1e-20
for i in range(len(data)):
    data[i][data[i]<=t] = data[i][data[i] > t ].min()

matplotlib.style.use( 'default' )
plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

N = 5
fig = plt.figure( figsize=(4*N, 5) )
dx = 1 / N
axs = [ fig.add_axes([i*dx, 1/10, dx, 4/5]) for i in range(N) ]
for i in range(N):
    axs.append( fig.add_axes( [i*dx+0.01, 1/10+0.01*N, 0.3*dx, 0.3*4/5 ] ) )

print( data[0].min(), data[0].max() )

for i in range(N):
    ax = axs[i]

    print( data[i].min(), data[i].max() )
    img = ax.imshow( data[i], norm = norms[i], cmap=cmaps[i] )
    #cbar = plt.colorbar( img, ax = ax, pad=0, fraction=0.0475, orientation='horizontal' )
    #cbar.ax.tick_params( direction='in', width=0.3, length=1.5, labelsize=10 )
    ax.set_title( Names[i] )
    ax.set_yticks( [] )
    ax.set_xticks( [] )

for i in range(N):

    ax = axs[i+N]
    img = ax.imshow( data[i], norm = norms[i], cmap=cmaps[i] )

    #cbar = plt.colorbar( img, ax = ax, pad=0, fraction=0.0475, orientation='horizontal' )
    #cbar.ax.tick_params( direction='in', width=0.3, length=1.5, labelsize=10 )
    ax.set_yticks( [] )
    ax.set_xticks( [] )

#fig.subplots_adjust(wspace=-1)
#fig.tight_layout()
plt.savefig( figs_dir + 'slices.pdf' )

