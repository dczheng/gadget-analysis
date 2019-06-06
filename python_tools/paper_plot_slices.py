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

cmaps = [
        cm.viridis, \
        #cm.gnuplot2, \
        cm.magma, \
        cm.viridis, \
        cm.viridis, \
        cm.viridis \
        ]
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
        data_dir + '/Density_0.00.dat', \
        data_dir + '/Density_0.00.dat', \
        data_dir + '/Density_0.00.dat'
                ]

norms = [
        mplc.LogNorm(),\
        #        None, \
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm()
        ]

#data = [ np.loadtxt( f ) for f in FileNames ]
data = [ np.loadtxt( FileNames[i] ) for i in range(2) ]
data.append( data[0] )
data.append( data[0] )
data.append( data[0] )

data = [ d[1:, :] for d in data ]

t = [ 0 ] * N
t[1] = 1
for i in range(len(data)):
    data[i][data[i]<=t[i]] = data[i][data[i] > t[i] ].min()

matplotlib.style.use( 'default' )
plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

N = 5
fig = plt.figure( figsize=(4*N, 5) )
dx = 1 / N
axs = [ fig.add_axes([i*dx, 1/10, dx, 4/5]) for i in range(N) ]

for i in range(N):
    axs.append( fig.add_axes( [i*dx, 2/30, dx, 1/30 ] ) )

for i in range(N):
    axs.append( fig.add_axes( [i*dx+0.01, 1/10+0.01*N, 0.4*dx, 0.4*4/5 ] ) )

for i in range(N):
    ax = axs[i]
    bar_ax = axs[i+N]

    print( data[i].min(), data[i].max() )
    img = ax.imshow( data[i], norm = norms[i], cmap=cmaps[i] )
    cbar = plt.colorbar( img, cax = bar_ax, orientation='horizontal' )
    cbar.ax.tick_params( direction='in', labelsize=15 )
    #help( cbar.ax.tick_params )
    cbar.ax.minorticks_off()
    ax.set_title( Names[i], fontsize=20 )
    ax.set_yticks( [] )
    ax.set_xticks( [] )

    #bar_ax.spines['top'].set_color( 'red' )
    #bar_ax.spines['bottom'].set_color( 'red' )
    #bar_ax.spines['left'].set_color( 'red' )
    #bar_ax.spines['right'].set_color( 'red' )

m, n = data[0].shape
wh = 50
rxy= [ n-wh-85,160 ]

for i in range(N):
    axs[i].add_patch( matplotlib.patches.Rectangle( rxy, wh, wh, color='w', fill=False ) )


for i in range(N):

    ax = axs[i+N*2]

    img = ax.imshow( data[i][rxy[1]:rxy[1]+wh:,rxy[0]:rxy[0]+wh], norm = norms[i], cmap=cmaps[i] )

    #cbar = plt.colorbar( img, ax = ax, pad=0, fraction=0.0475, orientation='horizontal' )
    #cbar.ax.tick_params( direction='in', width=0.3, length=1.5, labelsize=10 )
    ax.spines['top'].set_color( 'white' )
    ax.spines['bottom'].set_color( 'white' )
    ax.spines['left'].set_color( 'white' )
    ax.spines['right'].set_color( 'white' )
    ax.set_yticks( [] )
    ax.set_xticks( [] )


plt.savefig( figs_dir + 'slices.pdf' )

