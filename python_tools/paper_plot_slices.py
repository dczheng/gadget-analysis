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
Names = [
        r'$\rm Density \; [gcm^{-3}]$', \
        r'$\rm Mach \, Number$', \
        r'$\rm Magnetic \, Field \; [G]$', \
        r'$\rm CRE \,Number \, Density \; [cm^{-3}]$'
        ]

norms = [
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm()
        ]

FileNames1 = [
        data_dir + '/Density_131.dat', \
        data_dir + '/MachNumber_131.dat', \
        data_dir + '/MagneticField_131.dat', \
        data_dir + '/cre_n_131.dat', \
                ]
FileNames1_little = [
        data_dir + '/Density2_131.dat', \
        data_dir + '/MachNumber2_131.dat', \
        data_dir + '/MagneticField2_131.dat', \
        data_dir + '/cre_n2_131.dat'
                ]

FileNames2 = [
        data_dir + '/Density_141.dat', \
        data_dir + '/MachNumber_141.dat', \
        data_dir + '/MagneticField_141.dat', \
        data_dir + '/cre_n_141.dat', \
                ]
FileNames2_little = [
        data_dir + '/Density2_141.dat', \
        data_dir + '/MachNumber2_141.dat', \
        data_dir + '/MagneticField2_141.dat', \
        data_dir + '/cre_n2_141.dat'
                ]
#FileNames1 = [ data_dir + '/temp_141.dat' ]
#FileNames1_little = [ data_dir + '/temp_zoom_141.dat' ]
#FileNames2 = [ data_dir + '/cre_temp_141.dat' ]
#FileNames2_little = [ data_dir + '/cre_temp_zoom_141.dat' ]
N =  len( FileNames1 )


data1 = [ np.loadtxt( FileNames1[i] ) for i in range(N) ]
data1_little = [ np.loadtxt( FileNames1_little[i] ) for i in range(N) ]
data2 = [ np.loadtxt( FileNames2[i] ) for i in range(N) ]
data2_little = [ np.loadtxt( FileNames2_little[i] ) for i in range(N) ]

data_all = [ data1, data1_little, data2, data2_little ]
a = 1e-7
for i in range( 4 ):
    for j in range(N):
        data_all[i][j] = data_all[i][j][1:,:]
        v = data_all[i][j].max() * a
        if j == 1: #Mach Number
            v = 0.9
        data_all[i][j][ data_all[i][j]<=v ] = v

for i in range(4):
    for j in range(N):
        print( data_all[i][j].shape )

m, n = data1[0].shape

for i in range( N ):
    a = np.min( [ data_all[j][i].min() for j in range(4)] )
    b = np.max( [ data_all[j][i].max() for j in range(4)] )
    print( a, b )
    for j in range(4):
        data_all[j][i][0,0] = a
        data_all[j][i][0,1] = b

for j in range(N):
    for i in range(4):
        print( data_all[i][j][0,0], data_all[i][j][0,1] )


matplotlib.style.use( 'default' )
plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

fig = plt.figure( figsize=(4*N, 9) )
dx = 1 / N
ax_up = [ fig.add_axes([i*dx, 5/9, dx, 4/9]) for i in range(N) ]
t1 = 0.04
t2 = 0.4
ax_up_little = [ fig.add_axes([i*dx+t1/4, 5/9+t1*4/9, t2*dx, t2*4/9]) for i in range(N) ]
ax_middle = [ fig.add_axes([i*dx, 1/9, dx, 4/9]) for i in range(N) ]
ax_middle_little = [ fig.add_axes([i*dx+t1/4, 1/9+t1*4/9, t2*dx, t2*4/9]) for i in range(N) ]

t = 1/4
ax_bottom = [ fig.add_axes([i*dx, (1-t)/9, dx, t/9]) for i in range(N) ]

ax_all = [ax_up, ax_up_little, ax_middle, ax_middle_little, ax_bottom]
for ax in ax_all:
    for a in ax:
        a.set_xticks( [] )
        a.set_yticks( [] )

wh = int(5000 / 100000 * m)
rxy= [ (n-wh)//2,(m-wh)//2 ]

for i in range(2):
    for j in range(N):
        ax_all[i*2+1][j].spines['top'].set_color( 'white' )
        ax_all[i*2+1][j].spines['bottom'].set_color( 'white' )
        ax_all[i*2+1][j].spines['left'].set_color( 'white' )
        ax_all[i*2+1][j].spines['right'].set_color( 'white' )

for i in range(N):

    norm = norms[i]
    cmap = cmaps[i]
    Name = Names[i]

    #----------------------------------#
    ax = ax_up[i]
    data = data1[i]

    ax.imshow( data, norm = norm, cmap=cmap )
    ax.add_patch( matplotlib.patches.Rectangle( rxy, wh, wh, color='w', fill=False ) )

    ax = ax_up_little[i]
    data = data1_little[i]
    ax.imshow( data, norm = norm, cmap=cmap )

    #----------------------------------#
    ax = ax_middle[i]
    bar_ax = ax_bottom[i]
    data = data2[i]

    vmin = data.min()
    vmax = data.max()
    print( vmin, vmax )
    img = ax.imshow( data, norm = norm, cmap=cmap )
    ax.add_patch( matplotlib.patches.Rectangle( rxy, wh, wh, color='w', fill=False ) )
    cbar = plt.colorbar( img, cax = bar_ax, orientation='horizontal' )
    vmin = int(np.log10(vmin))
    vmax = int(np.log10(vmax))
    if i == 0:
        t = range( vmin, vmax )
    if i == 1:
        t = range( vmin, vmax+1 )
    if i == 2:
        t = range( vmin-1, vmax+1 )
    if i == 3:
        t = range( vmin-1, vmax-1 )
    tik = [ 10**v for v in t ]
    print( tik )
    cbar.set_ticks( tik )
    bar_ax.tick_params( axis='x', direction='in', labelsize=10, top=True, bottom=True )
    bar_ax.minorticks_off()
    bar_ax.set_xlabel( Name, fontsize=15 )

    ax = ax_middle_little[i]
    data = data2_little[i]
    img = ax.imshow( data, norm = norm, cmap=cmap )

    #----------------------------------#

plt.savefig( figs_dir + 'slices.pdf', dpi=300 )
#plt.show()

