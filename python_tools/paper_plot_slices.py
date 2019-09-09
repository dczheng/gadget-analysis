#!/usr/bin/env python3

from my_work_env import *

data_dir = sys.argv[1]
fn_out = sys.argv[2]
#matplotlib.style.use( 'default' )

cmaps = [
        cm.hot,\
        #cm.jet,\
        #plt.get_cmap( 'ds9b' ),\
        cm.nipy_spectral,\
        plt.get_cmap( 'ds9a' ),\
        cm.plasma
        ]
Names = [
        r'$\rm {\rho}/{\bar{\rho}}$', \
        r'$\rm Magnetic \, Field \; [\mu G]$', \
        r'$\rm Mach \, Number$', \
        r'$\rm {\epsilon}/{\epsilon_{\rm bar}}$'
        ]

norms = [
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm(),\
        mplc.LogNorm()
        ]

FileNames1 = [
        data_dir + '/Density_0.1-1.dat', \
        data_dir + '/MagneticField_0.1-1.dat', \
        data_dir + '/MachNumber_0.1-1.dat', \
        data_dir + '/cre_e_0.1-1.dat', \
                ]
FileNames1_little = [
        data_dir + '/Density_0.1-2.dat', \
        data_dir + '/MagneticField_0.1-2.dat', \
        data_dir + '/MachNumber_0.1-2.dat', \
        data_dir + '/cre_e_0.1-2.dat'
                ]

FileNames2 = [
        data_dir + '/Density_0-1.dat', \
        data_dir + '/MagneticField_0-1.dat', \
        data_dir + '/MachNumber_0-1.dat', \
        data_dir + '/cre_e_0-1.dat', \
                ]
FileNames2_little = [
        data_dir + '/Density_0-2.dat', \
        data_dir + '/MagneticField_0-2.dat', \
        data_dir + '/MachNumber_0-2.dat', \
        data_dir + '/cre_e_0-2.dat'
                ]

#d = np.loadtxt( "./Density_diffuse.dat" )[1:,:]
#print( d.shape )
##vmin = d[ d>0 ].min()
##d[ d<vmin ] = vmin
#plt.imshow( d, norm=mplc.LogNorm(), cmap=cm.hot )
#ax = plt.gca()
#ax.set_xticks( [] )
#ax.set_yticks( [] )
#ax.invert_yaxis()
#plt.savefig( fn_out, dpi=300 )
#exit()

N =  len( FileNames1 )

fs = 5 
t_cbar = 0.25
fig = plt.figure( figsize=(fs*N, fs*(2+t_cbar)) )
dx = 1 / N
dy = 1 / ( 2 + t_cbar )
r_cbar = 0.3
ax_up = [ fig.add_axes([i*dx, dy*(1+t_cbar), dx, dy]) for i in range(N) ]
ax_middle = [ fig.add_axes([i*dx, dy*t_cbar, dx, dy]) for i in range(N) ]
ax_bottom = [ fig.add_axes([i*dx, dy*t_cbar*(1-r_cbar), dx, dy*t_cbar*r_cbar])\
                for i in range(N) ]
r_little = 0.4
p_little = 0.05
ax_up_little = [ fig.add_axes([i*dx+dx*p_little, \
                               dy*(1+t_cbar) + dy*p_little,\
                               dx*r_little, \
                               dy*r_little ]) \
                    for i in range(N) ]
ax_middle_little = [ fig.add_axes([i*dx+dx*p_little, \
                                   dy*t_cbar + dy*p_little,\
                                   dx*r_little, \
                                   dy*r_little ]) \
                    for i in range(N) ]

ax_all = [ax_up, ax_up_little, ax_middle, ax_middle_little, ax_bottom]

data1 = [ np.loadtxt( FileNames1[i] ) for i in range(N) ]
data1_little = [ np.loadtxt( FileNames1_little[i] ) for i in range(N) ]
data2 = [ np.loadtxt( FileNames2[i] ) for i in range(N) ]
data2_little = [ np.loadtxt( FileNames2_little[i] ) for i in range(N) ]

data_all = [ data1, data1_little, data2, data2_little ]
for j in range(N):
    print( Names[j] )
    for i in range(4):
        data_all[i][j] = data_all[i][j][1:,:]
        print( data_all[i][j][data_all[i][j]>0].min(), data_all[i][j].max() )

for j in range(N):
    vmax = np.max( [ data_all[i][j].max() for i in range(4)] )
    vmin = np.min( [ data_all[i][j][data_all[i][j]>0].min() for i in range(4)] )
    if "Mach" in Names[j]:
        v = 0.9

    if "rho" in Names[j]:
        v = vmin * 0.1

    if "Magnetic" in Names[j]:
        v = vmax * 1e-12

    if "epsilon" in Names[j]:
        v = vmax * 1e-12

    for i in range( 4 ):
        data_all[i][j][ data_all[i][j]<=v ] = v
        data_all[i][j][0,0] = vmax

for i in range(4):
    for j in range(N):
        print( data_all[i][j].shape )

m, n = data1[0].shape
for j in range(N):
    print( Names[j] )
    for i in range(4):
        print( data_all[i][j].min(), data_all[i][j].max() )

for ax in ax_all:
    for a in ax:
        a.set_xticks( [] )
        a.set_yticks( [] )

for i in range(4):
    for j in range(N):
        if i % 2 == 1:
            ax_all[i][j].spines['top'].set_color( 'white' )
            ax_all[i][j].spines['bottom'].set_color( 'white' )
            ax_all[i][j].spines['left'].set_color( 'white' )
            ax_all[i][j].spines['right'].set_color( 'white' )
        else:
            ax_all[i][j].spines['top'].set_color( 'black' )
            ax_all[i][j].spines['bottom'].set_color( 'black' )
            ax_all[i][j].spines['left'].set_color( 'black' )
            ax_all[i][j].spines['right'].set_color( 'black' )

wh = int(5000 / 100000 * m)
rxy= [ (n-wh)//2,(m-wh)//2 ]

fx = 0.15 * n
fy = 0.8 * m

for i in range(N):

    print( Names[i] )
    norm = norms[i]
    cmap = cmaps[i]
    Name = Names[i]

    #----------------------------------#
    ax = ax_up[i]
    data = data1[i]

    ax.imshow( data, norm = norm, cmap=cmap )
    ax.add_patch( matplotlib.patches.Rectangle( rxy, wh, wh, color='w', fill=False ) )
    '''
    if "Mach" in Names[i]:
        ax.text( fx, fy, Names[i], color='white', fontsize=30, weight='bold' )
    if "rho" in Names[i]:
        ax.text( fx, fy, Names[i], color='white', fontsize=60, weight='bold' )
    if "Magnetic" in Names[i]:
        ax.text( fx, fy, Names[i], color='white', fontsize=30, weight='bold' )
    if "epsilon" in Names[i]:
        ax.text( fx, fy, Names[i], color='white', fontsize=60, weight='bold' )
    '''

    ax = ax_up_little[i]
    data = data1_little[i]
    ax.imshow( data, norm = norm, cmap=cmap )

    #----------------------------------#
    ax = ax_middle[i]
    bar_ax = ax_bottom[i]
    data = data2[i]

    vmin = data[data>0].min()
    vmax = data.max()
    print( vmin, vmax )
    img = ax.imshow( data, norm = norm, cmap=cmap )
    ax.add_patch( matplotlib.patches.Rectangle( rxy, wh, wh, color='w', fill=False ) )
    cbar = plt.colorbar( img, cax = bar_ax, orientation='horizontal' )
    vmin = int(np.log10(vmin))
    vmax = int(np.log10(vmax))
    if "rho" in Names[i]:
        t = range( vmin, vmax, 2 )
    if "Mach" in Names[i]:
        t = range( vmin, vmax+1 )
    if "Magnetic" in Names[i]:
        t = range( vmin, vmax+1, 2 )
    if "epsilon" in Names[i]:
        t = range( vmin, vmax+1, 2 )
    tik = [ 10**v for v in t ]
    print( tik )
    cbar.set_ticks( tik )
    set_tick_params( bar_ax, 25 )
    bar_ax.minorticks_off()
    bar_ax.set_xlabel( Name, fontsize=30 )

    ax = ax_middle_little[i]
    data = data2_little[i]
    img = ax.imshow( data, norm = norm, cmap=cmap )

    #----------------------------------#

for ax in ax_all:
    for a in ax:
        a.invert_yaxis()


plt.savefig( fn_out, dpi=300 )
#plt.show()

