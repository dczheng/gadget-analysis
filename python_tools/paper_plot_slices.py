#!/usr/bin/env python3

from my_work_env import *

data_dir = sys.argv[1]
fn_out = sys.argv[2]
#matplotlib.style.use( 'default' )

cmaps = [
        cm.jet,\
        cm.magma,\
        plt.get_cmap( 'ds9a' ),\
        cm.nipy_spectral,\
        #plt.get_cmap( 'ds9b' ),\
        cm.plasma,\
        cm.hot,\
        ]
Names = [
        r'$\rm {\rho}/{\bar{\rho}}$', \
        r'$\rm Magnetic \, Field \; [\mu G]$', \
        r'$\rm Mach \, Number$', \
        r'$\rm {\epsilon}/{\epsilon_{\rm bar}}$'
        ]

norms = [
        mplc.LogNorm,\
        mplc.LogNorm,\
        mplc.LogNorm,\
        mplc.LogNorm
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

vmaxs = []
vmins = []
for j in range(N):
    for i in range(4):
        d = data_all[j][i]
        d[d==0] = d[d>0].min() / 100

    vmaxs.append(np.max( [ data_all[i][j].max() for i in range(4)] ))
    vmins.append(np.min( [ data_all[i][j][data_all[i][j]>0].min() for i in range(4)] ))

    if "Mach" in Names[j]:
        print( "Mach" )
        vmins[j] = 1
        mmax =300 
        for i in range(4):
            d = data_all[i][j]
            print( d.min(), d.max() )
            n1 = len(d[d>mmax])
            n2 = len(d[d>vmins[j]])
            print( n1, n2, (n1/n2)*100 )
        vmaxs[j] = mmax

    if "rho" in Names[j]:
        pass

    if "Magnetic" in Names[j]:
        vmins[j] = 1e-9

    if "epsilon" in Names[j]:
        vmins[j] = vmaxs[j] * 1e-5

for i in range(4):
    for j in range(N):
        print( data_all[i][j].shape )

m, n = data1[0].shape

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

wh = int(10000.0 / 100000 * m)
rxy= [ (n-wh)//2,(m-wh)//2 ]

fx = 0.15 * n
fy = 0.8 * m

for i in range(N):

    print( '-' * 30 )
    print( Names[i] )
    norm = norms[i]
    cmap = cmaps[i]
    Name = Names[i]
    v1 = vmins[i]
    v2 = vmaxs[i]
    print( v1, v2 )

    #----------------------------------#
    ax = ax_up[i]
    data = data1[i]

    if norm:
        ax.imshow( data, norm = norm(vmin=v1, vmax=v2), cmap=cmap )
    else:
        ax.imshow( data, vmin=v1, vmax=v2, cmap=cmap )

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
    if norm:
        ax.imshow( data, norm = norm(vmin=v1, vmax=v2), cmap=cmap )
    else:
        ax.imshow( data, vmin=v1, vmax=v2, cmap=cmap )

    #----------------------------------#
    ax = ax_middle[i]
    bar_ax = ax_bottom[i]
    data = data2[i]

    if norm:
        img = ax.imshow( data, norm = norm(vmin=v1, vmax=v2), cmap=cmap )
    else:
        img = ax.imshow( data, vmin=v1, vmax=v2, cmap=cmap )
    ax.add_patch( matplotlib.patches.Rectangle( rxy, wh, wh, color='w', fill=False ) )
    cbar = plt.colorbar( img, cax = bar_ax, orientation='horizontal' )

    set_tick_params( bar_ax, 25 )
    bar_ax.set_xlabel( Name, fontsize=30 )

    if "Magnetic" in Names[i]:
        make_log_ticks( v1, v2, 2, a=2, axis=bar_ax.xaxis )

    if "Mach" in Names[i]:
        make_log_ticks( v1, v2, 2, a=1, axis=bar_ax.xaxis )

    if "rho" in Names[i]:
        make_log_ticks( v1, v2, 2, a=2, axis=bar_ax.xaxis )

    if "epsilon" in Names[i]:
        make_log_ticks( v1, v2, 2, a=1, axis=bar_ax.xaxis )

    ax = ax_middle_little[i]
    data = data2_little[i]

    if norm:
        ax.imshow( data, norm = norm(vmin=v1, vmax=v2), cmap=cmap )
    else:
        ax.imshow( data, vmin=v1, vmax=v2, cmap=cmap )

    #----------------------------------#

for ax in ax_all:
    for a in ax:
        a.invert_yaxis()


plt.savefig( fn_out, dpi=300 )
#plt.show()

