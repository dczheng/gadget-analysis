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

FileNames1 = [
        data_dir + '/g_0.dat', \
        data_dir + '/g_1.dat', \
        data_dir + '/g_2.dat', \
        data_dir + '/g_3.dat' \
                ]

FileNames2 = [
        data_dir + '/cre_g_0.dat', \
        data_dir + '/cre_g_1.dat', \
        data_dir + '/cre_g_2.dat', \
        data_dir + '/cre_g_3.dat' \
                ]

N = len( FileNames1 )

l = 5
tx = 0.0
ty = 0.1
fig = plt.figure( figsize=(l*N/(1-tx), l*2/(1-ty) ) )
dx = (1-tx)/N
dy = (1-ty)/2
ax_up = [ fig.add_axes([i*dx, dy+ty, dx, dy]) for i in range(N) ]
ax_bottom = [ fig.add_axes([i*dx, ty, dx, dy]) for i in range(N) ]


tt = 1/3
ax_cbar = [fig.add_axes( [i*dx, ty*(1-tt), dx, ty*tt]) for i in range(N)]

ty = 0.2
fig_residul = plt.figure( figsize=(l*N, l/(1-ty) ) )
dx = 1/N
dy = (1-ty)
ax_residul = [ fig_residul.add_axes([i*dx, ty, dx, dy]) for i in range(N) ]
tt = 1/3
ax_cbar_residul = [fig_residul.add_axes( [i*dx, ty*(1-tt), dx, ty*tt]) for i in range(N)]

ax_all = np.hstack([ax_up, ax_bottom,  ax_cbar, ax_residul, ax_cbar_residul])
for ax in ax_all:
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    ax.minorticks_off()

data1 = [ np.loadtxt( FileNames1[i] ) for i in range(N) ]
data2 = [ np.loadtxt( FileNames2[i] ) for i in range(N) ]

data_all = [ d for d in data1 ]
for d in data2:
    data_all.append( d  )
Ls = []
Tmin = 1e6
Tmax = 1e8
TTmin = 1e20
TTmax = -1
Ms = []
for i in range( 2*N ):
    Ls.append( data_all[i][0,2] )
    Ms.append( data_all[i][0, 15:15+6].sum() )
    #print( data_all[i][0,15:15+6] )
    data_all[i] = data_all[i][1:,:]
    print( data_all[i].shape )
    print( '<%f: '%Tmin, len( data_all[i][ (data_all[i]>0) * (data_all[i]<Tmin) ] ) )
    print( '>%f: '%Tmax, len( data_all[i][ data_all[i]>Tmax ] ) )
    #data_all[i][ data_all[i]<Tmin ] = np.nan

print( Ls )
#print( Ms )

Tmins = []
Tmaxs = []
for i in range(N):
    TTmin = np.min( [ data_all[i][data_all[i]>Tmin].min(),\
            data_all[i+N][ data_all[i+N]>Tmin ].min()] )
    TTmax = np.max( [ data_all[i].max(), \
                data_all[i+N].max()]  )
    print( 'TTmin: %f, TTmax: %f\n'%( TTmin, TTmax ) )

    data_all[i][data_all[i]<TTmin] = np.nan
    data_all[i+N][data_all[i+N]<TTmin] = np.nan

    Tmins.append( TTmin )
    Tmaxs.append( TTmax )
    #data_all[i][0,0] = TTmin
    #data_all[i+N][0,0] = TTmin

    #data_all[i][0,1] = TTmax
    #data_all[i+N][0,1] = TTmax

print( Tmins, Tmaxs )
m, n = data1[0].shape

Tunit = 1e7
Tunitstr = '10^7'
font = FontProperties()
font.set_size( 'xx-large' )
for i in range( 2*N ):
    data_all[i] /= Tunit
    if i < N:
        Tmins[i] /= Tunit
        Tmaxs[i] /= Tunit

for i in range(2*N):

    #cmap = cm.jet
    #cmap = cmaps[i%N]
    cmap = cmaps[1]

    ax = ax_all[i]
    data = data_all[i]
    img = ax.imshow( data,  cmap=cmap, vmin=Tmins[i%N], vmax=Tmaxs[i%N] )
    #img = ax.imshow( data,  cmap=cmap )

    t = data[ np.logical_not(np.isnan( data ))]
    print( "min: %g, max: %g"%(t.min(), t.max()) )

    t = "%e"%Ms[i]
    #print( t )
    t = t.split( 'e' )
    t = r'$\rm M={%.2f} \times 10^{%i} \, h^{-1}M_\odot$'%(float(t[0]),float(t[1]) )
    #print( t )
    ax.text( n*0.1,m*0.1, t, fontproperties=font )

    if i < N:
        cbar = ax_cbar[i]
        cbar = plt.colorbar( img, cax = cbar, orientation='horizontal' )
        cbar.ax.tick_params( axis='x', direction='in', labelsize=20 )
        cbar.ax.set_xlabel( r'$\rm Temperature \, [\times %s\,K]$'%Tunitstr, fontsize=20 )

    if i < N:
        #cmap = cm.jet
        ax = ax_residul[i]
        data = (data_all[i+N]-data_all[i]) / data_all[i]
        cbar = ax_cbar_residul[i]
        t = data[ np.logical_not( np.isnan( data ) ) ]
        print( t.min(), t.max() )
        img = ax.imshow( data,  cmap=cmap, norm=mplc.SymLogNorm( linthresh=1e-1 ) )
        cbar = plt.colorbar( img, cax = cbar, orientation='horizontal' )
        cbar.ax.tick_params( axis='x', direction='in', labelsize=20 )
        cbar.ax.set_xlabel( r'$\rm Relative \,difference$', fontsize=20 )


fig.savefig( figs_dir + 'group-temp.pdf', dpi=300 )
fig_residul.savefig( figs_dir + 'group-temp-residul.pdf', dpi=300 )
#plt.show()

