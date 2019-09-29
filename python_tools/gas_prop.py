#!/usr/bin/env python3

import h5py
from my_work_env import *

args = len( sys.argv )

sep_str =  '-' * 30 

if ( args == 1 ):
    print( "please give the command line arguments." )
    exit()

fn_pre = sys.argv[1]

fn_test = [ "%s.hdf5"%fn_pre, "%s.0.hdf5"%fn_pre ]
flag = -1
for i in range(2):
    fn = fn_test[i]
    if os.path.exists( fn ):
        flag = i
        f = h5py.File( fn, 'r' )
        fn_num = f[ '/Header' ].attrs[ 'NumFilesPerSnapshot' ]
        num_high = f[ '/Header' ].attrs[ 'NumPart_Total_HighWord' ]
        num_low = f[ '/Header' ].attrs[ 'NumPart_Total' ]
        num = num_high * 2**32 + num_low
        Redshift = f[ '/Header' ].attrs[ 'Redshift' ]
        f.close()

if flag == -1:
    print( "Error" )
    exit()

print( "NumPart: ", num )
print( "NumFilesPerSnapshot: %i"%fn_num )
print( "Redshift: %.2f"%Redshift )

gas_density = []
gas_m = []

read_temp = 1

if read_temp:
    gas_temp = []
else:
    gas_u = []

for i in range( fn_num ):

    if fn_num == 1:
        fn = "%s.hdf5"%fn_pre
    else:
        fn = "%s.%i.hdf5"%( fn_pre, i )

    print( 'read %s ...'%fn  )
    f = h5py.File( fn, 'r' )

    d = f[ '/PartType0/Density' ][()]
    gas_density.append( d )

    if read_temp:
        d = f[ '/PartType0/Temperature' ][()]
        gas_temp.append( d )
    else:
        d = f[ '/PartType0/InternalEnergy' ][()]
        gas_u.append( d )

    d = f[ '/PartType0/Masses' ][()]
    gas_m.append( d )

    f.close()

print( "stack data ..." )
gas_density = np.hstack( gas_density )
gas_m = np.hstack( gas_m )
if read_temp:
    gas_temp = np.hstack( gas_temp )
else:
    gas_u = np.hstack( gas_u )

print( sep_str )

print( "Density number: ", gas_density.shape[0] )

gas_density =  gas_density / mycc.rho_bar_in_gadget( Redshift )

if not read_temp:
    gas_u = gas_u * mycc.gadget_energy_in_erg / mycc.gadget_mass_in_g
    XH = mycc.Xh 
    yh = ( 1 - XH ) / ( 4*XH )
    ne = 1 + 2 * yh
    mu = ( 1 + 4 * yh ) / ( 1 + yh + ne )
    gas_temp = ( mycc.Gamma-1 ) / mycc.k_b * gas_u * mycc.m_p * mu

Tmin = gas_temp.min()
Tmax = gas_temp.max()
print( "Tmin: %g, Tmax: %g"%( Tmin, Tmax ) )
Dmin = gas_density.min()
Dmax = gas_density.max()
print( "Dmin: %g, Dmax: %g"%( Dmin, Dmax ) )

print( "compute gas ratio ..." )
index_dens = gas_density >= 1e3
index_t = gas_density < 1e3
index_cool = ( gas_temp < 1e5 ) * index_t
index_warm = ( gas_temp >= 1e5 ) * ( gas_temp < 1e7 ) * index_t
index_hot =  ( gas_temp >= 1e7 ) * index_t

m_tot = gas_m.sum()
m_dens = gas_m[ index_dens ].sum()
m_cool = gas_m[ index_cool ].sum()
m_warm = gas_m[ index_warm ].sum()
m_hot  = gas_m[ index_hot  ].sum()

r_dens = m_dens / m_tot
r_cool = m_cool / m_tot
r_warm = m_warm / m_tot
r_hot  = m_hot  / m_tot

print( 'cool: %g  warm: %g  hot: %g  dens: %g'\
        %( r_cool, r_warm, r_hot, r_dens ) )

exit()
fn_out = sys.argv[2]
print( "plot phase ..." )

N = 128
img = np.zeros( [N, N] )
dlogT = np.log10( Tmax/Tmin ) / ( N-1 )
dlogD = np.log10( Dmax/Dmin ) / ( N-1 )

temp_idx = (np.log10( gas_temp / Tmin ) / dlogT ).astype( 'int' )
d_idx    = (np.log10( gas_density / Dmin ) / dlogD ).astype( 'int' )

for i in range( len(temp_idx) ):
    img[ temp_idx[i], d_idx[i] ] += 1

img /= dlogT * dlogD

t_pad = 0.1
t_cbar = 0.05
fs = 4
fig = plt.figure( figsize=(fs/(1-2*t_pad-t_cbar), fs/(1-2*t_pad)) )
dx = 1 - 2*t_pad - t_cbar
dy = 1 - 2*t_pad

ax = fig.add_axes( [ t_pad, t_pad, dx, dy ] )
cax = fig.add_axes( [ t_pad+dx, t_pad, t_cbar, dy ] )

img = ax.imshow( img, norm=mplc.LogNorm(), cmap=cm.jet )
plt.colorbar( img, cax=cax )

make_log_ticks( Tmin, Tmax, N, axis=ax.yaxis )
make_log_ticks( Dmin, Dmax, N, axis=ax.xaxis )
ax.grid()
ax.invert_yaxis()
set_tick_params( ax, 15 )
set_tick_params( cax, 15 )
ax.set_xlabel( r'$\rho/\rho_{bar}$', fontsize=10 )
ax.set_ylabel( r'$T\,[K]$', fontsize=10 )

plt.savefig( fn_out )



