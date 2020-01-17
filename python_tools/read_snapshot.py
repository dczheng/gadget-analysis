#!/usr/bin/env python3

import h5py
from my_work_env import *
import pandas as pd

snap_prefix = sys.argv[1]
out_prefix = sys.argv[2]

sep_str =  '-' * 30 


print( sep_str )
fn_test = [ "%s.hdf5"%snap_prefix, "%s.0.hdf5"%snap_prefix ]
flag = -1
for i in range(2):
    fn = fn_test[i]
    if os.path.exists( fn ):
        flag = i
        f = h5py.File( fn, 'r' )
        fn_num = f[ '/Header' ].attrs[ 'NumFilesPerSnapshot' ]
        num_high = f[ '/Header' ].attrs[ 'NumPart_Total_HighWord' ]
        num_low = f[ '/Header' ].attrs[ 'NumPart_Total' ]
        BoxSize = f[ '/Header' ].attrs[ 'BoxSize' ]
        num = num_high * 2**32 + num_low
        Redshift = f[ '/Header' ].attrs[ 'Redshift' ]
        f.close()

if flag == -1:
    print( "Error" )
    exit()

print( "NumPart: ", num )
print( "NumFilesPerSnapshot: %i"%fn_num )
print( "Redshift: %.2f"%Redshift )

gas_data = { \
'Density': [],\
'Temperature': [],\
'InternalEnergy': [],\
'Masses': [],\
'Coordinates':[],\
'Velocities':[],\
}

dm_data = { \
'Coordinates':[],\
'Velocities':[],\
}

read_temp = 0

print( sep_str )
for i in range( fn_num ):

    if fn_num == 1:
        fn = "%s.hdf5"%snap_prefix
    else:
        fn = "%s.%i.hdf5"%( snap_prefix, i )

    print( 'read %s'%fn  )
    f = h5py.File( fn, 'r' )

    for k in gas_data.keys():
        if not read_temp and k == 'Temperature':
            continue
        print( 'load gas `%s`'%(k) )
        d = f[ '/PartType0/%s'%k ][()]
        gas_data[k].append( d )

    for k in dm_data.keys():
        print( 'load dark matter `%s`'%(k) )
        d = f[ '/PartType1/%s'%k ][()]
        dm_data[k].append( d )
    f.close()

print( sep_str )
for k in gas_data.keys():
    if not read_temp and k == 'Temperature':
        continue
    print( "stack gas `%s`"%k )
    if k == 'Coordinates':
        gas_data[k] = np.vstack( gas_data[k] )
        continue
    if k == 'Velocities':
        gas_data[k] = np.vstack( gas_data[k] )
        continue
    gas_data[k] = np.hstack( gas_data[k] )

for k in dm_data.keys():
    if not read_temp and k == 'Temperature':
        continue
    print( "stack dark matter `%s`"%k )
    if k == 'Coordinates':
        dm_data[k] = np.vstack( dm_data[k] )
        continue
    if k == 'Velocities':
        dm_data[k] = np.vstack( dm_data[k] )
        continue
    dm_data[k] = np.hstack( dm_data[k] )

print( sep_str )

if not read_temp:
    print( 'compute temperature ...' )
    gas_u = gas_data['InternalEnergy'] *\
                mycc.gadget_energy_in_erg / mycc.gadget_mass_in_g
    gas_density =  gas_data['Density'] / mycc.rho_bar_in_gadget( Redshift )
    XH = mycc.Xh 
    yh = ( 1 - XH ) / ( 4*XH )
    ne = 1 + 2 * yh
    mu = ( 1 + 4 * yh ) / ( 1 + yh + ne )
    gas_temp = ( mycc.Gamma-1 ) / mycc.k_b * gas_u * mycc.m_p * mu
    gas_data['Temperature'] = gas_temp

'''
N = 512
dL = BoxSize / N
x = gas_data['Coordinates'][:,0] / dL
y = gas_data['Coordinates'][:,1] / dL
x = x.astype( 'int' )
y = y.astype( 'int' )
print( x[ x<0 ] )
print( x[ x>N-1 ] )
print( y[ y<0 ] )
print( y[ y>N-1 ] )

x[ x<0 ] = 0
x[ x>N-1 ] = N-1
y[ y<0 ] = 0
y[ y>N-1 ] = N-1

img = np.zeros( [N, N] )
for i in range( len(x) ):
    img[x[i], y[i]] += 1

plt.imshow( img, norm=mplc.LogNorm() )
plt.grid()
plt.savefig( 'x.png' )
exit()
'''

print( sep_str )
print( 'save gas data ...' )
fn = out_prefix + '_gas.csv'
data = {\
    'X': gas_data['Coordinates'][:,0],\
    'Y': gas_data['Coordinates'][:,1],\
    'Z': gas_data['Coordinates'][:,2],\
    'Vx': gas_data['Velocities'][:,0],\
    'Vy': gas_data['Velocities'][:,1],\
    'Vz': gas_data['Velocities'][:,2],\
    'Temperature': gas_data['Temperature'],\
    'Density': gas_data['Density'],\
}

df = pd.DataFrame( data )
df.to_csv( fn, index=False )

print( sep_str )
print( 'save dark matter data ...' )
fn = out_prefix + '_dm.csv'
data = {\
    'X': dm_data['Coordinates'][:,0],\
    'Y': dm_data['Coordinates'][:,1],\
    'Z': dm_data['Coordinates'][:,2],\
    'Vx': dm_data['Velocities'][:,0],\
    'Vy': dm_data['Velocities'][:,1],\
    'Vz': dm_data['Velocities'][:,2],\
}

df = pd.DataFrame( data )
df.to_csv( fn, index=False )
