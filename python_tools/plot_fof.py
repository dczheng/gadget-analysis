#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import h5py
from matplotlib.patches import Circle

file = h5py.File( './fof.hdf5', 'r' )
fof_pos = file[ '/CenterOfMass' ].value
fof_radius = file[ '/VirialR200' ].value
file.close()

file = h5py.File( './test_snapshot_064.hdf5', 'r' )
particle_pos = file[ '/PartType1/Coordinates' ].value
file.close()

fig_size = 512
fig = plt.figure()
ax = fig.add_subplot( 111 )
#ax.plot( particle_pos[:,0], particle_pos[:,1], 'r.', markersize=0.005 );
xmin = np.min( particle_pos[:,0] )
xmax = np.max( particle_pos[:,0] )
ymin = np.min( particle_pos[:,1] )
ymax = np.max( particle_pos[:,1] )
delta_x = ( xmax-xmin ) / fig_size
delta_y = ( ymax-ymin ) / fig_size
fig_mat = np.zeros( (fig_size, fig_size ) )
x = particle_pos[:,0]
y = particle_pos[:,1]
x -= xmin
x = x / delta_x
x[ x>fig_size-1 ] = fig_size-1
x = x.astype( int )
#print( x )
y -= ymin
y = y / delta_y
y[ y>fig_size-1 ] = fig_size-1
y = y.astype( int )
#print( y )
for i in range( len( particle_pos ) ):
    fig_mat[ x[i], y[i] ] = fig_mat[ x[i], y[i] ] + 1
#print( fig_mat )
fig_mat[ fig_mat == 0 ] = 1e-10
fig_mat = np.log10( fig_mat )
fig_mat = np.transpose( fig_mat )
ax.imshow( fig_mat )

fof_x = fof_pos[ :,0 ]
fof_y = fof_pos[ :,1 ]
fof_r = fof_radius
print( fof_x.max(), fof_x.max() )
print( fof_y.max(), fof_y.min() )
print( fof_r.max(), fof_r.min() )
fof_x = ( fof_x - xmin ) / ( xmax-xmin ) * fig_size
fof_y = ( fof_y - ymin ) / ( ymax-ymin ) * fig_size
fof_r = ( fof_r - xmin ) / ( xmax-xmin ) * fig_size
print( fof_x.max(), fof_x.max() )
print( fof_y.max(), fof_y.min() )
print( fof_r.max(), fof_r.min() )
ax.plot( fof_x, fof_y, 'w.', markersize=1 );
for i in range( len( fof_pos ) ):
    cir = Circle( xy=(fof_x[i], fof_y[i]), radius=fof_r[i], alpha=0.5 );
    cir.set_edgecolor( [1,1,1,1] )
    ax.add_patch( cir )
'''
ax = fig.add_subplot( 111, projection='3d' )
ax.scatter( particle_pos[:,0], particle_pos[:,1], particle_pos[:,2], s=0.01 )
ax.scatter( fof_pos[:,0], fof_pos[:,1], fof_pos[:,2], s=250 )
'''
plt.axis( 'scaled' )
plt.axis( 'equal' )
plt.savefig( 'test.png' )

