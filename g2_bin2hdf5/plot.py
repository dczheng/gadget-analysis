#!/usr/bin/env python3
import yt
#ds = yt.load( './snapshot_005.hdf5' )
ds = yt.load( './out.hdf5' )
p = yt.ParticlePlot(ds,'particle_position_x','particle_position_y','particle_mass')
p.save()
