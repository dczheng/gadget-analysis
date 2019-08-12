#!/usr/bin/env python3

import h5py
from my_work_env import *

f = h5py.File( sys.argv[1], 'r' )
T = f[ '/PartType0/Temperature' ].value
dens = f[ '/PartType0/Density' ].value
m = f[ '/PartType0/Masses' ].value
z= f[ '/Header' ].attrs[ 'Redshift' ]

RhoBaryon = mycc.rho_bar_in_gadget( z )
dens  = dens / RhoBaryon

index_dens = dens >= 1e3 
index_t = dens < 1e3
index_cool = ( T < 1e5 ) * index_t
index_warm = ( 1e5 <= T ) * ( T < 1e7  ) * index_t
index_hot =  (T >= 1e7) * index_t

m_tot  = m.sum()
m_dens = m[ index_dens ].sum()
m_cool = m[ index_cool ].sum()
m_warm = m[ index_warm ].sum()
m_hot  = m[ index_hot  ].sum()

r_dens = m_dens / m_tot
r_cool = m_cool / m_tot
r_warm = m_warm / m_tot
r_hot  = m_hot  / m_tot

print( 'cool: %g, warm: %g, hot: %g, dens: %g'\
        %( r_cool, r_warm, r_hot, r_dens ) )
