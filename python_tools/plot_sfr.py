#!/usr/bin/env python3

from my_work_env import *

fn_sfr = sys.argv[1]
fn_sfr_cre = sys.argv[2]
fn_out = sys.argv[3]

fns = [ fn_sfr, fn_sfr_cre ]
ds = []
for f in fns:
    print( 'load %s ...'%f )
    ds.append( np.loadtxt(f) )

labels = [ 'SIM', 'SIM-CRE' ]
for i in range(2):
    d = ds[i]
    a = d[:,0]
    z = 1/a-1
    r = d[:,2]
    index = r>0
    z = z[index]
    r = r[index]
    plt.plot( z, r, label=labels[i] )

plt.yscale( 'log' )
plt.xlabel( r'$z$' )

plt.ylabel( r'$SFR({M_{\odot}}^{-1}{Mpc}^{-3})$' )
plt.savefig( fn_out )

