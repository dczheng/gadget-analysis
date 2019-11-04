#!/usr/bin/env python3

from my_work_env import *

fn_cr  = sys.argv[1]
fn_cre = sys.argv[2]
fn_o = sys.argv[3]

d_cr  = np.loadtxt( fn_cr )
d_cre = np.loadtxt( fn_cre )

x = d_cr[1:,0]
y = d_cr[1:,1]
#plt.step( x, y, label=r'$\frac{P_{\rm CRP}}{P_{\rm tot}}$' )
plt.plot( x, y, label=r'$\frac{P_{\rm CRP}}{P_{\rm tot}}$' )

x = d_cre[1:,0]
y = d_cre[1:,1]
#plt.step( x, y, label=r'$\frac{P_{\rm CRE}}{P_{\rm tot}}$' )
plt.plot( x, y, label=r'$\frac{P_{\rm CRE}}{P_{\rm tot}}$' )

plt.xscale( 'log' )
#plt.yscale( 'log' )
plt.legend( framealpha=0.1, fontsize=20 )
fs = 20
plt.xlabel( r'$\rm Pressure \, ratio$', fontsize=fs )
plt.ylabel( r'$\rm PDF$', fontsize=fs )
plt.xlim( [1e-8, 2] )

ax = plt.gca()
set_tick_params( ax, 15 )

plt.savefig( fn_o )
