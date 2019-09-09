#!/usr/bin/env python3

from my_work_env import *

print( "plot total radiation spectrum ..." )

rad_spec_file1 = sys.argv[1]
fn_out = sys.argv[2]

dat_1 = np.loadtxt( rad_spec_file1 )

fig = plt.figure()
ax = fig.add_subplot( 111 )

x = dat_1[:,0]
y = dat_1[:,1] / mycc.mJy

ax.loglog( x, y ) #label=r'$S2$' )
ax.set_xlabel( r'$\rm \nu \, [MHz]$', fontsize=20 )
ax.set_ylabel( r'$\rm I_{\nu}\; \rm [mJy\,sr^{-1}]$', fontsize=20 )
#ax.legend( framealpha=0.1 )
#ax.grid()
ax.tick_params( axis='both', direction='in', which='both', labelsize=20, pad=5 )
#ax.minorticks_off()

if ( len(sys.argv) < 2 ):
    file_pre = ''
else:
    file_pre = sys.argv[1] + '-'

fig.tight_layout()
fig.savefig( fn_out )
