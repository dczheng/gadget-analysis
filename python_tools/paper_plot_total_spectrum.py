#!/usr/bin/env python3

from my_work_env import *
from scipy.optimize import curve_fit

print( "plot total radiation spectrum ..." )

rad_spec_file1 = sys.argv[1]
fn_out = sys.argv[2]

dat_1 = np.loadtxt( rad_spec_file1 )

fig = plt.figure()
ax = fig.add_subplot( 111 )

def fit_f( x, a, b ):
    return x*a + b

x = dat_1[0,1:]
y = dat_1[-1,1:] / mycc.mJy
print( x )
print( y )

xx = np.log10( x )
yy = np.log10( y )
r = curve_fit( fit_f, xx, yy )

print( "alpha: %.2f"%r[0][0] )
ax.loglog( x, y , label=r'$\alpha=%.2f$'%(-r[0][0]) )
ax.set_xlabel( r'$\rm \nu \, [MHz]$', fontsize=20 )
ax.set_ylabel( r'$\rm I_{\nu}\; \rm [mJy\,sr^{-1}]$', fontsize=20 )
ax.legend( framealpha=0.1, prop={'size':20} )
#ax.grid()
ax.tick_params( axis='both', direction='in', which='both', labelsize=20, pad=5 )
#ax.minorticks_off()

if ( len(sys.argv) < 2 ):
    file_pre = ''
else:
    file_pre = sys.argv[1] + '-'

fig.tight_layout()
fig.savefig( fn_out )
