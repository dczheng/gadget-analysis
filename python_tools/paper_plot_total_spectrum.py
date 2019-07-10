#!/usr/bin/env python3

from my_work_env import *
import tools_and_constants as tc
from scipy.optimize import curve_fit

print( "plot total radiation spectrum ..." )

rad_spec_file1 = data_dir + './Spec_Tot.dat'

dat_1 = np.loadtxt( rad_spec_file1 )

fig = plt.figure()
ax = fig.add_subplot( 111 )

x = dat_1[:,0]
y = dat_1[:,1]

ax.loglog( x, y ) #label=r'$S2$' )
ax.set_xlabel( r'$\nu \, [\rm MHz]$', fontsize=20 )
ax.set_ylabel( r'$ I_{\nu}\; \rm [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1} ]$', fontsize=20 )
#ax.legend( framealpha=0.1 )
#ax.grid()
ax.tick_params( axis='both', direction='in', which='both', labelsize=20, pad=5 )
#ax.minorticks_off()

if ( len(sys.argv) < 2 ):
    file_pre = ''
else:
    file_pre = sys.argv[1] + '-'

fig.tight_layout()
fig.savefig( figs_dir + file_pre + 'rad_tot_spec.pdf' )
