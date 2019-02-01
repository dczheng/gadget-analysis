#!/usr/bin/env python3

from scipy.optimize import curve_fit

from my_work_env import *

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

def f( x, a, b ):
    return x*a+b

'''

fname = [
            './hge256_data/hge256_0.01_0.05_Spec_Tot.dat', \
            './hge256_data/hge256_0.005_0.05_Spec_Tot.dat', \
'./hge256_data/hge256_0.001_0.05_Spec_Tot.dat' \
        ]
label = [ 'CRE100\_01', 'CRE100\_005', 'CRE100\_01' ]
'''
ls = [ '-', '-.', '--' ]
fname = [ sys.argv[1] ]
label = [ 'x' ]

for i in range(len(label)):
    d = np.loadtxt( fname[i] )
    v = d[:,0]
    p = d[:, 1]
    logv = np.log10(v)
    logp = np.log10(p)
    #args, cov = curve_fit( f, logv, logp )
    #a = args[0]

    a = (logp[1]-logp[-1]) /(logv[1] - logv[-1])
    print( a )
    plt.plot( v, p, ls=ls[i], label=r'$\rm %s \quad \alpha=%.2f$'%(label[i], a) )


plt.xscale( 'log' )
plt.yscale( 'log' )

plt.xlabel( r'$\nu \; [\rm MHz]$', fontsize=20 )
#plt.ylabel( r'$I_{\nu} \; \rm [Jy/sr]$' )
plt.ylabel( r'$ I_{\nu}\; \rm [erg \, cm^{-2} \, sr^{-1} \, Hz^{-1} ]$', fontsize=20 )

#plt.title( "Radio Spectrum" )

ax = plt.gca()
ax.tick_params( axis='both', pad=5, direction='in', labelsize=20 )
ax.minorticks_off()

plt.legend(prop={'size':15}, framealpha=0.1)
plt.grid()

plt.tight_layout()

plt.savefig( output_dir + 'Tspec.pdf', figsize=(4,4) )

