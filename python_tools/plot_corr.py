#!/usr/bin/env python3

from  my_work_env import *

def myplot( ax, dcre, d ):
    ax.plot( dcre[:,0]/1000,dcre[:,1], '.', label='cre' )
    ax.plot( d[:,0]/1000, d[:,1], '.', label='no' )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    #ax.set_yscale( 'Symlog' )
    ax.set_xlabel( 'Mpc' )
    ax.set_ylabel( r'$\xi(r)$' )
    #ylim = list( ax.get_ylim() )
    #ylim[0] = 1e-3
    #ax.set_ylim( ylim )

    ax.legend()

fig, axs = plt.subplots( 2, 1, figsize=(6, 2*4) )
dmax = axs[0]
gax = axs[1]

dmcorr_cre = np.loadtxt( data_dir + 'DMCorr1d_000.dat' )
dmcorr = np.loadtxt( data_dir + 'DMCorr1d_001.dat' )

gcorr_cre = np.loadtxt( data_dir + 'GasCorr1d_000.dat' )
gcorr = np.loadtxt( data_dir + 'GasCorr1d_001.dat' )

myplot( dmax, dmcorr_cre, dmcorr )
myplot( gax, gcorr_cre, gcorr )

fig.tight_layout()
fig.savefig( figs_dir + 'corr.pdf' )

