#!/usr/bin/env python3

from my_work_env import *

fs = 5
t = 0.4
fig = plt.figure( figsize=(fs, fs/(1-t)) )
r = 0.7
s = 0.03
dax = fig.add_axes( [(1-r)/2, t+(1-t)*(1-r)/2, r, (1-t)*r] )
Tax = dax.twiny()
daxerr = fig.add_axes( [(1-r)/2, t*(1-r)/2+s, r, t*r] )
Taxerr = daxerr.twiny()

dpdf_cre = np.loadtxt( data_dir + './DensPdf_000.dat' )
dpdf = np.loadtxt( data_dir + './DensPdf_001.dat' )

Tpdf_cre = np.loadtxt( data_dir + './TPdf_000.dat' )
Tpdf = np.loadtxt( data_dir + './TPdf_001.dat' )

cs = [ 'b-','r--', 'g-', 'c--' ]
d = [ dpdf, dpdf_cre, Tpdf, Tpdf_cre ]
labels = [\
r'$\rho$', r'$\rho_{\rm CRE}$',\
r'$T$', r'$T_{\rm CRE}$'\
]
axs = [ dax, dax, Tax, Tax ]

for i in range(4):
    axs[i].plot( d[i][:,0], d[i][:,1], cs[i], label=labels[i] )

legend_loc = [ 'lower left', 'upper right' ]
xlabels = [ r'$\rho / \rho_{\rm bar}$', r'$T\,[K]$' ]
ylabels = [ r'$dN/d{\rm log_{10}}\rho \, [dN/d{\rm log_{10}}T]$' ]
axs = [ dax, Tax ]
for i in range(2):
    axs[i].set_xlabel( xlabels[i], fontsize=15 )
    axs[i].set_xscale( 'log' )
    axs[i].set_yscale( 'log' )
    axs[i].tick_params( axis='both', direction='in', pad=3 )
    axs[i].minorticks_off()
    axs[i].legend( loc=legend_loc[i], framealpha=0.1 )
    #axs[i].grid()
    if i == 0:
        axs[i].set_ylabel( ylabels[i] )

axs = [ daxerr, Taxerr ]
cs = [ 'b-','k-' ]
d = [ dpdf, dpdf_cre, Tpdf, Tpdf_cre ]
labels = [ r'$\rho$', r'$T$' ]
xlabels = [ r'$\rho / \rho_{\rm bar}$', r'$T\,[K]$' ]
legend_loc = [ 'lower right', 'upper left' ]
for i in range(2):
    axs[i].plot( d[i*2][:,0], (d[i*2+1][:,1]-d[i*2][:,1])/d[i*2][:,1],\
    cs[i], label=labels[i] )
    axs[i].minorticks_off()
    axs[i].set_xscale( 'log' )
    #axs[i].set_yscale( 'log' )
    axs[i].set_xlabel( xlabels[i], fontsize=15 )
    axs[i].tick_params( axis='both', direction='in', pad=3 )
    axs[i].minorticks_off()
    axs[i].legend( loc=legend_loc[i], framealpha=0.1 )
    if i == 0:
        axs[i].set_ylabel( r'$\rm Relative \,Diffeernce$' )

fig.savefig( figs_dir + './DensPdf_TPdf.pdf' )
