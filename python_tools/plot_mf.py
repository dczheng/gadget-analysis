#!/usr/bin/env python3
import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

axN = 3

axs = []
gs = gridspec.GridSpec( axN, 1, height_ratios=[1,1,1] )
for i in range( axN ):
    axs.append( plt.subplot( gs[i] ))

data_dir = './hge256_data/'
out_fn = './hge256_out/hge256_mf.pdf'
#out_fn = './hge256_out/hge256_mf.png'
simu_prefix = data_dir + 'hge256'
label_prefix = "CRE256_"
#zeta_list = [ 0.01, 0.001 ]
zeta_list = [ 0.01 ]
ls = [ '-', '-.' ]
ss = [ '*', '^' ]
colors = [ 'k', 'r' ]
chi = 0.05
ps_fn = simu_prefix + '_PS.dat'

ps = np.loadtxt( ps_fn )
axs[0].plot( ps[:,0], ps[:,1], label='PS' )

xlim_fac = 2

def my_hist( M, N ):
    MM = []
    NN = []
    for k in range( len( M ) ):

        if k == 0:
            NN.append( N[k] )
            MM.append( M[k] )
            continue
        MM.append( M[k] )
        MM.append( M[k] )
        if k == len(M)-1:
            NN.append( N[k-1] )
            NN.append( N.min() )
            continue
        NN.append( N[k-1] )
        NN.append( N[k] )
    MM.append( M.max() * xlim_fac )
    NN.append( N[-1] )

    return ( MM, NN )

simu_no_fn = simu_prefix + '_MF.dat'
simu_no = np.loadtxt( simu_no_fn )
N_no = simu_no[:, -1]
dndm_no = simu_no[:, 1]
M = simu_no[:, 0]
axs[0].plot( simu_no[:,0], simu_no[:,1], '.', label='CRE100')

MM, NN = my_hist( M, N_no )
axs[2].plot( MM, NN, '-.', label='CRE100')

for i in range( len(zeta_list) ):

    zeta = zeta_list[i]
    simu_fn = simu_prefix + '_%g_%g_MF.dat'%(zeta, chi)

    simu = np.loadtxt( simu_fn )

    str_zeta = str( zeta )
    str_zeta = str_zeta.split( '.' )
    dndm = simu[:,1]
    N = simu[:,-1]

    axs[0].plot( M, dndm, ss[i], \
            label='CRE100' + '\_%s'%str_zeta[1])# color=colors[i] )

    #axs[2].plot( simu[:,0], np.abs(simu[:,1]-simu_no[:,1]) / \
    #        simu_no[:,1], label='CRE100' + '\_%s'%str_zeta[1], \
    #        ls=ls[i] )

    MM, NN = my_hist( M, N )
    axs[2].plot( MM, NN, ls=ls[i], label='CRE100' + '\_%s'%str_zeta[1])# color=colors[i] )

    axs[1].plot( M, (dndm-dndm_no) / dndm_no, ls=ls[i])# color=colors[i] )



for i in range( axN ):
    axs[i].set_xscale( 'log' )
    axs[i].set_yscale( 'log' )
    axs[i].grid()
    axs[i].set_xlim( [M.min()/xlim_fac, M.max()*xlim_fac] )

    if i == 2:
        axs[i].set_ylabel( 'N', fontsize=20 )
        axs[i].legend( prop={'size':20}, framealpha=0.1)
        axs[i].set_xlabel( r'$M/M_{\odot}$', fontsize=20 )
    if i == 0:
        #axs[i].set_ylabel( r'$\frac{dn} {dlog_{10}(M/M_{\odot})} \,\rm [(h/Mpc)^3]$', fontsize=20 )
        axs[i].set_ylabel( r'$\rm MF\;[\frac{(h/Mpc)^3}{dlog_10(M/M_{\odot})}]$', fontsize=20 )
        axs[i].tick_params( axis='x', which='both',bottom=False, top=False,  labelbottom=False)
        axs[i].legend( prop={'size':20}, framealpha=0.1)
    if i == 1:
        axs[i].set_ylabel( r'$\frac{\Delta{\rm MF}}{\rm MF_{CRE100}}$', fontsize=20 )
        #axs[i].set_ylabel( r'${\Delta{\rm MF}}/{\rm MF_{CRE100}}$', fontsize=20 )
        axs[i].minorticks_off()
        #axs[i].yaxis.set_major_locator(plt.LinearLocator(numticks=8))
        axs[i].set_yscale( 'linear' )
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(0.02))
        axs[i].tick_params( axis='x', which='both',bottom=False, top=False,  labelbottom=False)
    axs[i].tick_params( axis='both', direction='in', pad=5, labelsize=20 )
    axs[i].minorticks_off()

fig = plt.gcf()
fig.set_size_inches(8,8)
fig.tight_layout()

fig.savefig( out_fn )
