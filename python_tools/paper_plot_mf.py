#!/usr/bin/env python3

from my_work_env import *
import matplotlib.gridspec as gridspec


ps_fn = sys.argv[1]
mf_fn = sys.argv[2]
cre_mf_fn =  sys.argv[3]
fn_out = sys.argv[4]

axN = 3
axs = []
gs = gridspec.GridSpec( axN, 1, height_ratios=[1,1,1] )
for i in range( axN ):
    axs.append( plt.subplot( gs[i] ))

data_dir = './'
ls = [ '-', '-.' ]
ss = [ '*', '^' ]
colors = [ 'k', 'r' ]

ps = np.loadtxt( ps_fn )
mf = np.loadtxt( mf_fn )
cre_mf = np.loadtxt( cre_mf_fn )

#axs[0].plot( ps[:,0], ps[:,1], label='PS' )

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

M = mf[:, 0]
dndm = mf[:, 1]
N = mf[:, -1]

cre_M = cre_mf[:, 0]
cre_dndm = cre_mf[:, 1]
cre_N = cre_mf[:, -1]

for i in range(len(M)):
    if M[i] != cre_M[i]:
        print( 'error' )
        exit()

index = (dndm>0) * (cre_dndm>0)
M = M[index]
dndm = dndm[index]
N = N[index]

cre_M = cre_M[index]
cre_dndm = cre_dndm[index]
cre_N = cre_N[index]

dndm_max = np.max( [ dndm.max(), cre_dndm.max() ] )
dndm_min = np.min( [ dndm.min(), cre_dndm.min() ] )

axs[0].plot( M, dndm, '-', label=r'$\rm SIM$' )
axs[0].plot( cre_M, cre_dndm, '--', label=r'$\rm SIM-CRE$' )
axs[0].plot( M, dndm, '*')
axs[0].plot( cre_M, cre_dndm, '.')

MM, NN = my_hist( M, N )
cre_MM, cre_NN = my_hist( cre_M, cre_N )

axs[2].plot( MM, NN, '-.', label=r'$\rm SIM$')
axs[2].plot( cre_MM, cre_NN, '-', label=r'$\rm SIM-CRE$')

axs[1].plot( M, (cre_dndm-dndm) / dndm )
axs[1].plot( M, (cre_dndm-dndm) / dndm, 'o' )
print( (cre_dndm-dndm) / dndm )

for i in range( axN ):
    axs[i].set_xscale( 'log' )
    axs[i].set_yscale( 'log' )
    #axs[i].grid()
    axs[i].set_xlim( [M.min()/xlim_fac, M.max()*xlim_fac] )

    if i == 0:
        #axs[i].set_ylabel( r'$\frac{dn} {dlog_{10}(M/M_{\odot})} \,\rm [(h/Mpc)^3]$', fontsize=20 )
        axs[i].set_ylabel( r'$\rm MF\;[\frac{(h/Mpc)^3}{dlog_{10}(M/M_{\odot})}]$', fontsize=20 )
        axs[i].legend( prop={'size':20}, framealpha=0.1)
        axs[i].set_ylim( [dndm_min/2, dndm_max*2] )
    if i == 1:
        axs[i].set_ylabel( r'$\rm Error$', fontsize=20 )
        #axs[i].yaxis.set_major_locator(plt.LinearLocator(numticks=8))
        axs[i].set_yscale( 'linear' )
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(0.02))

    if i == 2:
        axs[i].set_ylabel( r'$\rm Number$', fontsize=20 )
        axs[i].legend( prop={'size':20}, framealpha=0.1)
        axs[i].set_xlabel( r'$\rm M/M_{\odot}$', fontsize=20 )

    axs[i].tick_params( axis='both', direction='in', pad=5, labelsize=20 )
    if i == 1:
        axs[i].tick_params( axis='both', direction='in', pad=5, labelsize=15 )
    axs[i].minorticks_off()

fig = plt.gcf()
fig.set_size_inches(8,8)
fig.tight_layout()

fig.savefig( fn_out )
