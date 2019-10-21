#!/usr/bin/env python3

from my_work_env import *

ps_fn = sys.argv[1]
mf_fn = sys.argv[2]
cre_mf_fn =  sys.argv[3]
fn_out = sys.argv[4]

fs = 8
m = 3
n = 1
t_pad_x = 0.15
t_pad_y = 0.1
r = 0.2
dx = (1 - t_pad_x-t_pad_x*r) / n
dy = (1 - t_pad_y-t_pad_y*r) / m
fig = plt.figure( figsize=(fs, fs) )

axs = [ fig.add_axes( [t_pad_x, (m-1-i)*dy+t_pad_y, dx, dy] )\
        for i in range(m)]

ls = [ '-', '-.' ]
ss = [ '*', '^' ]
colors = [ 'k', 'r' ]

ps = np.loadtxt( ps_fn )
mf = np.loadtxt( mf_fn )
cre_mf = np.loadtxt( cre_mf_fn )

um = 1e2
#axs[0].plot( ps[:,0]/um, ps[:,1], label='PS' )

#index = (mf[:,-1]>100) * (cre_mf[:,-1]>100)
#mf = mf[index, :]
#cre_mf = cre_mf[index,:]

M = mf[:, 0] / um 
dndm = mf[:, 1]
N = mf[:, -1]

cre_M = cre_mf[:, 0] / um
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


axs[1].plot( M, (cre_dndm-dndm) / dndm )
axs[1].plot( M, (cre_dndm-dndm) / dndm, 'o' )
print( (cre_dndm-dndm) / dndm )

axs[2].step( M, N, '-.', label=r'$\rm SIM$')
axs[2].step( cre_M, cre_N, '-', label=r'$\rm SIM-CRE$')

for i in range( m ):
    axs[i].set_xscale( 'log' )
    axs[i].set_yscale( 'log' )

    if i == 0:
        axs[i].set_ylabel( r'$\rm MF \;[\frac{(h/Mpc)^3}{dlog_{10}(M/M_{\odot})}]$', fontsize=20 )
        axs[i].legend( prop={'size':20}, framealpha=0.1)
        axs[i].set_ylim( [dndm_min/2, dndm_max*2] )
        remove_tick_labels( axs[i], 'x' )
    if i == 1:
        axs[i].set_ylabel( r'$\rm Difference$', fontsize=20 )
        #axs[i].yaxis.set_major_locator(plt.LinearLocator(numticks=8))
        axs[i].set_yscale( 'linear' )
        #axs[i].yaxis.set_major_locator(plt.MultipleLocator(0.02))
        remove_tick_labels( axs[i], 'x' )

    if i == 2:
        axs[i].set_ylabel( r'$N$', fontsize=20 )
        axs[i].legend( prop={'size':20}, framealpha=0.1)
        axs[i].set_xlabel( r'$M [\times 10^{12}\,h^{-1}M_{\odot}]$', fontsize=20 )

    axs[i].tick_params( axis='both', direction='in', pad=5, labelsize=20 )
    if i == 1:
        axs[i].tick_params( axis='both', direction='in', pad=5, labelsize=15 )
    axs[i].minorticks_off()

fig.set_size_inches(8,8)

fig.savefig( fn_out )
