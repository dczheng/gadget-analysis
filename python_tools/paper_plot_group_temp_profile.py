#!/usr/bin/env python3

from my_work_env import *

FileNames1 = [
        data_dir + '/g_temp_profile_0.dat', \
        data_dir + '/g_temp_profile_1.dat', \
        data_dir + '/g_temp_profile_2.dat', \
        data_dir + '/g_temp_profile_3.dat' \
                ]

FileNames2 = [
        data_dir + '/cre_g_temp_profile_0.dat', \
        data_dir + '/cre_g_temp_profile_1.dat', \
        data_dir + '/cre_g_temp_profile_2.dat', \
        data_dir + '/cre_g_temp_profile_3.dat' \
                ]

cs = [ 'k', 'b', 'm', 'r', 'g', 'c']

N = len( FileNames1 )


Tunit = 1e7

fig, axs = plt.subplots( 1, N, figsize=( N*4, 4 ) )

for i in range( N ):
    d1 = np.loadtxt( FileNames1[i] )
    d2 = np.loadtxt( FileNames2[i] )

    c = cs[ i%len(cs) ]
    c1 = cs[ (i+1)%len(cs) ]
    c2 = cs[ (i+2)%len(cs) ]
    axs[i].errorbar(d1[:,0],d1[:,1]/Tunit,yerr=d1[:,2]/Tunit, \
            fmt='o',ecolor=c1, color=c, elinewidth=2,capsize=4, label='%i'%i)

    axs[i].errorbar(d2[:,0],d2[:,1]/Tunit,yerr=d2[:,2]/Tunit, \
            fmt='*',ecolor=c2, color=c, elinewidth=2,capsize=4, label='%i-CRE'%i)

    axs[i].set_xscale( 'log' )

    axs[i].set_xlabel( r'$\rm R\, [h^{-1}\, Mpc]$' )
    axs[i].set_ylabel( r'$\rm Temperature\, [10^7\,K]$' )
    axs[i].legend()



fig.savefig( figs_dir + 'group-temp-profile.pdf', dpi=300 )
#plt.show()

