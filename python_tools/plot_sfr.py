#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
UnitLength_in_cm = 3.08568321
UnitTime_in_s = 3.08568e16
UnitMass_in_g = 1.989e43
solar_mass = 1.989e33
sec_per_year = 3.155e7
for fn in [ 'sfr', 'sfr_supernova', 'sfr_shock', 'sfr_supernova_shock' ]:
    sfr = np.loadtxt( './'+fn+'.txt' )
    a = sfr[:,0]
    z = 1/a - 1
    r = sfr[:,2]
    zz = []
    rr = []
    for i in range(len( z )):
        if 1.6 < z[i]<6:
            zz.append( z[i] )
            rr.append( r[i] )
    zz = np.array( zz )
    rr = np.array( rr )
    rr = rr / (10/0.7)**3
    rr = rr * ( UnitMass_in_g / solar_mass ) / ( UnitTime_in_s / sec_per_year )
    plt.plot( zz, rr, label=fn )
    plt.yscale( 'log' )

plt.xlabel( r'$z$' )
plt.ylabel( r'$SFR({M_{\odot}}^{-1}{Mpc}^{-3})$' )
plt.legend()
plt.savefig("sfr.png")

